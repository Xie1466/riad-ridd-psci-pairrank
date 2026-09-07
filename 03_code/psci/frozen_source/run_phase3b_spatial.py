#!/usr/bin/env python3
"""Run restartable Phase 3B GROMACS distance and sparse-contact calculations."""

from __future__ import annotations
import argparse,csv,datetime as dt,importlib.util,json,os,re,shlex,subprocess,sys,tempfile,time
from pathlib import Path
import yaml

GMX=Path(os.environ.get('GMX_BIN','gmx'));VERSION='1.0.0'
FIELDS=['topology_id','system_id','canonical_system_id','trajectory','analysis_index','status','classification','started_utc','ended_utc','wall_seconds','output_bytes','gromacs_version','script_version','notes']
def args():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',default='config/analysis.yaml');p.add_argument('--simulation-map',default='metadata/simulation_file_map.csv');p.add_argument('--trajectory-manifest',default='metadata/analysis_trajectory_manifest.csv');p.add_argument('--component-evidence',default='metadata/component_mapping_evidence.csv');p.add_argument('--component-map',default='metadata/component_group_map.csv');p.add_argument('--manifest',default='metadata/spatial_run_manifest.csv');p.add_argument('--only',action='append');p.add_argument('--dry-run',action='store_true');return p.parse_args()
def read(p):
 with Path(p).open(newline='',encoding='utf-8') as h:return list(csv.DictReader(h))
def utc():return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace('+00:00','Z')
def atomic_csv(path,rows,fields):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);fd,tmp=tempfile.mkstemp(prefix=f'.{path.name}.',dir=path.parent,text=True)
 try:
  with os.fdopen(fd,'w',newline='',encoding='utf-8') as h:w=csv.DictWriter(h,fieldnames=fields,lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
  os.replace(tmp,path)
 finally:
  if os.path.exists(tmp):os.unlink(tmp)
def load_component_builder():
 p=Path('scripts/components/build_component_maps.py');s=importlib.util.spec_from_file_location('component_builder',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def parse_ndx(path):
 groups={};cur=None
 for line in Path(path).read_text().splitlines():
  if line.strip().startswith('['):cur=line.strip().strip('[] ').strip();groups[cur]=[]
  elif cur and line.strip():groups[cur].extend(int(x) for x in line.split())
 return groups
def write_ndx(path,groups):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);fd,tmp=tempfile.mkstemp(prefix=f'.{path.name}.',dir=path.parent,text=True)
 try:
  with os.fdopen(fd,'w',encoding='utf-8') as h:
   for name,ids in groups.items():
    h.write(f'[ {name} ]\n')
    for i in range(0,len(ids),15):h.write(' '.join(f'{x:6d}' for x in ids[i:i+15])+'\n')
  os.replace(tmp,path)
 finally:
  if os.path.exists(tmp):os.unlink(tmp)
def atom_metadata(sim,builder):
 info={};offset=0
 for chain,itp in enumerate(sim['itp_candidates'].split(';'),1):
  _,_,atoms,_=builder.parse_itp(itp)
  for a in atoms:
   # Re-read mass from the source atom line because the shared parser intentionally stores structural fields only.
   info[offset+a['atomnr']]={'chain':chain,'resnr':a['resnr'],'resname':a['resname'],'atomname':a['atomname'],'mass':None}
  section=''
  for line in Path(itp).read_text().splitlines():
   text=line.split(';',1)[0].strip()
   if text.startswith('['):section=re.sub(r'[\[\]\s]','',text)
   elif section=='atoms' and text:
    f=text.split();info[offset+int(f[0])]['mass']=float(f[7]) if len(f)>7 else 99.0
  offset+=len(atoms)
 return info
def run_logged(log,cmd,version,stdin=''):
 if log.exists():log=log.with_name(f'{log.stem}_retry_{int(time.time())}{log.suffix}')
 log.parent.mkdir(parents=True,exist_ok=True);start=time.monotonic();started=utc()
 with log.open('w',encoding='utf-8') as h:
  h.write(f'started_at_utc={started}\ncommand={shlex.join(cmd)}\nstdin={stdin!r}\n{version}\n--- output ---\n');h.flush()
  p=subprocess.run(cmd,input=stdin or None,text=True,stdout=h,stderr=subprocess.STDOUT,check=False);elapsed=int(round(time.monotonic()-start));h.write(f'\nended_at_utc={utc()}\nexit_status={p.returncode}\nwall_seconds={elapsed}\n')
 if elapsed>7200:raise RuntimeError(f'single task exceeded 2h: {log}')
 if p.returncode:raise RuntimeError(f'command failed: {log}')
 return elapsed
def log_succeeded(log):
 return any('exit_status=0' in p.read_text(encoding='utf-8',errors='replace') for p in log.parent.glob(log.stem+'*'+log.suffix))
def archive_if_interrupted(paths,log):
 if any(Path(p).exists() for p in paths) and not log_succeeded(log):
  stamp=int(time.time())
  for p in map(Path,paths):
   if p.exists():os.replace(p,p.with_name(f'{p.name}.interrupted_partial.{stamp}'))
def xvg(path):
 rows=[]
 for line in Path(path).read_text().splitlines():
  if line.strip() and not line.lstrip().startswith(('#','@')):rows.append([float(x) for x in line.split()])
 return rows
def contact_pairs(path):
 groups=parse_ndx(path);names=[n for n in groups if n.startswith('contacts_')]
 if len(names)!=1:raise ValueError(f'expected one contact group in {path}')
 ids=groups[names[0]]
 if len(ids)%2:raise ValueError('odd contact atom list')
 return list(zip(ids[::2],ids[1::2]))
def xpm_rows(path):
 text=Path(path).read_text();m=re.search(r'"(\d+) (\d+)\s+2 1"',text)
 if not m:raise ValueError('XPM header missing')
 width,height=map(int,m.groups());rows=re.findall(r'^"([ o]+)"[,]?$',text,re.M)[-height:]
 if len(rows)!=height or any(len(x)!=width for x in rows):raise ValueError('XPM matrix dimensions mismatch')
 return rows
def collapse_contacts(num_path,pairs_path,xpm_path,atom_info,group1,group2,out_series,out_pairs):
 times=xvg(num_path)
 if not Path(pairs_path).exists() and not Path(xpm_path).exists():
  if any(abs(r[1])>1e-12 for r in times):raise ValueError('missing sparse files but atom-contact count is nonzero')
  with Path(out_series).open('w',newline='',encoding='utf-8') as h:w=csv.writer(h,lineterminator='\n');w.writerow(['time_ns','unique_residue_contacts']);w.writerows((r[0]/1000.0,0) for r in times)
  fields=['chain_1','resnr_1','resname_1','chain_2','resnr_2','resname_2','occupancy','frames_present','total_frames','atom_pairs_observed']
  with Path(out_pairs).open('w',newline='',encoding='utf-8') as h:csv.DictWriter(h,fieldnames=fields,lineterminator='\n').writeheader()
  return 0
 times=xvg(num_path);pairs=contact_pairs(pairs_path);matrix=xpm_rows(xpm_path)
 if len(pairs)!=len(matrix):raise ValueError('pair/XPM row mismatch')
 set1=set(group1);by_res={}
 for i,(a,b) in enumerate(pairs):
  if a not in set1:a,b=b,a
  if a not in set1:raise ValueError('contact pair does not cross requested groups')
  ra=(atom_info[a]['chain'],atom_info[a]['resnr'],atom_info[a]['resname']);rb=(atom_info[b]['chain'],atom_info[b]['resnr'],atom_info[b]['resname'])
  by_res.setdefault((ra,rb),[]).append(i)
 counts=[0]*len(times);occupancies=[]
 for (ra,rb),indices in by_res.items():
  present=[any(matrix[i][frame]=='o' for i in indices) for frame in range(len(times))]
  for frame,value in enumerate(present):counts[frame]+=int(value)
  occupancies.append({'chain_1':ra[0],'resnr_1':ra[1],'resname_1':ra[2],'chain_2':rb[0],'resnr_2':rb[1],'resname_2':rb[2],'occupancy':sum(present)/len(present),'frames_present':sum(present),'total_frames':len(present),'atom_pairs_observed':len(indices)})
 if [int(r[1]) for r in times] != [sum(row[f]=='o' for row in matrix) for f in range(len(times))]:raise ValueError('GROMACS atom-contact count does not match XPM')
 with Path(out_series).open('w',newline='',encoding='utf-8') as h:w=csv.writer(h,lineterminator='\n');w.writerow(['time_ns','unique_residue_contacts']);w.writerows((r[0]/1000.0,c) for r,c in zip(times,counts))
 fields=['chain_1','resnr_1','resname_1','chain_2','resnr_2','resname_2','occupancy','frames_present','total_frames','atom_pairs_observed']
 with Path(out_pairs).open('w',newline='',encoding='utf-8') as h:w=csv.DictWriter(h,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(sorted(occupancies,key=lambda r:r['occupancy'],reverse=True))
 return len(by_res)
def main():
 a=args();root=Path(__file__).resolve().parents[2];os.chdir(root);cfg=yaml.safe_load(Path(a.config).read_text());sp=cfg['spatial_organization']
 if not ({'phase3B_spatial_organization','phase3C_DPA_T01_targeted_recovery'} & set(cfg['authorization_scope'])) or not cfg['project']['formal_analysis_authorized']:raise RuntimeError('Phase 3B/recovery spatial analysis not authorized')
 cutoffs=[float(x) for x in sp['sensitivity_contact_cutoffs_nm']];traj={r['topology_id']:r for r in read(a.trajectory_manifest)};evidence={r['topology_id']:r for r in read(a.component_evidence)};component_rows=read(a.component_map);all_maps=read(a.simulation_map)
 eligible=[r for r in all_maps if traj[r['topology_id']]['status']=='pass' and evidence[r['topology_id']].get('phase3A_gate')=='pass']
 if len(eligible) not in (15,16):raise RuntimeError(f'expected 15 or 16 evidence-eligible topologies, got {len(eligible)}')
 if a.only:eligible=[r for r in eligible if r['topology_id'] in set(a.only)]
 manifest_path=Path(a.manifest);existing=read(manifest_path) if manifest_path.exists() else [];by_id={r['topology_id']:r for r in existing};pending=[r for r in eligible if by_id.get(r['topology_id'],{}).get('status')!='pass']
 if a.dry_run:print(json.dumps({'eligible':len(eligible),'pending':len(pending),'cutoffs':cutoffs}));return 0
 vr=subprocess.run([str(GMX),'--version'],capture_output=True,text=True,check=False);version=vr.stdout+vr.stderr
 if vr.returncode:raise RuntimeError('gmx version failed')
 builder=load_component_builder()
 for sim in pending:
  tid=sim['topology_id'];started=utc();tic=time.monotonic();out=Path('results')/sim['canonical_system_id']/tid/'spatial_metrics';gmxout=out/'gromacs';out.mkdir(parents=True,exist_ok=True);gmxout.mkdir(exist_ok=True);logs=Path('logs/phase3B')/tid
  result={'topology_id':tid,'system_id':sim['system_id'],'canonical_system_id':sim['canonical_system_id'],'trajectory':traj[tid]['analysis_xtc'],'status':'failed','started_utc':started,'gromacs_version':'2025.2','script_version':VERSION}
  try:
   atom_info=atom_metadata(sim,builder);base_groups=parse_ndx(evidence[tid]['analysis_index']);groups=dict(base_groups)
   for name in ('enzyme_1','enzyme_2','RIAD','RIDD'):
    groups[name+'_heavy']=[i for i in base_groups[name] if atom_info[i]['mass']>2.0]
   for n in (1,2):
    ca=[i for i in base_groups[f'linker_{n}'] if atom_info[i]['atomname']=='CA']
    if len(ca)!=17:raise ValueError(f'linker_{n} does not contain 17 C-alpha atoms')
    groups[f'linker_{n}_ca']=ca
   idx=out/'analysis_index_phase3B.ndx'
   if not idx.exists():write_ndx(idx,groups)
   xtc=traj[tid]['analysis_xtc'];tpr=sim['primary_tpr']
   com_xvg=gmxout/'com_distances.xvg'
   archive_if_interrupted([com_xvg],logs/'com_distances.log')
   if not com_xvg.exists():run_logged(logs/'com_distances.log',[str(GMX),'distance','-f',xtc,'-s',tpr,'-n',str(idx),'-oall',str(com_xvg),'-tu','ns','-select','com of group "enzyme_1" plus com of group "enzyme_2" plus com of group "RIAD" plus com of group "RIDD"'],version)
   mindist_xvg=gmxout/'enzyme_min_distance.xvg';mindist_contacts=gmxout/'enzyme_atom_contacts_0p45.xvg'
   archive_if_interrupted([mindist_xvg,mindist_contacts],logs/'enzyme_min_distance.log')
   if not mindist_xvg.exists():run_logged(logs/'enzyme_min_distance.log',[str(GMX),'mindist','-f',xtc,'-s',tpr,'-n',str(idx),'-od',str(mindist_xvg),'-on',str(mindist_contacts),'-d','0.45','-tu','ns'],version,'enzyme_1_heavy\nenzyme_2_heavy\n')
   linker_xvg=gmxout/'linker_geometry.xvg'
   archive_if_interrupted([linker_xvg],logs/'linker_geometry.log')
   if not linker_xvg.exists():
    atoms=[]
    for n in (1,2):
     ca=groups[f'linker_{n}_ca'];atoms.extend([(ca[0],ca[-1])]+list(zip(ca[:-1],ca[1:])))
    selection=' plus '.join(f'atomnr {x}' for pair in atoms for x in pair)
    run_logged(logs/'linker_geometry.log',[str(GMX),'distance','-f',xtc,'-s',tpr,'-n',str(idx),'-oall',str(linker_xvg),'-tu','ns','-select',selection],version)
   for interface,g1,g2 in (('enzyme',groups['enzyme_1_heavy'],groups['enzyme_2_heavy']),('riad_ridd',groups['RIAD_heavy'],groups['RIDD_heavy'])):
    name1='enzyme_1_heavy' if interface=='enzyme' else 'RIAD_heavy';name2='enzyme_2_heavy' if interface=='enzyme' else 'RIDD_heavy'
    for cutoff in cutoffs:
     tag=f'{cutoff:.2f}'.replace('.','p');series=out/f'{interface}_contacts_{tag}_timeseries.csv';pairs_csv=out/f'{interface}_contacts_{tag}_pair_occupancy.csv'
     if series.exists() and pairs_csv.exists():continue
     with tempfile.TemporaryDirectory(prefix=f'phase3B_{tid}_{interface}_{tag}_',dir='/tmp') as tmp:
      num=Path(tmp)/'num.xvg';pairs=Path(tmp)/'pairs.ndx';matrix=Path(tmp)/'matrix.xpm'
      cmd=[str(GMX),'hbond-legacy','-contact','-f',xtc,'-s',tpr,'-n',str(idx),'-num',str(num),'-hbn',str(pairs),'-hbm',str(matrix),'-r',f'{cutoff:.2f}']
      run_logged(logs/f'{interface}_contacts_{tag}.log',cmd,version,f'{name1}\n{name2}\n')
      collapse_contacts(num,pairs,matrix,atom_info,g1,g2,series,pairs_csv)
   analyze=[sys.executable,'scripts/spatial/analyze_spatial_metrics.py','--topology-id',tid,'--system-id',sim['system_id'],'--canonical-system-id',sim['canonical_system_id'],'--input-dir',str(out),'--component-map',a.component_map,'--config',a.config]
   run_logged(logs/'analyze.log',analyze,version)
   qc=json.loads((out/'spatial_qc.json').read_text());result.update({'analysis_index':str(idx),'status':'pass' if qc['qc_pass'] else 'failed','classification':'validated_spatial_metrics' if qc['qc_pass'] else 'technical_failure','notes':'phase3B_complete' if qc['qc_pass'] else 'spatial_qc_failed'})
  except Exception as exc:result['notes']=f'{type(exc).__name__}: {exc}'
  result['ended_utc']=utc();result['wall_seconds']=int(round(time.monotonic()-tic));result['output_bytes']=sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
  by_id[tid]=result;ordered=[by_id[x['topology_id']] for x in all_maps if x['topology_id'] in by_id];atomic_csv(manifest_path,ordered,FIELDS);print(json.dumps({'topology_id':tid,'status':result['status'],'wall_seconds':result['wall_seconds'],'output_bytes':result['output_bytes']},sort_keys=True),flush=True)
  total=sum(int(r.get('output_bytes') or 0) for r in by_id.values())
  if total>20_000_000_000:raise RuntimeError('Phase 3B output exceeds 20 GB gate')
 return 0 if all(by_id.get(r['topology_id'],{}).get('status')=='pass' for r in eligible) else 2
if __name__=='__main__':
 try:raise SystemExit(main())
 except Exception as e:print(f'run_phase3b_spatial: {type(e).__name__}: {e}',file=sys.stderr);raise SystemExit(2)
