#!/usr/bin/env python3
"""Summarize Phase 3B distance, contact, occupancy, and modeled-linker metrics."""

from __future__ import annotations
import argparse,csv,json,math,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/matplotlib-protein-scaffold')
import matplotlib.pyplot as plt
import numpy as np
import yaml

def args():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--topology-id',required=True);p.add_argument('--system-id',required=True);p.add_argument('--canonical-system-id',required=True);p.add_argument('--input-dir',required=True);p.add_argument('--component-map',default='metadata/component_group_map.csv');p.add_argument('--config',default='config/analysis.yaml');return p.parse_args()
def xvg(path):
 rows=[]
 for line in Path(path).read_text().splitlines():
  if line.strip() and not line.lstrip().startswith(('#','@')):rows.append([float(x) for x in line.split()])
 x=np.asarray(rows,float)
 if x.ndim!=2 or len(x)!=10001:raise ValueError(f'unexpected XVG shape {path}: {x.shape}')
 return x
def table(path):
 with Path(path).open(newline='',encoding='utf-8') as h:return list(csv.DictReader(h))
def write(path,fields,rows):
 with Path(path).open('w',newline='',encoding='utf-8') as h:w=csv.DictWriter(h,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
def stats(y):
 return {'n':len(y),'mean':float(np.mean(y)),'sd':float(np.std(y,ddof=1)),'median':float(np.median(y)),'min':float(np.min(y)),'max':float(np.max(y)),'q05':float(np.quantile(y,.05)),'q25':float(np.quantile(y,.25)),'q75':float(np.quantile(y,.75)),'q95':float(np.quantile(y,.95))}
def plot_line(path,t,series,ylabel,title,late=True):
 fig,ax=plt.subplots(figsize=(7.2,4.2))
 for name,y in series.items():ax.plot(t,y,lw=.7,label=name)
 if late:ax.axvspan(80,100,color='0.85',alpha=.45,label='exploratory 80–100 ns')
 ax.set(xlabel='Time (ns)',ylabel=ylabel,title=title)
 if len(series)>1:ax.legend(frameon=False,fontsize=8)
 fig.tight_layout();fig.savefig(path.with_suffix('.png'),dpi=220);fig.savefig(path.with_suffix('.svg'));plt.close(fig)
def main():
 a=args();d=Path(a.input_dir);cfg=yaml.safe_load(Path(a.config).read_text());sp=cfg['spatial_organization'];cutoffs=[float(x) for x in sp['sensitivity_contact_cutoffs_nm']]
 rows=[r for r in table(a.component_map) if r['topology_id']==a.topology_id];cres={r['component']:int(r['residue_count']) for r in rows};norm=math.sqrt(cres['enzyme_1']*cres['enzyme_2'])
 com=xvg(d/'gromacs/com_distances.xvg');mind=xvg(d/'gromacs/enzyme_min_distance.xvg');link=xvg(d/'gromacs/linker_geometry.xvg')
 if com.shape[1]!=3 or mind.shape[1]<2 or link.shape[1]!=35:raise ValueError(f'unexpected metric columns: {com.shape} {mind.shape} {link.shape}')
 t=com[:,0]
 if not (np.allclose(t,mind[:,0]) and np.allclose(t,link[:,0]) and np.all(np.diff(t)>0) and abs(t[0])<1e-6 and abs(t[-1]-100)<1e-6):raise ValueError('time alignment/window QC failed')
 metrics={'enzyme_com_distance_nm':com[:,1],'riad_ridd_com_distance_nm':com[:,2],'enzyme_min_distance_nm':mind[:,1],
          'linker_1_end_to_end_nm':link[:,1],'linker_1_extension_ratio':link[:,1]/np.sum(link[:,2:18],axis=1),
          'linker_2_end_to_end_nm':link[:,18],'linker_2_extension_ratio':link[:,18]/np.sum(link[:,19:35],axis=1)}
 pair_occupancy_summary=[]
 for interface in ('enzyme','riad_ridd'):
  for cutoff in cutoffs:
   tag=f'{cutoff:.2f}'.replace('.','p');cr=table(d/f'{interface}_contacts_{tag}_timeseries.csv');ct=np.array([float(r['time_ns']) for r in cr]);y=np.array([float(r['unique_residue_contacts']) for r in cr])
   if not np.allclose(t,ct):raise ValueError(f'{interface} {cutoff}: contact time mismatch')
   key=f'{interface}_contacts_{tag}';metrics[key]=y
   if interface=='enzyme':metrics[f'enzyme_contact_density_{tag}']=y/norm
   occ=table(d/f'{interface}_contacts_{tag}_pair_occupancy.csv');ov=np.array([float(r['occupancy']) for r in occ]) if occ else np.array([])
   pair_occupancy_summary.append({'topology_id':a.topology_id,'interface':interface,'cutoff_nm':cutoff,'observed_residue_pairs':len(ov),'mean_pair_occupancy':float(ov.mean()) if len(ov) else 0.0,'median_pair_occupancy':float(np.median(ov)) if len(ov) else 0.0,'max_pair_occupancy':float(ov.max()) if len(ov) else 0.0,'interface_occupancy':float(np.mean(y>0))})
 write(d/'contact_occupancy_summary.csv',list(pair_occupancy_summary[0]),pair_occupancy_summary)
 series_rows=[]
 for i,time_ns in enumerate(t):
  row={'time_ns':f'{time_ns:.6f}'};row.update({k:f'{v[i]:.9g}' for k,v in metrics.items()});series_rows.append(row)
 write(d/'spatial_metrics_timeseries.csv',list(series_rows[0]),series_rows)
 summary=[];blocks=[];hist=[];stability=[]
 for name,y in metrics.items():
  for window,start,end,mask in [('full_0_100',0,100,np.ones(len(t),bool)),('exploratory_80_100',80,100,t>=80)]:
   s=stats(y[mask]);s.update({'topology_id':a.topology_id,'system_id':a.system_id,'canonical_system_id':a.canonical_system_id,'metric':name,'window':window,'start_ns':start,'end_ns':end,'future_confirmed_window':''})
   if '_contacts_' in name:
    mean=s['mean'];defined=mean>float(sp['near_zero_contact_mean']);cv=s['sd']/mean if defined else None;s['contact_cv']=cv if defined else '';s['contact_stability']=1/(1+cv) if defined else '';s['contact_stability_status']='defined' if defined else 'undefined_near_zero_mean'
   else:s['contact_cv']='';s['contact_stability']='';s['contact_stability_status']='not_applicable'
   summary.append(s)
  for start in np.arange(0,100,10):
   mask=(t>=start)&((t<start+10) if start<90 else (t<=100));s=stats(y[mask]);s.update({'topology_id':a.topology_id,'metric':name,'start_ns':start,'end_ns':start+10})
   if '_contacts_' in name:
    defined=s['mean']>float(sp['near_zero_contact_mean']);cv=s['sd']/s['mean'] if defined else None;s['contact_cv']=cv if defined else '';s['contact_stability']=1/(1+cv) if defined else '';s['contact_stability_status']='defined' if defined else 'undefined_near_zero_mean'
    stability.append({'time_ns':start+5,'metric':name,'contact_stability':s['contact_stability']})
   else:s['contact_cv']='';s['contact_stability']='';s['contact_stability_status']='not_applicable'
   blocks.append(s)
  counts,edges=np.histogram(y,bins=30)
  for left,right,count in zip(edges[:-1],edges[1:],counts):hist.append({'topology_id':a.topology_id,'metric':name,'bin_left':left,'bin_right':right,'count':int(count),'fraction':count/len(y)})
 write(d/'spatial_metrics_summary.csv',list(summary[0]),summary);write(d/'spatial_time_block_summary.csv',list(blocks[0]),blocks);write(d/'spatial_metrics_distribution.csv',list(hist[0]),hist);write(d/'contact_stability_block_timeseries.csv',list(stability[0]),stability)
 for name in ('enzyme_com_distance_nm','riad_ridd_com_distance_nm','enzyme_min_distance_nm','linker_1_end_to_end_nm','linker_1_extension_ratio','linker_2_end_to_end_nm','linker_2_extension_ratio'):
  plot_line(d/name,t,{name:metrics[name]},name,a.topology_id+' '+name)
 for interface in ('enzyme','riad_ridd'):
  s={f'{x:.2f} nm':metrics[f'{interface}_contacts_{str(f"{x:.2f}").replace(".","p")}'] for x in cutoffs};plot_line(d/f'{interface}_contacts',t,s,'Unique residue pairs',a.topology_id+' '+interface+' contacts')
  stab={}
  for x in cutoffs:
   key=f'{interface}_contacts_{str(f"{x:.2f}").replace(".","p")}'
   stab[f'{x:.2f} nm']=np.array([float(r['contact_stability']) if r['contact_stability']!='' else np.nan for r in stability if r['metric']==key],float)
  plot_line(d/f'{interface}_contact_stability',np.arange(5,100,10),stab,'Contact stability',a.topology_id+' '+interface+' contact stability',late=False)
  occ_rows=[r for r in pair_occupancy_summary if r['interface']==interface]
  fig,ax=plt.subplots(figsize=(6.2,4.0));xx=np.arange(len(occ_rows));ax.bar(xx,[r['interface_occupancy'] for r in occ_rows],color='0.35')
  ax.set_xticks(xx,[f"{r['cutoff_nm']:.2f}" for r in occ_rows]);ax.set(xlabel='Contact cutoff (nm)',ylabel='Interface occupancy',title=a.topology_id+' '+interface+' interface occupancy');ax.set_ylim(0,1.05)
  fig.tight_layout();fig.savefig(d/f'{interface}_interface_occupancy.png',dpi=220);fig.savefig(d/f'{interface}_interface_occupancy.svg');plt.close(fig)
 if any(k.startswith('enzyme_contact_density_') for k in metrics):
  s={f'{x:.2f} nm':metrics[f'enzyme_contact_density_{str(f"{x:.2f}").replace(".","p")}'] for x in cutoffs};plot_line(d/'enzyme_contact_density',t,s,'Normalized contact density',a.topology_id+' contact density')
 qc={'topology_id':a.topology_id,'frames':len(t),'start_ns':float(t[0]),'end_ns':float(t[-1]),'strictly_monotonic_time':bool(np.all(np.diff(t)>0)),'finite':bool(all(np.isfinite(y).all() for y in metrics.values())),'nonnegative':bool(all((y>=0).all() for y in metrics.values())),'linker_extension_ratio_within_geometry':bool((metrics['linker_1_extension_ratio']<=1.0001).all() and (metrics['linker_2_extension_ratio']<=1.0001).all()),'cutoffs_nm':cutoffs,'primary_cutoff_nm':float(sp['primary_contact_cutoff_nm']),'equilibrium_confirmed':False,'exploratory_window':'80-100 ns','qc_pass':True}
 qc['qc_pass']=all(qc[x] for x in ('strictly_monotonic_time','finite','nonnegative','linker_extension_ratio_within_geometry'))
 (d/'spatial_qc.json').write_text(json.dumps(qc,indent=2,sort_keys=True)+'\n');print(json.dumps(qc,sort_keys=True));return 0 if qc['qc_pass'] else 2
if __name__=='__main__':raise SystemExit(main())
