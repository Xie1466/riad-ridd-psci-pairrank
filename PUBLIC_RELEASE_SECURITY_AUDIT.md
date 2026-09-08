# Public release security/path audit

Status: **PASS for the final v1.0.2 public package**

The complete release tree was recursively scanned for machine-specific absolute server paths, user-home paths, private-key headers, credential markers, unintended email addresses, caches and excluded large binary classes. All 74 frozen scientific/configuration/code/test/reference files are byte-identical to v1.0.1.

All obsolete public-deposit wording was removed. Runtime-local temporary directories used by plotting libraries are not server dependencies. No credential, private key, server IP, unintended email address, third-party article/SI file, or proprietary figure was found.
