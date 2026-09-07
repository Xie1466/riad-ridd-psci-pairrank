# Public release security/path audit

Status: **PASS after v1.0.1 patch**

The complete release tree was recursively scanned after the public-release patch for machine-specific absolute server paths, user-home paths, private-key headers, password/token/API-key markers, unintended email addresses, caches and excluded large binary classes. Historical absolute server paths found in one stored reference PSCI output were replaced with logical `repos/truefit/...` provenance labels. No model value or scientific field was changed.

Invalid DOI/GitHub placeholder identifiers were removed for public staging. Runtime-local temporary directories used by plotting libraries are not server dependencies. No credential, private key, server IP, unintended email address, third-party article/SI file, or proprietary figure was found.
