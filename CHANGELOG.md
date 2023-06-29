## Unreleased
## 2.3.2 - 2023-05-30
- Assembly-specific code in olive moved inside the wdl, changes to output file name
## 2.3.1 - 2022-08-31
- Changed workflow name to `delly_matched_by_tumor_group` for consistency
## 2.3.0 - 2022-01-25
- Delly upgrade to 0.9.1
- [GP-2874](https://jira.oicr.on.ca/browse/GP-2874) - Making regression test script more robust
## 2.2.0  - 2021-06-01
- Migrating to vidarr
## 2.1   - 2021-01-22
- Added filtering of PASS variants, additional output
## 2.0.2 - 2020-10-14
- Added "set -euo pipefail" to the runDelly task
## 2.0.1 - 2020-04-23
- Adding timeout setting to the first step (marking duplicate reads with picard)
## 2.0   - 2019-08-22
- [GP-2058](https://jira.oicr.on.ca/browse/GP-2058) - Converting to Cromwell
## 1.2.1 - 2019-02-05
- [GP-1913](https://jira.oicr.on.ca/browse/GP-1913) - Update to pipedev 2.4.7 (SeqWare 2.0.2) to fix provisioning of empty files
- [GP-1916](https://jira.oicr.on.ca/browse/GP-1916) - Rename germline mode to unmatched mode
## 1.2 - 2016-01-26
- Delly upgrade to 0.7.2
## 1.1 - 2015-03-19
- Upgrade to SeqWare 1.1.0, common-utilities 1.6 and workflow-utilities 1.6.
