# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.6.1] - 2025-04-03
### Added
- runtime attribute io_slots

## [2.6.0] - 2024-08-25
### Added
- exposed parameters which are helpful for getting through difficult data:
  * insertSizeCutoff
  * translocationQuality
  * minClip
  * minCliqueSize
  * minRefSeparation
  * maxReadSeparation
- additionalPrameters for even more fine-tuning

## [2.5.0] - 2024-06-25
### Added
[GRD-797](https://jira.oicr.on.ca/browse/GRD-797)] - add vidarr labels to outputs (changes to medata only)

## [2.4.2] - 2024-03-15
### Changed
- Workflow definitions for delly_normal_only were made after the last tag and thus are not available to shesmu. Retagging to be able to use delly_normal_only

## [2.4.1] - 2024-03-14
### Added
- Adding normal_only alias to use delly files downstream for dellyGermline

## [2.4.0] - 2023-05-30
### Changed
- Same code as 2.4.0, use major version to mark upgrades.

## [2.3.2] - 2023-05-30
### Changed
- Assembly-specific code in olive moved inside the wdl, changes to output file name

## [2.3.1] - 2022-08-31
### Changed
- Changed workflow name to `delly_matched_by_tumor_group` for consistency

## [2.3.0] - 2022-01-25
### Changed
- Delly upgrade to 0.9.1
- [GP-2874](https://jira.oicr.on.ca/browse/GP-2874)] - Making regression test script more robust

## [2.2.0] - 2021-06-01
### Changed
- Migrating to vidarr

## [2.1.0] - 2021-01-22
### Changed
- Added filtering of PASS variants, additional output

## [2.0.2] - 2020-10-14
### Added
- Added "set -euo pipefail" to the runDelly task

## [2.0.1] - 2020-04-23
### Added
- Adding timeout setting to the first step (marking duplicate reads with picard)

## [2.0.0] - 2019-08-22
### Changed
- [GP-2058](https://jira.oicr.on.ca/browse/GP-2058)] - Converting to Cromwell

## [1.2.1] - 2019-02-05
### Changed
- [GP-1913](https://jira.oicr.on.ca/browse/GP-1913)] - Update to pipedev 2.4.7 (SeqWare 2.0.2) to fix provisioning of empty files
- [GP-1916](https://jira.oicr.on.ca/browse/GP-1916)] - Rename germline mode to unmatched mode

## [1.2.0] - 2016-01-26
### Changed
- Delly upgrade to 0.7.2

## [1.1.0] - 2015-03-19
### Changed
- Upgrade to SeqWare 1.1.0, common-utilities 1.6 and workflow-utilities 1.6.
