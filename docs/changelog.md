# Change Log

This document follows guidelines from [Keep a Changelog](http://keepachangelog.com/en/0.3.0/) and  adheres to [semantic versioning](http://semver.org/).

## [0.8.0] - 2023-10-03
### Fixed
- Updated Python version
- Cleand and lint the code

## [0.7.0] - 2019-08-15
### Changed
- Use command callability checker from `ubiquerg`.
- Moved to logmuse, removing the logging functions from pararead

### Fixed
- Fixed bug with example that prevented it from working on multiple cores

## [0.6.0] - 2019-03-25
- Made compatible with python 3

## [0.5.0] - 2018-09-20
### Added
- Added `retain_temp` argument that is useful for debugging.
- Added modified version of read counting example to demonstrate logging.

### Changed
- Updated the CLI for logging options

## [0.4.0] - 2018-02-27
### Changed
- Updated packaging and documentation for release on pypi.


## [0.3.0] - 2017-06-09
### Changed
- A `ParaReadProcessor`'s `__call__` function is no longer applied to any 
chromosome that's empty (i.e., to which no reads mapped). Instead, an action 
(`empty_action`, which may be overridden) is taken for such a chromosome. 
This frees an implementor of `ParaReadProcessor` from the need to worry about 
empty chromosomes while allowing the flexibility for custom action in that case.

## [0.2.0] - 2017-06-09
### Added
- The read processor's constructor now establishes a root logger if that's 
not already been done by a client application.
- API functions to the `logs` module that assist clients in providing and 
parsing logging options, and with integrating and using a logging system
in their applications.
### Changed
- Previously, an attempt was made to optimize distribution of read chunks 
(e.g., by chromosome) by interleaving the read chunk keys (e.g., chromosome 
names) according to the genomic size (in base pairs) of the chromosome. This 
behavior is now turned off by default, as doing so seems to accelerate real 
runtime.
- The name for the function with which to fetch a file that's been registered 
with `pararead` changes from `fetch` to the clearer `fetch_file`, to 
distinguish this from a function that, say, fetches a block of reads or 
genomic positions/regions from a file.
- Simplified use of more information-rich development mode in logging system.
- The name for the function to create and configure a logger changes from 
`setup_pararead_logger` to the more general `setup_logger` to reflect the 
intention for this to serve as a common infrastructure for other repositories 
and projects.
