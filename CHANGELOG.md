# Changelog

## [v1.5.0](https://github.com/jsdoublel/DISCO/releases/tag/v1.5.0) - 2026-11-02

DISCO is now an installable Python package (`pip install disco-gdl`).

- `disco.py` and `ca_disco.py` are replaced by a single `disco` command with `decomp` and `concat` subcommands
- Flags renamed for consistency: `--single_tree` is now `-s`/`--single-tree` and `--remove_in_paralogs` is now `--remove-in-paralogs`
- Fixed `concat` rerooting at the best-scoring node instead of on its incident edge, which occasionally produced a non-optimal tagging/decomposition (and hence different alignments than `decomp` would imply). It also no longer spams TreeSwift rerooting warnings
- Fixed default `decomp` output/outgroup file names for input paths without an extension or starting with `./`
- Added a test suite

## [v1.4.1](https://github.com/jsdoublel/DISCO/releases/tag/v1.4.1) - 2025-03-30 

Suppressed unecessary warning added by treeswift (explained [here](https://github.com/jsdoublel/DISCO/issues/6)) that spammed output.
 
## [v1.4.0](https://github.com/jsdoublel/DISCO/releases/tag/v1.4) - 2024-04-21

Added ability to make a partition file to CA-DISCO

## [v1.3.1](https://github.com/jsdoublel/DISCO/releases/tag/v1.3.1) - 2022-08-31

Fixes issue where extra newline character at the end of the file causes a crash 

## [v1.3.0](https://github.com/jsdoublel/DISCO/releases/tag/v1.3) - 2022-07-08

Allows for more flexibility with leaf labels

## [v1.2.0](https://github.com/jsdoublel/DISCO/releases/tag/v1.2) - 2022-06-08

Tree rerooting fix (addresses issue with treeswift rooting with no edge lengths)

## [v1.1.0](https://github.com/jsdoublel/DISCO/releases/tag/v1.1) - 2021-07-19

Added `ca_disco.py` and adjustable filter (instead of just filtering all trees with less than 4 taxa)

## [v1.0.0](https://github.com/jsdoublel/DISCO/releases/tag/v1.0) - 2021-05-21

Initial Release

