# DISCO

Decomposition Into Single-COpy gene trees ([DISCO](https://doi.org/10.1093/sysbio/syab070)) is a method for decomposing multi-copy gene-family trees while attempting to preserve orthologs and discard paralogs. These single-copy gene trees can be subsequently used by methods that can estimate species trees from single-copy gene trees such as [ASTRAL](https://github.com/smirarab/ASTRAL) or [ASTRID](https://github.com/pranjalv123/ASTRID) in order to obtain an accurate estimation of the species tree. Additionally, DISCO can be paired with concatenation analysis using the script `ca_disco.py`. 

**NOTE:** For species tree estimation default settings are recommended; however, for orthology detection using `-m 2` is recommended so small groups are retrieved. It is also highly recommended that you use the most recent version of DISCO, as it deals with some limitations of TreeSwift.

## Versions

- v1.0 (May 21st 2021): Initial Release
- v1.1 (July 19th, 2021): Added `ca_disco.py` and adjustable filter (instead of just filtering all trees with less than 4 taxa)
- v1.2 (June 8th, 2022): Tree rerooting fix (addresses issue with treeswift rooting with no edge lengths)
- v1.3 (July 8th, 2022): Allows for more flexibility with leaf labels
- v1.3.1 (October 31st, 2022): Fixes issue where extra newline character at the end of the file causes a crash 
- v1.4.0 (April 21st, 2024): Added ability to make a partition file to CA-DISCO
- v1.4.1 (March 30th, 2025): Suppressed unecessary warning added by treeswift (explained [here](https://github.com/jsdoublel/DISCO/issues/6)) that spammed output.

**Note:** At present, it is recommended to use the latest version if you want to run the DISCO algorithm as described in Willson et al. 2021, as using earlier versions can causes results inconsistent with the DISCO algorithm in rare cases. 

## Citation

If you use DISCO, please cite:
```
@article{willson2022disco,
  title={DISCO: Species tree inference using multicopy gene family tree decomposition},
  author={Willson, James and Roddur, Mrinmoy Saha and Liu, Baqiao and Zaharias, Paul and Warnow, Tandy},
  journal={Systematic Biology},
  volume={71},
  number={3},
  pages={610--629},
  year={2022},
  publisher={Oxford University Press}
}
```

## Algorithm

Given a list of multi-copy gene trees, DISCO does the following for each tree:

1. Root the tree and tag each internal vertex as either a duplication event or a speciation event in such a way that minimizes the total number of duplications and losses. We do this with the ASTRAL-Pro rooting and tagging algorithm ([Zhang et. al. 2020](https://doi.org/10.1093/molbev/msaa139)).
2. Decompose gene tree by splitting off the smallest subtree under every vertex tagged as a duplication from the leaves to the root until all duplication events are resolved; it returns the set of single-copy trees produced.

## Dependencies

- Python 3
- [TreeSwift](https://github.com/niemasd/TreeSwift)
- [Biopython](https://biopython.org/) (for `ca_disco.py)

Treeswift can be installed with `pip install treeswift` and Biopython can be installed with `pip install biopython`.

## Usage

### disco.py

**Input**: File containing list of multi-copy trees in newick format

**Output**: File containing resulting list of single-copy trees after decomposition in newick format

```
python3 disco.py -i <input_file> -o <ouput_file> -d <delimiter>
```

#### Arguments

##### Required

```
-i, --input           Input newick tree file
```

##### Optional

```
-o, --output          Output newick tree file
-d, --delimiter       Delimiter separating species name from rest of leaf label
-s, --single-tree     Output only single large tree per gene tree 
-n, --nth-delimiter   Split on Nth delimiter (only works with -d)
-m, --minimum         Minimum number of taxa required for tree to be outputted
-v, --verbose         Enable verbose output
--keep-labels         Keep original leaf labels instead of using species name
--no-decomp           Outputs rooted trees without decomposition
--outgroups           Write outgroups (including ties) to txt file
--remove-in-paralogs  Remove in-paralogs before rooting/scoring
```

#### Example

```bash
python3 disco.py -i example/gtrees-mult.trees
```

### ca_disco.py

**Input**: File containing list of multi-copy trees in newick format and set of alignment files corresponding to the gene families.

**Output**: Concatenated alignment file

```
python3 ca_disco.py -i <input_trees> -a <aln_list> -o <output> -d <delimiter> -m <number> 
```

`disco.py` must be present in the same directory as `ca_disco.py` in order for it to run. Also, unlike `disco.py`, it is necessary for the input newick trees given to `ca_disco.py` to have unique leaf labels where the taxon name comes first and is separated from the rest of the name by some delimiter. 

The `-a` argument should be given the path to an "alignment list" file containing the path to each alignment file you want to concatinate separated by a new line. The order of the alignment files is important---for each tree in the input newick tree file, there should be a corresponding alignment on the same line in the respective alignment list file.  

#### Arguments

##### Required

```
-i, --input           Input newick tree file
-a, --alignment       Text file containing paths to alignment files
-f, --format          Format of alignment file (either "fasta" or "phylip")
-o, --output          Output concatenated alignment file
```

##### Optional
```
-m, --filter          Minimum number of taxa required sequence group to be included
-d, --delimiter       Delimiter separating species name from rest of leaf label
-p, --partition       Create partition file
```

#### Example

```bash 
python3 ca_disco.py -i example/g_100.trees -o example.phy -a example/seq_list.csv -f phylip
```

#### Partition File

`ca_disco.py` now has the option to create a partition file; this can be done with the optional `-p` argument. If this argument is specified, it is necessary to provide relavent information in the alignment list file. For example, see `example/seq_list.csv`. This would generate a partition file like so:
```
GTR+G, 0001=1-200
GTR+G, 0002=201-700
GTR+G, 0004=701-1300
...
```
**Note:** Gene 0003 is missing. This is not a mistake; this is because it does not have enough informative information to generate large enough DISCO subtrees and thus is not included in the alignment.

