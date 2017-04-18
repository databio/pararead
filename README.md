# Pararead: parallel processing of sequencing reads

This project aims to simplify parallel processing of sequencing reads (BAM or SAM files), 
using the chromosome as the unit on which to group/partition. The main model provided is 
`ParaReadProcessor`, for which concrete children are created by implementing a `__call__`
method. This makes creates a callable instance that is then mapped over chromosomes.

## Install
```
pip install --user --update https://github.com/databio/pararead/zipball/master
```
