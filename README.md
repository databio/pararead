# Pararead: parallel processing of sequencing reads

Pararead is a python package that simplifies parallel processing of sequencing reads (BAM or SAM files), 
using the chromosome as the unit on which to group/partition. Pararead is build for developers of python scripts that process data read-by-read. It enables you to quickly and easily parallelize your script.

## Install
User:
```
pip install --user --upgrade https://github.com/databio/pararead/zipball/master
```
Within specific active environment:
```
pip install --upgrade https://github.com/databio/pararead/zipball/master
```

## Developing tools that use pararead

The main model provided is an abstract class called`ParaReadProcessor`, for which concrete children are created by implementing a `__call__`
method. This creates a callable instance that is then mapped over chromosomes.

The concept is generally described in this early [blog post](http://databio.org/posts/tabix_files.html), which initiated the project that eventually became `pararead`. More details will be forthcoming.
