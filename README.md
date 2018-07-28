# Pararead: parallel processing of sequencing reads

Pararead is a python package that simplifies parallel processing of DNA sequencing reads (BAM or SAM files), 
by parallelizing across chromosomes. Pararead is built for developers of python scripts that process data read-by-read. It enables you to quickly and easily parallelize your script.

## Install

`pararead` is hosted on pypi. Install with:

```
pip install --user pararead
```

Or, within active environment:
```
pip install --upgrade pararead
```

## Minimum working example

In the `examples` folder you can find [examples/count_reads.py](examples/count_reads.py), which will count the number of reads in a sam/bam file in parallel.

Run this on your bam file like this:

```
python count_reads.py file.bam -O output.txt --cores 2
```

Look at the code to see how this is implemented.

## Developing tools that use pararead

The main model provided is an abstract class called`ParaReadProcessor`, for which concrete children are created by implementing a `__call__` method. This creates a callable instance that is then mapped over chromosomes.

The concept is generally described in this early [blog post](http://databio.org/posts/tabix_files.html), which initiated the project that eventually became `pararead`. More details will be forthcoming.
