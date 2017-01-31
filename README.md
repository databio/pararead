ParaRead is a module for processing sequencing reads (bam or sam files) in
parallel. You use it by writing a child class that extends the included
ParaReadProcessor class. It should implement the __call__() function at a minimum,
defining what to do with each chromosome.