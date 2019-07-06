# blast
simple plotter able to handle big chromosomes,  
works with python3,  
use option -h or --help to get further instructions  
  
## MANUAL  

    blast.py					Version 1.2	23/04/2019
        --f1 <fasta/flat>
        --f2 <fasta/flat> *opt*
        --first <number> *opt, only take first n sequences*
        --contrast *opt,reduces markersize*
        --equal *opt,axis ratios equal*
        -p image_name *opt,default=noimage*
        -w width *opt,sliding window size, default=20*
        -l *opt,create txt with axis legend*


## EXAMPLEs

Example blast of 10MB reference chromosome against 10MB assembly contig
![10_MB_chr plot](/example/example_10MB_chr.png)
Example blast of two assemblies of the same Bacteria species, using first 20 sequences 
![bacteria flatfile plot](/example/example_bacteria.png)
