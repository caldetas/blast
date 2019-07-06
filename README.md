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


Example of legend.txt from the Bacteria plot
        contigs x-axis: bacterium01.fa
        0 >NODE_118_length_352008_cov_29.228
        1 >NODE_238_length_228980_cov_29.6623
        2 >NODE_256_length_215098_cov_31.518
        3 >NODE_292_length_195792_cov_29.5205
        4 >NODE_325_length_178895_cov_29.2128
        5 >NODE_344_length_171183_cov_30.6578
        6 >NODE_357_length_165959_cov_29.4406
        7 >NODE_364_length_163596_cov_28.86
        8 >NODE_414_length_147097_cov_29.793
        9 >NODE_453_length_135363_cov_29.3416
        10 >NODE_475_length_129330_cov_29.8447
        11 >NODE_492_length_125801_cov_28.2851
        12 >NODE_514_length_120819_cov_30.6622
        13 >NODE_640_length_95778_cov_28.4116
        14 >NODE_769_length_82167_cov_29.0989
        15 >NODE_773_length_82016_cov_29.5208
        16 >NODE_883_length_73596_cov_31.1848
        17 >NODE_918_length_71453_cov_30.2534
        18 >NODE_949_length_69467_cov_34.0411
        19 >NODE_966_length_68402_cov_28.6601


        contigs y-axis: bacterium02.fa
        0 >NODE_48_length_304421_cov_25.982110
        1 >NODE_62_length_275972_cov_26.232827
        2 >NODE_109_length_216134_cov_25.791569
        3 >NODE_133_length_197743_cov_25.800871
        4 >NODE_135_length_196203_cov_28.442474
        5 >NODE_213_length_148308_cov_25.846244
        6 >NODE_232_length_139918_cov_27.682598
        7 >NODE_250_length_136235_cov_26.766950
        8 >NODE_264_length_133629_cov_27.257244
        9 >NODE_343_length_115119_cov_26.436307
        10 >NODE_346_length_114847_cov_26.356133
        11 >NODE_357_length_113688_cov_27.432666
        12 >NODE_460_length_95666_cov_26.299177
        13 >NODE_510_length_89767_cov_26.162297
        14 >NODE_538_length_87057_cov_26.139913
        15 >NODE_580_length_82780_cov_26.816486
        16 >NODE_581_length_82750_cov_27.807276
        17 >NODE_599_length_80902_cov_27.474114
        18 >NODE_702_length_72549_cov_26.684391
        19 >NODE_712_length_72043_cov_27.677281
