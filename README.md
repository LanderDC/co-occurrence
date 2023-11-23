# Co-occurrence analysis
Scripts to correlate the presence of contigs in samples. Useful to find all segments of segmented viruses.

The script takes an abundance table (contigs should be rows, samples columns) and by default filters out all contigs that are present in less than 10% of the samples. The output is a correlation matrix,  which is then filtered to the preferred correlation coefficient (by default 0.3) and put into a pairwise dataframe for all contigs above the correlation threshold.

The abundance table is either transformed to a presence/absence table or, if you provide a file with the contig lenths (`-l/--lengths`), the read count is divided by the contig length. This would take also the abundance of ecah contig in account without bias towards the length of the contig (larger contigs have more reads, although they might not be that abundant).

The script also allows you to correlate specific contigs of interest with all other contigs in your study. It needs the names of your contigs of interest (`-s/--segments`). The output is a correlation matrix of all contigs above the correlation threshold with your contigs of interest.

`co-occurrence.py` and `co-occurrence.R` are basically the same (apart from the segment specific and contig length corrected options), the python script only runs much faster (eg. for ~9000 contigs 3m in python vs 1h3m for the R script with 42 threads). So use the R script only if you have no other option.

**Caveats:**
1. You need a decent amount of samples with your virus a) present *and* b) absent, otherwise there will be a) no correlations or b) too many correlations with non-related contigs (eg. host sequences).
2. The correlations should always be checked, the output is not perfect and is more a tool to help you recover unknown segments. You could for example check for approximately equal coverage of your related contigs, check the open reading frames of all your contigs, etc.

**Examples of input:**
1. Abundance table (`-i/--input`)
```
         Sample1  Sample2  Sample3  Sample4
Contig1  1839     0        868      0
Contig2  0        729      0        0
Contig3  1303     0        69       0
Contig4  0        0        0        90
```

2. Segments file (`-s/--segments`)
```
Contig1 
Contig3
etc.
```

3. Lengths file (`-l/--lengths`)
```
Contig1 7493
Contig2 2923
Contig3 3092
Contig4 1490
```

**Dependencies:**
1. `pandas`
2. `numpy`
