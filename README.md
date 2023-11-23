# segment_co-occurrence
Scripts to correlate the presence of contigs in samples. Useful to find all segments of segmented viruses.

`co-occurrence.py` and `co-occurrence.R` are basically the same, the python script only runs much faster (eg. for ~9000 contigs 3m vs 1h3m for the R script with 42 threads). So use the R script only if you have no other option.
The script takes an abundance table (contigs should be rows, samples columns) and by default filters out all contigs that are present in less than 10% of the samples. The output is a correlation matrix,  which is then filtered to the preferred correlation coefficient (by default 0.3) and put into a pairwise dataframe for all contigs above the correlation threshold.

`specific_co-occurrence.py` allows you to correlate specific contigs of interest with all other contigs in your study. It needs an abundance table (contigs should be rows, samples columns) and the names of your contigs of interest (`-s/--segments`). The output is a correlation matrix of all contigs above the correlation threshold with your contigs of interest.

In both scripts, the abundance table is either transformed to a presence/absence table or, if you provide a file with the contig lenths (`-l/--lengths`), the read count is divided by the contig length. This would take also the abundance of ecah contig in account without bias towards the length of the contig (larger contigs have more reads, although they might not be that abundant).

