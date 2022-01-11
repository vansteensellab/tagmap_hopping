# tagmap_hopping
pipeline for tagmentation on local transposon "hopping"

This version was used to create the data for Zuin *et al.*:

Before aligning paired-end sequencing reads, reads were filtered using an adaptation of cutadapt, processing each read pair in multiple steps. Sequence patterns originating from Tn5 and each ITR were removed. The paired end reads coming from both ITRs were treated the same. First the presence of the unique part of the 5′ITR and 3′ITR sequence were detected at the start of the second read of the pair, and if present this sequence was trimmed. Next, sequence up to and including the TTAA site that were found on both the 5′ITR and 3′ITR was trimmed off. This sequence only partly contained the respective primers used for each ITR, and was used to filter reads that contained the sequence expected for a correct PCR product starting at the transposon. The sequence up to, but not including the TTAA was removed. Next all other sequence patterns coming from either Tn5 or the ITR were removed from the 5′ of the first read in the pair and the 3′ of both reads.

After filtering and trimming of the reads, reads were aligned to a reference genome with an *in silico* insertion of the split-GFP construct, but with a single TTAA motif instead of the PiggyBac transposon. This was done by aligning the homology arms found in the plasmid against mm10 reference genome. The complete sequence on the reference matching both arms was replaced by the plasmid sequence inserted.

Alignment was done using Bowtie2 with fragment length set to a minimum of 0 and max of 2000bp and the very-sensitive option was used. After reads were aligned to the genome, sambamba was used to remove duplicates and samtools was used to filter out read pairs that were not properly paired. We then designated, for each read-pair, the position of the first 4 nucleotides of the second read as a putative insertion site. To calculate fraction of reads originating from the non-mobilised position, the number of read pairs that overlapped the non-mobilised position (the TTAA replacing the PiggyBac of the in silico insert) was divided over the total number of reads originating from putative insertion sites supported by at least one read pair with a mapping quality higher than 2. 



The pipeline was run in a semi-automated fashion:
A manual edit was made to the chain file created in the pipeline.
This chain file can be found in the config folder together with the configuration files used.

The pipeline considers the sequence including the homology arms as a replacement.
The result is that the location of the piggybac insertions on the arms is lost during liftOver (insertions are still in the list, but without location on mm10).

A quick fix was used for this:
A new chain file was made with the right coordinates to replace the automaticaly generated chain file.
