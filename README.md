## Overview
Copy Number Variation pipeline is being replaced by a set of stand-alone tool-specific workflows and deciders, so the elements of CNV pipeline
are described elsewhere.

Delly2 workflow produces a set of vcf files with different types of structural variant calls: Translocation, Deletion, Inversion and Duplications
It uses .bam files as input. The below graph describes the process:

![delly flowchart](workflow-structural-variation/docs/delly-wf.png)

## Preprocessing

The expected inputs for the DELLY tool are library-level BAMs with distinct insert size and median. In most cases, this means that the BAM files will not need to be merged prior to processing. However, the DELLY website recommends marking duplicate reads. We may also [optionally] realign around indels and perform base recalibration.

### Mark duplicates

Picard Tools MarkDuplicates is used to flag reads as PCR or optical duplicates.
```
 java -jar MarkDuplicates.jar
 INPUT=sample.bam
 OUTPUT=sample.dedup.bam
 METRICS_FILE=sample.metrics
```

### Detect deletions
```
delly
-t DEL
-x excludeList.tsv
-o sample.DEL.bam
-q 0
-g hn19.fa
sample.bam
```
 
### Detect tandem duplications
```
delly
-t DUP
-x excludeList.tsv
-o sample.DUP.bam
-q 0
-g hn19.fa
sample.bam
```
 
### Detect inversions
```
delly
-t INV
-x excludeList.tsv
-o sample.INV.bam
-q 0
-g hn19.fa
sample.bam
Detect translocations
```
### Detecting translocations
```
delly
-t TRA
-x excludeList.tsv
-o sample.TRA.bam
-q 0
-g hn19.fa
sample.bam
```
### Post-processing
Each DELLY tool produces several files, which will all need to be merged together after the chromosomes are finished processing. The output format is described on the DELLY webpage. The merging script may require a small parser to combine the output from multiple runs in together.
Merge DELLY results with vcftools.

## Workflow Options
This workflow runs Structural Variation discovery tool DELLY using its four modes, merging the results into a single file at the end.

Parameter|Value|Description
---|---|---
ref_fasta | hg19.fa | Reference fasta file - Genome sequence, default is hg19.fa
exclude_list | human.hg19h.excl.tsv | To save runtime it is advisable to exclude telomere and centromere regions. For human, DELLY ships with such an exclude list
sample_name | TestSample | Sample name - should be derived automatically by the Decider
mapping_quality | 30 | customizable
input_bams | | comma-separated list of input .bam files
output_dir | seqware-results | A standard SeqWare parameter specifying the sub-directory where the output files will be provisioned
output_prefix | ./ | A standard SeqWare parameter specifying the root directory where the output files will be provisioned
manual_output | false | Whether or not to use manual output. When false, a random integer will be inserted into the path of the file in order to ensure uniqueness. When true, the output files will be moved to the location of output_prefix/output_dir
mode | | Mode: somatic OR germline

## Decider Options
This decider prepares .ini files for Structural Variation workflow runs, it extracts .bam files with assumption that they were not sorted or indexed (the former is conditioned on a flag). It will also pick the latest file if duplicates are present

Parameter|Value|Description
---|---|---
--ref-fasta | hg19.fa | Reference fasta file - Genome sequence, default is hg19.fa
--exclude-list | human.hg19h.excl.tsv | To save run time it is advisable to exclude telomere and centromere regions. For human, DELLY ships with such an exclude list.
--mapping-quality | | Not set by default, but customizable
--picard-memory | 6000 | memory (in Megabytes) allocated to Picard
--delly-memory | 8000 | memory (in Megabytes) allocated to Delly
--queue | | cluster queue
--output-dir | seqware-results | A standard SeqWare parameter specifying the sub-directory where the output files will be provisioned
--output-prefix | ./ | A standard SeqWare parameter specifying the root directory where the output files will be provisioned
--mode | | Mode: somatic OR germline
