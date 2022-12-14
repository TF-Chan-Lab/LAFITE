LAFITE
======

Low-abundance Aware Full-length Isoform clusTEr

Overview
--------
LAFITE is designated to identify high-consensus full-length isoforms from Nanopore Direct RNA-seq data. LAFITE combines multiple features from reference annotation and DRS reads (TSS, TES, splicing junction, and read polyadenylation event) and is more sensitive to Low-abundance transcripts.



Prerequisites
-------------
* [bedtools](https://github.com/arq5x/bedtools2)
* [Minimap2](https://github.com/lh3/minimap2)
* [nanopolish](https://github.com/jts/nanopolish)
* [samtools](http://www.htslib.org)
* Python 3.11

Installation
------------
To avoid potential conflicts, we recommend running LAFITE in a conda environment.
```
conda create -n LAFITE_env -c conda-forge python=3.11
conda activate LAFITE_env

# stable release
pip install LAFITE

# or the latest development version 
pip install git+https://github.com/TF-Chan-Lab/LAFITE
```

Usage
-----
1. **Run minimap2 and samtools to generate alignment file in bam format**
	```
	minimap2 -ax splice -u f -k 14 -G 500000 --secondary=no REFERENCE_FA FASTQ > ALIGNMENT_SAM
	samtools view -bS ALIGNMENT_SAM|samtools sort - > ALIGNMENT_BAM
	```
	LAFITE also supports other splicing-aware long read alignment tools.
2. **Run Nanopolish polya to generate read polyadenylation result (optional but recommend)**  
Current long-read sequencing technologies (Nanopore cDNA/DRS or PacBio Iso-Seq) are all designed to capture RNA molecules with poly(A) tail. However, RNA fragmentation and pore blocking may bring a considerable part of truncated reads which will interfere downstream analysis. Therefore, LAFITE utilizes the read polyadenylation status reported by Nanopolish to filter reads that have completed the sequencing process.  

   ```
   nanopolish index -d PATH_TO_FAST5 -s GUPPY_SEQUENCING_SUMMARY FASTQ
   nanopolish polya -t NUM_OF_THREADS -r FASTQ -b ALIGNMENT_BAM -g REFERENCE_FA > Nanopolish_PolyA_RES
   ```
	LAFITE also provides an alternative approach to estimate read polyadenylation status by scanning any poly(A) motifs that existed at the read 3'-end.  

1. **Run LAFITE**  
	```
	usage: lafite [-h] -b BAM [-B BEDTOOLS] -g GTF -f GENOME -o OUTPUT
              [-n MIN_COUNT_TSS_TES] [-i MIS_INTRON_LENGTH]
              [-c MIN_NOVEL_TRANS_COUNT] [-s MIN_SINGLE_EXON_COVERAGE]
              [-l MIN_SINGLE_EXON_LEN] [-L LABEL] [-p POLYA]
              [-m POLYA_MOTIF_FILE] [-r RELATIVE_ABUNDANCE_THRESHOLD]
              [-j SHORT_SJ_TAB] [-w SJ_CORRECTION_WINDOW] [--no_full_cleanup]
              [-t THREAD] [-T TSS_PEAK] [-d TSS_CUTOFF]

	Low-abundance Aware Full-length Isoform clusTEr

	optional arguments:
	  -h, --help            show this help message and exit
	  -b BAM                path to the alignment file in bam format
	  -B BEDTOOLS           path to the executable bedtools
	  -g GTF                path to the reference gene annotation in GTF format
	  -f GENOME             path to the reference genome fasta
	  -o OUTPUT             path to the output file
	  -n MIN_COUNT_TSS_TES  minimum number of reads supporting a alternative TSS or TES, default: 3
	  -i MIS_INTRON_LENGTH  length cutoff for correcting unexpected small intron, default: 150
	  -c MIN_NOVEL_TRANS_COUNT
	                        minimum occurrences required for a isoform from novel loci, default: 3
	  -s MIN_SINGLE_EXON_COVERAGE
	                        minimum read coverage required for a novel single-exon transcript, default: 4
	  -l MIN_SINGLE_EXON_LEN
	                        minimum length for single-exon transcript, default: 100
	  -L LABEL              name prefix for output transcripts, default: LAFT
	  -p POLYA              path to the file contains read Polyadenylation event
	  -m POLYA_MOTIF_FILE   path to the polya motif file
	  -r RELATIVE_ABUNDANCE_THRESHOLD
	                        minimum abundance of the predicted multi-exon transcripts as a fraction of the
							total transcript assembled at a given locus, default: 0.01
	  -j SHORT_SJ_TAB       path to the short read splice junction file
	  -w SJ_CORRECTION_WINDOW
	                        edit distance to reference splicing site for splicing correction, default: 40
	  --no_full_cleanup     keep all intermediate files
	  -t THREAD             number of the threads, default: 4
	  -T TSS_PEAK           path to the TSS peak file
	  -d TSS_CUTOFF         minimum TSS distance for a transcript to be considered as a novel transcript
	```
- LAFITE can run with the following arguments:
   ```
   lafite -b ALIGNMENT_BAM -g REFERENCE_GTF -f REFERENCE_FA -o OUTPUT_GTF -t NUM_OF_THREADS -p Nanopolish_PolyA_RES
   ```
- LAFITE can also run without the result from *nanoplish polya*. Then, a Poly(A) motif list must be provided for the corresponding species.  
   We have provided the Poly(A) motif list for human and mouse retrieved from [*Tian* et al.](https://academic.oup.com/nar/article/33/1/201/2401035) .
   
   ```
   lafite -b ALIGNMENT_BAM -g REFERENCE_GTF -f REFERENCE_FA -o OUTPUT_GTF -t NUM_OF_THREADS -m POLYA_MOTIFS_OF_SPECIES
   ```
- LAFITE accepts the TSS peaks from 5'-end CAGE data for identifying high-confidence TSSs. Users can prepare the TSS data in the following format where:
  - The first column is the chromosome name
  - The second column is the 0-based start position of the TSS peak
  - The third column is the 1-based end position of the TSS peak
  - The fourth column is the strand information  
- LAFITE also accepts the splicing junctions from Illumina short read RNA-seq data to proof the long reads. LAFITE supports the SJ.out.tab from STAR aligner. Users can also prepare the splicing junctions in the following format where:
  - The first column is the chromosome name
  - The second column is the 0-based start position of the splicing junction
  - The third column is the 1-based end position of the splicing junction
  - The fourth column is the strand information  

Development
-----------

LAFITE was developed following the [fastai/nbdev](https://github.com/fastai/nbdev) framework.
