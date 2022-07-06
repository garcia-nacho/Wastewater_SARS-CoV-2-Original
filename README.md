# Wastewater SARS-CoV-2 Surveillance Pipeline
[![Linux](https://svgshare.com/i/Zhy.svg)](https://svgshare.com/i/Zhy.svg)   [![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://https://docker.com/)

:warning: This repository is no longer maintained here. For updates visit [FHI's github](https://github.com/folkehelseinstituttet/Wastewater_SARS-CoV-2)

SARS-CoV-2 wastewater surveillance using long reads (Tested with Nanopore).

## Introduction
This pipeline identifies and quantifies SARS-CoV-2 lineages at the read-level.   
<img src="/Overview.png" width="700">     

## Installation   
<code>git clone https://github.com/garcia-nacho/Wastewater_SARS-CoV-2/ </code>  
<code> docker build -t wastewater Wastewater_SARS-CoV-2 </code>
   
## Run   
Basic run using default settings:   
<code>docker run -it --rm -v $(pwd):/Data wastewater </code>  
   
To filter reads by quality add the flag *-e qual*. E.g:   
<code>docker run -it --rm -e qual=10 -v $(pwd):/Data wastewater </code>   
(Use *-e qual=-1* to skip the filtering step)   
   
To change the default noise cut-off (0.15) use the flag *-e noise*. E.g:       
<code>docker run -it --rm -e noise=0.3 -v $(pwd):/Data wastewater </code>

To change the region to analyze (default 1250-2250 by default) use the flags *-e start* and *-e end*. E.g:    
<code>docker run -it --rm -e start=1000 -e end=2000 -v $(pwd):/Data wastewater </code>

The script must be run in a folder with the following structure:

<pre>
./ExpXX         
  |-Sample1     
      |-File_XXXXX_1.fastq.gz       
      |-File_XXXXX_2.fastq.gz
      |-File_XXXXX_3.fastq.gz
      |-...
  |-Sample2      
      |-File_XXXXX_1.fastq.gz       
      |-File_XXXXX_2.fastq.gz
      |-File_XXXXX_3.fastq.gz
      |-... 
  |-...   

</pre>

The filename of the *.fastq.gz* files are irrelevant and the samples are named using the folder that containes them as name    

Alternatively, you can use the prebuilt docker image stored at [dockerhub](https://hub.docker.com/repository/docker/garcianacho/wastewater)

<code>docker pull garcianacho/wastewate && docker run -it --rm -v $(pwd):/Data wastewater</code>
   
Note that older versions of docker might require the flag <code>--privileged</code> to run properly.       
## Output   
The pipeline generates four folders: analysis, bam, QC, sequences   
   
**analysis**    
Here you will find:   
-Sankey plots showing the most abundant combinations of mutations at nucleotide and amino acid levels   
-Barplots showing the relative abundance of the different combinations. Differente granularity levels are included (amino acid sequence level, nextclade clade level, pangolin level)   
-Excel files with the raw results    
-The nextclade output for the different variants identified   

**bam**   
Here you will find the *bam* files containing the reads aligned against the Spike gene of the *wuhan-hu-1* strain   
   
**QC**   
Here you will find:
-Plots showing the coverage of the samples over the spike gene   
-Noise at the different positions (Noise is defined here as ratio of bases not included in the consensus)   
-A noise.tsv file that contains the raw data regarding coverage and noise   
-consensus_qual.txt. Quality of the different bases called in the consensus   
   
**sequences**   
Here you find the fasta sequences of the consensus and variants found in each sample    
     
## Under the hood
The pipeline filters the reads according a quality cut-off using *[seqkit](https://bioinf.shenwei.me/seqkit/)*. Reads are mapped against the reference Spike using *[minimap2](https://github.com/lh3/minimap2)* and filtered and sorted using *[samtools](http://www.htslib.org/)*. The resulting *bam* file is indexed using *[samtools](http://www.htslib.org/)*. The positions showing a mix of bases are identified using *[noisefinder](https://github.com/garcia-nacho/NoisExtractor)* and the bases connected with their read-ids are retreived using *[bbasereader](https://github.com/garcia-nacho/bbasereader)*. All positions are merged using the read-id as pivot column and the different variants are indenfied and saved in a fasta file. The resulting fasta files are analyzed with *[Nextclade](https://clades.nextstrain.org/)* to find the mutations at the amino acid level, the clades and *Pangolin lineages*. The plots and analyses are generated in *[R](https://www.r-project.org/)* 
