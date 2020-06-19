# Comparison of different tools for Influenza genome sequencing data analysis

## Project purpose 
The project purpose is to compare different tools designed for Influenza virus whole-genome NGS data analysis (consensus sequence assembly + SNP calling).

## Project objectives  
1. Available tools analysis:
- searching for available tools
- local building
2. Processing of Illumina & Nanopore data with different tools:
- running tools with test data
- running tools with real data (in process)
3. Output comparison:
- output investigation
- coding script for analysis
- visualization & comparison  
4. Create recommendations for BWAcycle improvement (facultative, in process)

## Tools to compare

[1. BWAcycle – Smorodintsev Research Institute of Influenza](https://github.com/Molecular-virology-lab/bwacycle)  
Usage: run from directory with .fq.gz files listed in sample_list.txt:  
- For illumina: ```python3 path/to/map_all.py -list path/to/sample_list.txt```    
- For nanopore: ```python3 path/to/map_all.py -list path/to/sample_list.txt -data_type nanopore```

[2. IRMA – CDC Atlanta](https://wonder.cdc.gov/amd/flu/irma/)  
Usage:  
- For illumina: ```path/to/IRMA FLU path/to/R1.fastq.gz path/to/R2.fastq.gz <sample_lable>```    
- For nanopore: ```path/to/IRMA FLU-minion path/to/fastq.gz <sample_lable>```  

[3. INSaFLU – Instituto Nacional de Saude (INSA) Doutor Ricardo Jorge](https://github.com/INSaFLU/INSaFLU)  
Usage: online version.

[4. FluLINE – WHO CC Melbourne](https://github.com/UmaSangumathi/FluLINE)  
Excluded from comparison after exploratory analysis.  

[5. FluSeq – Singapore](https://github.com/hkailee/FluSeq)  
Excluded from comparison after exploratory analysis.  

## Script for comparison  
**analyzer.py** is a script that processes the BWAcycle, IRMA, and INSaFLU (this list may expand in the near future) outputs  and draws graphs comparing coverage and found SNPs for each of 8 genome segments. Now it supports only INSaFLU web-version output for Illumina reads (If the output for other technologies is the same, then it also supports them, but we didn't succeded in processing Nanopore reads with INSaFLU).

### Usage

Input: paths to BWAcycle, IRMA and INSaFLU sample output directory. The order is important.  

```python3 analyzer.py path/to/bwacycle/output/ path/to/IRMA/output/ path/to/INSaFLU/output/```  

Since this is the first very crude version, there are problems with file paths, so you need to run the analysis from the directory with the script. "/" at the end of the path is required.

Examples:
- for Illumina reads:  
```python3 analyzer.py bwacycle/Samples/36/ IRMA/36/ INSaFLU/```  
- for Nanopore reads:  
```python3 analyzer.py bwacycle/Samples/8_100_76/ IRMA/8_100_76_NP/```  

### Output description

Examples of output for 3 samples are avaliable in the analysis_results folder (36 - Illumina reads, other - Nanopore).  

Script draws graphs comparing coverage and found SNPs for each of 8 Influenza genome segments.  

#### Coverage comparison

Each graph reflects the coverage obtained for each segment of Influenza virus genome. BWAcycle and IRMA create statistics per position, INSaFLU shows only mean coverage for whole segment.  
![](https://github.com/MrDunn0/Influenza_2020/blob/master/analysis_results/36_analysis/36_images/coverage/36_HA.png)  

#### SNP comparison

Each graph displays the position of the variant and its meaning - the type of replacement determines its color.
![](https://github.com/MrDunn0/Influenza_2020/blob/master/analysis_results/36_analysis/36_images/variants/36_HA.png)  

The most important point in the analysis of whole-genome sequencing data is SNP calling. The new variants in the Influenza virus genome determine the characteristics of the circulating strains, and therefore the specifics of seasonal vaccines. The accuracy and reliability of the obtained SNPs are based on two components: 
- correct assembly and identification of the consensus sequence
- SNP calling actually  

Created plots allow us to identify gross errors in the assembly of the consensus sequence and / or its identification, without resorting to an eyeball analysis of dozens of output files. For example, in this way we found out IRMA cannot adequately analyze our test sample while BWAcycle didn’t found SNPs at all :( The most likely reason is the reference that was incorrectly chosen by the tool for alignment.

### Conclusions (*in process*)  
To answer which tools is the best it is not enough to look at a couple (of dozens) of graphs. We need to compare the data obtained in the analysis of a large number of “reference” samples, i.e. samples for which we know the strain and other metadata exactly. In addition to the data themselves, such analysis requires significant computing power.  
Since all human and computational resources of Smorodintsev Research Institute of Influenza are thrown into the study and fight against the COVID-19 epidemic, we were forced to suspend the analysis and draw intermediate conclusions.  
Based on the capabilities and limitations of the tools being compared, as well as the results of the analysis of the test sample, we came to the conclusion that the most suitable tools are **BWAcycle** and **IRMA**.   
*We hope to confirm this conclusion statistically soon, as well as develop a strategy for improving BWAcycle based on the difficulties and limitations discovered during the analysis.*  

### Dependencies

Python 3 libs:

- pandas 1.0.1

- matplotlib 3.1.3

- plotnine 0.6.0

- biopython 1.76

Out of Python:

MAFFT v7.450 (already in the repo)
