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

[2. IRMA – CDC Atlanta](https://wonder.cdc.gov/amd/flu/irma/)

[3. INSaFLU – Instituto Nacional de Saude (INSA) Doutor Ricardo Jorge](https://github.com/INSaFLU/INSaFLU)

[4. FluLINE – WHO CC Melbourne](https://github.com/UmaSangumathi/FluLINE)

[5. FluSeq – Singapore](https://github.com/hkailee/FluSeq)

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

### Dependencies

Python 3 libs:

- pandas 1.0.1

- matplotlib 3.1.3

- plotnine 0.6.0

- biopython 1.76

Out of Python:

MAFFT v7.450 (already in the repo)



