# Comparison of different tools for Influenza genome sequencing data analysis

The project purpose is:
to compare different tools designed for Influenza virus whole-genome NGS data analysis (consensus sequence assembly + SNP calling).

Tools to compare:

1. BWAcycle – Research Institute of Influenza (https://github.com/Molecular-virology-lab/bwacycle)
2. IRMA – CDC Atlanta (https://wonder.cdc.gov/amd/flu/irma/)
3. InsaFLU – Instituto Nacional de Saude (INSA) Doutor Ricardo Jorge (https://github.com/INSaFLU/INSaFLU)
4. FluLINE – WHO CC Melbourne (https://github.com/UmaSangumathi/FluLINE)
5. FluSeq – Singapore (https://github.com/hkailee/FluSeq


analyzer.py is a script that processes the BWAcycle, IRMA, and INSaFLU (this list may expand in the near future) outputs  and draws graphs comparing coverage and found SNPs for each gene. Now it supports only INSaFLU web-version output for Illumina reads (If the output for other technologies is the same, then it also supports them, but we didn't succeded in processing nanopore reads with INSaFLU).

Input: paths to BWAcycle, IRMA and INSaFLU sample output directory. The order is important.

Examples of output for 3 samples are avaliable in the analysis_results folder (36 - Illumina reads, other - nanopore).

Usage example:

Illumina:

python3 analyzer.py bwacycle/Samples/36/ IRMA/36/ INSaFLU/

Nanopore:

python3 analyzer.py bwacycle/Samples/8_100_76/ IRMA/8_100_76_NP/

Since this is the first very crude version, there are problems with file paths, so you need to run the analysis from the directory with the script. "/" at the end of the path are required.

Dependencies:

Python 3 libs:

pandas 1.0.1
matplotlib 3.1.3
plotnine 0.6.0
biopython 1.76

Out of Python:

MAFFT v7.450 (already in the repo)



