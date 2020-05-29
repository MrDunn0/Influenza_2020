# This is a script for comparing outputs of tools analyzing the influenza virus
# -*- coding: utf-8 -*-
# argv[1] - always bwacycle sample dir
# argv[2] - always irma sample dir
# argv[3] - insaFLU sample dir(optional)
import glob
from sys import argv
import os
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
from Bio import SeqIO
from plotnine import *

GENE_NAMES = ("HA", "M", "NA", "NP", "NS", "PA", "PB1", "PB2")
ORDERED_GENES = ("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")


class EightError(Exception):
    """" Raised when number of something is not equal to eight"""
    pass


def get_bwacycle_coverage_file_paths(path_to_sample_dir):
    return glob.glob(f"{path_to_sample_dir}/*/*.table")


def get_irma_coverage_file_paths(path_to_sample_dir):
    return glob.glob(f"{path_to_sample_dir}/tables/*coverage.txt")


def get_irma_fasta_paths(path_to_sample_dir):
    return glob.glob(f"{path_to_sample_dir}/*.fasta")


def get_bwacycle_gene_coverage(path_to_coverage_file):  # coverage_file - SampleName_GeneName.table file
    with open(path_to_coverage_file, "r") as file:
        file.readline()
        coverage = []
        for line in file:
            temp_line = line.strip().split()
            if len(temp_line) <= 9:
                continue
            coverage.append(temp_line[9])
        return coverage


def get_irma_gene_coverage(path_to_coverage_file):
    with open(path_to_coverage_file, "r") as file:
        file.readline()
        return [line.split()[2] for line in file]


def get_insaflu_gene_coverage(path_to_insaflu_dir):
    with open(f"{path_to_insaflu_dir}/coverage.tsv", "r") as coverage_file:
        for line in coverage_file:
            if line.startswith("\"Mean"):
                genes_coverage = [float(coverage.replace('"', "")) for coverage in
                                  coverage_file.readline().split()[1:9]]
                return [genes_coverage[3], genes_coverage[6], genes_coverage[5], genes_coverage[4], genes_coverage[7],
                        genes_coverage[2], genes_coverage[1], genes_coverage[0]]


def get_insa_flu_fasta(path_to_sample_dir):
    return [file for file in os.listdir(f"{path_to_sample_dir}") if file.endswith("fasta")][0]


def separate_genes(bwacycle_fasta, irma_fasta, insaflu_fasta=None, sample_name=""):
    # writes appropriate genes from n files to separate fasta file.
    subprocess.run("cat {} > {}_IRMA.fasta".format(" ".join([fasta for fasta in irma_fasta]), sample_name), shell=True)
    with open("{}_IRMA.fasta".format(sample_name), "r") as irma_file, open(bwacycle_fasta, "r") as bwacycle_file:
        if insaflu_fasta:
            insaflu_file = open(insaflu_fasta, "r")
        bwacycle_reads = [(seq_record.id, seq_record.seq) for seq_record in SeqIO.parse(bwacycle_file, format="fasta")]
        irma_reads = [(seq_record.id, seq_record.seq) for seq_record in SeqIO.parse(irma_file, format="fasta")]
        if insaflu_fasta:
            insaflu_reads = [(seq_record.id, seq_record.seq) for seq_record in
                             SeqIO.parse(insaflu_file, format="fasta")]

        if not os.path.exists("separate_genes"):
            os.mkdir("separate_genes")
        if not os.path.exists(f"separate_genes/{sample_name}"):
            os.mkdir(f"separate_genes/{sample_name}")
        for i in range(8):
            with open(f"separate_genes/{sample_name}/{sample_name}_{ORDERED_GENES[i]}.fasta", "w") as file:
                file.write(">" + str(bwacycle_reads[i][0]) + "\n")
                file.write(str(bwacycle_reads[i][1]) + "\n")
                file.write(">" + str(irma_reads[i][0]) + "\n")
                file.write(str(irma_reads[i][1]) + "\n")
                if insaflu_fasta:
                    file.write(">" + str(insaflu_reads[i][0]) + "\n")
                    file.write(str(insaflu_reads[i][1]) + "\n")


def align_all_genes(sample_name):
    if not os.path.exists(f"separate_genes/{sample_name}/alignments"):
        os.mkdir(f"separate_genes/{sample_name}/alignments")
    for gene in ORDERED_GENES:
        print(f"Aligning {gene}")
        subprocess.run(
            f"mafft-linux64/mafft.bat --auto separate_genes/{sample_name}/{sample_name}_{gene}.fasta 1>separate_genes/{sample_name}/alignments/{sample_name}_{gene}_aln.fa 2>/dev/null",
            shell=True)
    print(f"Alignments are saved to separate_genes/{sample_name}")
    if len(glob.glob(f"separate_genes/{sample_name}/alignments/*aln.fa")) != 8:
        raise EightError(
            f"The directory separate_genes/{sample_name}/alignments does not contain 8 files with alignment. Please, make sure that you deserve it."
        )


def build_gap_matrix(sample_name, insaflu=False):
    gaps = [[], [], []]  # order: bwacycle, irma, insaflu
    for i in range(8):
        with open(f"separate_genes/{sample_name}/alignments/{sample_name}_{ORDERED_GENES[i]}_aln.fa",
                  "r") as file:
            fasta_parser = SeqIO.parse(file, format="fasta")
            for tool_i in range(3 if insaflu else 2):
                seq_record = next(fasta_parser)
                gap_counter = 0
                for char in str(seq_record.seq):
                    if char != "-":
                        break
                    gap_counter += 1
                gaps[tool_i].append(gap_counter)
    with open(f"{sample_name}_gaps.tsv", "w") as gap_file:
        if insaflu:
            gap_file.write("Gene\tBWAcycle\tIRMA\tInsaFLU\n")
        else:
            gap_file.write("Gene\tBWAcycle\tIRMA\n")
        for i in range(8):
            if insaflu:
                gap_file.write(f"{ORDERED_GENES[i]}\t{gaps[0][i]}\t{gaps[1][i]}\t{gaps[2][i]}\n")
            else:
                gap_file.write(f"{ORDERED_GENES[i]}\t{gaps[0][i]}\t{gaps[1][i]}\n")
    print(f"Counted gaps are saved to {sample_name}_gaps.tsv")
    return gaps


def coverage_plot(gene_name, bwacycle_coverage, irma_coverage, insaflu_coverage=None, output_dir=".", sample_name=""):
    bwacycle_coverage = [int(i) for i in bwacycle_coverage]
    irma_coverage = [int(i) for i in irma_coverage]
    length_diff = len(bwacycle_coverage) - len(irma_coverage)
    if length_diff > 0:
        irma_coverage.extend([0 for i in range(abs(length_diff))])
    else:
        bwacycle_coverage.extend([0 for i in range(abs(length_diff))])
    # print(len(bwacycle_coverage) == len(irma_coverage))
    gene_coverage = {"Position": [i for i in range(1, len(irma_coverage) + 1)],
                     "bwacycle_coverage": bwacycle_coverage, "irma_coverage": irma_coverage}

    df = pd.DataFrame(gene_coverage)
    plt.plot("Position", "bwacycle_coverage", data=df)
    plt.plot("Position", "irma_coverage", data=df)
    if insaflu_coverage:
        plt.hlines(y=insaflu_coverage, label="INSaFLU mean", xmin=0, xmax=len(df.Position))
    plt.legend()
    plt.xlabel("Position")
    plt.ylabel("Coverage")
    plt.title(f"{gene_name} gene coverage")
    plt.savefig(f"{output_dir}/{sample_name}_{gene_name}.png")
    plt.close()


def get_bwacycle_gene_names(path_to_bwacycle_output):
    return [file for file in os.listdir(path_to_bwacycle_output) if os.path.isdir(f"{path_to_bwacycle_output}{file}")]


def get_irma_gene_names(path_to_irma_output):
    return [file.split(".")[0] for file in os.listdir(path_to_irma_output) if file.endswith("fasta")]


def make_irma_variants_table(path_to_irma_output, sample_name=""):
    variant_files = [file for file in os.listdir(f"{path_to_irma_output}/tables") if file.endswith("variants.txt")]
    variant_files.sort()
    ordered_filenames = [variant_files[7], variant_files[6], variant_files[5], variant_files[0],
                         variant_files[3], variant_files[2], variant_files[1], variant_files[4]]
    with open(f"{sample_name}_irma_variants.tsv", "w") as variants_out:
        for file in ordered_filenames:
            with open(f"{path_to_irma_output}/tables/{file}", "r") as gene_variants:
                gene_variants.readline()
                for line in gene_variants:
                    temp_line = line.strip().split()
                    variants_out.write(
                        "{}\t{}\t{}\t{}\n".format(temp_line[0], temp_line[1], temp_line[3], temp_line[4]))


def get_insaflu_variants_file_name(path_to_insaflu_dir):
    return [file for file in os.listdir(path_to_insaflu_dir) if file.endswith("variants.tsv")][0]


def check_isaflu_gene_names(path_to_insaflu_dir):
    variants_table = get_insaflu_variants_file_name(path_to_insaflu_dir)
    genes = set(GENE_NAMES)
    with open(f"{path_to_insaflu_dir}{variants_table}", "r") as insaflu_variants:
        insaflu_variants.readline()
        insaflu_genes = set()
        for line in insaflu_variants:
            insaflu_genes.add(line.strip().split()[16].split('"')[1])
    return len(insaflu_genes & genes) == len(insaflu_genes)


def build_variant_tables(bwacycle_dir, irma_dir, gaps_matrix, insaflu_dir=None, sample_name=""):
    gene_variants_tables = [[], [], [], [], [], [], [], []]
    # Write IRMA variants to the variants table
    irma_gene_names = sorted(get_irma_gene_names(irma_dir))
    ordered_irma_genes = [irma_gene_names[7], irma_gene_names[6], irma_gene_names[5], irma_gene_names[0],
                          irma_gene_names[3], irma_gene_names[2], irma_gene_names[1], irma_gene_names[4]]
    with open(f"{sample_name}_irma_variants.tsv", "r") as irma_variants:
        for line in irma_variants:
            temp_line = line.strip().split()
            gene_index = ordered_irma_genes.index(temp_line[0])
            gene_variants_tables[gene_index].append(
                "{}\t{}\t{}\t{}\tIRMA\t{}\t{}\n".format(temp_line[0], temp_line[1], temp_line[2], temp_line[3],
                                                        gaps_matrix[1][gene_index],
                                                        int(temp_line[1]) + gaps_matrix[1][gene_index]))

    # Write bwacycle variants to the variants_tables
    bwacycle_gene_names = sorted(get_bwacycle_gene_names(bwacycle_dir))
    ordered_bwacycle_genes = [bwacycle_gene_names[7], bwacycle_gene_names[6], bwacycle_gene_names[5],
                              bwacycle_gene_names[0],
                              bwacycle_gene_names[3], bwacycle_gene_names[2], bwacycle_gene_names[1],
                              bwacycle_gene_names[4]]
    with open(f"{bwacycle_dir}{sample_name}.vcf", "r") as bwacycle_variants:
        for line in bwacycle_variants:
            if line.startswith("#"):
                continue
            temp_line = line.strip().split()
            if len(temp_line[3] + temp_line[4]) != 2:
                continue
            splitted_chr = temp_line[0].split("_")
            gene_name = splitted_chr[-2] + "_" + splitted_chr[-1]
            gene_index = ordered_bwacycle_genes.index(gene_name)
            gene_variants_tables[gene_index].append(
                "{}\t{}\t{}\t{}\tBWAcycle\t{}\t{}\n".format(gene_name, temp_line[1], temp_line[3], temp_line[4],
                                                            gaps_matrix[1][gene_index],
                                                            int(temp_line[1]) + gaps_matrix[1][gene_index]))

    # Write insaflu variants to the variants table
    if insaflu_dir:
        if not check_isaflu_gene_names(insaflu_dir):
            print("WARNING! INSaFLU gene names do not match these:{}. The output will be incorrect".format(
                " ".join(ORDERED_GENES)))
        with open(f"{insaflu_dir}{get_insaflu_variants_file_name(insaflu_dir)}", "r") as insaflu_variants:
            insaflu_variants.readline()
            for line in insaflu_variants:
                temp_line = line.strip().replace('"', "").split()
                if len(temp_line[4] + temp_line[5]) != 2:
                    continue
                gene_index = ORDERED_GENES.index(temp_line[16])
                gene_variants_tables[gene_index].append(
                    "{}\t{}\t{}\t{}\tINSaFLU\t{}\t{}\n".format(temp_line[16], temp_line[2], temp_line[4], temp_line[5],
                                                               gaps_matrix[2][gene_index],
                                                               int(temp_line[2]) + gaps_matrix[2][gene_index]))

    if not os.path.exists("separate_genes/variants"):
        os.mkdir("separate_genes/variants")

    for i, gene_name in enumerate(ORDERED_GENES):
        with open(f"separate_genes/variants/{sample_name}_{gene_name}.variants.tsv", "w") as variants_table:
            variants_table.write("gene_name\tpos\tref\talt\ttool\tgap\ttrue_pos\n")
            for record in gene_variants_tables[i]:
                variants_table.write(record)
    print(f"Variant tables are saved to separate_genes/variants/")


def plot_variants(sample_name="", output_dir="."):
    variant_tables = sorted(os.listdir("separate_genes/variants"))
    if len(variant_tables) != 8:
        raise EightError("WARNING! There are not 8 tables in the separate_genes/variants directory")
    for i, table in enumerate(variant_tables):
        with open(f"separate_genes/variants/{table}", "r") as variant_table:
            df = pd.read_csv(variant_table, sep="\t")
            (ggplot()
             + geom_segment(df, aes("true_pos", xend="true_pos", y=0, yend=1, color="alt"), show_legend=False)
             + ggtitle(f"Validated SNPs for {GENE_NAMES[i]} gene")
             + facet_wrap("tool", nrow=3, ncol=1)
             + theme_bw()
             + theme(axis_ticks=element_blank(), axis_text_y=element_blank())
             + xlab("Position")
             + ylab("")).save(filename=f"{output_dir}/{sample_name}_{GENE_NAMES[i]}.png", format="png",
                              height=12, width=9, verbose=False)


def parse_args(argv):
    insaflu = len(argv) >= 4
    insaflu_fasta_path = f"{argv[3]}{get_insa_flu_fasta(argv[3])}" if insaflu else None
    bwacycle_sample_dir = argv[1]
    irma_sample_dir = argv[2]
    insaflu_sample_dir = argv[3] if insaflu else None

    bwacycle_coverage_paths = sorted(get_bwacycle_coverage_file_paths(argv[1]))
    irma_coverage_paths = sorted((get_irma_coverage_file_paths(argv[2])))
    irma_fasta_paths = sorted(get_irma_fasta_paths(argv[2]))
    bwacycle_fasta_path = f"{argv[1]}/{argv[1].split('/')[-2]}.fasta"
    ordered_irma_fasta_paths = [irma_fasta_paths[7], irma_fasta_paths[6], irma_fasta_paths[5], irma_fasta_paths[0],
                                irma_fasta_paths[3], irma_fasta_paths[2], irma_fasta_paths[1], irma_fasta_paths[4]]

    if len(bwacycle_coverage_paths) + len(irma_coverage_paths) != 16:
        raise EightError("Number of genes is not equal to 8")

    sample_name = argv[1].split('/')[-2]

    return (sample_name, insaflu, insaflu_fasta_path, insaflu_sample_dir,
            bwacycle_sample_dir, irma_sample_dir, bwacycle_coverage_paths,
            irma_coverage_paths, bwacycle_fasta_path, ordered_irma_fasta_paths)


def main(sample_name, insaflu, insaflu_fasta_path, insaflu_sample_dir,
         bwacycle_sample_dir, irma_sample_dir, bwacycle_coverage_paths,
         irma_coverage_paths, bwacycle_fasta_path, ordered_irma_fasta_paths):
    coverage_images_dir = f"{sample_name}_images/coverage"
    if not os.path.exists(coverage_images_dir):
        os.makedirs(coverage_images_dir)

    variants_dir = f"{sample_name}_images/variants"
    if not os.path.exists(variants_dir):
        os.mkdir(variants_dir)

    insaflu_coverage = get_insaflu_gene_coverage(insaflu_sample_dir) if insaflu else None

    for i in range(8):
        print(f"Plotting {GENE_NAMES[i]} gene coverage")
        insaflu_gene_coverage = insaflu_coverage[i] if insaflu else None
        bwacycle_coverage = get_bwacycle_gene_coverage(bwacycle_coverage_paths[i])
        irma_coverage = get_irma_gene_coverage(irma_coverage_paths[i])
        coverage_plot(GENE_NAMES[i], get_bwacycle_gene_coverage(bwacycle_coverage_paths[i]),
                      get_irma_gene_coverage(irma_coverage_paths[i]), output_dir=coverage_images_dir,
                      sample_name=sample_name, insaflu_coverage=insaflu_gene_coverage)

    print(f"Images are saved to {coverage_images_dir}")

    # Processing .fasta files
    print("Processing .fasta files")

    # Arrange in the accepted order (PB2, PB1, PA, HA, NP, N, M, NS)

    print("Writing genes to separate .fasta files")
    separate_genes(bwacycle_fasta_path, ordered_irma_fasta_paths, sample_name=sample_name,
                   insaflu_fasta=insaflu_fasta_path)

    align_all_genes(sample_name)

    print("Counting gaps")
    gap_matrix = build_gap_matrix(sample_name, insaflu=insaflu)

    print("Building variant tables")
    make_irma_variants_table(irma_sample_dir, sample_name=sample_name)
    build_variant_tables(bwacycle_sample_dir, irma_sample_dir, gaps_matrix=gap_matrix, sample_name=sample_name,
                         insaflu_dir=insaflu_sample_dir)

    print("Drawing variants")
    plot_variants(sample_name=sample_name, output_dir=variants_dir)
    print(f"Plots are saved to {variants_dir}")
    print("Finished!")


if __name__ == "__main__":
    main(*parse_args(argv))
