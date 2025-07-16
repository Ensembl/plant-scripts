#!/usr/bin/env python
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Merge Genes script.

This script aims to merge the genes from Nipponbare 3 annotation sources (RAPDB, MSU and Gramene).

Below are some of the criteria to consider:
1. (plus) and (minus) strand genes are considered separately 
2. The genes are considered as overlapping if the overlap length of the sequence is >=50% of the length of the smallest line.

"""
import argparse
import logging
from pathlib import Path

import pandas

from ensembl.utils import StrPath
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def update_gene_info(
    update_gene: dict,
    id_list: list,
    merged_chromosome: str,
    merged_gene_start: str,
    merged_gene_end: str,
    merged_strand: str,
) -> dict:
    """
    Takes the id_list, other parameters from main, update the gene information and returns to the main

    Args:
        update_gene (dict): column names as keys and corresponding parameters are values in teh form of list
        id_list (list): A list of ids from genes to merge
        merged_chromosome: Chromosome number from the merged gene
        merged_gene_start: New gene start coordinate after merging the gene
        merged_gene_end: New gene end coordinate after merging the gene
        merged_strand: Strand from the merged genes

    Returns:
        The updated gene.
    """
    update_gene["chr"].append(merged_chromosome)
    update_gene["source"].append("panoryza")
    update_gene["gene"].append("gene")
    update_gene["start"].append(merged_gene_start)
    update_gene["end"].append(merged_gene_end)
    update_gene["score"].append(".")
    update_gene["strand"].append(merged_strand)
    update_gene["frame"].append(".")
    merge_ids = ",".join(id_list)
    update_gene["attribute"].append(merge_ids)
    return update_gene


def sort_genes(output_file: dict) -> pandas.DataFrame:
    """Sort genes.

    Args:
        output_file: gene information in a dictionary
    Returns:
        A sorted gene matrix.

    """
    gene_matrix = pandas.DataFrame.from_dict(output_file)
    gene_matrix["chr_no"] = gene_matrix["chr"].str.extract("(\d+)", expand=False).astype(int)
    gene_matrix.sort_values(by=["chr_no", "start"], inplace=True)
    gene_matrix.drop("chr_no", axis=1, inplace=True)
    return gene_matrix


def merge_genes(input_gff_file: StrPath) -> dict:
    """Merge genes from the input GFF file.

    Args:
        input_gff_file: Path to the GFF file.

    Returns:
        The merged genes in a dictionary.
    """

    list_id = []
    dict_merged_gene = {
        "chr": [],
        "source": [],
        "gene": [],
        "start": [],
        "end": [],
        "score": [],
        "strand": [],
        "frame": [],
        "attribute": [],
    }

    logging.info("The input file is " + str(input_gff_file))

    with open(input_gff_file) as gff_file:
        lines = gff_file.readlines()
        ## if line is starting with 'Chr' then split it and assign to different variables
        for line in lines:
            if line.startswith("#"):
                continue
            chromo, _, _, start, end, _, strand, _, description = line.split("\t")
            gene_id = description.split(";Name=")[0].split(";biotype=")[0]
            start = int(start)
            end = int(end)

            ### Assign the values for the first gene in the gff
            if line == lines[0]:
                gene_start = start
                gene_end = end
                current_chromosome = str(line.split("\t")[0])
                current_strand = str(line.split("\t")[6])

            ## Checking if the strand and chromosome are same between the previous and current line; otherwise print and start with the new gene
            if strand != current_strand or chromo != current_chromosome:
                dict_merged_gene = update_gene_info(
                    dict_merged_gene, list_id, current_chromosome, gene_start, gene_end, current_strand
                )
                current_strand = strand
                gene_start = start
                gene_end = end
                current_chromosome = chromo
                list_id = []

            ## if start of the next gene is lesser than the current gene, raise exception
            if start < gene_start:
                raise Exception("Error: File is not sorted properly")

            ## if start coordinate of a gene is equal or greater than the previous AND if start coordinate is lesser than the gene end then start merging the genes
            elif (start == gene_start or start > gene_start) and start < gene_end:
                ## if end coord of current gene is lesser than the previous, then merge the gene into the previous one
                if end <= gene_end:
                    list_id.append(gene_id)

                ## if end coord of current gene is greater than the previous, then calculate the length of both genes and check overlap criteria
                elif end > gene_end:
                    gene_length = gene_end - gene_start
                    length = end - start
                    min_length = min(length, gene_length)
                    overlap_length = gene_end - start

                    ## Checking the overlap criteria
                    if int(overlap_length) * 2 >= int(min_length):
                        gene_end = end
                        list_id.append(gene_id)
                        logging.info("overlap\n")
                    else:
                        dict_merged_gene = update_gene_info(
                            dict_merged_gene,
                            list_id,
                            current_chromosome,
                            gene_start,
                            gene_end,
                            current_strand,
                        )
                        list_id = []
                        gene_start = start
                        gene_end = end
                        list_id.append(gene_id)

            ## if start coord of a gene is greater than the end coord of the previous gene then print the previous gene (current gene will not merge into that)
            if start >= gene_end:
                dict_merged_gene = update_gene_info(
                    dict_merged_gene, list_id, current_chromosome, gene_start, gene_end, current_strand
                )
                list_id = []
                gene_start = start
                gene_end = end
                list_id.append(gene_id)

            ## If last line in the file, print the gene(s)
            if line == lines[-1]:
                dict_merged_gene = update_gene_info(
                    dict_merged_gene, list_id, current_chromosome, gene_start, gene_end, current_strand
                )
    return dict_merged_gene


def parse_args(arg_list: list[str] | None) -> argparse.Namespace:
    """Parse the aruments from the command line.

    Args:
        arg_list: List of arguments to parse. If `None`, grab them from the command line.

    Returns:
        Parsed arguments.

    """
    parser = ArgumentParser(description="Merge genes script")
    parser.add_argument_src_path(
        "--input_gff_file",
        default="sorted_all_line.gff",
        help="Input file with gene info from RAPDB, MSU and Gramene",
    )
    parser.add_argument_dst_path(
        "--output_dir", default=".", help="Output directory for the merged genes file"
    )
    parser.add_log_arguments()
    return parser.parse_args(arg_list)


def main(arg_list: list[str] | None = None) -> None:
    """Main script entry point.

    Args:
        arg_list: Arguments to parse passing list to `parse_args()`.
    """
    args = parse_args(arg_list)
    init_logging_with_args(args)

    merge_genes_metadata = merge_genes(args.input_gff_file)
    sorted_genes = sort_genes(merge_genes_metadata)
    sorted_genes.to_csv(args.output_dir / "merged_genes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
