## source /isiseqruns/jfreeman_tmp_home/GO_GENE/GO_gene_profile/.venv/bin/activate

import cmdlogtime
import sys
import pdb
from pathlib import Path
import numpy as np
import pandas as pd


GO_SEQ_OUT = Path(
    "input/GO_termsWith5orMoreGenes_EnrichedInUpRegBamgfpVsTkvqd_cystoblast_BP.tsv"
)
GENE_LENGTHS = Path("input/gene_flydb_lengths.txt")
GO_TERMS = Path("input/drosophila_genesets_flydb_at_least_5_genes_forGOseq.txt")
COMMAND_LINE_DEF_FILE = "./go_gene.txt"


# There is also an example python file that uses all this stuff in this directory.
def main():
    # Call cmdlogtime.begin() at the beginning of main()
    (start_time_secs, pretty_start_time, my_args, logfile) = cmdlogtime.begin(
        COMMAND_LINE_DEF_FILE, sys.argv[0]
    )

    go_gene(my_args)
    # the command line arguments will be in the my_args dictionary returned, so you can access them like this:
    # get_treats_vectors = my_args["get_treats_vectors"]

    #  Then you put all of your code here.....

    # if you want to add stuff to the logfile:
    # logfile.write("Run Skim here for each of the B terms files in B_TERMS_DIR")

    # The only functions you probably will ever need in cmdlogtime are begin(), end(), and maybe make_dir. Here's an example:
    # cmdlogtime.make_dir(intermediate_out_dir)

    # Call cmdlogtime.end() at the end of main()
    cmdlogtime.end(logfile, start_time_secs)


def make_len_dic(gene_lengths):
    length_dict = {}
    with open(gene_lengths, "r") as f:
        for line in f:
            k, v = line.strip().split()
            length_dict[k] = int(v)
    return length_dict


def make_go_dic(go_terms):

    go_dict = {}
    with open(go_terms, "r") as f:
        for line in f:
            k = line.split()[0]
            v = [el for el in line.strip().split()[1:] if el.startswith("FBgn")]
        go_dict[k] = v
    return go_dict


def make_merge_dict(go_dict, length_dict):
    go_dict = make_go_dic(GO_TERMS)
    length_dict = make_len_dic(GENE_LENGTHS)
    merged_dict = {}
    for key in go_dict.keys():
        merged_dict[key] = list(
            map(
                lambda x: length_dict[x] if x in length_dict.keys() else np.NaN,
                go_dict[key],
            )
        )

    return merged_dict


def get_mean(merged_dict):
    means = {key: np.nanmean(value) for key, value in merged_dict.items()}
    means = pd.DataFrame.from_dict(means, orient="index")
    means.reset_index(inplace=True)
    means = means.rename(columns={"index": "category", 0: "means"})
    return means


def get_stdv(merged_dict):
    stdv = {key: np.nanstd(value) for key, value in merged_dict.items()}
    stdv = pd.DataFrame.from_dict(stdv, orient="index")
    stdv.reset_index(inplace=True)
    stdv = stdv.rename(columns={"index": "category", 0: "stdv"})
    return stdv


def append_go_seq_out(go_seq_out, means, stdv):
    main = pd.read_csv(go_seq_out, sep="\t")
    merged = means.merge(
        stdv,
        how="left",
        on="category",
    )

    final = main.merge(merged, how="left", on="category")
    final_rounded = final.copy()
    for col, fmt in [("over_represented_pvalue", "%.5e"), ("Mean_Gene_Len", "%.1f")]:
        final_rounded[col] = final_rounded[col].map(lambda x: fmt % x)
    return final_rounded


def go_gene(my_args):
    go_dict = make_go_dic(GO_TERMS)
    length_dict = make_len_dic(GENE_LENGTHS)
    merged_dict = make_merge_dict(go_dict, length_dict)
    means = get_mean(merged_dict=merged_dict)
    stdv = get_stdv(merged_dict=merged_dict)
    append_go_seq_out(GO_SEQ_OUT, means=means, stdv=stdv)

    # output = calc_output(my_args["GO_SEQ_OUT"], my_args["GO_SEQ_OUT"], my_args["go_term_genes"])
    # print(output)
    # with open(my_args["out_dir"] + "/go_gene.out", "w") as f:
    # f.write(output)


if __name__ == "__main__":
    main()
