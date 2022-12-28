## source /isiseqruns/jfreeman_tmp_home/GO_GENE/GO_gene_profile/.venv/bin/activate
## python /isiseqruns/jfreeman_tmp_home/GO_GENE/GO_gene_profile/go_gene.py /isiseqruns/jfreeman_tmp_home/GO_GENE/GO_gene_profile/ /isiseqruns/jfreeman_tmp_home/GO_GENE/input/GO_termsWith5orMoreGenes_EnrichedInUpRegBamgfpVsTkvqd_cystoblast_BP.tsv /isiseqruns/jfreeman_tmp_home/GO_GENE/input/gene_flydb_lengths.txt /isiseqruns/jfreeman_tmp_home/GO_GENE/input/drosophila_genesets_flydb_at_least_5_genes_forGOseq.txt
import cmdlogtime
import sys
import pdb
import numpy as np
import pandas as pd
import os


COMMAND_LINE_DEF_FILE = (
    "/isiseqruns/jfreeman_tmp_home/GO_GENE/GO_gene_profile/go_gene.txt"
)


def main():

    (start_time_secs, pretty_start_time, my_args, logfile) = cmdlogtime.begin(
        COMMAND_LINE_DEF_FILE, sys.argv[0]
    )
    append_length_stats_for_go_gene(my_args)

    cmdlogtime.end(logfile, start_time_secs)


def post_process_david_output(my_args):
    pass


def make_merge_dict(go_terms, gene_lengths):
    go_dict = {}
    with open(go_terms, "r") as f:
        for line in f:
            k = line.split()[0]
            v = [el for el in line.strip().split()[1:] if el.startswith("FBgn")]
            go_dict[k] = v
    length_dict = {}
    with open(gene_lengths, "r") as f:
        for line in f:
            k, v = line.strip().split()
            length_dict[k] = int(v)
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
    means = means.rename(columns={"index": "category", 0: "Mean_Gene_Len"})
    return means


def get_stdv(merged_dict):
    stdv = {key: np.nanstd(value) for key, value in merged_dict.items()}
    stdv = pd.DataFrame.from_dict(stdv, orient="index")
    stdv.reset_index(inplace=True)
    stdv = stdv.rename(columns={"index": "category", 0: "Stdv_Gene_Len"})
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
    for col, fmt in [
        ("over_represented_pvalue", "%.5e"),
        ("under_represented_pvalue", "%.5e"),
        ("fdr", "%.5e"),
        ("Mean_Gene_Len", "%.1f"),
        ("Stdv_Gene_Len", "%.1f"),
    ]:
        final_rounded[col] = final_rounded[col].map(lambda x: fmt % x)
    return final_rounded


def append_length_stats_for_go_gene(my_args):

    merged_dict = make_merge_dict(
        go_terms=my_args["go_term_genes"], gene_lengths=my_args["gene_lengths"]
    )
    means = get_mean(merged_dict=merged_dict)
    stdv = get_stdv(merged_dict=merged_dict)
    final = append_go_seq_out(
        go_seq_out=my_args["go_seq_output"], means=means, stdv=stdv
    )

    ofile = os.path.join(my_args["out_dir"] + "/go_seq_post_processed.tsv")
    final.to_csv(ofile, sep="\t")


if __name__ == "__main__":
    main()
