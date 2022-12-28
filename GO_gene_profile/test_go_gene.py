import pdb
from pathlib import Path
import numpy as np
import pandas as pd


go_seq_path = Path(
    "/w5home/jfreeman/GO_GENE/jira_stuff/GO_termsWith5orMoreGenes_EnrichedInUpRegBamgfpVsTkvqd_cystoblast_BP.tsv"
)
gene_length_path = Path("/w5home/jfreeman/GO_GENE/jira_stuff/gene_flydb_lengths.txt")
go_term_path = Path(
    "/w5home/jfreeman/GO_GENE/jira_stuff/drosophila_genesets_flydb_at_least_5_genes_forGOseq.txt"
)


def test_make_len_dic():
    f = open(gene_length_path, "r")
    length_dict = {}
    for line in f:
        k, v = line.strip().split()
        length_dict[k.strip()] = v.strip()
    assert length_dict["FBgn0031081"] == "3870"
    f.close()
    for k in length_dict:
        length_dict[k] = int(length_dict[k])
    assert length_dict["FBgn0031081"] == 3870
    return length_dict


def test_make_go_dic():
    f = open(go_term_path, "r")
    go_dict = {}
    for line in f:
        k = line.split()[0]
        v = line.strip().split()[1:]
        v = [x for x in v if x.startswith("FBgn") is True]
        go_dict[k] = v
    f.close()
    assert go_dict["GO:0000274"] == [
        "FBgn0016691",
        "FBgn0016120",
        "FBgn0285943",
        "FBgn0019644",
        "FBgn0036345",
    ]
    return go_dict


go_dict = test_make_go_dic()
length_dict = test_make_len_dic()


def test_merge_dict():
    a = go_dict
    b = length_dict
    c = {}
    for key in a.keys():
        c[key] = list(map(lambda x: b[x] if x in b.keys() else np.NaN, a[key]))
    assert c["GO:0000445"] == [1063, 2311, 2035, 927, 5402]
    return c


c = test_merge_dict()


def test_mean():
    means = {key: np.nanmean(value) for key, value in c.items()}
    assert means["GO:0035332"] == 4645.371428571429
    means = pd.DataFrame.from_dict(means, orient="index")
    means.reset_index(inplace=True)
    means = means.rename(columns={"index": "category", 0: "Mean_Gene_Len"})
    return means


def test_stdv():
    stdv = {key: np.nanstd(value) for key, value in c.items()}
    assert stdv["GO:0000978"] == 1806.7573407450259
    stdv = pd.DataFrame.from_dict(stdv, orient="index")
    stdv.reset_index(inplace=True)
    stdv = stdv.rename(columns={"index": "category", 0: "Stdv_Gene_Len"})
    return stdv


means = test_mean()
stdv = test_stdv()


def test_append_go_seq_out():
    pdb.set_trace()
    main = pd.read_csv(go_seq_path, sep="\t")
    merged = means.merge(stdv, how="left", on="category")
    final = main.merge(merged, how="left", on="category")
    final.to_csv("Rounded_go_seq_out_with_mean_and_stdv.tsv", sep="\t", float_format= '%.5e')
    final_rounded = final.copy()
    for col, fmt in [
        ('over_represented_pvalue', '%.5e'),
        ('Mean_Gene_Len', '%.1f')
    ]:
        final_rounded[col] = final_rounded[col].map(lambda x: fmt % x)
    '''
        'over_represented_pvalue', 'under_represented_pvalue',
       'numDEInCat', 'numInCat', 'term', 'ontology', 'fdr', 'Mean_Gene_Len',
       'Stdv_Gene_Len'
    ]#def make_len_dic(gene_lengths):
    #length_dict = {}
    #with open(gene_lengths, "r") as f:
        #for line in f:
           # k, v = line.strip().split()
           # length_dict[k] = int(v)
   # return length_dict


#def make_go_dic(go_terms):

    go_dict = {}
    with open(go_terms, "r") as f:
        for line in f:
            k = line.split()[0]
            v = [el for el in line.strip().split()[1:] if el.startswith("FBgn")]
        go_dict[k] = v
    return go_dict


#def make_merge_dict(go_dict, length_dict):
    go_dict = make_go_dic(g["go_term_genes"])
    length_dict = make_len_dic(my_args["gene_lengths"])
    merged_dict = {}
    for key in go_dict.keys():
        merged_dict[key] = list(
            map(
                lambda x: length_dict[x] if x in length_dict.keys() else np.NaN,
                go_dict[key],
            )
        )

    return merged_dict
    '''
