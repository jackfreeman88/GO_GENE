import pdb
from pathlib import Path
import numpy as np
import pandas as pd


WandA_categoriesEnriched_GSCLike_GFP_highBMP = Path(
    "/w5home/jfreeman/GO_GENE/jira_stuff/WandA_categoriesEnriched_GSCLike_GFP_highBMP.txt"
)
WandA_Cystoblast_enrichedGOTerms = Path(
    "/w5home/jfreeman/GO_GENE/jira_stuff/WandA_Cystoblast_enrichedGOTerms.txt"
)
WandA_enrichedGO_terms_GSClike_lowBMP_upRegBamKD = Path(
    "/w5home/jfreeman/GO_GENE/jira_stuff/WandA_enrichedGO_terms_GSClike_lowBMP_upRegBamKD.txt"
)
WandA_GSCLikeHIghBMPsignalingDownRegBamKD_enrichedGO_terms = Path(
    "/w5home/jfreeman/GO_GENE/jira_stuff/WandA_GSCLikeHIghBMPsignalingDownRegBamKD_enrichedGO_terms.txt"
)

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
    pdb.set_trace()
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


def test_WandA_merge():
    keys1 = [
        "GO:0045944",
        "GO:0006325",
        "GO:0098609",
        "GO:0008406",
        "GO:0016477",
        "GO:0043005",
    ]
    keys2 = ["GO:0032543", "GO:0006397", "GO:0006457"]
    keys3 = ["GO:0006457", "GO:0035080", "GO:0070482"]
    keys4 = ["GO:0007399", "GO:0030154", "GO:0006357", "GO:0050808"]
    enriched_GSC_like_high_bmp_dict = {x: c[x] for x in keys1}
    cytoblast_enriched_dict = {x: c[x] for x in keys2}
    enriched_GSC_like_low_bmp_dict = {x: c[x] for x in keys3}
    gsc_high_bmp_down_reg = {x: c[x] for x in keys4}

    return enriched_GSC_like_high_bmp_dict


WANDA = test_WandA_merge()


def test_wanda_mean():
    means = {key: np.nanmean(value) for key, value in WANDA.items()}
    means = pd.DataFrame.from_dict(means, orient="index")
    means.reset_index(inplace=True)
    wanda_means = means.rename(columns={"index": "category", 0: "Mean_Gene_Len"})
    return wanda_means


wanda_means = test_wanda_mean()


def test_wanda_stdv():
    stdv = {key: np.nanstd(value) for key, value in WANDA.items()}
    stdv = pd.DataFrame.from_dict(stdv, orient="index")
    stdv.reset_index(inplace=True)
    wanda_stdv = stdv.rename(columns={"index": "category", 0: "Stdv_Gene_Len"})
    return wanda_stdv


wanda_stdv = test_wanda_stdv()


def test_append_wanda_go_seq_out():
    main = pd.read_table(
        WandA_categoriesEnriched_GSCLike_GFP_highBMP,
        sep="\t",
        names=["category", "description"],
    )
    merged = wanda_means.merge(wanda_stdv, how="left", on="category")
    final = main.merge(merged, how="left", on="category")
    final_rounded = final.round(5)
    final_rounded.to_csv("FRRWandA_categoriesEnriched_GSCLike_GFP_highBMP", sep="\t")
