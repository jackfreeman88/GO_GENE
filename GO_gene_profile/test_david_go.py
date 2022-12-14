import pdb
from pathlib import Path
import numpy as np
import pandas as pd

## Capitalize
DavidAnalysis_GO_BP_and_Kegg_Tkvqd_GSC_likeHighBMP_downRegBamGFP = Path(
    "/w5home/jfreeman/GO_GENE/jira_stuff/DavidAnalysis_GO_BP_and_Kegg_Tkvqd_GSC_likeHighBMP_downRegBamGFP.txt"
)
DavidAnalysis_GO_BP_and_Kegg_UpRegBamKD_Tkvqd_GSClikeLowBMP = Path(
    "/w5home/jfreeman/GO_GENE/jira_stuff/DavidAnalysis_GO_BP_and_Kegg_UpRegBamKD_Tkvqd_GSClikeLowBMP.txt"
)
DavidAnalysis_WandA_cystoblast_keggAndGO_BP = Path(
    "/w5home/jfreeman/GO_GENE/jira_stuff/DavidAnalysis_WandA_cystoblast_keggAndGO_BP.txt"
)
DAVIDanalysisDownRegBamKD_GSClikeHighBMP_GO_BP_and_Kegg = Path(
    "/w5home/jfreeman/GO_GENE/jira_stuff/DAVIDanalysisDownRegBamKD_GSClikeHighBMP_GO_BP_and_Kegg.txt"
)
go_term_path = Path(
    "/w5home/jfreeman/GO_GENE/jira_stuff/drosophila_genesets_flydb_at_least_5_genes_forGOseq.txt"
)
gene_flydb_lengths = Path("/w5home/jfreeman/GO_GENE/jira_stuff/gene_flydb_lengths.txt")


def test_make_len_dic():
    f = open(gene_flydb_lengths, "r")
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


length_dict = test_make_len_dic()


def test_trim_david():
    f = open(DAVIDanalysisDownRegBamKD_GSClikeHighBMP_GO_BP_and_Kegg, "r")
    d_go_dict = {}
    next(f)
    for line in f:
        d_go = line.split()[1]
        k = d_go.split("~")[0]
        v = line.strip().split()[2:]
        for item in v:
            if item.endswith(",") is True:
                item = item[:-1]
                v.append(item)
        v = [x for x in v if x.startswith("FBGN") is True and x.endswith(",") is False]
        v_filt = []
        for item in v:
            item = item.replace("GN", "gn")
            v_filt.append(item)
        d_go_dict[k] = v_filt
    f.close()

    return d_go_dict


d_go_dict = test_trim_david()


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


def test_updated_dict():
    updated_dict = dict()
    for key in d_go_dict:
        try:
            updated_dict[key] = go_dict[key]
        except KeyError:
            updated_dict[key] = np.NaN


    pdb.set_trace()
    return updated_dict


def test_d_merge_dict():
    a = d_go_dict
    b = length_dict
    d_merge = {}
    for key in a.keys():
        d_merge[key] = list(map(lambda x: b[x] if x in b.keys() else np.NaN, a[key]))

    return d_merge


d_merge = test_d_merge_dict()


def test_d_mean():
    means = {key: np.nanmean(value) for key, value in d_merge.items()}

    means = pd.DataFrame.from_dict(means, orient="index")
    means.reset_index(inplace=True)
    d_means = means.rename(columns={"index": "Term", 0: "Mean_Gene_Len"})
    return d_means


def test_d_stdv():
    stdv = {key: np.nanstd(value) for key, value in d_merge.items()}
    # assert stdv['GO:0000978'] == 1806.7573407450259
    stdv = pd.DataFrame.from_dict(stdv, orient="index")
    stdv.reset_index(inplace=True)
    d_stdv = stdv.rename(columns={"index": "Term", 0: "Stdv_Gene_Len"})
    return d_stdv


d_means = test_d_mean()
d_stdv = test_d_stdv()


def test_append_d_go_seq_out():
    main = pd.read_table(
        DAVIDanalysisDownRegBamKD_GSClikeHighBMP_GO_BP_and_Kegg, sep="\t"
    )
    main["Term"] = main["Term"].str.split("~")
    main[["Term", "Description"]] = pd.DataFrame(main.Term.tolist(), index=main.index)
    main.insert(2, "Description", main.pop("Description"))
    merged = d_means.merge(d_stdv, how="left", on="Term")
    final = main.merge(merged, how="left", on="Term")
    final_rounded = final.round(5)
    final_rounded.to_csv(
        "DAVIDanalysisDownRegBamKD_GSClikeHighBMP_GO_BP_and_Kegg.txt",
        sep="\t",
    )


"""
    kegg = final.copy()
    final = final.loc[final["Mean_Gene_Len"].notnull()]
    kegg = kegg.loc[kegg["Mean_Gene_Len"].isnull()]
    kegg["Term"] = kegg["Term"].str.split()

    kegg[["Term", "Description", "f", "g", "i"]] = pd.DataFrame(
        kegg.Term.tolist(), index=kegg.index
    )
    fkegg = kegg.merge(merged, how="left", on="Term")
    fkegg.pop("Mean_Gene_Len_x")
    fkegg.pop("f")
    fkegg.pop("g")
    fkegg.pop("i")
    fkegg.pop("Stdv_Gene_Len_x")

    kegg = fkegg.rename(
        columns={"Mean_Gene_Len_y": "Mean_Gene_Len", "Stdv_Gene_Len_y": "Stdv_Gene_Len"}
    )
    final = final.append(kegg, ignore_index=True)"""
