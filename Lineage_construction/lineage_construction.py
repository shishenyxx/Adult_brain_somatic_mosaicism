### This python script is used to conduct lineage construction based on snMPAS and MPAS quantification in single nuclei and bulk samples
### The script is written by Xin Xu
#----import-----Python 3.7.6-----------------------------------------------------------------
import pandas as pd #version 1.0.1
import numpy as np #version 1.18.1
import seaborn as sns #version 0.10.0
import matplotlib.pyplot as plt #version 3.1.3
import networkx as nx #version 2.4
from networkx.drawing.nx_agraph import graphviz_layout
import scipy.stats #version 1.4.1
from sklearn import linear_model #version 0.22.1

#----Single cell binary profiles-------------------------------------------------------------
dbsm_data = pd.read_csv("./ipynbs/20200618_4dbsm_mos_259.csv")
#all the single cell samples
cells = ['NeuN_C2', 'NeuN_F1', 'NeuN_E3', 'NeuN_D5', 'DAPI_E12', 'DAPI_E9',
       'NeuN_C1', 'NeuN_F3', 'DAPI_H11', 'NeuN_H2', 'DAPI_H7', 'NeuN_D1',
       'NeuN_A2', 'NeuN_B5', 'NeuN_B6', 'DAPI_D7', 'DAPI_B8', 'DAPI_H9',
       'DAPI_C12', 'NeuN_E4', 'DAPI_D10', 'DAPI_C8', 'DAPI_H10',
       'NeuN_D4', 'NeuN_D2', 'NeuN_H4', 'DAPI_B12', 'NeuN_F5', 'NeuN_G5',
       'NeuN_G3', 'NeuN_E5', 'DAPI_A11', 'DAPI_A7', 'DAPI_E11', 'DAPI_G8',
       'DAPI_A9', 'NeuN_G1', 'DAPI_B10', 'NeuN_A4', 'NeuN_G2', 'DAPI_G9',
       'DAPI_D8', 'DAPI_F10', 'NeuN_E1', 'NeuN_F6', 'DAPI_C7', 'DAPI_E8',
       'NeuN_A6', 'DAPI_G7', 'DAPI_F8', 'NeuN_B4', 'NeuN_A5', 'DAPI_C9',
       'DAPI_D11', 'DAPI_A10', 'DAPI_H8', 'NeuN_F2', 'DAPI_A12',
       'DAPI_G10', 'DAPI_F9', 'DAPI_G12', 'NeuN_G6', 'DAPI_E7',
       'DAPI_B11', 'DAPI_C11', 'DAPI_B7', 'DAPI_F7', 'DAPI_D9', 'NeuN_B3',
       'NeuN_C3', 'DAPI_F12', 'NeuN_E6', 'NeuN_D3', 'NeuN_F4', 'NeuN_H6',
       'DAPI_E10']
cells = ["7614-Ctx-T-L-Sml-" + "-".join(a.split('_')) for a in cells]

single_cell_df = dbsm_data.loc[(dbsm_data.TISSUE_CELLS == "cells") &
                                (dbsm_data.SET_SN_FIGURE == 1) &
                                (dbsm_data.ID.isin(cells)) & 
                                (dbsm_data.CLADE_FIG != "ZONK"),
                                 ["CHR_POS_REF_ALT","ID", "VARIANT_ORDER_FIG"]].drop_duplicates()

variants = single_cell_df.CHR_POS_REF_ALT.unique()

var_map = {}
lineage_dict = {}
for var, num in zip(single_cell_df.CHR_POS_REF_ALT, single_cell_df.VARIANT_ORDER_FIG):
    var_map[var] = "group" + str(int(num))
    lineage_dict["group" + str(int(num))] = var

#get the single cell binary profile for each variant
single_cell_df = single_cell_df[["ID"]].drop_duplicates()
one_hot_dict = {}
for group in lineage_dict.keys():
    var = lineage_dict[group]
    data = dbsm_data.loc[dbsm_data.CHR_POS_REF_ALT == var, ["ID", "SET_CP_TL_PRESENT_05"]]
    data = pd.merge(single_cell_df, data, how="left", on="ID")
    one_hot_dict[group] = np.array(data["SET_CP_TL_PRESENT_05"])

level_hash = {}
level_dict = {}
for group in one_hot_dict.keys():
    vector = one_hot_dict[group]
    level = int(sum(vector))
    level_hash[group] = level
    if level not in level_dict.keys():
        level_dict[level] = []
    level_dict[level].append(group)

levels = sorted(list(level_dict.keys()), reverse=True)

#-----Bulk samples & sorted population maf profiles-----------------------------------------
dbsm_data = dbsm_data.loc[(~np.isnan(dbsm_data["VARIANT_ORDER_FIG"])) & 
                          (dbsm_data.TISSUE_CELLS == "tissues") &
                          (dbsm_data.DEPTH>=100) & 
                          (~dbsm_data.ID.isin(["JGG_cells", "7614_Ctx_PF_L_PU1", "7614_Ctx_F_L_PU1", 
                                          "7614_Ctx_F_R_PU1", "7614_Ctx_F_R_NeuN","7614_Ctx_P_L_Lhx2", 
                                              "7614_Ctx_P_L_PU1", "7614_Ctx_P_R_PU1", "7614_Ctx_O_L_PU1", 
                                            "7614_Ctx_O_R_PU1", "7614_Ctx_O_R_Lhx2", "7614_Ctx_T_L_Olig2", 
                                              "7614_Ctx_T_L_PU1", "7614_Ctx_T_R_Lhx2", "7614_Ctx_PF_L_TBR1"]))
                          ,]
bulk_samples = []
subsample_samples = []
sorted_samples = []
single_samples = []
for sample in dbsm_data.ID.unique():
    if sample.split("_")[-1] == "Lrg" or sample.split("_")[-1] == "Sml":
        bulk_samples.append(sample)
    elif sample.split("_")[-1].startswith("Sml"):
        subsample_samples.append(sample)
    elif ("Olig2" in sample) or ("NeuN" in sample) or ("Lhx2" in sample) or ("PU1" in sample) or ("TBR1" in sample):
        sorted_samples.append(sample)
    else:
        single_samples.append(sample)    

bulk_right_samples = []
bulk_left_samples = []
for sample in bulk_samples:
    if sample.split("_")[-2] == "L":
        bulk_left_samples.append(sample)
    elif sample.split("_")[-2] == "R":
        bulk_right_samples.append(sample)

subsample_right_samples = []
subsample_left_samples = []
for sample in subsample_samples:
    if sample.split("_")[-2] == "L":
        subsample_left_samples.append(sample)
    elif sample.split("_")[-2] == "R":
        subsample_right_samples.append(sample)

sorted_right_samples = []
sorted_left_samples = []
for sample in sorted_samples:
    if sample.split("_")[-2] == "L":
        sorted_left_samples.append(sample)
    elif sample.split("_")[-2] == "R":
        sorted_right_samples.append(sample)

def get_mean_maf(dbsm_data, sample_df):
    maf_dict = {}
    for group in lineage_dict.keys():
        variant = lineage_dict[group]
        final_data = dbsm_data.loc[dbsm_data.CHR_POS_REF_ALT == variant, ["ID", "LOWER_CI", "MAF"]]
        final_data.loc[final_data.LOWER_CI < 0.0014, "MAF"] = 0
        final_data = final_data[["ID", "MAF"]]
        final_data = pd.merge(sample_df, final_data, how="left", on="ID")
        values = final_data.MAF.astype(float)
        maf_dict[group] = final_data.MAF.values.astype(float)
    return maf_dict

large_left_df = pd.DataFrame(np.array(bulk_left_samples + subsample_left_samples).reshape((-1,1)), columns=["ID"])
large_right_df = pd.DataFrame(np.array(bulk_right_samples + subsample_right_samples).reshape((-1,1)), columns=["ID"])
sorted_left_df = pd.DataFrame(np.array(sorted_left_samples).reshape((-1,1)), columns=["ID"])
sorted_right_df = pd.DataFrame(np.array(sorted_right_samples).reshape((-1,1)), columns=["ID"])
single_df = pd.DataFrame(np.array(single_samples[:-1]).reshape((-1,1)), columns=["ID"])

large_left_dict = get_mean_maf(dbsm_data, large_left_df)
large_right_dict = get_mean_maf(dbsm_data, large_right_df)
sorted_left_dict = get_mean_maf(dbsm_data, sorted_left_df)
sorted_right_dict = get_mean_maf(dbsm_data, sorted_right_df)
single_maf_dict = get_mean_maf(dbsm_data, single_df)

lineage_maf_dict = {}
for key in large_left_dict.keys():
    mafs = np.array(large_left_dict[key])
    maf = np.mean(mafs[~np.isnan(mafs)])
    lineage_maf_dict[key] = [maf]
    
for key in large_right_dict.keys():
    mafs = np.array(large_right_dict[key])
    maf = np.mean(mafs[~np.isnan(mafs)])
    lineage_maf_dict[key].append(maf)

for key in sorted_left_dict.keys():
    mafs = np.array(sorted_left_dict[key])
    maf = np.mean(mafs[~np.isnan(mafs)])
    lineage_maf_dict[key].append(maf)
    
for key in sorted_right_dict.keys():
    mafs = np.array(sorted_right_dict[key])
    maf = np.mean(mafs[~np.isnan(mafs)])
    lineage_maf_dict[key].append(maf)

for key in single_maf_dict.keys():
    lineage_maf_dict[key] = lineage_maf_dict[key] + list(single_maf_dict[key])



#----Correlation between variants-----------------------------------
bulk_dbsm_data = dbsm_data[(dbsm_data.TISSUE_CELLS=="tissues") & (dbsm_data.ID!="JGG_Tissue") &
           (dbsm_data.DEPTH>=100) & (dbsm_data.EXPERIMENT=="bulk")]
sorted_dbsm_data = dbsm_data[(dbsm_data.TISSUE_CELLS=="tissues") & (~dbsm_data.ID.isin(["JGG_Tissue",
            "7614_Ctx_PF_L_PU1", "7614_Ctx_F_L_PU1", "7614_Ctx_F_R_PU1", "7614_Ctx_F_R_NeuN",
            "7614_Ctx_P_L_Lhx2", "7614_Ctx_P_L_PU1", "7614_Ctx_P_R_PU1", "7614_Ctx_O_L_PU1", 
            "7614_Ctx_O_R_PU1", "7614_Ctx_O_R_Lhx2", "7614_Ctx_T_L_Olig2", "7614_Ctx_T_L_PU1", 
            "7614_Ctx_T_R_Lhx2", "7614_Ctx_PF_L_TBR1"])) & (dbsm_data.DEPTH>=100) & 
                    (~dbsm_data.EXPERIMENT.isin(["bulk", "subsample"]))]

corr_dbsm_data = pd.concat([bulk_dbsm_data, sorted_dbsm_data])

def condense_maf_data(dbsm_data, variants):
    final_data = dbsm_data.loc[dbsm_data.CHR_POS_REF_ALT == variants[0], ["ID", "MAF"]]
    final_data.columns = ["ID", variants[0]]
    for variant in variants[1:len(variants)]:
        data = dbsm_data.loc[dbsm_data.CHR_POS_REF_ALT == variant, ["ID", "MAF"]]
        data.columns = ["ID", variant]
        final_data = pd.merge(final_data, data, on="ID", how="outer")
    final_data = final_data.set_index(final_data.ID)
    final_data = final_data[final_data.columns[1:]]
    return final_data 

def compute_pairwise_corr(maf_df, variants):
    coef_matrix = np.zeros([len(variants), len(variants)])
    pval_matrix = np.zeros([len(variants), len(variants)])
    for i in range(len(variants)):
        var_i = variants[i]
        mafs_i = maf_df[var_i].values
        for j in range(len(variants)):
            var_j = variants[j]
            mafs_j = maf_df[var_j].values
            data = pd.DataFrame(np.vstack([mafs_i, mafs_j]).T).dropna()
            coef, pval = scipy.stats.pearsonr(data[0], data[1])
            coef_matrix[i,j] = coef
            pval_matrix[i,j] = pval
    coef_df = pd.DataFrame(coef_matrix, columns = variants)
    coef_df = coef_df.set_index(coef_df.columns)
    pval_df = pd.DataFrame(pval_matrix, columns = variants)
    pval_df = pval_df.set_index(pval_df.columns)
    return coef_df, pval_df

corr_maf_df = condense_maf_data(corr_dbsm_data, variants)
corr_coef_df, corr_pval_df = compute_pairwise_corr(corr_maf_df, variants)


#----Lineage constraint network construction------------------------------------

def check_form_edge(group_i, group_j, maf_dict, onehot_dict, corr_coef_df, corr_pval_df,
    maf_threshold, onehot_threshold):
    var_i = lineage_dict[group_i]
    var_j = lineage_dict[group_j]
    corr_ij = corr_coef_df.loc[var_i, var_j]
    pval_ij = corr_pval_df.loc[var_i, var_j]
    if corr_ij < 0 and pval_ij < 0.05:
        return False
    else:
        vafs_i = maf_dict[group_i]
        vafs_j = maf_dict[group_j]
        onehots_i = onehot_dict[group_i]
        onehots_j = onehot_dict[group_j]
        maf_count = 0
        for vaf_i, vaf_j in zip(vafs_i, vafs_j):
            if vaf_i == np.nan or vaf_j == np.nan:
                continue
            if vaf_i < vaf_j:
                maf_count += 1
        onehot_count = 0
        for onehot_i, onehot_j in zip(onehots_i, onehots_j):
            if onehot_j == 1 and onehot_i == 0:
                onehot_count += 1
        if maf_count <= maf_threshold and onehot_count <= onehot_threshold:
            return True
        else:
            return False

def construct_lineage(lineage_dict, maf_dict, onehot_dict, 
                    corr_coef_df, corr_pval_df):
    lineage_network = nx.DiGraph()
    #add all variants to the network
    all_nodes = list(lineage_dict.keys())
    #first round of construction; check all node pairs
    for group_i in maf_dict.keys():
        for group_j in maf_dict.keys():
            if group_i == group_j:
                continue
            if check_form_edge(group_i, group_j, maf_dict, onehot_dict, 
                corr_coef_df, corr_pval_df,0, 0):
                lineage_network.add_edge(group_i, group_j)
    #remove redundunt edges for nodes that are connected to multiple parents
    parents_dict = {}
    for edge in lineage_network.edges:
        if edge[1] not in parents_dict.keys():
            parents_dict[edge[1]] = []
        parents_dict[edge[1]].append(edge[0])
    for node in parents_dict.keys():
        parent_levels = []
        if len(parents_dict[node]) > 1:
            parents = parents_dict[node]
            for parent in parents:
                parent_level= level_hash[parent]
                parent_levels.append([parent_level, parent])
        for item in sorted(parent_levels)[:-1]:
            lineage_network.remove_edge(item[1], node)

    #check the nodes that are not connected to any parent (loosing up the threshold)
    all_children = [e[1] for e in lineage_network.edges]
    for node in all_nodes:
        if (node not in all_children):
            level = level_hash[node]
            parent_levels = sorted(np.array(levels)[np.array(levels) > level])
            for parent_level in parent_levels:
                parent_nodes = level_dict[parent_level]
                connected= False
                for parent_node in parent_nodes:
                    if check_form_edge(parent_node, node, maf_dict, onehot_dict, 
                        corr_coef_df, corr_pval_df, 1, 0):
                        lineage_network.add_edge(parent_node, node)
                        connected = True
                if connected:
                    break
    #check the nodes that are not connected to any parent (loosing up the threshold further)
    all_children = [e[1] for e in lineage_network.edges]
    for node in all_nodes:
        if (node not in all_children):
            level = level_hash[node]
            parent_levels = sorted(np.array(levels)[np.array(levels) > level])
            for parent_level in parent_levels:
                parent_nodes = level_dict[parent_level]
                connected= False
                for parent_node in parent_nodes:
                    if check_form_edge(parent_node, node, maf_dict, onehot_dict, 
                        corr_coef_df, corr_pval_df, 1, 1):
                        lineage_network.add_edge(parent_node, node)
                        connected = True
                if connected:
                    break
    #connect every subtree to "root"
    all_children = [e[1] for e in lineage_network.edges]
    for node in all_nodes:
        if (node not in all_children):
            lineage_network.add_edge("root", node)

    return lineage_network

def draw_network(network, path=None):
    new_network = nx.DiGraph()
    for edge in network.edges:
        node_1 = edge[0]
        if node_1 != "root":
            node_1 = int(node_1[5:])
        node_2 = int(edge[1][5:])
        new_network.add_edge(node_1, node_2)
    pos =graphviz_layout(new_network, prog='dot')
    nx.draw(new_network, pos, with_labels=True, arrows=True)
    if path!= None:
        plt.savefig(path)
    plt.show()


lineage_network = construct_lineage(lineage_dict, lineage_maf_dict, one_hot_dict, 
                                    corr_coef_df, corr_pval_df)

draw_network(lineage_network, path="./lineage_network.pdf")

#----Compute contribution from each clade--------------------------------------

#map nodes names to number for easy indexing
new_lineage_dict = {}
new_lineage_dict[0] = "root"
new_var_dict = {}
new_var_dict["root"] = 0
i = 1
for node in lineage_network.nodes:
    if node == "root":
        continue
    var = lineage_dict[node]
    new_lineage_dict[i] = var
    new_var_dict[var] = i
    i += 1

#replace node names in the tree with the indices
final_tree = nx.DiGraph()
for edge in lineage_network.edges():
    if edge[0] == "root":
        node_i = 0
    else:
        node_i = new_var_dict[lineage_dict[edge[0]]]
    node_j = new_var_dict[lineage_dict[edge[1]]]
    final_tree.add_edge(node_i, node_j)

#function to code 0/1 for lineage
def traverse_lineage(tree, node):
    one_hot_coding = lineage_one_hot_dict[node]
    if len(tree.edges(node)) == 0:
        return
    for edge in tree.edges(node):
        child = edge[1]
        this_one_hot = one_hot_coding.copy()
        this_one_hot[child] = 1
        lineage_one_hot_dict[child] = this_one_hot
        traverse_lineage(tree, child)
    return

#code lineage in by traversing the tree
lineage_one_hot_dict = {}
root = 0
one_hot_coding = [0]*30
lineage_one_hot_dict[root] = one_hot_coding
traverse_lineage(final_tree, root)

#concatenate lineage_one_hot_dict to numpy array for regression
lineage_one_hot = []
for i in range(30):
    lineage_one_hot.append(lineage_one_hot_dict[i])
lineage_one_hot = np.vstack(lineage_one_hot) 
lineage_one_hot = lineage_one_hot[1:, 1:]


def compute_sample_means(samples, dbsm_data):
    sample_means = []
    for sample in samples:
        group_means = []
        for i in range(1,30):
            variant = new_lineage_dict[i]
            data = dbsm_data.loc[(dbsm_data.CHR_POS_REF_ALT == variant )&
                                        (dbsm_data.ID == sample), 
                                        ["LOWER_CI", "MAF"]]
            data.loc[data.LOWER_CI < 0.0014, "MAF"] = 0
            group_mean = np.mean(data.MAF) * 2
            group_means.append(group_mean)
        sample_means.append(group_means)
    sample_means = np.vstack(sample_means)
    return sample_means


sorted_dbsm_data = dbsm_data[(dbsm_data.TISSUE_CELLS=="tissues") & (~dbsm_data.ID.isin(["JGG_Tissue",
            "7614_Ctx_PF_L_PU1", "7614_Ctx_F_L_PU1", "7614_Ctx_F_R_PU1", "7614_Ctx_F_R_NeuN",
            "7614_Ctx_P_L_Lhx2", "7614_Ctx_P_L_PU1", "7614_Ctx_P_R_PU1", "7614_Ctx_O_L_PU1", 
            "7614_Ctx_O_R_PU1", "7614_Ctx_O_R_Lhx2", "7614_Ctx_T_L_Olig2", "7614_Ctx_T_L_PU1", 
            "7614_Ctx_T_R_Lhx2", "7614_Ctx_PF_L_TBR1"])) & (dbsm_data.DEPTH>=100) & 
                    (~dbsm_data.EXPERIMENT.isin(["bulk", "subsample"]))]
sorted_samples = sorted_dbsm_data.ID.unique()
sorted_sample_df = pd.DataFrame(sorted_samples, columns=["ID"])               


large_dbsm_data = dbsm_data[(dbsm_data.TISSUE_CELLS=="tissues") & (dbsm_data.ID!="JGG_Tissue") &
          (dbsm_data.DEPTH>=100) & (dbsm_data.EXPERIMENT.isin(["bulk", "subsample"]))]
large_samples = large_dbsm_data.ID.unique()
large_sample_df =  pd.DataFrame(large_samples, columns=["ID"])

sorted_sample_means = compute_sample_means(sorted_samples, sorted_dbsm_data)
large_sample_means = compute_sample_means(large_samples, large_dbsm_data)


def regression_compute_coefs(sample_means, samples):
    coef_dict = {}
    for sample_values, sample in zip(sample_means, samples):
        if sum(np.isnan(sample_values)) > 0:
            not_nan = ~np.isnan(sample_values)
            sample_values = sample_values[not_nan]
            compute_lineage_one_hot = lineage_one_hot[:,not_nan].copy()   
        else:
            compute_lineage_one_hot = lineage_one_hot.copy()
        clf = linear_model.Lasso(alpha=0, fit_intercept=False, max_iter=10000, positive=True)
        clf.fit(compute_lineage_one_hot.T, sample_values)
        sample_coef = clf.coef_
        diff = np.matmul(compute_lineage_one_hot.T, sample_coef.reshape(-1,1)).squeeze() - sample_values
        sample_mse = sum(diff*diff)
        coef_dict[sample] = [sample_coef, sample_mse]
    return coef_dict

def coef_dict_to_df(samples, coef_dict):
    coef_np = []
    for i, sample in enumerate(samples):
        coefs = coef_dict[sample][0]
        res = 1 - sum(coefs)
        coef_np.append([res] + list(coef_dict[sample][0]))
    coef_np = np.array(coef_np)
    variants = ["germline"] + [new_lineage_dict[a] for a in range(1,30)]
    coef_df = pd.DataFrame(coef_np, columns=variants)
    coef_df["SAMPLE"] = samples
    coef_df = coef_df[["SAMPLE"] + variants]
    return coef_df

sorted_coef_dict=regression_compute_coefs(sorted_sample_means, sorted_samples)
large_coef_dict=regression_compute_coefs(large_sample_means, large_samples)

sorted_coef_df = coef_dict_to_df(sorted_samples, sorted_coef_dict)
large_coef_df = coef_dict_to_df(large_samples, large_coef_dict)

sorted_coef_df.to_csv("sorted_group_29_coef_07_31_2020.csv", index=None)
large_coef_df.to_csv("large_group_29_coef_07_31_2020.csv", index=None)

def compute_composition(target_to_sorted, coef_dict):
    comp_samples = list(target_to_sorted.keys())
    suffix_df = pd.DataFrame(np.array(["NeuN", "Olig2", "Lhx2", "PU1", "TBR1"]).reshape((-1,1)), columns=["suffix"])
    composition_data = []
    composition_dict = {}
    for sample in comp_samples:
        coefs = coef_dict[sample][0]
        to_samples = target_to_sorted[sample]
        sorted_coefs = []
        for sorted_sample in to_samples:
            sorted_coefs.append(sorted_coef_dict[sorted_sample][0])
        sorted_coefs = np.vstack(sorted_coefs)
        clf = linear_model.Lasso(alpha=0, fit_intercept=True, max_iter=10000, positive=True)
        clf.fit(sorted_coefs.T, coefs)
        diff = np.matmul(sorted_coefs.T, clf.coef_) - coefs
        sample_mse = sum(diff*diff)
        composition_dict[sample] = [to_samples, sorted_coefs, clf.coef_, sample_mse]

        samples = [a.split("_")[-1] for a in to_samples]
        composition_df = pd.DataFrame(np.vstack([samples, clf.coef_]).T, columns=["suffix","coef"])
        composition_df = pd.merge(suffix_df, composition_df, on="suffix",how="left")
        composition_data.append(list(composition_df.coef.values))
    composition_data = np.vstack(composition_data)
    composition_df = pd.DataFrame(composition_data, columns = ["NeuN", "Olig2", "Lhx2", "PU1", "TBR1"])
    composition_df["SAMPLE"] = comp_samples
    composition_df = composition_df[["SAMPLE","NeuN", "Olig2", "Lhx2", "PU1", "TBR1"]]
    return composition_df, composition_dict

large_to_sorted = {}
for large_sample in large_samples:
    if "Ctx" in large_sample:
        large_sorted_pop = []
        for sorted_sample in sorted_samples:
            if sorted_sample.startswith("_".join(large_sample.split("_")[:-1])):
                large_sorted_pop.append(sorted_sample)
        large_to_sorted[large_sample] = large_sorted_pop

large_composition_df, large_composition_dict = compute_composition(large_to_sorted, large_coef_dict)

large_composition_df.to_csv("./bulk_sorted_composition_07_28_2020.csv")
