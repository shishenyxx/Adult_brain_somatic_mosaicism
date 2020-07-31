# -*- coding: utf-8 -*-
"""
Created on May 20th 2020
Last modified on July 17th 2020
​
@author: Xin Xu
"""


#Import modules (Python version 3.7.6)
import pandas as pd #version 1.0.1
import numpy as np #version 1.18.1
import scipy.stats #version 1.4.1 
import matplotlib.pyplot as plt #version 3.1.3
import seaborn as sns #version 0.10.0

plt.rcParams['svg.fonttype'] = 'none'
plt.ioff()
​
#Import modules done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Functions

#summarize maf
def condense_maf_data(dbsm_data):
    final_data = dbsm_data.loc[dbsm_data.CHR_POS_REF_ALT == variants[0], ["ID", "MAF"]]
    final_data.columns = ["ID", variants[0]]
    for variant in variants[1:len(variants)]:
        data = dbsm_data.loc[dbsm_data.CHR_POS_REF_ALT == variant, ["ID", "MAF"]]
        data.columns = ["ID", variant]
        final_data = pd.merge(final_data, data, on="ID", how="outer")
    final_data = final_data.set_index(final_data.ID)
    final_data = final_data[final_data.columns[1:]]
    return final_data 

#compute correlations among all variant pairs
def compute_pairwise_corr(maf_df):
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
    coef_df = coef_df.set_index(variants)
    pval_df = pd.DataFrame(pval_matrix, columns = variants)
    pval_df = pval_df.set_index(variants)
    return coef_df, pval_df

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define color coding
pf_color_dict = {True: "#5C3763", False: "#FFFFFF"}
f_color_dict = {True: "#88527B", False: "#FFFFFF"}
p_color_dict = {True: "#AF738D", False: "#FFFFFF"}
o_color_dict = {True: "#CF9DA4", False: "#FFFFFF"}
t_color_dict = {True: "#E8CCC7", False: "#FFFFFF"}

cerebellum_color_dict = {True: "#173E20", False: "#FFFFFF"}
heart_color_dict = {True: "#3D7247", False: "#FFFFFF"}
liver_color_dict = {True: "#679D72", False: "#FFFFFF"}
kidney_color_dict = {True: "#A0C7A7", False: "#FFFFFF"}


laterization_color_dict = {'Both': "#C7C7C7",
                           'Left_only': "#2472A9",
                           'Right_only': "#EE7C1B",
                           'Left_enriched': "#A6CEE3",
                           'Right_enriched': "#FDBF6F",
                           'Not_lateral_bulk': "#919191"}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Load Data
dbsm_data = pd.read_csv("./20200618_4dbsm_mos_259.csv")

#value mapping between variants and some columns
pf_dict = dict(zip(dbsm_data.CHR_POS_REF_ALT, dbsm_data.IN_PF))
f_dict = dict(zip(dbsm_data.CHR_POS_REF_ALT, dbsm_data.IN_F))
p_dict = dict(zip(dbsm_data.CHR_POS_REF_ALT, dbsm_data.IN_P))
o_dict = dict(zip(dbsm_data.CHR_POS_REF_ALT, dbsm_data.IN_O))
t_dict = dict(zip(dbsm_data.CHR_POS_REF_ALT, dbsm_data.IN_T))
cerebellum_dict = dict(zip(dbsm_data.CHR_POS_REF_ALT, dbsm_data.SET_IN_CEREBELLUM))
heart_dict = dict(zip(dbsm_data.CHR_POS_REF_ALT, dbsm_data.SET_IN_HEART))
liver_dict = dict(zip(dbsm_data.CHR_POS_REF_ALT, dbsm_data.SET_IN_LIVER))
kidney_dict = dict(zip(dbsm_data.CHR_POS_REF_ALT, dbsm_data.SET_IN_KIDNEY))
laterization_dict =  dict(zip(dbsm_data.CHR_POS_REF_ALT, dbsm_data.LATERALIZATION))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#parse data
bulk_dbsm_data = dbsm_data[(dbsm_data.TISSUE_CELLS=="tissues") & (dbsm_data.ID!="JGG_Tissue") &
           (dbsm_data.DEPTH>=100) & (dbsm_data.EXPERIMENT=="bulk")]

sorted_dbsm_data = dbsm_data[(dbsm_data.TISSUE_CELLS=="tissues") & (~dbsm_data.ID.isin(["JGG_Tissue",
            "7614_Ctx_PF_L_PU1", "7614_Ctx_F_L_PU1", "7614_Ctx_F_R_PU1", "7614_Ctx_F_R_NeuN",
            "7614_Ctx_P_L_Lhx2", "7614_Ctx_P_L_PU1", "7614_Ctx_P_R_PU1", "7614_Ctx_O_L_PU1", 
            "7614_Ctx_O_R_PU1", "7614_Ctx_O_R_Lhx2", "7614_Ctx_T_L_Olig2", "7614_Ctx_T_L_PU1", 
            "7614_Ctx_T_R_Lhx2", "7614_Ctx_PF_L_TBR1"])) & (dbsm_data.DEPTH>=100) & 
                    (~dbsm_data.EXPERIMENT.isin(["bulk", "subsample"]))]

variants = sorted_dbsm_data.CHR_POS_REF_ALT.unique()
sorted_samples = sorted_dbsm_data.ID.unique()
bulk_samples = bulk_dbsm_data.ID.unique()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#summarize mafs for sorted population data and bulk data
sorted_maf_df = condense_maf_data(sorted_dbsm_data)
bulk_maf_df = condense_maf_data(bulk_dbsm_data)

#compute pairwise correlation coef for sorted population data and bulk data
sorted_coef_df, sorted_pval_df = compute_pairwise_corr(sorted_maf_df)
bulk_coef_df, bulk_pval_df = compute_pairwise_corr(bulk_maf_df)

#remove NA rows and columns
sorted_coef_df = sorted_coef_df.dropna()
sorted_variants = sorted_coef_df.index
sorted_coef_df = sorted_coef_df[sorted_variants]

bulk_coef_df = bulk_coef_df.dropna()
bulk_variants = bulk_coef_df.index
bulk_coef_df = bulk_coef_df[bulk_variants]

#get variants shared by both dataset
intersect_variants = set(sorted_variants).intersection(set(bulk_variants))
sorted_coef_df = sorted_coef_df[intersect_variants]
sorted_coef_df = sorted_coef_df.loc[intersect_variants, ]
bulk_coef_df = bulk_coef_df[intersect_variants]
bulk_coef_df = bulk_coef_df.loc[intersect_variants, ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#plot cluster map

# bulk data 
bulk_variants = bulk_coef_df.columns
pf_colors = bulk_variants.map(pf_dict).map(pf_color_dict)
f_colors = bulk_variants.map(f_dict).map(f_color_dict)
p_colors = bulk_variants.map(p_dict).map(p_color_dict)
o_colors = bulk_variants.map(o_dict).map(o_color_dict)
t_colors = bulk_variants.map(t_dict).map(t_color_dict)
cerebellum_colors = bulk_variants.map(cerebellum_dict).map(cerebellum_color_dict)
heart_colors = bulk_variants.map(heart_dict).map(heart_color_dict)
liver_colors = bulk_variants.map(liver_dict).map(liver_color_dict)
kidney_colors = bulk_variants.map(kidney_dict).map(kidney_color_dict)
laterization_colors = bulk_variants.map(laterization_dict).map(laterization_color_dict)

clustermap_bulk = sns.clustermap(bulk_coef_df, figsize=(50,50), cmap="RdBu_r",
                row_colors=[cerebellum_colors, heart_colors, liver_colors, kidney_colors,
                            pf_colors, f_colors, p_colors, o_colors, t_colors, laterization_colors],
              colors_ratio=0.005, dendrogram_ratio=0.1)
plt.savefig('bulk_clustermap.svg')   


#sorted population clustermap
sorted_variants = sorted_coef_df.columns
pf_colors = sorted_variants.map(pf_dict).map(pf_color_dict)
f_colors = sorted_variants.map(f_dict).map(f_color_dict)
p_colors = sorted_variants.map(p_dict).map(p_color_dict)
o_colors = sorted_variants.map(o_dict).map(o_color_dict)
t_colors = sorted_variants.map(t_dict).map(t_color_dict)
cerebellum_colors = sorted_variants.map(cerebellum_dict).map(cerebellum_color_dict)
heart_colors = sorted_variants.map(heart_dict).map(heart_color_dict)
liver_colors = sorted_variants.map(liver_dict).map(liver_color_dict)
kidney_colors = sorted_variants.map(kidney_dict).map(kidney_color_dict)
laterization_colors = sorted_variants.map(laterization_dict).map(laterization_color_dict)

clustermap_sorted = sns.clustermap(sorted_coef_df, figsize=(50,50), cmap="RdBu_r",
                row_colors=[cerebellum_colors, heart_colors, liver_colors, kidney_colors,
                            pf_colors, f_colors, p_colors, o_colors, t_colors, laterization_colors],
              colors_ratio=0.005, dendrogram_ratio=0.1)
plt.savefig('sorted_population_clustermap.svg')     


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#place bulk and sorted together 

secondary_hybrid_sorted = bulk_coef_df.loc[list(clustermap_sorted.data2d.columns),]
secondary_hybrid_bulk = sorted_coef_df.loc[list(clustermap_bulk.data2d.columns),]

#upper triange is sorted population
hybrid_sorted = np.zeros([249, 249])
for i in range(249):
    for j in range(249):
        if i <= j:
            hybrid_sorted[i,j]= clustermap_bulk.data2d.values[i,j]
        else:
            hybrid_sorted[i,j]= secondary_hybrid_sorted.values[i,j]
            
hybrid_sorted = pd.DataFrame(hybrid_sorted, columns = list(clustermap_sorted.data2d.columns))
hybrid_sorted = hybrid_sorted.set_index(clustermap_sorted.data2d.columns)

fig, ax = plt.subplots(figsize=(50, 50))
sns.heatmap(hybrid_sorted, cmap="RdBu_r", cbar=False)
plt.savefig('hybrid_sorted_heatmap.svg')

#upper triange is bulk
hybrid_bulk = np.zeros([249, 249])
for i in range(249):
    for j in range(249):
        if i <= j:
            hybrid_bulk[i,j]= clustermap_bulk.data2d.values[i,j]
        else:
            hybrid_bulk[i,j]= secondary_hybrid_bulk.values[i,j]
            
hybrid_bulk = pd.DataFrame(hybrid_bulk, columns = list(clustermap_bulk.data2d.columns))
hybrid_bulk = hybrid_bulk.set_index(clustermap_bulk.data2d.columns)

fig, ax = plt.subplots(figsize=(50, 50))
sns.heatmap(hybrid_bulk, cmap="RdBu_r", cbar=False)
plt.savefig('hybrid_bulk_heatmap.svg')

