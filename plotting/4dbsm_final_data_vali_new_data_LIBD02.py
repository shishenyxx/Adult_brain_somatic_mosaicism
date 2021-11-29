# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 12:49:43 2021

@author: Martin
"""

#Import modules

import sys, os, gzip
import pandas as pd
import numpy as np
import math
import scipy as sp
import seaborn as sns
import matplotlib.pyplot as plt
import pingouin as pg

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle, Polygon, Patch
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D as Line

plt.rcParams['svg.fonttype'] = 'none'
plt.ioff()

#Import modules done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Functions

#note that the raw input is coming from the 7614 split; see below
#------------------------------------------------------------------------------
#change column name for CHROM-POS...; harmonize ID column
    
def change_cpra(df):
    
    df.rename({'CHROM-POS-REF-ALT': 'CHR_POS_REF_ALT'}, axis=1, inplace=True)

def change_ID(df):
    
    df['ID'] = df.apply(lambda row: '_'.join(row.ID.split('-')), axis=1)

#------------------------------------------------------------------------------
#separate 7614 and LIBD

def separate_7614_LIBD(df):
    '''split the two experiments with the appropriate controls. Note that after
    this, only code for 7614 is in this .py file. LIBD data will be processed
    separately.'''
    
    _7614 = df[(df.ID.str.contains('7614')) | (df.ID == 'JGG_Cells')]
    LIBD = df[(df.ID.str.contains('830')) | (df.ID == 'JGG_LIBD')]
    
    return _7614, LIBD

#------------------------------------------------------------------------------
#annotate categories

def annotate_experiment_tissues(df):
    '''use the ID tags to define different features, such as experiment type,
    orientation etc. This has been changed to accomodate the new labels.'''
    
    df['LIST'] = df.apply(extract_organ_lat_experiment, axis=1)
    
    df[['INDIVIDUAL', 'ORGAN', 'BRAIN_REGION', 'L_R', 'EXPERIMENT']] = \
        pd.DataFrame(df.LIST.to_list(), index=df.index)
        
    df.drop('LIST', axis=1, inplace=True)

def extract_organ_lat_experiment(row):
    
    items = row.ID.split('_')
    
    if items[0] == 'JGG':
        return ['JGG', 'JGG', 'no_region', 'no_lat', 'ctrl']
    
    elif len(items) == 2:
        return [items[0], items[1], 'no_region', 'no_lat', 'bulk']
    
    elif len(items) == 5:
        return [items[0], items[1], items[2], items[3], 'bulk']
    
    elif len(items) == 6:
        return [items[0], items[1], items[2], items[3], items[5]]
    
    else:
        return ['zonk', 'zonk', 'zonk', 'zonk', 'zonk']

#------------------------------------------------------------------------------

def add_indel_info(df):
    '''add info whether a variant is an indel or a long indel (more than 1bp
    difference).'''
    
    df['INDEL'] = (df.REF.str.len() != df.ALT.str.len())
    df['LONG_INDEL'] = ((df.REF.str.len() > 2) | (df.ALT.str.len() > 2))


def add_depth(df):
    '''add depth'''
    
    df['DEPTH'] = df.ALT_COUNT + df.REF_COUNT
    df['DEPTH'] = df.DEPTH.fillna(0)


def add_depth30_flag(df):
    '''use depth <30 as a flag for coverage issues. might need additional
    flags for other coverage types.'''
    
    value_bulk = {'bulk': 1}
    
    df['FLAG_30DEPTH'] = df.apply(lambda row: row.DEPTH < 30, axis=1)
    df['FLAG_30DEPTH_5'] = df.apply(lambda row: row.FLAG_30DEPTH *
                                             value_bulk.get(row.EXPERIMENT, 0),
                                     axis=1)

def remove_dups(df):
    '''remove duplicated lines if present.'''
    
    df.drop_duplicates(keep=False, inplace=True)

#------------------------------------------------------------------------------
#annotate hets and mosaics

def annotate_noise_het(df):
    '''annotate noise with determined threshold for lower to not be considered
    noise, in addition to being higher than the control upper. het is assigned
    for tissues if upper is higher than 0.40. These
    values are based on empirical analysis from het variants and hom variants
    by XY. also uses a 30x depth filter for both noise and het.'''
    
    df['FLAG_NOISE'] = df.apply(flag_noise, axis=1)
    df['FLAG_HET'] = df.apply(flag_het, axis=1)


def flag_noise(row):
    
    '''based on the 95% CI:
    data[(data.CATEGORY == 'MOSAIC') & (data.FLAG_30DEPTH == False) &
         (data.EXPERIMENT == 'bulk') &
         (data.INDIVIDUAL_DETECTED.astype('str') != data.INDIVIDUAL) &
         (data.UPPER_CI <= 0.4)].LOWER_CI
     
    Out[26]: 
        count    1.202100e+04
        mean     1.218374e-03
        std      9.254369e-03
        min     -3.470000e-18
        1%      -2.440000e-20
        50%      7.930000e-05
        90%      9.234610e-04
        95%      3.524257e-03
        99%      2.572567e-02
        max      3.900199e-01'''
    
    above_thresh = row.LOWER_CI > 3.524257e-03
    above_ctrl = row.LOWER_CI > row.NORMAL_UPPER_CI
    depth = (row.DEPTH >= 30)
    alt = (row.ALT_COUNT >= 3)
    not_JGG = (row.ORGAN != 'JGG')
    
    return not bool(above_thresh and above_ctrl and depth and alt and not_JGG)

def flag_het(row):
    
    return bool((row.UPPER_CI >= 0.4) & (row.DEPTH >= 30))


#------------------------------------------------------------------------------

def variant_is_present(df):
    '''essentially only the reverse of the noise flag.'''
    
    value_bulk = {'bulk': 1}
    
    df['IN_SAMPLE'] = df.apply(lambda row: not row['FLAG_NOISE'], axis=1)
    df['IN_5'] = df.apply(lambda row: row.IN_SAMPLE *
                                      value_bulk.get(row.EXPERIMENT, 0) *
                              (row.INDIVIDUAL == str(row.INDIVIDUAL_DETECTED)),
                           axis=1)
    df['IN_WRONG_5'] = df.apply(lambda row: row.IN_SAMPLE *
                                            value_bulk.get(row.EXPERIMENT, 0) *
                              (row.INDIVIDUAL != str(row.INDIVIDUAL_DETECTED)),
                                axis=1)
    
    df['FLAG_HET_5'] = df.apply(lambda row: row.FLAG_HET *
                                            value_bulk.get(row.EXPERIMENT, 0) *
                              (row.INDIVIDUAL == str(row.INDIVIDUAL_DETECTED)),
                                 axis=1)


#------------------------------------------------------------------------------
#make sorter

def make_sorter(df):
    '''make a sorter column.'''
    
    organ = {'Ctx': 'A', 'Cbl': 'B', 'JGG': 'C'}
    region = {'PF': 'A', 'T': 'B', 'JGG': 'C'}
    lr = {'L': 'A', 'R': 'B', 'no_lat': 'C'}
    exp = {'bulk': 'A', 'NeuN': 'B', 'Olig2': 'C', 'ctrl': 'D'}
    
    df['SORTER_STR'] = df.apply(lambda row: row.INDIVIDUAL+
                                            organ[row.ORGAN] +
                                            lr[row.L_R] +
                                            region.get(row.BRAIN_REGION, 'X') +
                                            exp.get(row.EXPERIMENT, 'X'),
                                axis=1)
    
    sorters = df.SORTER_STR.sort_values().unique()
    
    sorter_dic = {sorter: i for sorter, i
                            in zip(sorters, range(1, (len(sorters) + 1)))}
    
    df['SORTER'] = df.apply(lambda row: sorter_dic[row.SORTER_STR], axis=1)
    
    df['CHROM_SORTER'] = df.apply(get_chrom_str, axis=1)
    
    df.sort_values(['CHROM_SORTER', 'CHR_POS_REF_ALT', 'SORTER'], inplace=True)
    df.reset_index(drop=True, inplace=True)

def get_chrom_str(row):
    '''use normal string manipulation to get chrN with N bein 01-22 or X.'''
    
    chrom = str(row['CHROM'])
    
    if len(chrom) != 2 and chrom != 'X':
        chrom = '0' + chrom
    
    return ('chr' + chrom)
#------------------------------------------------------------------------------

def add_seed_info(df):
    '''add a recurrence flag, so analysis on the variant level that do not
    need organ information can be done.'''
    
    df['SEED'] = ~df.duplicated('CHR_POS_REF_ALT', keep='first')


#------------------------------------------------------------------------------
#annotate signal, mosaic_support, real_mosaic

def annotate_mos(df):
    '''use the flags above to determine whether a variant is a real mosaic or
    not.'''
    
    df.sort_values(['CHROM_SORTER', 'CHR_POS_REF_ALT', 'SORTER'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    grp = df.groupby('CHR_POS_REF_ALT')
    
    df['FLAG_PRESENT_ALL'] = grp['IN_SAMPLE'].transform(flag_signal)
    df['FLAG_PRESENT_5'] = grp['IN_5'].transform(flag_signal)
    df['FLAG_PRESENT_WRONG_5'] = grp['IN_WRONG_5'].transform(flag_signal)
    
    df['SUM_HET_5'] = grp['FLAG_HET_5'].transform('sum')
    
    df['SUM_DEPTH_FLAG'] = grp['FLAG_30DEPTH'].transform('sum')
    df['SUM_DEPTH_FLAG_5'] = grp['FLAG_30DEPTH_5'].transform('sum')
    
    df['SET_MOSAIC'] = (
                        (df.CATEGORY == 'MOSAIC') &
                        (df.FLAG_PRESENT_5 == True) &
                        (df.FLAG_PRESENT_WRONG_5 == False) &
                        (df.SUM_HET_5 < 3) &
                        (df.NORMAL_LOWER_CI < 3.524257e-03) &
                        (df.SUM_DEPTH_FLAG_5 < 3)
                        )
    
def flag_signal(col):
    
    return sum(col) > 0

#------------------------------------------------------------------------------

def make_mos(df):
    '''clean up df and only take mosaics.'''
    
    df = df[df.SET_MOSAIC == True]
    
    df.drop(['GNOMAD_FREQ', 'REPEAT_MASKER', 'SEGDUP', 'HOMOPOLYMER',
             'DINUCLEOTIDE', 'NEAR_INDEL', 'UCSC_RPMSK', 'SORTER_STR'],
            axis=1, inplace=True)
    
    return df
    

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#all functions to annotate following mosaic variant determination (size!)

#------------------------------------------------------------------------------
#mutsigs for plotting

def add_mut_cats(df):
    '''extract the mutational category from the sequence snippet.'''
    
    df['TRI_REF'] = df.REF_SEQ.str[7:10]
    
    df['TRI_CAT'] = df.apply(lambda row: categorize_mut(row)[0], axis=1)
    df['CAT'] = df.apply(lambda row: categorize_mut(row)[1], axis=1)


def categorize_mut(row):
    '''use standard approach to make mutationl signature categories.'''
    
    if row['INDEL'] == True:
        
        return ('INDEL', 'INDEL')
    
    elif row['REF'] in set('CT'):
        
        tri_cat = row['TRI_REF']
        cat = row['REF'] + '>' + row['ALT']
        
        return (tri_cat, cat)
    
    else:
        
        tri_cat = rev_comp(row['TRI_REF'])
        cat = rev_comp(row['REF']) + '>' + rev_comp(row['ALT'])
        
        return (tri_cat, cat)

    
def rev_comp(seq):
    '''make reverse complement of sequence. simple, assumes that all are upper
    case and regular nucleotides. left without failsaves to make sure that the
    data passed on is as expected.'''
    reverse_comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    rev_seq = ''
    
    for n in seq[::-1]:
        rev_seq += reverse_comp[n]
    
    return rev_seq

#------------------------------------------------------------------------------
#loose functions
    
def sorter(df):
    '''convenience function to sort by variant and reset index.'''
    
    df.sort_values(['CHROM_SORTER', 'CHR_POS_REF_ALT', 'SORTER'], inplace=True)
    df.reset_index(drop=True, inplace=True)

def count_positive_tissues(df):
    '''used to assess number of positive calls in 25/tissues/cells'''
    
    df['SUM_POS_ALL'] = df.groupby('CHR_POS_REF_ALT')['IN_SAMPLE']\
                                        .transform('sum')
    
    df['SUM_POS_5'] = df.groupby('CHR_POS_REF_ALT')['IN_5']\
                                        .transform('sum')

#------------------------------------------------------------------------------
#make a tissue sorter
def tissue_sorter(df):
    '''sort for the 25 regions.'''
    
    df['INTERMEDIATE'] = df.apply(lambda row: row.ORGAN + '_' +
                                              row.BRAIN_REGION,
                                  axis=1)
    dictionary = {'Ctx_PF': 1, 'Ctx_T': 2, 'Cbl_no_region': 3,
                  'JGG_no_region':4}
    
    df['ORGAN_SORTER'] = df.apply(lambda row:
                                          dictionary.get(row.INTERMEDIATE, 10),
                                  axis=1)
    
    df.drop('INTERMEDIATE', axis=1, inplace=True)

#------------------------------------------------------------------------------
#get max is complicated by cells and by tissues that have very low coverage

def get_max_af(df):
    '''get the maximum AF for each variant and the lower/upper coordinates.'''
    
    df['CONSIDERED_MAF'] = df.apply(lambda row: row.MAF * (row.IN_SAMPLE) *
                              (row.INDIVIDUAL == str(row.INDIVIDUAL_DETECTED)),
                                    axis=1)
    
    df['CONSIDERED_MAF_5'] = df.apply(lambda row: row.MAF * (row.IN_5),
                                      axis=1)
    
    grp = df.groupby('CHR_POS_REF_ALT') #generates groupby object
    
    df['MAX_MAF'] = grp['CONSIDERED_MAF'].transform('max')
    df['MAX_MAF_BOOL'] = df.apply(lambda row:
                                       row['CONSIDERED_MAF'] == row['MAX_MAF'],
                                  axis=1)
    
    df['MAX_MAF_5'] = grp['CONSIDERED_MAF_5'].transform('max')
    df['MAX_MAF_5_BOOL_int'] = df.apply(lambda row:
                                   row['CONSIDERED_MAF_5'] == row['MAX_MAF_5'],
                                        axis=1)
    df['MAX_MAF_5_BOOL'] = (df.MAX_MAF_5_BOOL_int) & (df.FLAG_PRESENT_5)
    
    df.drop('MAX_MAF_5_BOOL_int', axis=1, inplace=True)

#+-----------------------------------------------------------------------------
#functions to annotate categories such as lateralized, brain specific etc.
def annotate_expression_pattern(df):
    '''use transform to annotate patterns as bool. note that some of the groups
    require a clean mosaic table, i.e. variants must be determined as positive
    in at least one tissue.'''
    
    sorter(df)
    
    df['SET_ONE_TISSUE'] = one_tissue(df)
    df['SET_CORTEX_ONLY'] = cortex_only(df)
    df['SET_LEFT_ONLY'] = left_only(df)
    df['SET_RIGHT_ONLY'] = right_only(df)
    df['SET_IN_CORTEX'] = in_cortex(df)
    df['SET_IN_CEREBELLUM'] = in_cerebellum(df)
    #laterlized in brain, but can be anywhere in other tissues
    df['SET_LEFT_ONLY_CTX'] = left_only_ctx(df)
    df['SET_RIGHT_ONLY_CTX'] = right_only_ctx(df)

def one_tissue(df):
    '''could be done easier with df.apply(lambda row: row.SUM_POS_25 ==1,
    axis=1); but used as template.'''
    
    df_bulk = df[df.EXPERIMENT == 'bulk']
    
    df_grouped = df_bulk.groupby(['CHR_POS_REF_ALT'])['IN_5'].sum()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_5 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    
    return df_out.OUTPUT

def cortex_only(df):
    '''all only calls rely on a variant being mosaic, but absent in all other 
    tissues. e.g. for cortex it has to be not present in any non-ctx tissues.
    logic is used for all tyypes of _only.'''
    
    df_nonctx = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN != 'Ctx') &
                   (df.SUM_POS_5 > 0)]
    
    df_grouped = df_nonctx.groupby(['CHR_POS_REF_ALT'])['IN_5'].max()\
                                                               .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_5 == 0, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def left_only(df):
    
    df_r = df[(df.EXPERIMENT == 'bulk') & (df.L_R != 'L') &
              (df.SUM_POS_5 > 0)]
    
    df_grouped = df_r.groupby(['CHR_POS_REF_ALT'])['IN_5'].max()\
                                                          .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_5 == 0, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def right_only(df):
    
    df_l = df[(df.EXPERIMENT == 'bulk') & (df.L_R != 'R') &
              (df.SUM_POS_5 > 0)]
    
    df_grouped = df_l.groupby(['CHR_POS_REF_ALT'])['IN_5'].max()\
                                                          .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_5 == 0, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_cortex(df):
    
    df_ctx = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN == 'Ctx')]
    
    df_grouped = df_ctx.groupby(['CHR_POS_REF_ALT'])['IN_5'].max()\
                                                            .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_5 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_cerebellum(df):
    
    df_cbl = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN =='Cbl')]
    
    df_grouped = df_cbl.groupby(['CHR_POS_REF_ALT'])['IN_5'].max()\
                                                            .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_5 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def left_only_ctx(df):
    
    df_cr = df[(df.EXPERIMENT == 'bulk') & (df.SET_IN_CORTEX == True) &
                ((df.L_R == 'R') & (df.ORGAN == 'Ctx'))]
    
    df_grouped = df_cr.groupby(['CHR_POS_REF_ALT'])['IN_5'].max()\
                                                           .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_5 == 0, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def right_only_ctx(df):
    
    df_cl = df[(df.EXPERIMENT == 'bulk') & (df.SET_IN_CORTEX == True) &
                ((df.L_R == 'L') & (df.ORGAN == 'Ctx'))]
    
    df_grouped = df_cl.groupby(['CHR_POS_REF_ALT'])['IN_5'].max()\
                                                           .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_5 == 0, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#Functions to define a category label
def define_cat_labels(df):
    '''make several labels that are mutually exclusive. Organ-specific (HLK),
    Organ+brain, Across-brain, lateralized brain, single tissue only.'''
    
    df['CAT_LABEL'] = df.apply(categorize_label, axis=1)

def categorize_label(row):
    '''use bools to define label.'''
    
    if row.SUM_POS_5 == 0:
        
        return 'not_in_5'
    
    if row.SET_IN_CEREBELLUM == True:
        if row.SET_ONE_TISSUE == True:
            return 'cerebellum only'
        else:
            return 'across brain'
        
    else: #i.e. variant is only in cortex
            
        if row.SET_LEFT_ONLY == False and row.SET_RIGHT_ONLY == False:
            return 'across hemispheres'
        elif row.SET_ONE_TISSUE == True:
            return 'single brain tissue'
        elif row.SET_LEFT_ONLY == True:
            return 'left only'
        elif row.SET_RIGHT_ONLY == True:
            return 'right only'
        else:
            return 'ZONK'
    
    return 'ZONK2'
            
#------------------------------------------------------------------------------

###############################################################################
###############################################################################
###############################################################################

#Plotting~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Base Plotting
#------------------------------------------------------------------------------

def plot_bars_for_categories(df):
    '''plot several categories as boxes side by side with the same
    visualization.'''
    
    df = df[df.MAX_MAF_5_BOOL == True]
    
    fig, axs = plt.subplots(nrows=4, ncols=4)
    
    output_data = []
    
    for i, ID in enumerate(df.INDIVIDUAL.sort_values().unique()):
        
        df_ = df[df.INDIVIDUAL == ID]
        
        all_mos = df_
        cortex_only = df_[df_.SET_CORTEX_ONLY]
        brain_side_restricted = df_[((df_.SET_LEFT_ONLY) |
                                     (df_.SET_RIGHT_ONLY)) &
                                   (df_.SET_CORTEX_ONLY)]
        brain_one_tissue_only = df_[(df_.SET_ONE_TISSUE) &
                                    (df_.SET_CORTEX_ONLY)]
        
        plots = [all_mos, cortex_only, brain_side_restricted,
                 brain_one_tissue_only]
        
        colors = ['xkcd:dark blue', 'xkcd:blue', 'xkcd:sky blue',
                  'xkcd:pale blue']
        
        for plot, color, j in zip(plots, colors, range(4)):
            
            sns.countplot(x='MAX_MAF_5_BOOL', data=plot, color=color,
                          edgecolor='k', ax=axs[i][j])
            
            axs[i][j].set_ylim(-5,300)
            axs[i][j].set_yticks([0, 150, 300])
            axs[i][j].set_xlabel('')
            axs[i][j].set_xticks([])
            
            if j != 0:
                sns.despine(left=True, bottom=True, offset=5, trim=True,
                            ax=axs[i][j])
                axs[i][j].set_yticks([])
                axs[i][j].set_ylabel('')
            else:
                sns.despine(bottom=True, offset=5, trim=True, ax=axs[i][j])
                axs[i][j].set_ylabel('Number of Variants')
        
        output_data.append([len(pl) for pl in plots])
    
    all_mos = df
    cortex_only = df[df.SET_CORTEX_ONLY]
    brain_side_restricted = df[((df.SET_LEFT_ONLY) | (df.SET_RIGHT_ONLY)) &
                                (df.SET_CORTEX_ONLY)]
    brain_one_tissue_only = df[(df.SET_ONE_TISSUE) & (df.SET_CORTEX_ONLY)]
    
    plots = [all_mos, cortex_only, brain_side_restricted,
             brain_one_tissue_only]
    
    colors = ['xkcd:dark blue', 'xkcd:blue', 'xkcd:sky blue',
              'xkcd:pale blue']
    
    names = ['All Variants', 'Cortex Only', 'Brain Lateralized',
             'Brain Single Sample']
    
    for plot, color, name, i in zip(plots, colors, names, range(4)):
        
        sns.countplot(x='MAX_MAF_5_BOOL', data=plot, color=color,
                      edgecolor='k', ax=axs[3][i])
        
        axs[3][i].set_ylim(-5,500)
        axs[3][j].set_yticks([0, 250, 500])
        axs[3][i].set_xlabel(name, rotation=45, ha='right')
        axs[3][i].set_xticks([])
        
        if i != 0:
            sns.despine(left=True, bottom=True, offset=5, trim=True,
                        ax=axs[3][i])
            axs[3][i].set_yticks([])
            axs[3][i].set_ylabel('')
        else:
            sns.despine(bottom=True, offset=5, trim=True, ax=axs[3][i])
            axs[3][i].set_ylabel('Number of Variants')
    
    output_data.append([len(pl) for pl in plots])
    
    plt.show()
    
    return output_data


def plot_boxes_for_categories_split(df):
    '''plot several categories as boxes side by side with the same
    visualization..'''
    
    df = df[df.MAX_MAF_5_BOOL == True]
    
    df['MAF'] = df.MAF**0.5
    
    fig, axs = plt.subplots(nrows=2, ncols=4,
                            gridspec_kw={'height_ratios':[1,5]})

    all_mos = df
    cortex_only = df[df.SET_CORTEX_ONLY]
    brain_side_restricted = df[((df.SET_LEFT_ONLY) | (df.SET_RIGHT_ONLY)) &
                                (df.SET_CORTEX_ONLY)]
    brain_one_tissue_only = df[(df.SET_ONE_TISSUE) & (df.SET_CORTEX_ONLY)]
    
    plots = [all_mos, cortex_only, brain_side_restricted,
             brain_one_tissue_only]
    
    colors = ['xkcd:dark blue', 'xkcd:blue', 'xkcd:sky blue',
              'xkcd:pale blue']
    
    names = ['All Variants', 'Cortex Only', 'Brain Lateralized',
             'Brain Single Sample']
    
    for plot, color, name, i in zip(plots, colors, names, range(4)):
        
        sns.boxplot(x='MAX_MAF_5_BOOL', y='MAF', data=plot, color=color,
                    ax=axs[0,i])
        sns.boxplot(x='MAX_MAF_5_BOOL', y='MAF', data=plot, color=color,
                    ax=axs[1,i])
        
        axs[0,i].set_ylim(0.2**0.5, 1)
        axs[0,i].set_yticks([0.2**0.5, 1])
        axs[0,i].set_yticklabels(['0.2', '1.00'])
        
        axs[1,i].set_ylim(0, 0.2)
        axs[1,i].set_yticks([0, 0.05**0.5, 0.1**0.5, 0.2**0.5])
        axs[1,i].set_yticklabels(['0.00', '0.05', '0.10', '0.20'])
        
        axs[0,i].set_xlabel('')
        axs[0,i].set_xticks([])
        axs[1,i].set_xlabel(name, rotation=45, ha='right')
        axs[1,i].set_xticks([])
        
        if i != 0:
            sns.despine(left=True, bottom=True, offset=5, trim=True,
                        ax=axs[0,i])
            axs[0,i].set_yticks([])
            axs[0,i].set_ylabel('')
            sns.despine(left=True, bottom=True, offset=5, trim=True,
                        ax=axs[1,i])
            axs[1,i].set_yticks([])
            axs[1,i].set_ylabel('')
        else:
            sns.despine(bottom=True, offset=5, trim=True, ax=axs[0,i])
            sns.despine(bottom=True, offset=5, trim=True, ax=axs[1,i])
            axs[1,i].set_ylabel('Max. Allelic Fraction\n(sqrt-t)')
    
    plt.show()


def plot_sml_1tissues(df):
    '''Plot bars counting the number of variants that are only present in one
    small region.'''
    
    df = df[(df.MAX_MAF_5_BOOL) & (df.SET_CORTEX_ONLY) &
            (df.SET_ONE_TISSUE)].sort_values(by=['ORGAN_SORTER', 'L_R'])
    
    fig, axs = plt.subplots(nrows=4, ncols=1)
    
    for i, ID in enumerate(df.INDIVIDUAL.sort_values().unique()):
        
        df_ = df[df.INDIVIDUAL == ID]
        
        if ID == '8307':
            #to fix the issue with 0 bars
            fake = pd.DataFrame(data={'BRAIN_REGION' : ['T'], 'L_R': ['R']})
            df_ = df_.append(fake, sort=False).reset_index(drop=True)
            
        sns.countplot(x='BRAIN_REGION', hue='L_R', data=df_, order=['PF', 'T'],
                      edgecolor='k', ax=axs[i])
        
        axs[i].set_xticklabels([])
        axs[i].set_xlabel('')
        axs[i].set_ylim(0, 120)
        axs[i].set_yticks([0, 40, 80, 120])
        axs[i].set_ylabel('')
        
    sns.countplot(x='BRAIN_REGION', hue='L_R', data=df, order=['PF', 'T'],
                  edgecolor='k', ax=axs[3])
    
    axs[3].set_ylim(0, 120)
    axs[3].set_yticks([0, 40, 80, 120])
    axs[3].set_ylabel('Number of Variants')
    axs[3].set_xlabel('Cortical Lobes')
    
    axs[1].get_legend().remove()
    axs[2].get_legend().remove()
    axs[3].get_legend().remove()
    
    plt.title('Remove R-T 8307!')
    
    sns.despine(bottom=True, offset=5, trim=True)
    
    plt.show()


def plot_number5_vsAF(df):
    '''note that OR lines are drawn at 5.5 and 0.05 manually)'''
    
    df = df.copy()
    
    df.fillna(0, inplace=True)
    
    df = df[(df.EXPERIMENT == 'bulk') & (df.CONSIDERED_MAF_5 != 0)]
    
    df['MAF_'] = df.CONSIDERED_MAF_5**0.5
    
    sns.regplot(x='MAF_', y='SUM_POS_5', data=df, color='k',
                scatter_kws={'s':4, 'alpha':0.1}, fit_reg=False)
                #, line_kws={'color':'r'})
    #sns.jointplot('MAF_', 'SUM_POS_25', kind='kde', data=df,
     #             color='r')#.plot_joint(sns.regplot, marker='o', color='k')
    
    plt.xlabel('Allelic Fraction (sqrt-t)')
    plt.xlim(0, 1)
    plt.xticks(ticks=[0., 0.01**0.5, 0.05**0.5, 0.1**0.5, 0.2**0.5, 0.5**0.5],
               labels=['0.00', '0.01', '0.05', '0.10', '0.20', '0.50'])
    
    plt.ylabel('Number of Samples')
    plt.ylim(0, 5.1)
    plt.xlim(0, 1)
    plt.yticks(ticks=[0, 1, 2, 3, 4, 5])
    
    sns.despine(offset=5, trim=True)
    
    plt.show()
    #return(sp.stats.spearmanr(df.MAF_, df.SUM_POS_25))


#---
def lr_ap_plot(df, depth=1000):
    '''when doing individuals, for ID03, set
    sns.distributions._has_statsmodels = False.'''
    
    df = df[(df.EXPERIMENT == 'bulk') & (df.DEPTH >= depth) &
            (df.SET_IN_CORTEX) & (df.ORGAN != 'Cbl')]
    
    df['SUM_SORTS_POS'] = df.groupby('CHR_POS_REF_ALT')['IN_5']\
                            .transform('sum')
    df = df[df.SUM_SORTS_POS > 0].reset_index()
    
    df['MAX_CONSIDERED'] = df.groupby('CHR_POS_REF_ALT')\
                                          ['CONSIDERED_MAF_5'].transform('max')
    df['LEFT_AVG_CMAF'] = l_r_acmaf(df, l_r='L')
    df['RIGHT_AVG_CMAF'] = l_r_acmaf(df, l_r='R')
    df['NORM_DELTA_LR'] = df.apply(norm_delta, axis=1)
    
    df['ANT_AVG_CMAF'] = a_p_acmaf(df, a_p='A')
    df['POS_AVG_CMAF'] = a_p_acmaf(df, a_p='P')
    df['NORM_DELTA_AP'] = df.apply(norm_delta_ap, axis=1)
    
    plot_df = df[~(df.CHR_POS_REF_ALT.duplicated())][['INDIVIDUAL_DETECTED',
                                                      'NORM_DELTA_LR',
                                                      'NORM_DELTA_AP']]
    plot_df = plot_df.dropna()
    
    #only takes one color as argument
    #col_dic = {'8305' : 'b', '8306' : 'g', '8307': 'xkcd:magenta'}
    #colors = [col_dic[str(ID)] for ID in plot_df.INDIVIDUAL_DETECTED]
    
    #plt.scatter(plot_df.NORM_DELTA_LR, plot_df.NORM_DELTA_AP, alpha=0.7)
    sns.jointplot('NORM_DELTA_LR', 'NORM_DELTA_AP', kind='kde', data=plot_df,
                  color='k', joint_kws={'levels' : 30}).plot_joint(
                                        sns.regplot, marker='+', fit_reg=False,
                                        color='b')
    
    plt.xlabel('Normalized d (R-L)')
    plt.ylabel('Normalized d (A-P)')
    
    plt.show()
    
def l_r_acmaf(df, l_r):
    
    df_lat = df[df.L_R == l_r]
    df_grp = df_lat.groupby('CHR_POS_REF_ALT').mean().reset_index()
    
    df_tomerge = df_grp[['CHR_POS_REF_ALT', 'CONSIDERED_MAF_5']]
    df_tomerge.rename({'CONSIDERED_MAF_5': 'OUTPUT'}, axis=1,
                      inplace=True)
    
    df_out = pd.merge(df, df_tomerge, how='left')
    
    return df_out.OUTPUT

def norm_delta(row):
    
    l = row.LEFT_AVG_CMAF
    r = row.RIGHT_AVG_CMAF
    maxi = max([l, r])
    
    return (r-l)/maxi

def a_p_acmaf(df, a_p):
    
    if a_p == 'A':
        areas = ['PF']
    elif a_p == 'P':
        areas = ['T']
    
    df_ap = df[df.BRAIN_REGION.isin(areas)]
    df_grp = df_ap.groupby('CHR_POS_REF_ALT').mean().reset_index()
    
    df_tomerge = df_grp[['CHR_POS_REF_ALT', 'CONSIDERED_MAF_5']]
    df_tomerge.rename({'CONSIDERED_MAF_5': 'OUTPUT'}, axis=1,
                      inplace=True)
    
    df_out = pd.merge(df, df_tomerge, how='left')
    
    return df_out.OUTPUT
    
def norm_delta_ap(row):
    
    a = row.ANT_AVG_CMAF
    p = row.POS_AVG_CMAF
    maxi = max([a, p])
    
    return (a-p)/maxi
#---

#---
def lr_ap_plot_sortID02(df, depth=1000):
    
    df = df[(df.EXPERIMENT.isin(['NeuN', 'Olig2'])) & (df.DEPTH >= depth) &
            (df.INDIVIDUAL == '8305') & (df.INDIVIDUAL_DETECTED == 8305)]
    
    df['SUM_SORTS_POS'] = df.groupby('CHR_POS_REF_ALT')['IN_SAMPLE']\
                            .transform('sum')
    df = df[df.SUM_SORTS_POS > 0].reset_index()
    
    df['MAX_SORTED_CONSIDERED'] = df.groupby('CHR_POS_REF_ALT')\
                                            ['CONSIDERED_MAF'].transform('max')
    df['LEFT_AVG_CMAF'] = l_r_acmaf_sort(df, l_r='L')
    df['RIGHT_AVG_CMAF'] = l_r_acmaf_sort(df, l_r='R')
    df['NORM_DELTA_LR'] = df.apply(norm_delta, axis=1)
    
    df['ANT_AVG_CMAF'] = a_p_acmaf_sort(df, a_p='A')
    df['POS_AVG_CMAF'] = a_p_acmaf_sort(df, a_p='P')
    df['NORM_DELTA_AP'] = df.apply(norm_delta_ap, axis=1)
    
    plot_df = df[~(df.CHR_POS_REF_ALT.duplicated())][['NORM_DELTA_LR',
                                                     'NORM_DELTA_AP']]
    plot_df = plot_df.dropna()
    
    #plt.scatter(plot_df.NORM_DELTA_LR, plot_df.NORM_DELTA_AP, alpha=0.7)
    sns.jointplot('NORM_DELTA_LR', 'NORM_DELTA_AP', kind='kde', data=plot_df,
                  color='k').plot_joint(sns.regplot, marker='+', fit_reg=False,
                                        color='r')
    
    plt.xlabel('Normalized d (R-L)')
    plt.ylabel('Normalized d (A-P)')
    
    plt.show()
    
    
def l_r_acmaf_sort(df, l_r):
    
    df_lat = df[df.L_R == l_r]
    df_grp = df_lat.groupby('CHR_POS_REF_ALT').mean().reset_index()
    
    df_tomerge = df_grp[['CHR_POS_REF_ALT', 'CONSIDERED_MAF']]
    df_tomerge.rename({'CONSIDERED_MAF': 'OUTPUT'}, axis=1,
                      inplace=True)
    
    df_out = pd.merge(df, df_tomerge, how='left')
    
    return df_out.OUTPUT

def a_p_acmaf_sort(df, a_p):
    
    if a_p == 'A':
        areas = ['PF', 'F']
    elif a_p == 'P':
        areas = ['P', 'O', 'T']
    
    df_ap = df[df.BRAIN_REGION.isin(areas)]
    df_grp = df_ap.groupby('CHR_POS_REF_ALT').mean().reset_index()
    
    df_tomerge = df_grp[['CHR_POS_REF_ALT', 'CONSIDERED_MAF']]
    df_tomerge.rename({'CONSIDERED_MAF': 'OUTPUT'}, axis=1,
                      inplace=True)
    
    df_out = pd.merge(df, df_tomerge, how='left')
    
    return df_out.OUTPUT
#---

#---
def plot_sort_overview(df, POS):
    '''plot both sorted populations in an overview plot.'''
    
    markers = ['NeuN', 'Olig2']
    colors = ['xkcd:brick red', 'xkcd:orange']
    
    df = df[(df.EXPERIMENT.isin(markers)) & (df.POS == POS) &
            (df.INDIVIDUAL == '8305')]
    sorter(df)
    
    df['MAF_'] = (df.CONSIDERED_MAF)**0.5
    
    maxi, semi, maxi_l, semi_l = get_ticks(df)
    
    fig, axs = plt.subplots(ncols=4, sharey=False)
    
    for i, lr in enumerate(['L', 'R']):
        for j, marker, color in zip(range(2), markers, colors):
            k = i * 2 + j
            
            df_m = df[(df.L_R == lr) & (df.EXPERIMENT == marker)]
            
            colors_m = make_colors_markers(df_m, color)
            axs[k].vlines(range(len(df_m)), ymin=[0 for i in df_m.MAF_],
                          ymax=df_m.MAF_)
            axs[k].scatter(range(len(df_m)), df_m.MAF_, color=colors_m, s=100,
                           edgecolors='k', zorder=100)
            
            axs[k].set_xlim(-1, 5)
            axs[k].set_ylim(-((maxi/10)**0.5), maxi**0.5)
            
            if k != 0:
                sns.despine(left=True, bottom=True, ax=axs[k])
                axs[k].set_xticks([])
                axs[k].set_yticks([])
    
    line = Line([0.5125, 0.5125], [0.25, 0.88], transform=fig.transFigure,
                figure=fig, color='k', linestyle='--')
    fig.lines.extend([line])
    
    axs[0].set_yticks([0., semi**0.5, maxi**0.5])
    axs[0].set_yticklabels(['0.0', semi_l, maxi_l])
    sns.despine(bottom=True, trim=True, ax=axs[0])
    axs[0].set_xticks([])
    axs[0].set_ylabel('AF (sqrt-t)')
    axs[1].set_title(df.CHR_POS_REF_ALT.unique()[0])
                
    plt.show()

def make_colors_markers(df, color, depth=100):
    
    colors = []
    
    for index, row in df.iterrows():
        if row.DEPTH < depth:
            colors.append('w')
        elif row.DEPTH < 1000:
            colors.append('0.8')
        else:
            colors.append(color)
    
    return colors

def get_ticks(df):
    
    maxAF = max(df.MAF_)
    
    if maxAF < 0.05**0.5:
        return 0.05, 0.02, '0.05', '0.02'
    elif maxAF < 0.1**0.5:
        return 0.1, 0.02, '0.10', '0.02'
    elif maxAF < 0.2**0.5:
        return 0.2, 0.1, '0.2', '0.1'
    elif maxAF < 0.5**0.5:
        return 0.5, 0.2, '0.5', '0.2'
    else:
        return 1.0, 0.2, '1.0', '0.2'
#---

def make_abspres_plot(df, POS):
    '''just use the IN_5 to make a presence/absence plot for the entire phage
    or just the brain regions dependin on flag brain_only.'''
    
    #initialize df
    df = df[(df.EXPERIMENT == 'bulk') & (df.POS == POS) &
            (df.INDIVIDUAL == '8305')]
    sorter(df)
    
    #make color list
    colors = ['xkcd:pale pink' if in5 == 1 else '0.9' for in5 in df.IN_5]
    colors.insert(0, '0.5')
    for i in range(7):
        colors.insert(2, '0.5')
    colors.insert(10, '0.5')
    for i in range(7):
        colors.insert(12, '0.5')
    
    #make boxes
    fig,ax = plt.subplots(1)

    boxes = make_boxes_abspres()
        
    pc = PatchCollection(boxes, facecolor=colors, edgecolor='w')
    ax.add_collection(pc)
    
    ax.plot()
    ax.set(adjustable='box', aspect='equal')
    sns.despine(left=True, bottom=True)
    plt.tick_params(bottom=False, left=False,
                    labelbottom=False, labelleft=False)
    plt.xlim(0,20)
    plt.ylim(-1, 59)
    
    plt.title(df.CHR_POS_REF_ALT.unique()[0])

    plt.show()

def make_boxes_abspres():
    
    boxes = [Rectangle((0, 40), 10, 5), Rectangle((2, 41), 6, 3),
             Rectangle((0, 34.5), 10, 5), Rectangle((2, 35.5), 6, 3),
             Rectangle((0, 29), 10, 5), Rectangle((2 ,30), 6, 3),
             Rectangle((0, 23.5), 10, 5), Rectangle((2, 24.5), 6, 3),
             Rectangle((0, 18), 10, 5), Rectangle((2, 19), 6, 3),
             Rectangle((10.5, 40), 10, 5), Rectangle((12.5, 41), 6, 3),
             Rectangle((10.5, 34.5), 10, 5), Rectangle((12.5, 35.5), 6, 3),
             Rectangle((10.5, 29), 10, 5), Rectangle((12.5, 30), 6, 3),
             Rectangle((10.5, 23.5), 10, 5), Rectangle((12.5, 24.5), 6, 3),
             Rectangle((10.5, 18), 10, 5), Rectangle((12.5, 19), 6, 3),
             Rectangle((5,12.5), 10.5, 5)]
    
    return boxes
#---

def combined_sort_abspres(df_, POS):
    '''combine both in one plot.'''
    
    markers = ['NeuN', 'Olig2']
    colors = ['xkcd:brick red', 'xkcd:orange']
    
    df = df_[(df_.EXPERIMENT.isin(markers)) & (df_.POS == POS) &
            (df_.INDIVIDUAL == '8305')]
    sorter(df)
    
    df['MAF_'] = (df.CONSIDERED_MAF)**0.5
    
    maxi, semi, maxi_l, semi_l = get_ticks(df)
    
    fig, axs = plt.subplots(ncols=5, sharey=False)
    
    for i, lr in enumerate(['L', 'R']):
        for j, marker, color in zip(range(2), markers, colors):
            k = i * 2 + j
            
            df_m = df[(df.L_R == lr) & (df.EXPERIMENT == marker)]
            
            colors_m = make_colors_markers(df_m, color)
            axs[k].vlines(range(len(df_m)), ymin=[0 for i in df_m.MAF_],
                          ymax=df_m.MAF_)
            axs[k].scatter(range(len(df_m)), df_m.MAF_, color=colors_m, s=100,
                           edgecolors='k', zorder=100)
            
            axs[k].set_xlim(-1, 5)
            axs[k].set_ylim(-((maxi/10)**0.5), maxi**0.5)
            
            if k != 0:
                sns.despine(left=True, bottom=True, ax=axs[k])
                axs[k].set_xticks([])
                axs[k].set_yticks([])
    
    line = Line([0.4, 0.4], [0.25, 0.88], transform=fig.transFigure,
                figure=fig, color='k', linestyle='--')
    fig.lines.extend([line])
    
    axs[0].set_yticks([0., semi**0.5, maxi**0.5])
    axs[0].set_yticklabels(['0.0', semi_l, maxi_l])
    sns.despine(bottom=True, trim=True, ax=axs[0])
    axs[0].set_xticks([])
    axs[0].set_ylabel('AF (sqrt-t)')
    axs[1].set_title(df.CHR_POS_REF_ALT.unique()[0])
    
    df = df_[(df_.EXPERIMENT == 'bulk') & (df_.POS == POS) &
            (df_.INDIVIDUAL == '8305')]
    sorter(df)
    
    #make color list
    colors = ['xkcd:pale pink' if in5 == 1 else '0.9' for in5 in df.IN_5]
    colors.insert(0, '0.5')
    for i in range(7):
        colors.insert(2, '0.5')
    colors.insert(10, '0.5')
    for i in range(7):
        colors.insert(12, '0.5')
    
    boxes = make_boxes_abspres()
        
    pc = PatchCollection(boxes, facecolor=colors, edgecolor='w')
    axs[4].add_collection(pc)
    
    axs[4].plot()
    axs[4].set(adjustable='box', aspect='equal')
    sns.despine(left=True, bottom=True, ax=axs[4])
    axs[4].set_xticks([])
    axs[4].set_yticks([])
    axs[4].set_xlim(0,20)
    axs[4].set_ylim(-1, 59)
    
    plt.show()
#---

def correlate_markers(df, m1='NeuN', m2='Olig2', depth=1000):
    '''markers: NeuN, Olig2'''
    
    df = df[((df.EXPERIMENT == m1) | (df.EXPERIMENT == m2)) &
            (df.INDIVIDUAL=='8305') & (df.DEPTH >= depth) &
            (df.INDIVIDUAL_DETECTED == 8305)]

    sorter(df)
    
    df['SUM_PRESENCE'] = df.groupby(['POS', 'L_R', 'BRAIN_REGION'])\
                            ['IN_SAMPLE'].transform('sum')
    
    df = df[df.SUM_PRESENCE != 0]
    df['SUM_COUNT'] = df.groupby(['POS', 'L_R', 'BRAIN_REGION'])\
                                               ['IN_SAMPLE'].transform('count')
    df= df[df.SUM_COUNT == 2]
    
    
    x = df[df.EXPERIMENT == m1].reset_index().CONSIDERED_MAF
    y = df[df.EXPERIMENT == m2].reset_index().CONSIDERED_MAF
    
    plot_df = df[df.EXPERIMENT == m1].reset_index()
    plot_df['X'] = x**0.5
    plot_df['Y'] = y**0.5

    sns.regplot(x='X', y='Y', data=plot_df,
               scatter_kws={'edgecolor': 'k', 'color': '0.8', 'alpha': 0.5},
               line_kws={'color': 'k'})
    
    plt.xlim(-0.015, 1.)
    plt.ylim(-0.015, 1.)
    
    plt.xticks(ticks=[0., 0.2**0.5, 1.0], labels=['0.0', '0.2', '1.0'])
    plt.yticks(ticks=[0., 0.2**0.5, 1.0], labels=['0.0', '0.2', '1.0'])
    
    name = {'Olig2': 'OLIG2'}
    
    plt.xlabel('AF (sqrt-t) ' + name.get(m1, m1))
    plt.ylabel('AF (sqrt-t) ' + name.get(m2, m2))
    
    sns.despine(trim=True, offset=5)
    
    ax = plt.gca()
    ax.set(adjustable='box', aspect='equal')
    
    plt.show()
    
    return(sp.stats.spearmanr(x, y))





