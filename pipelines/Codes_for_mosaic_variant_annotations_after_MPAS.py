##Codes for varint annotation and mosaic variant filter after MPAS and snMPAS
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 04:33:04 2019

@author: Martin W. Breuss
"""

#Import modules

import sys, os, gzip
import pandas as pd
import numpy as np
import math
import scipy as sp
import pingouin as pg
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle, Polygon, Patch
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D as Line

plt.rcParams['svg.fonttype'] = 'none'
plt.ioff()

#Import modules done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Functions

#First level - Mosaic Output
#import------------------------------------------------------------------------

def import_tables(path_to_tissues, path_to_nuclei):
    
    tis = pd.read_csv(path_to_tissues, sep='\t')
    nuc = pd.read_csv(path_to_nuclei, sep='\t')
    
    df = tis.append(nuc, ignore_index=True)
    
    return df

#------------------------------------------------------------------------------
#change column name for CHROM-POS...
    
def change_cpra(df):
    
    df.rename({'CHROM-POS-REF-ALT': 'CHR_POS_REF_ALT'}, axis=1, inplace=True)

#annotation--------------------------------------------------------------------

#------------------------------------------------------------------------------
#annotate categories
def annotate_experiment_tissues(df):
    '''use the ID tags to define different features, such as experiment type,
    orientation etc. This assumes that ID descriptors stay the same moving
    forward.'''
    
    df['TISSUE_CELLS'] = df.apply(tissue_or_cell, axis=1)
    
    df['LIST'] = df.apply(extract_organ_lat_experiment, axis=1)
    
    df[['ORGAN', 'BRAIN_REGION', 'L_R', 'LG_SM', 'EXPERIMENT']] = \
        pd.DataFrame(df.LIST.to_list(), index=df.index)
        
    df.drop('LIST', axis=1, inplace=True)

def tissue_or_cell(row):
    
    if '-' in row.ID:
        return 'cells'
    else:
        return 'tissues'

def extract_organ_lat_experiment(row):
    
    if row.TISSUE_CELLS == 'cells':
        items = row.ID.split('-')
    elif row. TISSUE_CELLS == 'tissues':
        items = row.ID.split('_')
    else:
        return ['zonk1', 'zonk1', 'zonk1', 'zonk1', 'zonk1']
    
    if items[0] == 'JGG':
        return ['JGG', 'no_region', 'no_lat', 'no_lrg_sml', 'ctrl']
    
    elif len(items) == 2:
        return [items[1], 'no_region', 'no_lat', 'no_lrg_sml', 'bulk']
    
    elif len(items) == 3:
        return [items[1], 'no_region', items[2], 'no_lrg_sml', 'bulk']
    
    elif len(items) == 5:
        if items[-1] in ['Sml', 'Lrg']:
            return [items[1], items[2], items[3], items[4], 'bulk']
        elif 'Sml' in items[-1]:
            return [items[1], items[2], items[3], items[4], 'subsample']
        else:
            return [items[1], items[2], items[3], 'Lrg', items[4]]
    
    elif len(items) == 7:
        return [items[1], items[2], items[3], items[4], '_'.join(items[5:])]
    
    else:
        return ['zonk2', 'zonk2', 'zonk2', 'zonk2', 'zonk2']

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


#------------------------------------------------------------------------------
#annotate hets and mosaics

def annotate_noise_het(df):
    '''annotate noise with threshold of 0.005 for lower to not be considered
    noise, in addition to being higher than the control upper. het is assigned
    for cells if not noise, for tissues if upper is higher than 0.40. These
    values are based on empirical analysis from het variants and hom variants
    by XY. also uses a 30x depth filter for both noise and het.'''
    
    df['FLAG_NOISE'] = df.apply(flag_noise, axis=1)
    df['FLAG_NOISE01'] = df.apply(flag_noise01, axis=1)
    df['FLAG_NOISE95'] = df.apply(flag_noise95, axis=1)
    df['FLAG_NOISE95rn'] = df.apply(flag_noise95rn, axis=1)
    df['FLAG_HET'] = df.apply(flag_het, axis=1)


def flag_noise(row):
    
    #if additional depth filter is needed
    #value_tc = {'cells': 0, 'tissues': 1}
    #depth = (row.DEPTH >= (30 + value_tc[row.TISSUE_CELLS * 170]))
    
    above_thresh = row.LOWER_CI > 0.005
    above_ctrl = row.LOWER_CI > row.NORMAL_UPPER_CI
    depth = (row.DEPTH >= 30)
    not_JGG = (row.ORGAN != 'JGG')
    
    return not bool(above_thresh and above_ctrl and depth and not_JGG)

def flag_noise01(row):
    
    #if additional depth filter is needed
    #value_tc = {'cells': 0, 'tissues': 1}
    #depth = (row.DEPTH >= (30 + value_tc[row.TISSUE_CELLS * 170]))
    
    above_thresh = row.LOWER_CI > 0.001
    above_ctrl = row.LOWER_CI > row.NORMAL_UPPER_CI
    depth = (row.DEPTH >= 30)
    not_JGG = (row.ORGAN != 'JGG')
    
    return not bool(above_thresh and above_ctrl and depth and not_JGG)

def flag_noise95(row):
    
    '''based on the 95% CI:
    data[(data.CATEGORY == 'REF_HOMO') & (data.FLAG_30DEPTH == False) &
    (data.TISSUE_CELLS == 'tissues') & (data.ORGAN != 'JGG')]\
    .LOWER_CI.describe(percentiles=[0.01, 0.5, 0.90, 0.95, 0.99])
    Out[331]: 
        count    2.077000e+03
        mean     8.897367e-04
        std      9.281449e-03
        min     -3.470000e-18
        1%      -2.170000e-19
        50%      5.430000e-05
        90%      8.373440e-04
        95%      1.436920e-03
        99%      8.406548e-03
        max      3.201055e-01'''
    
    #if additional depth filter is needed
    #value_tc = {'cells': 0, 'tissues': 1}
    #depth = (row.DEPTH >= (30 + value_tc[row.TISSUE_CELLS * 170]))
    
    above_thresh = row.LOWER_CI > 0.00143692
    above_ctrl = row.LOWER_CI > row.NORMAL_UPPER_CI
    depth = (row.DEPTH >= 30)
    not_JGG = (row.ORGAN != 'JGG')
    
    return not bool(above_thresh and above_ctrl and depth and not_JGG)

def flag_noise95rn(row):
    
    '''based on the 95% CI:0.1
    data[(data.CATEGORY == 'REF_HOMO') & (data.FLAG_30DEPTH == False) &
    (data.TISSUE_CELLS == 'tissues') & (data.ORGAN != 'JGG')]\
    .LOWER_CI.describe(percentiles=[0.01, 0.5, 0.90, 0.95, 0.99])
     Out[26]: 
         count    2.089000e+03
         mean     8.811356e-04
         std      9.254298e-03
         min     -3.470000e-18
         1%      -2.170000e-19
         50%      5.410000e-05
         90%      8.266642e-04
         95%      1.397682e-03
         99%      8.367806e-03
         max      3.201055e-01'''
         
    #if additional depth filter is needed
    #value_tc = {'cells': 0, 'tissues': 1}
    #depth = (row.DEPTH >= (30 + value_tc[row.TISSUE_CELLS * 170]))
    
    above_thresh = row.LOWER_CI > 0.001397682
    above_ctrl = row.LOWER_CI > row.NORMAL_UPPER_CI
    depth = (row.DEPTH >= 30)
    alt = (row.ALT_COUNT >= 3)
    not_JGG = (row.ORGAN != 'JGG')
    
    return not bool(above_thresh and above_ctrl and depth and alt and not_JGG)

def flag_het(row):
    
    if row.TISSUE_CELLS == 'cells':
        return not (row.FLAG_NOISE)
    
    elif row.TISSUE_CELLS == 'tissues':
        return bool((row.UPPER_CI >= 0.4) & (row.DEPTH >= 30))

#------------------------------------------------------------------------------

def variant_is_present(df):
    '''essentially only the reverse of the noise flag.'''
    
    value_tc = {'cells': 0, 'tissues': 1}
    value_bulk = {'bulk': 1}
    
    df['IN_SAMPLE'] = df.apply(lambda row: not row['FLAG_NOISE'], axis=1)
    df['IN_TISSUE'] = df.apply(lambda row: row.IN_SAMPLE *
                                           value_tc[row.TISSUE_CELLS],
                               axis=1)
    df['IN_25'] = df.apply(lambda row: row.IN_SAMPLE *
                                       value_bulk.get(row.EXPERIMENT, 0),
                           axis=1)
    df['IN_CELLS'] = df.apply(lambda row: row.IN_SAMPLE *
                                          (1 - value_tc[row.TISSUE_CELLS]),
                              axis=1)
    
    df['IN_SAMPLE01'] = df.apply(lambda row: not row['FLAG_NOISE01'], axis=1)
    df['IN_TISSUE01'] = df.apply(lambda row: row.IN_SAMPLE01*
                                           value_tc[row.TISSUE_CELLS],
                               axis=1)
    df['IN_2501'] = df.apply(lambda row: row.IN_SAMPLE01*
                                       value_bulk.get(row.EXPERIMENT, 0),
                           axis=1)
    
    df['IN_SAMPLE95'] = df.apply(lambda row: not row['FLAG_NOISE95'], axis=1)
    df['IN_TISSUE95'] = df.apply(lambda row: row.IN_SAMPLE95 *
                                           value_tc[row.TISSUE_CELLS],
                               axis=1)
    df['IN_2595'] = df.apply(lambda row: row.IN_SAMPLE95 *
                                       value_bulk.get(row.EXPERIMENT, 0),
                           axis=1)
    
    df['IN_SAMPLE95rn'] = df.apply(lambda row: not row['FLAG_NOISE95rn'], axis=1)
    df['IN_TISSUE95rn'] = df.apply(lambda row: row.IN_SAMPLE95rn *
                                           value_tc[row.TISSUE_CELLS],
                               axis=1)
    df['IN_2595rn'] = df.apply(lambda row: row.IN_SAMPLE95rn *
                                       value_bulk.get(row.EXPERIMENT, 0),
                           axis=1)
    
    df['FLAG_HET_TISSUE'] = df.apply(lambda row: row.FLAG_HET *
                                                 value_tc[row.TISSUE_CELLS],
                                     axis=1)
    df['FLAG_HET_25'] = df.apply(lambda row: row.FLAG_HET *
                                             value_bulk.get(row.EXPERIMENT, 0),
                                 axis=1)
    df['FLAG_HET_CELLS'] = df.apply(lambda row: row.FLAG_HET *
                                             (1 - value_tc[row.TISSUE_CELLS]),
                                    axis=1)
    
#------------------------------------------------------------------------------
#make sorter and fix sn order

def make_sorter(df):
    '''make a sorter column.'''
    
    tc = {'cells': 'B', 'tissues': 'A'}
    organ = {'Ctx': 'A', 'Kidney': 'E', 'Heart': 'C', 'JGG': 'F',
             'Liver': 'D', 'Cbl': 'B'}
    region = {'T': 'E', 'no_region': 'F', 'PF': 'A', 'F': 'B', 'P': 'C',
              'O': 'D'}
    lr = {'L': 'A', 'R': 'B', 'no_lat': 'C'}
    ls = {'Sml': 'B', 'no_lrg_sml': 'X', 'Lrg': 'A', 'Sml11': 'L', 'Sml3': 'D',
          'Sml12': 'M', 'Sml7': 'H', 'Sml2': 'C', 'Sml8': 'I','Sml6': 'G',
          'Sml10': 'K', 'Sml9': 'J', 'Sml4': 'E', 'Sml5': 'F', 'Sml13': 'N'}
    exp = {'bulk': 'A', 'subsample': 'B', 'NeuN': 'C', 'TBR1': 'D',
           'Olig2': 'E', 'Lhx2': 'F', 'PU1': 'G', 'ctrl': 'X'}
    
    df['EXPERIMENT_SORTABLE'] = df.apply(fix_sn, axis=1)
    
    df['SORTER_STR'] = df.apply(lambda row: tc[row.TISSUE_CELLS] +
                                            organ[row.ORGAN] +
                                            lr[row.L_R] +
                                            region[row.BRAIN_REGION] +
                                            exp.get(row.EXPERIMENT_SORTABLE,
                                             ('W' + row.EXPERIMENT_SORTABLE)) +
                                            ls[row.LG_SM],
                                axis=1)
    
    sorters = df.SORTER_STR.sort_values().unique()
    
    sorter_dic = {sorter: i for sorter, i
                            in zip(sorters, range(1, (len(sorters) + 1)))}
    
    df['SORTER'] = df.apply(lambda row: sorter_dic[row.SORTER_STR], axis=1)
    
    df['CHROM_SORTER'] = df.apply(get_chrom_str, axis=1)
    
    df.sort_values(['CHROM_SORTER', 'CHR_POS_REF_ALT', 'SORTER'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    
def fix_sn(row):
    
    if '_' in row.EXPERIMENT:
        
        lst = row.EXPERIMENT.split('_')
        
        if len(lst[-1]) == 2:
            
            lst = [lst[0], ('0' + lst[-1][1]) + lst[-1][0]]
        
        else:
            lst = [lst[0], (lst[-1][1:]) + lst[-1][0]]

        
        return '_'.join(lst[::-1])
        
    else:
        return row.EXPERIMENT

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
    
def add_depth30_flag(df):
    '''use depth <30 as a flag for coverage issues. might need additional
    flags for other coverage types.'''
    
    value_bulk = {'bulk': 1}
    
    df['FLAG_30DEPTH'] = df.apply(lambda row: row.DEPTH < 30, axis=1)
    df['FLAG_30DEPTH_25'] = df.apply(lambda row: row.FLAG_30DEPTH *
                                             value_bulk.get(row.EXPERIMENT, 0),
                                     axis=1)
#------------------------------------------------------------------------------
#annotate signal, mosaic_support, real_mosaic

def annotate_mos(df):
    '''use the flags above to determine whether a variant is a real mosaic or
    not.'''
    
    df.sort_values(['CHROM_SORTER', 'CHR_POS_REF_ALT', 'SORTER'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    grp = df.groupby('CHR_POS_REF_ALT')
    
    df['FLAG_PRESENT_ALL'] = grp['IN_SAMPLE'].transform(flag_signal)
    df['FLAG_PRESENT_TISSUE'] = grp['IN_TISSUE'].transform(flag_signal)
    df['FLAG_PRESENT_CELLS'] = grp['IN_CELLS'].transform(flag_signal)
    df['FLAG_PRESENT_25'] = grp['IN_25'].transform(flag_signal)
    
    df['FLAG_PRESENT_ALL01'] = grp['IN_SAMPLE01'].transform(flag_signal)
    df['FLAG_PRESENT_TISSUE01'] = grp['IN_TISSUE01'].transform(flag_signal)
    df['FLAG_PRESENT_2501'] = grp['IN_2501'].transform(flag_signal)
    
    df['FLAG_PRESENT_ALL95'] = grp['IN_SAMPLE95'].transform(flag_signal)
    df['FLAG_PRESENT_TISSUE95'] = grp['IN_TISSUE95'].transform(flag_signal)
    df['FLAG_PRESENT_2595'] = grp['IN_2595'].transform(flag_signal)
    
    df['FLAG_PRESENT_ALL95rn'] = grp['IN_SAMPLE95rn'].transform(flag_signal)
    df['FLAG_PRESENT_TISSUE95rn'] = grp['IN_TISSUE95rn'].transform(flag_signal)
    df['FLAG_PRESENT_2595rn'] = grp['IN_2595rn'].transform(flag_signal)
    
    df['SUM_HET_TISSUES'] = grp['FLAG_HET_TISSUE'].transform('sum')
    df['SUM_HET_CELLS'] = grp['FLAG_HET_CELLS'].transform('sum')
    df['SUM_HET_25'] = grp['FLAG_HET_25'].transform('sum')
    
    df['SUM_DEPTH_FLAG'] = grp['FLAG_30DEPTH'].transform('sum')
    df['SUM_DEPTH_FLAG_25'] = grp['FLAG_30DEPTH_25'].transform('sum')
    
    df['SET_MOSAIC'] = (
                        (df.CATEGORY == 'MOSAIC') &
                        (df.FLAG_PRESENT_TISSUE95rn == True) &
                        (df.SUM_HET_25 < 13) &
                        (df.NORMAL_LOWER_CI < 0.00143692) &
                        (df.SUM_DEPTH_FLAG_25 < 13)
                        )
    
def flag_signal(col):
    
    return sum(col) > 0

#------------------------------------------------------------------------------

def make_mos(df):
    '''clean up df and only take mosaics.'''
    
    df = df[df.SET_MOSAIC == True]
    
    df.drop(['GNOMAD_FREQ', 'REPEAT_MASKER', 'SEGDUP', 'HOMOPOLYMER',
             'DINUCLEOTIDE', 'NEAR_INDEL', 'UCSC_RPMSK', 'FLAG_NOISE',
             'FLAG_NOISE01', 'FLAG_NOISE95','IN_SAMPLE', 'IN_TISSUE', 'IN_25',
             'IN_CELLS', 'IN_SAMPLE01', 'IN_TISSUE01', 'IN_2501',
             'IN_SAMPLE95', 'IN_TISSUE95', 'IN_2595', 'SORTER_STR',
             'FLAG_PRESENT_ALL', 'FLAG_PRESENT_TISSUE', 'FLAG_PRESENT_CELLS',
             'FLAG_PRESENT_25', 'FLAG_PRESENT_ALL01', 'FLAG_PRESENT_TISSUE01',
             'FLAG_PRESENT_2501', 'FLAG_PRESENT_ALL95',
             'FLAG_PRESENT_TISSUE95', 'FLAG_PRESENT_2595'],
            axis=1, inplace=True)
    
    df.rename({'IN_SAMPLE95rn': 'IN_SAMPLE', 'IN_TISSUE95rn': 'IN_TISSUE',
               'IN_2595rn': 'IN_25',
               'FLAG_PRESENT_ALL95rn': 'FLAG_PRESENT_ALL',
               'FLAG_PRESENT_TISSUE95rn': 'FLAG_PRESENT_TISSUE',
               'FLAG_PRESENT_2595rn': 'FLAG_PRESENT_25',
               'FLAG_NOISE95rn': 'FLAG_NOISE'}, axis=1, inplace=True)
        
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
                                        
    df['SUM_POS_TISSUE'] = df.groupby('CHR_POS_REF_ALT')['IN_TISSUE']\
                                        .transform('sum')
    
    df['SUM_POS_CELLS'] = df.SUM_POS_ALL - df.SUM_POS_TISSUE
    
    df['SUM_POS_25'] = df.groupby('CHR_POS_REF_ALT')['IN_25']\
                                        .transform('sum')


#------------------------------------------------------------------------------
#make a tissue sorter
def tissue_sorter(df):
    '''sort for the 25 regions.'''
    
    df['INTERMEDIATE'] = df.apply(lambda row: row.ORGAN + '_' +
                                              row.BRAIN_REGION,
                                  axis=1)
    dictionary = {'Ctx_PF': 1, 'Ctx_F': 2, 'Ctx_P': 3, 'Ctx_O': 4, 'Ctx_T': 5,
                  'Cbl_no_region': 6, 'Heart_no_region': 7,
                  'Liver_no_region': 8, 'Kidney_no_region': 9,
                  'JGG_no_region':10}
    
    df['ORGAN_SORTER'] = df.apply(lambda row: dictionary[row.INTERMEDIATE],
                                  axis=1)
    
    df.drop('INTERMEDIATE', axis=1, inplace=True)
    

#------------------------------------------------------------------------------
#get max is complicated by cells and by tissues that have very low coverage

def get_max_af(df):
    '''get the maximum AF for each variant and the lower/upper coordinates.'''
    
    df['CONSIDERED_MAF_TISSUE'] = df.apply(lambda row:
                                            row.MAF * (row.IN_TISSUE),
                                           axis=1)
    
    df['CONSIDERED_MAF_25'] = df.apply(lambda row: row.MAF * (row.IN_25),
                                           axis=1)
    
    grp = df.groupby('CHR_POS_REF_ALT') #generates groupby object
    
    df['MAX_MAF_TISSUE'] = grp['CONSIDERED_MAF_TISSUE'].transform('max')
    df['MAX_MAF_TISSUE_BOOL'] = df.apply(lambda row:
                                row['CONSIDERED_MAF_TISSUE'] ==
                                row['MAX_MAF_TISSUE'],
                                         axis=1)
    
    df['MAX_MAF_25'] = grp['CONSIDERED_MAF_25'].transform('max')
    df['MAX_MAF_25_BOOL_int'] = df.apply(lambda row: row['CONSIDERED_MAF_25']
                                                     == row['MAX_MAF_25'],
                                         axis=1)
    df['MAX_MAF_25_BOOL'] = (df.MAX_MAF_25_BOOL_int) & (df.FLAG_PRESENT_25)
    
    df.drop('MAX_MAF_25_BOOL_int', axis=1, inplace=True)
    
    
#+-----------------------------------------------------------------------------
#functions to annotate categories such as lateralized, brain specific etc.
def annotate_expression_pattern(df):
    '''use transform to annotate patterns as bool. note that some of the groups
    require a clean mosaic table, i.e. variants must be determined as positive
    in at least one tissue.'''
    
    sorter(df)
    
    df['SET_ONE_TISSUE'] = one_tissue(df)
    df['SET_ONE_REGION'] = one_region(df) #cortical regions
    df['SET_CORTEX_ONLY'] = cortex_only(df)
    df['SET_BRAIN_ONLY'] = brain_only(df)
    df['SET_KIDNEY_ONLY'] = kidney_only(df)
    df['SET_LEFT_ONLY'] = left_only(df)
    df['SET_RIGHT_ONLY'] = right_only(df)
    df['SET_IN_CORTEX'] = in_cortex(df)
    df['SET_IN_BRAIN'] = in_brain(df)
    df['SET_IN_CEREBELLUM'] = in_cerebellum(df)
    df['SET_IN_HEART'] = in_heart(df)
    df['SET_IN_LIVER'] = in_liver(df)
    df['SET_IN_KIDNEY'] = in_kidney(df)
    #laterlized in brain, but can be anywhere in other tissues
    df['SET_LEFT_ONLY_CTX'] = left_only_ctx(df)
    df['SET_RIGHT_ONLY_CTX'] = right_only_ctx(df)

def one_tissue(df):
    '''could be done easier with df.apply(lambda row: row.SUM_POS_25 ==1,
    axis=1); but used as template.'''
    
    df_bulk = df[df.EXPERIMENT == 'bulk']
    
    df_grouped = df_bulk.groupby(['CHR_POS_REF_ALT'])['IN_25'].sum()\
                                                              .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    
    return df_out.OUTPUT

def one_region(df):
    
    df_region = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN == 'Ctx') &
                   (df.SUM_POS_25 == 2)]
    
    df_grouped = df_region.groupby(['CHR_POS_REF_ALT', 'BRAIN_REGION', 'L_R'])\
                          ['IN_25'].sum().reset_index()
    df_grouped = df_grouped.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                                 .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 2, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def cortex_only(df):
    '''all only calls rely on a variant being mosaic, but absent in all other 
    tissues. e.g. for cortex it has to be not present in any non-ctx tissues.
    logic is used for all tyypes of _only.'''
    
    df_nonctx = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN != 'Ctx') &
                   (df.SUM_POS_25 > 0)]
    
    df_grouped = df_nonctx.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                                 .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 0, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def brain_only(df):
    
    df_nonbr = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN != 'Ctx') &
                  (df.ORGAN != 'Cbl') & (df.SUM_POS_25 > 0)]
    
    df_grouped = df_nonbr.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                               .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 0, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def kidney_only(df):
    
    df_nokid = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN != 'Kidney') &
                  (df.SUM_POS_25 > 0)]
    
    df_grouped = df_nokid.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                               .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 0, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def left_only(df):
    
    df_r = df[(df.EXPERIMENT == 'bulk') & (df.L_R != 'L') &
              (df.SUM_POS_25 > 0)]
    
    df_grouped = df_r.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                            .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 0, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def right_only(df):
    
    df_l = df[(df.EXPERIMENT == 'bulk') & (df.L_R != 'R') &
              (df.SUM_POS_25 > 0)]
    
    df_grouped = df_l.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                           .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 0, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_cortex(df):
    
    df_ctx = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN == 'Ctx')]
    
    df_grouped = df_ctx.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_brain(df):
    
    df_br = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN.isin(['Ctx', 'Cbl']))]
    
    df_grouped = df_br.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                            .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_cerebellum(df):
    
    df_cbl = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN =='Cbl')]
    
    df_grouped = df_cbl.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_heart(df):
    
    df_hea = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN =='Heart')]
    
    df_grouped = df_hea.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_liver(df):
    
    df_liv = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN =='Liver')]
    
    df_grouped = df_liv.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_kidney(df):
    
    df_kid = df[(df.EXPERIMENT == 'bulk') & (df.ORGAN == 'Kidney')]
    
    df_grouped = df_kid.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def left_only_ctx(df):
    
    df_cr = df[(df.EXPERIMENT == 'bulk') & (df.SET_IN_CORTEX == True) &
                ((df.L_R == 'R') & (df.ORGAN == 'Ctx'))]
    
    df_grouped = df_cr.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 0, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def right_only_ctx(df):
    
    df_cl = df[(df.EXPERIMENT == 'bulk') & (df.SET_IN_CORTEX == True) &
                ((df.L_R == 'L') & (df.ORGAN == 'Ctx'))]
    
    df_grouped = df_cl.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 0, axis=1)
    
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
    
    if row.SUM_POS_25 == 0:
        
        return 'not_in_25'
    
    if row.SET_IN_BRAIN == False: #look at organ only calls
        
        if row.SET_ONE_TISSUE == False:
            
            if row.SET_KIDNEY_ONLY == True:
                return 'both kidneys' #both kidneys still are not one_tissue
            else:
                return 'organs'
        
        elif row.SET_IN_HEART == True:
            return 'heart'
        elif row.SET_IN_LIVER == True:
            return 'liver'
        elif row.SET_IN_KIDNEY == True:
            return 'one kidney'
        else:
            return 'ZONK1' #control, should not appear
        
    else: #i.e. variant is in brain
        
        if row.SET_BRAIN_ONLY == False: #i.e. not ONLY in brain
            return 'body'
        
        else: #i.e. variant is only in brain
            
            if row.SET_LEFT_ONLY == False and row.SET_RIGHT_ONLY == False:
                return 'across hemispheres'
            elif row.SET_ONE_TISSUE == True:
                return 'single brain tissue'
            elif row.SET_LEFT_ONLY == True:
                return 'left only'
            elif row.SET_RIGHT_ONLY == True:
                return 'right only'
            else:
                return 'ZONK2'
    
    return 'ZONK3'
            
#------------------------------------------------------------------------------
    
#------------------------------------------------------------------------------
#extra annotations for plotting

def annotate_locations(df,
    path_to_key='additional_files/20200531_location_key.csv'):
    '''use the annotations provided by XY to annotate location (mainly for the
    small subsamples). path is the one given when working in the current
    folder. might have to be adjusted accordingly.'''
    
    key = pd.read_csv(path_to_key)
    
    dictionary = {}
    
    for index, row in key.iterrows():
        dictionary[row.ID] = row.LOCATION
    
    df['LOCATION'] = df.apply(lambda row: dictionary[row.ID], axis=1)
    
#---
def annotate_lateralization_presence(df):
    '''replace the LEFT_RIGHT column from XY with a lateralization column.'''
    
    df['NUMBER_LEFT'] = number_left(df)
    df['NUMBER_RIGHT'] = number_right(df)
    df['LATERALIZATION'] = df.apply(lateralization, axis=1)
    
    df['IN_PF'] = in_pf(df)
    df['IN_F'] = in_f(df)
    df['IN_P'] = in_p(df)
    df['IN_O'] = in_o(df)
    df['IN_T'] = in_t(df)
    
    
def number_left(df):
    
    df_bulk = df[(df.EXPERIMENT == 'bulk') & (df.L_R == 'L')]
    
    df_grouped = df_bulk.groupby(['CHR_POS_REF_ALT'])['IN_25'].sum()\
                                                           .reset_index()
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'IN_25']]
    df_tomerge.rename({'IN_25': 'OUTPUT'}, axis=1, inplace=True)
    df_out = pd.merge(df, df_tomerge, how='left')
    
    return df_out.OUTPUT

def number_right(df):
    
    df_bulk = df[(df.EXPERIMENT == 'bulk') & (df.L_R == 'R')]
    
    df_grouped = df_bulk.groupby(['CHR_POS_REF_ALT'])['IN_25'].sum()\
                                                           .reset_index()
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'IN_25']]
    df_tomerge.rename({'IN_25': 'OUTPUT'}, axis=1, inplace=True)
    df_out = pd.merge(df, df_tomerge, how='left')
    
    return df_out.OUTPUT
    
def lateralization(row):
    
    if row.SET_LEFT_ONLY == True:
        return 'Left_only'
    
    elif row.SET_RIGHT_ONLY == True:
        return 'Right_only'
    
    elif (row.NUMBER_LEFT == 0) and (row.NUMBER_RIGHT == 0):
        return 'Not_lateral_bulk'
    
    elif (row.SET_IN_CEREBELLUM or row.SET_IN_HEART or row.SET_IN_LIVER) and\
         (row.NUMBER_LEFT == 0):
        return 'Right_enriched'
    elif (row.SET_IN_CEREBELLUM or row.SET_IN_HEART or row.SET_IN_LIVER) and\
         (row.NUMBER_RIGHT == 0):
        return 'Right_enriched'
    
    elif row.NUMBER_LEFT/row.NUMBER_RIGHT >= 1.5:
        return 'Left_enriched'
    
    elif row.NUMBER_RIGHT/row.NUMBER_LEFT >= 1.5:
        return 'Right_enriched'
    
    elif (row.NUMBER_LEFT/row.NUMBER_RIGHT < 1.5 and
          row.NUMBER_RIGHT/row.NUMBER_LEFT < 1.5):
        return 'Both'
    
    else:
        return 'Zonk'
    
def in_pf(df):
    
    df_ctx = df[(df.EXPERIMENT == 'bulk') & (df.BRAIN_REGION == 'PF')]
    
    df_grouped = df_ctx.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_f(df):
    
    df_ctx = df[(df.EXPERIMENT == 'bulk') & (df.BRAIN_REGION == 'F')]
    
    df_grouped = df_ctx.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_p(df):
    
    df_ctx = df[(df.EXPERIMENT == 'bulk') & (df.BRAIN_REGION == 'P')]
    
    df_grouped = df_ctx.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_o(df):
    
    df_ctx = df[(df.EXPERIMENT == 'bulk') & (df.BRAIN_REGION == 'O')]
    
    df_grouped = df_ctx.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_t(df):
    
    df_ctx = df[(df.EXPERIMENT == 'bulk') & (df.BRAIN_REGION == 'T')]
    
    df_grouped = df_ctx.groupby(['CHR_POS_REF_ALT'])['IN_25'].max()\
                                                             .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_25 == 1, axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT
#---
#---
def annotate_presence_pf_tl(df):
    '''function to annotate presence in pf l and r or t l, for the purpose
    of the subsample plotting.'''
    
    df['IN_PF_L'] = in_area(df, 'L', 'PF')
    df['IN_PF_R'] = in_area(df, 'R', 'PF')
    df['IN_T_L'] = in_area(df, 'L', 'T')
    df['IN_ALL_PF_L'] = in_all_area(df, 'L', 'PF')
    df['IN_ALL_PF_R'] = in_all_area(df, 'R', 'PF')
    df['IN_ALL_T_L'] = in_all_area(df, 'L', 'T')
    
    
def in_area(df, side, brreg):
    
    df_area = df[(df.EXPERIMENT.isin(['bulk', 'subsample'])) &
                (df.L_R == side) & (df.BRAIN_REGION == brreg)]
    
    df_grouped = df_area.groupby(['CHR_POS_REF_ALT'])['IN_TISSUE'].max()\
                                                                 .reset_index()
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_TISSUE == 1,
                                            axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT

def in_all_area(df, side, brreg):
    
    df_area = df[(df.EXPERIMENT.isin(['bulk', 'subsample'])) &
                (df.L_R == side) & (df.BRAIN_REGION == brreg)]
    
    df_grouped = df_area.groupby(['CHR_POS_REF_ALT'])['IN_TISSUE'].sum()\
                                                                 .reset_index()
    df_grouped['COUNT'] = df_area.groupby(['CHR_POS_REF_ALT'])\
                            ['IN_TISSUE'].count().values
    df_grouped['OUTPUT'] = df_grouped.apply(lambda row: row.IN_TISSUE ==
                                                        row.COUNT,
                                            axis=1)
    
    df_tomerge = df_grouped[['CHR_POS_REF_ALT', 'OUTPUT']]
    df_out = pd.merge(df, df_tomerge, how='left')
    df_out.fillna(False, inplace=True)
    
    return df_out.OUTPUT
    
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Sorted Populations Annotations

def annotate_exclusion_sorts(df):
    '''based on coverage the following dictionary was defined with XY. Only
    include those sorted populations that are not listed here.'''
    
    bad_sorts = {'7614_Ctx_F_L_PU1': 0, '7614_Ctx_F_R_NeuN': 0,
                 '7614_Ctx_F_R_PU1' : 0, '7614_Ctx_O_L_PU1': 0,
                 '7614_Ctx_O_R_Lhx2': 0, '7614_Ctx_O_R_PU1': 0,
                 '7614_Ctx_PF_L_TBR1': 0, '7614_Ctx_P_L_PU1': 0,
                 '7614_Ctx_T_L_PU1': 0, '7614_Ctx_T_R_Lhx2': 0,
                 '7614_Ctx_PF_L_PU1': 0, '7614_Ctx_P_R_PU1': 0,
                 '7614_Ctx_T_L_Olig2': 0, '7614_Ctx_P_L_Lhx2':0}
    
    df['SET_SORTED'] = df.apply(lambda row: int(row.EXPERIMENT in
                                                ['NeuN', 'TBR1', 'Olig2',
                                                 'Lhx2', 'PU1']),
                                axis=1)
    
    df['SET_SORT_INCLUDED'] = df.apply(lambda row: row.SET_SORTED *
                                                   bad_sorts.get(row.ID, 1),
                                       axis=1)


#------------------------------------------------------------------------------
#sn annotations

def annotate_set_sn_positive(df):
    ''' Using a cutoff of 5.425525e-04 for single nuclei based on similar logic
    that was used to define the lower threshold for mosaicism. using the
    complete data set, cutoff was found using
    
    data[(data.CATEGORY == 'REF_HOMO') & (data.FLAG_30DEPTH == False) &
    (data.TISSUE_CELLS == 'cells') & (data.ORGAN != 'JGG')]\
    .LOWER_CI.describe(percentiles=[0.01, 0.5, 0.90, 0.95, 0.99])
    Out[443]: 
        count    2.135000e+03
        mean     1.452234e-04
        std      1.410194e-03
        min     -3.469447e-18
        1%      -8.673617e-19
        50%      0.000000e+00
        90%      3.037908e-04
        95%      5.425525e-04
        99%      2.039662e-03
        max      6.087990e-02
        Name: LOWER_CI, dtype: float64
        
        Because of noise issues a _05 category was added for everything. This
        might be the better category to use in any case.
    '''
    
    df['SET_CELL_POSITIVE'] = (
                               (df.TISSUE_CELLS == 'cells') &
                               (df.LOWER_CI > 5.425525e-04) &
                               (df.LOWER_CI > df.NORMAL_UPPER_CI) &
                               (df.DEPTH >= 30) & (df.ORGAN != 'JGG')
                              ).astype(int)
    
    df['SET_CP_TL_PRESENT'] = (
                               (df.SET_CELL_POSITIVE == True) &
                               (df.IN_T_L == True)
                              ).astype(int)
    
    df['NUMBER_CP'] = df.groupby('CHR_POS_REF_ALT')['SET_CELL_POSITIVE']\
                        .transform('sum')
    df['NUMBER_CP_TLP'] = df.groupby('CHR_POS_REF_ALT')['SET_CP_TL_PRESENT']\
                            .transform('sum')
    df['NUMBER_CP_CELL'] = df.groupby('EXPERIMENT')['SET_CELL_POSITIVE']\
                             .transform('sum')
    df['NUMBER_CP_TLP_CELL'] = df.groupby('EXPERIMENT')['SET_CP_TL_PRESENT']\
                                 .transform('sum')
    
    df['SET_CELL_POSITIVE_05'] = (
                                  (df.TISSUE_CELLS == 'cells') &
                                  (df.LOWER_CI > 0.005) &
                                  (df.LOWER_CI > df.NORMAL_UPPER_CI) &
                                  (df.DEPTH >= 30) & (df.ORGAN != 'JGG')
                                 ).astype(int)
    
    df['SET_CP_TL_PRESENT_05'] = (
                                  (df.SET_CELL_POSITIVE_05 == True) &
                                  (df.IN_T_L == True)
                                 ).astype(int)
    
    df['NUMBER_CP_05'] = df.groupby('CHR_POS_REF_ALT')['SET_CELL_POSITIVE_05']\
                           .transform('sum')
    df['NUMBER_CP_TLP_05'] = df.groupby('CHR_POS_REF_ALT')['SET_CP_TL_PRESENT_05']\
                               .transform('sum')
    df['NUMBER_CP_CELL_05'] = df.groupby('EXPERIMENT')['SET_CELL_POSITIVE_05']\
                                .transform('sum')
    df['NUMBER_CP_TLP_CELL_05'] = df.groupby('EXPERIMENT')['SET_CP_TL_PRESENT_05']\
                                    .transform('sum')

#---
def annotation_sn_usedforfig(df):
    
    cp = 'NUMBER_CP_TLP_05'
    dictionary_sn_clades = make_sn_dict()
    
    df['SET_SN_FIGURE'] = ((df[cp] > 1) & (df[cp] < 21)).astype(int)
    df['VARIANT_ORDER_FIG'] = sn_fig_variant_order(df)
    df['CLADE_FIG'] = df.apply(lambda row:
                                          dictionary_sn_clades.get(
                                           row.VARIANT_ORDER_FIG, ['ZONK'])[0],
                               axis=1)
    df['SUB_CLADE_FIG'] = df.apply(lambda row:
                                              dictionary_sn_clades.get(
                                        row.VARIANT_ORDER_FIG, [0, 'ZONK'])[1],
                                   axis=1)

def sn_fig_variant_order(df):
    
    cp = 'NUMBER_CP_TLP_05'
    present = 'SET_CP_TL_PRESENT_05'
    number_cell = 'NUMBER_CP_TLP_CELL_05'
    
    df_ = df[(df[cp] > 0) & (df[cp] < 21) & (df.TISSUE_CELLS == 'cells')]
    
    df_ = df_[['CHR_POS_REF_ALT', 'EXPERIMENT', cp, present, number_cell]]
    
    df_ = df_[(df_[cp] > 1) & (df_[number_cell] > 0)]
    
    df_.sort_values(by=[cp, 'CHR_POS_REF_ALT'], ascending=False, inplace=True)
    df_.reset_index(drop=True, inplace=True)
    df_.reset_index(inplace=True)
    
    variants = pd.Series(df_.CHR_POS_REF_ALT.unique()).reset_index()
    variants.rename({'index':'VARIANT_ORDER_FIG', 0: 'CHR_POS_REF_ALT'},
                    axis=1, inplace=True)
    
    df_out = pd.merge(df, variants, how='left')
    
    return df_out.VARIANT_ORDER_FIG

def make_sn_dict():
    
    '''manually assigned after generating the df used for plotting with above
    sections'''
    
    dictionary = {#Variant 0 Clade: I
                  0: ('I', 'Founder'), 8: ('I', 'I-1'), 10: ('I', 'I-1'),
                  13: ('I', 'I-2'), 31: ('I', 'I-1'),
                  #Variant 1 Clade: II
                  1: ('II', 'Founder'), 6: ('II', 'II-1'), 9: ('II', 'II-1'),
                  15: ('II', 'II-1'), 16: ('II', 'II-2'), 19: ('II', 'II-2'),
                  25: ('II', 'II-1'), 26: ('II', 'No_Sub'),
                  28: ('II', 'No_Sub'), 32: ('II', 'II-1'), 
                      #II-1 by imputation
                  23: ('II', 'II-1'),
                  #Variant 2 Clade: III
                  2: ('III', 'Founder'), 5: ('III', 'No_Sub'),
                  7: ('III', 'No_Sub'), 11: ('III', 'No_Sub'),
                  18: ('III', 'No_Sub'), 30: ('III', 'No_Sub'),
                  #Variant 12 Clade: IV
                  12: ('IV', 'Founder'), 17: ('IV', 'No_Sub'),
                  21: ('IV', 'No_Sub'),
                  #Variant 20 Clade: V
                  20: ('V', 'Founder'),
                  #Variant 22 Clade: VI
                  22: ('VI', 'Founder'), 29: ('VI', 'No_Sub'),
                  #Variant 27 Clade: VII
                  27: ('VII', 'Founder')
                  }
    
    return dictionary
#---

def make_normAF(df):
    
    '''make a normalized AF for each variant (only considering positive tissues
    and excluding 0s).'''
    
    df_no0 = df.replace(0, np.NaN)
    
    df['CON_MAF_T_MEAN'] = df_no0.groupby('CHR_POS_REF_ALT')\
                               .CONSIDERED_MAF_TISSUE.transform('mean')
    df['NORM_AF'] = df.CONSIDERED_MAF_TISSUE/df.CON_MAF_T_MEAN
    
    
