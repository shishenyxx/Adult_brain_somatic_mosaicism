##Codes for the plotting of all main figure panels based on MPAS and snMPAS annotation

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
    
    
    
###############################################################################
###############################################################################
###############################################################################

#Plotting~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Base Plotting
#------------------------------------------------------------------------------

def plot_bars_for_categories(df):
    '''plot several categories as boxes side by side with the same
    visualization. swarmflag indicates whether swarmplot is included or not.'''
    
    df = df[df.MAX_MAF_25_BOOL == True]
    
    fig, axs = plt.subplots(nrows=1, ncols=6)

    all_mos = df
    across_body = df[(df.SET_IN_BRAIN) & ~(df.SET_BRAIN_ONLY)]
    brain_only = df[df.SET_BRAIN_ONLY]
    cortex_only = df[df.SET_CORTEX_ONLY]
    brain_side_restricted = df[((df.SET_LEFT_ONLY) | (df.SET_RIGHT_ONLY)) &
                                (df.SET_CORTEX_ONLY)]
    brain_one_tissue_only = df[(df.SET_ONE_TISSUE) & (df.SET_BRAIN_ONLY)]
    
    plots = [all_mos, across_body, brain_only, cortex_only,
             brain_side_restricted, brain_one_tissue_only]
    
    colors = ['0.5', '0.8', 'xkcd:dark blue', 'xkcd:blue', 'xkcd:sky blue',
              'xkcd:pale blue']
    
    names = ['All Variants', 'Brain and Organs', 'Brain Only', 'Cortex Only',
             'Brain Lateralized', 'Brain Single Sample']
    
    for plot, color, name, i in zip(plots, colors, names, range(7)):
        
        sns.countplot(x='MAX_MAF_25_BOOL', data=plot, color=color,
                      edgecolor='k', ax=axs[i])
        
        axs[i].set_ylim(-5,300)
        axs[i].set_xlabel(name, rotation=45, ha='right')
        axs[i].set_xticks([])
        
        if i != 0:
            sns.despine(left=True, bottom=True, offset=5, trim=True, ax=axs[i])
            axs[i].set_yticks([])
            axs[i].set_ylabel('')
        else:
            sns.despine(bottom=True, offset=5, trim=True, ax=axs[i])
            axs[i].set_ylabel('Number of Variants')
    
    plt.show()
    
    return [len(pl) for pl in plots]
    
    
def plot_boxes_for_categories(df, swarm_flag=False):
    '''plot several categories as boxes side by side with the same
    visualization. swarmflag indicates whether swarmplot is included or not.'''
    
    df = df[df.MAX_MAF_25_BOOL == True]
    
    df['MAF'] = df.MAF**0.5
    
    fig, axs = plt.subplots(nrows=1, ncols=6)

    all_mos = df
    across_body = df[(df.SET_IN_BRAIN) & ~(df.SET_BRAIN_ONLY)]
    brain_only = df[df.SET_BRAIN_ONLY]
    cortex_only = df[df.SET_CORTEX_ONLY]
    brain_side_restricted = df[((df.SET_LEFT_ONLY) | (df.SET_RIGHT_ONLY)) &
                                (df.SET_CORTEX_ONLY)]
    brain_one_tissue_only = df[(df.SET_ONE_TISSUE) & (df.SET_BRAIN_ONLY)]
    
    plots = [all_mos, across_body, brain_only, cortex_only,
             brain_side_restricted, brain_one_tissue_only]
    
    colors = ['0.5', '0.8', 'xkcd:dark blue', 'xkcd:blue', 'xkcd:sky blue',
              'xkcd:pale blue']
    
    names = ['All Variants', 'Brain and Organs', 'Brain Only', 'Cortex Only',
             'Brain Lateralized', 'Brain Single Sample']
    
    for plot, color, name, i in zip(plots, colors, names, range(7)):
        
        sns.boxplot(x='MAX_MAF_25_BOOL', y='MAF', data=plot, color=color,
                    ax=axs[i])
        if swarm_flag == True:
            sns.swarmplot(x='MAX_MAF_25_BOOL', y='MAF', data=plot, color='w',
                          edgecolor='k', linewidth=1, alpha=0.5, ax=axs[i])
        
        axs[i].set_ylim(0, 1)
        axs[i].set_yticks([0, 0.05**0.5, 0.2**0.5, 0.5**0.5, 1])
        axs[i].set_yticklabels(['0.00', '0.05', '0.20', '0.50', '1.00'])
        axs[i].set_xlabel(name, rotation=45, ha='right')
        axs[i].set_xticks([])
        
        if i != 0:
            sns.despine(left=True, bottom=True, offset=5, trim=True, ax=axs[i])
            axs[i].set_yticks([])
            axs[i].set_ylabel('')
        else:
            sns.despine(bottom=True, offset=5, trim=True, ax=axs[i])
            axs[i].set_ylabel('Allelic Fraction (sqrt-t)')
    
    plt.show()
    
    
def plot_with_custum_ci_SP_highnumber(df, x=300, y=0.51):
    '''use the reduced dataframe as input (e.g. use only brain variants).'''
    
    df = df[df.MAX_MAF_25_BOOL == True]
    df = df.sort_values(by=['MAF'], ascending=False)
    df = df.reset_index(drop=True).reset_index()
    plt.plot((df.iloc[:,0] + 1), df['MAF'], color='g')
    plt.plot((df.iloc[:,0] + 1), df['UPPER_CI'],
             color='0.5')
    plt.plot((df.iloc[:,0] + 1), df['LOWER_CI'],
             color='0.5')
    plt.xlabel('Ranked Mosaic Variants')
    plt.ylabel('Allelic Fraction')
    plt.xlim(0,x)
    plt.ylim(0,y)
    sns.despine(offset=5, trim=True)
    plt.show()
    
    
def plot_line(df, x=300, y=0.51):
    '''use the reduced dataframe as input (e.g. use only brain variants).'''
    
    df = df[df.MAX_MAF_25_BOOL == True]
    df = df.sort_values(by=['MAF'], ascending=False)
    df = df.reset_index(drop=True).reset_index()
    plt.fill_between((df.iloc[:,0] + 1), df['UPPER_CI'], df['LOWER_CI'],
                     color='0.8')
    plt.plot((df.iloc[:,0] + 1), df['MAF'], color='g')
    plt.xlabel('Ranked Mosaic Variants')
    plt.ylabel('Allelic Fraction')
    plt.xlim(0,x)
    plt.ylim(0,y)
    sns.despine(offset=5, trim=True)
    plt.show()
    
    
def plot_violin(df, y=0.61):
    '''use the reduced dataframe as input (e.g. use only brain variants).'''
    
    df = df[df.MAX_MAF_25_BOOL == True]
    
    sns.violinplot(x='SET_MOSAIC', y='MAF', data=df, inner=None, color='0.4',
                   cut=0)
    
    plt.xlabel('Mosaic Variants')
    plt.xticks(ticks=[], labels=[])
    plt.yticks(ticks=[0, 0.3, 0.6])
    plt.ylabel('')
    plt.ylim(0,y)
    sns.despine(bottom=True, offset=5, trim=True)
    plt.show()

    
def plot_sml_1tissues(df, y=50):
    '''Plot bars counting the number of variants that are only present in one
    small region.'''
    
    df = df[(df.MAX_MAF_25_BOOL) & (df.SET_CORTEX_ONLY) & (df.SET_ONE_TISSUE) &
            (df.LG_SM == 'Sml')].sort_values(by=['ORGAN_SORTER', 'L_R'])
    
    #to fix the issue with 0bars
    fake = pd.DataFrame(data={'BRAIN_REGION': ['F', 'O'],
                              'L_R': ['R', 'L']})
    
    df = df.append(fake, sort=False).reset_index(drop=True)
    
    sns.countplot(x='BRAIN_REGION', hue='L_R', data=df,
                  order=['PF', 'F', 'P', 'O', 'T'], edgecolor='k')
    
    plt.ylim(0, y)
    plt.ylabel('Number of Variants')
    plt.xlabel('Cortical Lobes')
    plt.title('Remove R-F and L-O!')
    sns.despine(bottom=True, offset=5, trim=True)
    
    plt.show()
    
    
#base function phage af plot---------------------------------------------------

def make_AF_plot_variant(df, POS, depth=100, return_flag=False):
    '''take the position and plot AF in all tissues for this position in boxes
    that are graded by AF level.'''
    
    #initialize data
    pos = df[(df.POS == POS) & (df.EXPERIMENT == 'bulk')][['CHR_POS_REF_ALT',
                                                           'MAF', 'SORTER',
                                                           'DEPTH',
                                                           'FLAG_NOISE']]
    pos.sort_values(by='SORTER', inplace=True)
    
    #make color list/hatches
    colors = make_colors(pos)
    hatches = make_hatches(pos, depth)
    
    #make boxes
    fig, ax = plt.subplots(1)

    boxes = make_boxes()
    
    for box, color, hatch in zip(boxes, colors, hatches):
        box.set_facecolor(color)
        box.set_edgecolor('k')
        box.set_hatch(hatch)
        ax.add_patch(box)
        
    ax.plot()
    sns.despine(left=True, bottom=True)
    plt.tick_params(bottom=False, left=False,
                    labelbottom=False, labelleft=False)
    plt.xlim(-14.75, 35.25)
    plt.ylim(0, 50)
    
    plt.title(pos.CHR_POS_REF_ALT.unique()[0])

    if return_flag == True:
        return ax
    
    else:
        plt.show()
    

#support functions+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def make_colors(pos):
    
    colors = []
    
    for index, row in pos.iterrows():
        
        if row.FLAG_NOISE == 'ZONK':
            
            colors.append('k')
        
        elif row.FLAG_NOISE == True:
            
            colors.append('w')
        
        else:
            colors.append((1. ,0. , 0.,row.MAF**(1/2)))
    
    for index in list(range(1, 29, 3)):
        
        colors.insert(index, 'w')
    
    return colors

def make_hatches(pos, depth):
    
    hatches = []
    
    for index, row in pos.iterrows():
        
        if row.DEPTH >= depth:
            hatches.append('')
            
        else:
            hatches.append('\\\\')
    
    for index in list(range(1, 29, 3)):
        hatches.insert(index, '')
    
    return hatches


def make_boxes():
    
    boxes = [Rectangle((0, 40), 10, 5), Rectangle((2, 41), 6, 3),
             Rectangle((2, 41), 6, 3),
             Rectangle((0, 34.5), 10, 5), Rectangle((2, 35.5), 6, 3),
             Rectangle((2, 35.5), 6, 3),
             Rectangle((0, 29), 10, 5), Rectangle((2 ,30), 6, 3),
             Rectangle((2 ,30), 6, 3),
             Rectangle((0, 23.5), 10, 5), Rectangle((2, 24.5), 6, 3),
             Rectangle((2, 24.5), 6, 3),
             Rectangle((0, 18), 10, 5), Rectangle((2, 19), 6, 3),
             Rectangle((2, 19), 6, 3),
             Rectangle((10.5, 40), 10, 5), Rectangle((12.5, 41), 6, 3),
             Rectangle((12.5, 41), 6, 3),
             Rectangle((10.5, 34.5), 10, 5), Rectangle((12.5, 35.5), 6, 3),
             Rectangle((12.5, 35.5), 6, 3),
             Rectangle((10.5, 29), 10, 5), Rectangle((12.5, 30), 6, 3),
             Rectangle((12.5, 30), 6, 3),
             Rectangle((10.5, 23.5), 10, 5), Rectangle((12.5, 24.5), 6, 3),
             Rectangle((12.5, 24.5), 6, 3),
             Rectangle((10.5, 18), 10, 5), Rectangle((12.5, 19), 6, 3),
             Rectangle((12.5, 19), 6, 3),
             Rectangle((5,12.5), 10.5, 5),
             Rectangle((5, 9.5), 10.5, 2.5),
             Rectangle((5, 6.5), 10.5, 2.5),
             Rectangle((0, 3.5), 10, 2.5), Rectangle((10.5, 3.5), 10, 2.5)]
    
    return boxes

#multiple plot of phage--------------------------------------------------------

def multi_phage_plot(df, cols=4):
    '''use phage plot function for all variants making master figure.'''
    
    number = len(df.CHR_POS_REF_ALT.unique())
    rows = math.ceil(number/cols) # rounds up
    
    fig, axs = plt.subplots(nrows=rows, ncols=cols)
    
    indexes = []
    
    for i in range(rows):
        for j in range(cols):
            indexes.append((i,j))
    
    indexes
    
    for pos, ij in zip(df.POS.unique(), indexes[:number]):
        
        make_AF_plot_variant_onax(df, pos, axs[ij[0], ij[1]])
        axs[ij[0], ij[1]].set_xticks([], [])
        axs[ij[0], ij[1]].set_yticks([], [])
    
    for i, j in indexes[number:]:
        axs[i, j].set_xticks([], [])
        axs[i, j].set_yticks([], [])
    
    plt.show()


def multi_phage_plot_save(df, name, cols=4):
    '''use phage plot function for all variants making master figure. this is
    the save version of the function that takes name. name, in the easiest
    version should just be the file's name and then it will be saved locally.
    Note that savefig can be used with arguments.'''
    
    number = len(df.CHR_POS_REF_ALT.unique())
    rows = math.ceil(number/cols) # rounds up
    
    fig, axs = plt.subplots(nrows=rows, ncols=cols)
    
    indexes = []
    
    for i in range(rows):
        for j in range(cols):
            indexes.append((i,j))
    
    indexes
    
    for pos, ij in zip(df.POS.unique(), indexes[:number]):
        
        make_AF_plot_variant_onax(df, pos, axs[ij[0], ij[1]])
        axs[ij[0], ij[1]].set_xticks([], [])
        axs[ij[0], ij[1]].set_yticks([], [])
    
    for i, j in indexes[number:]:
        axs[i, j].set_xticks([], [])
        axs[i, j].set_yticks([], [])
    
    plt.savefig(name, orientation='landscape')
    
    
def make_AF_plot_variant_onax(df, POS, ax, depth=100):
    '''take the position and plot AF in all tissues for this position in boxes
    that are graded by AF level.'''
    
    #initialize data
    
        
    pos = df[(df.POS == POS) & (df.EXPERIMENT == 'bulk')][['CHR_POS_REF_ALT',
                                                           'MAF', 'SORTER',
                                                           'DEPTH',
                                                           'FLAG_NOISE']]
    pos.sort_values(by='SORTER', inplace=True)
    
    #make color list/hatches
    colors = make_colors(pos)
    hatches = make_hatches(pos, depth)
    
    #make boxes
    boxes = make_boxes()
    
    for box, color, hatch in zip(boxes, colors, hatches):
        box.set_facecolor(color)
        box.set_edgecolor('k')
        box.set_hatch(hatch)
        ax.add_patch(box)
        
    ax.plot()
    sns.despine(left=True, bottom=True)
    ax.set_xlim(-14.75, 35.25)
    ax.set_ylim(0, 50)
    
    ax.set_title(pos.CHR_POS_REF_ALT.unique()[0], fontsize=6) #fontsize=6

    return ax
    

#make tissues instead of phage-------------------------------------------------

def make_tissues():
    
    tissues = [
               #right brain
               Polygon([(0.882, 24.605), (2.254, 29.331), (7.527, 34.549),
                        (14.173, 24.546), (12.289, 19.854), (4.018, 19.173)]),
               Circle((6.822, 24.921), 2), Circle((6.822, 24.921), 2),
               Polygon([(7.527, 34.549), (19.817, 38.959), (20.921, 25.51),
                        (14.173, 24.546)]),
               Circle((15.415, 30.382), 2), Circle((15.415, 30.382), 2),
               Polygon([(19.817, 38.959), (29.125, 37.856), (34.955, 32.497),
                        (28.258, 26.558), (20.921, 25.51)]),
               Circle((25.49, 32.08), 2), Circle((25.49, 32.08), 2),
               Polygon([(28.258, 26.558), (34.955, 32.497), (38.08, 26.938),
                        (37.121, 20.218), (29.496, 17.889), (29.968, 23.639)]),
               Circle((33.729, 24.921), 2), Circle((33.729, 24.921), 2),
               Polygon([(12.289, 19.854), (14.173, 24.546), (28.258, 26.558),
                        (29.968, 23.639), (29.496, 17.889), (15.411, 15.877)]),
               Circle((21.689, 21.217), 2), Circle((21.689, 21.217), 2),
            
               #left brain
               Polygon([(-0.882, 24.605), (-2.254, 29.331), (-7.527, 34.549),
                        (-14.173, 24.546), (-12.289, 19.854),
                        (-4.018, 19.173)]),
               Circle((-6.822, 24.921), 2), Circle((-6.822, 24.921), 2),
               Polygon([(-7.527, 34.549), (-19.817, 38.959), (-20.921, 25.51),
                        (-14.173, 24.546)]),
               Circle((-15.415, 30.382), 2), Circle((-15.415, 30.382), 2),
               Polygon([(-19.817, 38.959), (-29.125, 37.856),
                        (-34.955, 32.497), (-28.258, 26.558),
                        (-20.921, 25.51)]),
               Circle((-25.49, 32.08), 2), Circle((-25.49, 32.08), 2),
               Polygon([(-28.258, 26.558), (-34.955, 32.497), (-38.08, 26.938),
                        (-37.121, 20.218), (-29.496, 17.889),
                        (-29.968, 23.639)]),
               Circle((-33.729, 24.921), 2), Circle((-33.729, 24.921), 2),
               Polygon([(-12.289, 19.854), (-14.173, 24.546),
                        (-28.258, 26.558), (-29.968, 23.639),
                        (-29.496, 17.889), (-15.411, 15.877)]),
               Circle((-21.689, 21.217), 2), Circle((-21.689, 21.217), 2),

               #cerebellum, heart, liver, kidney-r, kidney-l
               Polygon([(-38.08, 18.152), (-37.121, 20.218), (-29.496, 17.889),
                        (-30.723, 14.447), (-34.09, 12.857),
                        (-38.08, 14.486)]),
               Polygon([(-3.137, 12.612), (-1.569, 15.329), (0, 12.668),
                        (1.569, 15.329), (3.137, 12.612), (2.258, 10.816),
                        (0, 8.509), (-2.258, 10.816)]),
               Polygon([(-5.25, 4.003), (-3.871, 6.392), (4.321, 6.392),
                        (5.251, 4.869), (1.186, 3.246), (-1.625, 2.966),
                        (-2.186, 0), (-4.945, 0)]),
               Polygon([(13.312, 7.371), (11.722, 10.125), (8.541, 10.125),
                        (10.126, 7.371), (8.531, 4.616), (11.722, 4.616)]),
               Polygon([(-13.312, 7.371), (-11.722, 10.125), (-8.541, 10.125),
                        (-10.126, 7.371), (-8.531, 4.616), (-11.722, 4.616)]),
            ]
               
    return tissues


def make_AF_plot_variant_tissues(df, POS):
    '''take the position and plot AF in all tissues for this position in boxes
    that are graded by AF level.'''
    
    #initialize data
    
    #needed b/c of incomplete list of tissues
    if len(df.SORTER.unique()) != 25:
    
        pos = df[df.POS == POS][['MAF', 'SORTER']]
        flag = df[df.POS == POS][['CHR_POS_REF_ALT', 'FLAG_NOISE', 'SORTER']]
    
        pos.sort_values(by='SORTER', inplace=True)
    
        dummy = pd.DataFrame(data=[i for i in range(25)], columns=['SORTER'])
        dummy['MAF'] = 0
        pos = pos.append(dummy, sort=True).groupby('SORTER').sum()\
                                                            .reset_index()
        
        pos = pd.merge(pos, flag, how='left')
        pos.fillna('ZONK', inplace=True)
    
    else:
        
        pos = df[df.POS == POS][['CHR_POS_REF_ALT', 'MAF', 'SORTER',
                                 'FLAG_NOISE']]
        pos.sort_values(by='SORTER', inplace=True)
    
    #make color list
    
    colors = make_colors(pos)
    
    #make boxes
    
    fig,ax = plt.subplots(1)

    boxes = make_tissues()
        
    pc = PatchCollection(boxes, facecolor=colors, edgecolor='k')
    ax.add_collection(pc)
    
    ax.plot()
    sns.despine(left=True, bottom=True)
    plt.tick_params(bottom=False, left=False,
                    labelbottom=False, labelleft=False)
    plt.xlim(-40, 40)
    plt.ylim(-1, 59)
    
    plt.title(pos.CHR_POS_REF_ALT.unique()[0])

    plt.show()

#------------------------------------------------------------------------------
#subsample plot
    
def make_subplot_variant(df, POS, side='L', brreg='T', depth=100):
    '''take the position and plot AF in subtissues for this position in boxes
    that are graded by AF level.'''
    
    #initialize data
    df = df[(df.POS == POS) & (df.EXPERIMENT.isin(['bulk', 'subsample'])) &
            (df.L_R == side) & (df.BRAIN_REGION == brreg)][['CHR_POS_REF_ALT',
                                                           'LOCATION', 'DEPTH',
                                                      'CONSIDERED_MAF_TISSUE']]
    
    if (brreg != 'PF') and not (brreg == 'T' and side == 'L'):
        return 'NO SUBSAMPLES AVAILABLE'
    
    all_sites = ['LRG', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
                 'L', 'M']
    available = df.LOCATION.tolist()
    missing = [loc for loc in all_sites if loc not in available]
    
    #make colors list for sml
    colors = make_colors_sub(df, missing)
    hatches = make_hatches_sub(df, missing, depth)
    
    #make boxes
    fig, ax = plt.subplots(1)

    boxes = make_boxes_sub()
    
    for box, color, hatch in zip(boxes, colors, hatches):
        box.set_facecolor(color)
        box.set_edgecolor('k')
        box.set_hatch(hatch)
        ax.add_patch(box)
        

    #make poly for lrg
    poly_lrg = Polygon([(8, 41), (8, 62), (29, 62), (29, 52), ( 18, 52),
                        (18, 41)])
    
    poly_color = (1. ,0. , 0.,
                  df[df.LOCATION == 'LRG'].iloc[0, 3]**(1/2))
    
    if df[df.LOCATION == 'LRG'].iloc[0, 2] >= depth:
        hatch_lrg = ''
    else:
        hatch_lrg = '\\\\'
    
    pc2 = PatchCollection([poly_lrg], facecolor=[poly_color], edgecolor='k',
                          hatch=hatch_lrg)
    ax.add_collection(pc2)
    
    #plot and prepare
    ax.plot()
    ax.set(adjustable='box', aspect='equal')
    sns.despine(left=True, bottom=True)
    plt.tick_params(bottom=False, left=False,
                    labelbottom=False, labelleft=False)
    plt.xlim(0, 70)
    plt.ylim(0, 70)
    
    plt.title(df.CHR_POS_REF_ALT.unique()[0])

    plt.show()

def make_subplot_variant_onax(df, POS, ax, side, brreg, depth=100):
    '''take the position and plot AF in subtissues for this position in boxes
    that are graded by AF level. return ax after taking ax as argument.'''
    
    #initialize data
    df = df[(df.POS == POS) & (df.EXPERIMENT.isin(['bulk', 'subsample'])) &
            (df.L_R == side) & (df.BRAIN_REGION == brreg)][['CHR_POS_REF_ALT',
                                                           'LOCATION', 'DEPTH',
                                                      'CONSIDERED_MAF_TISSUE']]
    
    if (brreg != 'PF') and not (brreg == 'T' and side == 'L'):
        return 'NO SUBSAMPLES AVAILABLE'
    
    all_sites = ['LRG', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
                 'L', 'M']
    available = df.LOCATION.tolist()
    missing = [loc for loc in all_sites if loc not in available]
    
    #make colors list for sml
    colors = make_colors_sub(df, missing)
    hatches = make_hatches_sub(df, missing, depth)
    
    #make boxes
    boxes = make_boxes_sub()
    
    for box, color, hatch in zip(boxes, colors, hatches):
        box.set_facecolor(color)
        box.set_edgecolor('k')
        box.set_hatch(hatch)
        ax.add_patch(box)

    #make poly for lrg
    poly_lrg = Polygon([(8, 41), (8, 62), (29, 62), (29, 52), ( 18, 52),
                        (18, 41)])
    
    poly_color = (1. ,0. , 0.,
                  df[df.LOCATION == 'LRG'].iloc[0, 3]**(1/2))
    
    if df[df.LOCATION == 'LRG'].iloc[0, 2] >= depth:
        hatch_lrg = ''
    else:
        hatch_lrg = '\\\\'
    
    pc2 = PatchCollection([poly_lrg], facecolor=[poly_color], edgecolor='k',
                          hatch=hatch_lrg)
    ax.add_collection(pc2)
    
    #plot and prepare
    ax.plot()
    ax.set(adjustable='box', aspect='equal')
    sns.despine(left=True, bottom=True)
    ax.tick_params(bottom=False, left=False,
                   labelbottom=False, labelleft=False)
    ax.set_xlim(0, 70)
    ax.set_ylim(0, 70)
    
    ax.set_title(df.CHR_POS_REF_ALT.unique()[0], fontsize=6)

    return ax


def multi_subplot(df, side='L', brreg='T', cols=10):
    '''use subplot function for all variants making master figure.'''
    
    number = len(df.CHR_POS_REF_ALT.unique())
    rows = math.ceil(number/cols) # rounds up
    
    fig, axs = plt.subplots(nrows=rows, ncols=cols)
    
    indexes = []
    
    for i in range(rows):
        for j in range(cols):
            indexes.append((i,j))
    
    indexes
    
    for pos, ij in zip(df.POS.unique(), indexes[:number]):
        
        make_subplot_variant_onax(df, pos, axs[ij[0], ij[1]], side, brreg)
        axs[ij[0], ij[1]].set_xticks([], [])
        axs[ij[0], ij[1]].set_yticks([], [])
    
    for i, j in indexes[number:]:
        axs[i, j].set_xticks([], [])
        axs[i, j].set_yticks([], [])
    
    plt.show()

######Combination Plot for Phage and all three Regions######
    
def combination_plot(df, POS):
    '''use the four onax plots to obtain a full view of a single variant.'''
    
    ax1 = plt.subplot2grid((2, 3), (0, 0))
    ax2 = plt.subplot2grid((2, 3), (1, 0))
    ax3 = plt.subplot2grid((2, 3), (0, 1), rowspan=2)
    ax4 = plt.subplot2grid((2, 3), (0, 2))
    
    make_subplot_variant_onax(df, POS, ax1, 'L', 'PF')
    make_subplot_variant_onax(df, POS, ax2, 'L', 'T')
    make_AF_plot_variant_onax(df[df.EXPERIMENT == 'bulk'], POS, ax3)
    make_subplot_variant_onax(df, POS, ax4, 'R', 'PF')
    
    ax3.set_xlim(-1, 21.5)
    ax3.set_ylim(2.5, 46)
    ax3.tick_params(bottom=False, left=False,
                    labelbottom=False, labelleft=False)
    
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_title('')
    
    ax1.set_title('L-PF')
    ax2.set_title('L-T')
    ax3.set_title(df[df.POS == POS].CHR_POS_REF_ALT.unique()[0])
    ax4.set_title('R-PF')
    
    plt.show()
    
    
def combination_plot_pdf(df, path):
    '''use the four onax plots to obtain a full view of a single variant.'''
    
    pdf_pages = PdfPages(path)
    
    df = df.fillna(0)
    
    for POS in df.POS.unique():
        
        fig = plt.figure(figsize=(10, 5), dpi=100)
    
        ax1 = plt.subplot2grid((2, 3), (0, 0))
        ax2 = plt.subplot2grid((2, 3), (1, 0))
        ax3 = plt.subplot2grid((2, 3), (0, 1), rowspan=2)
        ax4 = plt.subplot2grid((2, 3), (0, 2))
        
        make_subplot_variant_onax(df, POS, ax1, 'L', 'PF')
        make_subplot_variant_onax(df, POS, ax2, 'L', 'T')
        make_AF_plot_variant_onax(df[df.EXPERIMENT == 'bulk'], POS, ax3)
        make_subplot_variant_onax(df, POS, ax4, 'R', 'PF')
        
        ax3.set_xlim(-1, 21.5)
        ax3.set_ylim(2.5, 46)
        ax3.tick_params(bottom=False, left=False,
                        labelbottom=False, labelleft=False)
    
        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_title('')
    
        ax1.set_title('L-PF')
        ax2.set_title('L-T')
        ax3.set_title(df[df.POS == POS].CHR_POS_REF_ALT.unique()[0])
        ax4.set_title('R-PF')
        
        pdf_pages.savefig(fig)
        plt.close()
    
    pdf_pages.close()
    
#---
def make_boxes_sub():
    
    boxes = [Rectangle((30, 30), 10, 10), Rectangle((30, 41), 10, 10),
             Rectangle((41, 30), 10, 10), Rectangle((30, 19), 10, 10),
             Rectangle((19, 30), 10, 10), Rectangle((30, 52), 10, 10),
             Rectangle((52, 30), 10, 10), Rectangle((30, 8), 10, 10),
             Rectangle((8, 30), 10, 10), Rectangle((41, 41), 10, 10),
             Rectangle((41, 19), 10, 10), Rectangle((19, 19), 10, 10),
             Rectangle((19, 41), 10, 10)]
    
    return boxes

def make_colors_sub(df, missing):
    
    pairs = df[df.LOCATION != 'LRG'][['LOCATION', 'CONSIDERED_MAF_TISSUE']]
    
    missing_pairs = [(miss, 'no_data') for miss in missing]
    
    miss_df = pd.DataFrame(data=missing_pairs, columns=['LOCATION',
                                                      'CONSIDERED_MAF_TISSUE'])
        
    complete = pairs.append(miss_df)
    complete.sort_values('LOCATION', inplace=True)
    
    colors = []
    
    for index, row in complete.iterrows():
        
        if row.CONSIDERED_MAF_TISSUE == 'no_data':
            
            colors.append('0.8')
        
        else:
            colors.append((1. ,0. , 0.,row.CONSIDERED_MAF_TISSUE**(1/2)))
    
    return colors

def make_hatches_sub(df, missing, depth):
    
    trips = df[df.LOCATION != 'LRG'][['LOCATION', 'CONSIDERED_MAF_TISSUE',
                                      'DEPTH']]
    
    missing_trips = [(miss, 'no_data', 'no_depth') for miss in missing]
    
    miss_df = pd.DataFrame(data=missing_trips, columns=['LOCATION',
                                                       'CONSIDERED_MAF_TISSUE',
                                                       'DEPTH'])
        
    complete = trips.append(miss_df)
    complete.sort_values('LOCATION', inplace=True)
    
    hatches = []
    
    for index, row in complete.iterrows():
        
        if row.DEPTH == 'no_depth':
            hatches.append('')
        
        elif row.DEPTH >= depth:
            hatches.append('')
            
        else:
            hatches.append('\\\\')
    
    return hatches
#---
#------------------------------------------------------------------------------

#make clustermap---------------------------------------------------------------
    
def make_clustermap(df, category='CHR_POS_REF_ALT', method='weighted'):
    '''use pandas to generate a table for clustering. category is 'ID' or
    it is 'CHR_POS_REF_ALT'.'''
    
    if category == 'ID':
        a = 'ID'
        b = 'CHR_POS_REF_ALT'
    
    elif category == 'CHR_POS_REF_ALT':
        a = 'CHR_POS_REF_ALT'
        b = 'ID'
    
    df_ = df.groupby([a, b]).mean()['MAF'].unstack()
    
    df_ = df_**0.5
    
    sns.clustermap(df_.fillna(0), method=method)
    
    plt.show()
    
def make_clustermap_sub(df, area = 'tl', category='CHR_POS_REF_ALT',
                        method='weighted', cmap='Blues'):
    '''use pandas to generate a table for clustering. category is 'ID' or
    it is 'LOCATION'. use tl, pfl, or pfr to determine the region.'''
    
    if area == 'tl':
        df = df[(df.BRAIN_REGION == 'T') & (df.L_R == 'L') &
                (df.IN_T_L == True) & (df.LOCATION != 'LRG') &
                (df.EXPERIMENT.isin(['bulk', 'subsample']))]
    
    elif area == 'pfl':
        df = df[(df.BRAIN_REGION == 'PF') & (df.L_R == 'L') &
                (df.IN_PF_L == True) & (df.LOCATION != 'LRG') &
                (df.EXPERIMENT.isin(['bulk', 'subsample']))]
        
    elif area == 'pfr':
        df = df[(df.BRAIN_REGION == 'PF') & (df.L_R == 'R') &
                (df.IN_PF_R == True) & (df.LOCATION != 'LRG') &
                (df.EXPERIMENT.isin(['bulk', 'subsample']))]
    
    else:
        return 'NO SUBSAMPLING AVAILABLE!'
    
    if category == 'LOCATION':
        a = 'LOCATION'
        b = 'CHR_POS_REF_ALT'
    
    elif category == 'CHR_POS_REF_ALT':
        a = 'CHR_POS_REF_ALT'
        b = 'LOCATION'
    
    df_ = df.groupby([a, b]).mean()['CONSIDERED_MAF_TISSUE'].unstack()
    
    df_ = df_**0.5
    
    sns.clustermap(df_.fillna(0), method=method, cmap=cmap)
    
    plt.show()
    
#------------------------------------------------------------------------------
#mut signatures

def mutsig_overview(df):
    
    df = df[(df.SEED == True) & (df.INDEL == False)]
    
    all_mos = df.copy()
    across_body = df[(df.SET_IN_BRAIN) & ~(df.SET_BRAIN_ONLY)]
    brain_only = df[df.SET_BRAIN_ONLY]
    cortex_only = df[df.SET_CORTEX_ONLY]
    brain_side_restricted = df[((df.SET_LEFT_ONLY) | (df.SET_RIGHT_ONLY)) &
                                (df.SET_CORTEX_ONLY)]
    brain_one_tissue_only = df[(df.SET_ONE_TISSUE) & (df.SET_BRAIN_ONLY)]
    
    all_mos['PLOT_CLASS'] = '6_all'
    across_body['PLOT_CLASS'] = '5_across'
    brain_only['PLOT_CLASS'] = '4_brain'
    cortex_only['PLOT_CLASS'] = '3_cortex'
    brain_side_restricted['PLOT_CLASS'] = '2_lat'
    brain_one_tissue_only['PLOT_CLASS'] = '1_single'
    
    df = pd.concat([all_mos, across_body, brain_only, cortex_only,
                    brain_side_restricted, brain_one_tissue_only])
    
    df0 = make_0_cat_indplot(df)
    
    df = df.groupby(['PLOT_CLASS', 'CAT']).count()
    df.reset_index(inplace=True)
    df = df.append(df0, ignore_index=True, sort=False)
    df = df.groupby(['PLOT_CLASS', 'CAT']).sum()
    df.reset_index(inplace=True)
    
    df['SUM'] = df.groupby(['PLOT_CLASS'])['POS'].transform('sum')
    df['REL_COUNT'] = df.POS / df.SUM
    df = df.reset_index()
    
    df.sort_values(by=['PLOT_CLASS', 'CAT'], ascending=False, inplace=True)
    
    pos = [0, 1, 2, 3, 4, 5]
    width = 1
    names = ['All Variants', 'Brain and Organs', 'Brain Only', 'Cortex Only',
             'Brain Lateralized', 'Brain Single Sample']
    
    bars_lst = []
    bottoms_lst = []
    bottoms = [0, 0, 0, 0, 0, 0]
    
    for cat in df.CAT.unique():
        
        bottoms_lst.append(bottoms)
        bars_lst.append(df[df.CAT == cat].REL_COUNT)
        
        bottoms = np.add(bottoms, df[df.CAT == cat].REL_COUNT).tolist()
        
    colors = ['xkcd:light blue', 'xkcd:black', 'xkcd:bright red',
              'xkcd:light gray', 'xkcd:slime green', 'xkcd:soft pink'][::-1]
    
    for i in range(6):
        
        plt.bar(pos, bars_lst[i], bottom=bottoms_lst[i], color=colors[i],
                edgecolor='white', width=width)
    
    plt.xlabel('Labels')
    plt.ylabel('Relative Contribution')
    
    sns.despine(bottom=True, offset=5, trim=True)
    plt.xticks(pos, names, rotation=45)
    plt.show()
    
def make_0_cat_indplot(df):
    '''makes all 6 cats for Individual/Age etc. after plot_cat is added'''
    
    df = df[['PLOT_CLASS']].drop_duplicates()
    
    refs = ['C', 'T']
    nts = ['A', 'T', 'G', 'C']
    
    cats = [ref + '>' + alt for ref in refs for alt in nts]
    cats.remove('T>T')
    cats.remove('C>C')
    
    columns = []
    
    for cat in cats:
        for pltclass in df.PLOT_CLASS:
            columns.append((pltclass, cat))
    
    df = pd.DataFrame(data=columns, columns=['PLOT_CLASS', 'CAT'])
    df['POS'] = 0
    
    return df

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#sorted pop plots

#---
def plot_sort_overview(df, POS):
    '''plot all sorted populations except for TBR1 in an overview plot.'''
    
    markers = ['NeuN', 'Olig2', 'Lhx2', 'PU1']
    colors = ['xkcd:brick red', 'xkcd:orange', 'xkcd:dark blue', 'xkcd:teal']
    
    df = df[(df.EXPERIMENT.isin(markers)) & (df.POS == POS)]
    sorter(df)
    
    df['MAF_'] = (df.CONSIDERED_MAF_TISSUE * df.SET_SORT_INCLUDED)**0.5
    
    maxi, semi, maxi_l, semi_l = get_ticks(df)
    
    
    fig, axs = plt.subplots(ncols=8, sharey=False)
    
    for i, lr in enumerate(['L', 'R']):
        for j, marker, color in zip(range(4), markers, colors):
            k = i * 4 + j
            
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
    axs[3].set_title(df.CHR_POS_REF_ALT.unique()[0])
                
    plt.show()

def plot_sort_overview_onax(df, POS, axs):
    '''as above, but with axes passed.'''
    
    markers = ['NeuN', 'Olig2', 'Lhx2', 'PU1']
    colors = ['xkcd:brick red', 'xkcd:orange', 'xkcd:dark blue', 'xkcd:teal']
    
    df = df[(df.EXPERIMENT.isin(markers)) & (df.POS == POS)]
    sorter(df)
    
    df['MAF_'] = (df.CONSIDERED_MAF_TISSUE * df.SET_SORT_INCLUDED)**0.5
    
    maxi, semi, maxi_l, semi_l = get_ticks(df)
    
    
    for i, lr in enumerate(['L', 'R']):
        for j, marker, color in zip(range(4), markers, colors):
            k = i * 4 + j
            
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
    
    axs[0].set_yticks([0., semi**0.5, maxi**0.5])
    axs[0].set_yticklabels(['0.0', semi_l, maxi_l])
    sns.despine(bottom=True, trim=True, ax=axs[0])
    axs[0].set_xticks([])
                
    return axs

def make_colors_markers(df, color, depth=100):
    
    colors = []
    
    for index, row in df.iterrows():
        if row.SET_SORT_INCLUDED == 0:
            colors.append('0.4')
        elif row.DEPTH < depth:
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
#more sort plots under combination and digital plot
        
######Combination Plot for Phage and all three Regions and Sorted######
    
def combination_sort_plot(df, POS):
    '''use the four onax plots to obtain a full view of a single variant.'''
    
    ax1 = plt.subplot2grid((2, 3), (0, 0))
    ax2 = plt.subplot2grid((2, 3), (1, 0))
    ax3 = plt.subplot2grid((2, 3), (0, 1), rowspan=2)
    ax4 = plt.subplot2grid((2, 3), (0, 2))
    ax5 = plt.subplot2grid((2, 24), (1, 16))
    ax6 = plt.subplot2grid((2, 24), (1, 17))
    ax7 = plt.subplot2grid((2, 24), (1, 18))
    ax8 = plt.subplot2grid((2, 24), (1, 19))
    ax9 = plt.subplot2grid((2, 24), (1, 20))
    ax10 = plt.subplot2grid((2, 24), (1, 21))
    ax11 = plt.subplot2grid((2, 24), (1, 22))
    ax12 = plt.subplot2grid((2, 24), (1, 23))
    
    make_subplot_variant_onax(df, POS, ax1, 'L', 'PF')
    make_subplot_variant_onax(df, POS, ax2, 'L', 'T')
    make_AF_plot_variant_onax(df[df.EXPERIMENT == 'bulk'], POS, ax3)
    make_subplot_variant_onax(df, POS, ax4, 'R', 'PF')
    plot_sort_overview_onax(df, POS, [ax5, ax6, ax7, ax8, ax9, ax10, ax11, 
                                      ax12])
    ax3.set_xlim(-1, 21.5)
    ax3.set_ylim(2.5, 46)
    ax3.tick_params(bottom=False, left=False,
                    labelbottom=False, labelleft=False)
    
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_title('')
    
    ax1.set_title('L-PF')
    ax2.set_title('L-T')
    ax3.set_title(df[df.POS == POS].CHR_POS_REF_ALT.unique()[0])
    ax4.set_title('R-PF')
    ax8.set_title('AF (sqrt-t) Sorted Populations')
    
    plt.show()
    
    
def combination_sort_plot_pdf(df, path):
    '''use the four onax plots to obtain a full view of a single variant.'''
    
    pdf_pages = PdfPages(path)
    
    df = df.fillna(0)
    
    for POS in df.POS.unique():
        
        fig = plt.figure(figsize=(16, 9), dpi=100)
    
        ax1 = plt.subplot2grid((2, 3), (0, 0))
        ax2 = plt.subplot2grid((2, 3), (1, 0))
        ax3 = plt.subplot2grid((2, 3), (0, 1), rowspan=2)
        ax4 = plt.subplot2grid((2, 3), (0, 2))
        ax5 = plt.subplot2grid((2, 24), (1, 16))
        ax6 = plt.subplot2grid((2, 24), (1, 17))
        ax7 = plt.subplot2grid((2, 24), (1, 18))
        ax8 = plt.subplot2grid((2, 24), (1, 19))
        ax9 = plt.subplot2grid((2, 24), (1, 20))
        ax10 = plt.subplot2grid((2, 24), (1, 21))
        ax11 = plt.subplot2grid((2, 24), (1, 22))
        ax12 = plt.subplot2grid((2, 24), (1, 23))
        
        make_subplot_variant_onax(df, POS, ax1, 'L', 'PF')
        make_subplot_variant_onax(df, POS, ax2, 'L', 'T')
        make_AF_plot_variant_onax(df[df.EXPERIMENT == 'bulk'], POS, ax3)
        make_subplot_variant_onax(df, POS, ax4, 'R', 'PF')
        plot_sort_overview_onax(df, POS, [ax5, ax6, ax7, ax8, ax9, ax10, ax11, 
                                          ax12])
        
        ax3.set_xlim(-1, 21.5)
        ax3.set_ylim(2.5, 46)
        ax3.tick_params(bottom=False, left=False,
                        labelbottom=False, labelleft=False)
    
        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_title('')
    
        ax1.set_title('L-PF')
        ax2.set_title('L-T')
        ax3.set_title(df[df.POS == POS].CHR_POS_REF_ALT.unique()[0])
        ax4.set_title('R-PF')
        ax8.set_title('AF (sqrt-t) Sorted Populations')
        
        pdf_pages.savefig(fig)
        plt.close()
    
    pdf_pages.close()
    

#------------------------------------------------------------------------------
#digital presence/absence plot

def make_abspres_plot(df, POS, brain_only=False):
    '''just use the IN_25 to make a presence/absence plot for the entire phage
    or just the brain regions dependin on flag brain_only.'''
    
    #initialize df
    df = df[(df.EXPERIMENT == 'bulk') & (df.POS == POS)]
    sorter(df)
    
    #make color list
    colors = ['xkcd:pale pink' if in25 == 1 else '0.9' for in25 in df.IN_25]
    
    #make boxes
    fig,ax = plt.subplots(1)

    boxes = make_boxes_abspres(brain_only)
        
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
    

def make_boxes_abspres(brain_only):
    
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
             Rectangle((5,12.5), 10.5, 5),
             Rectangle((5, 9.5), 10.5, 2.5),
             Rectangle((5, 6.5), 10.5, 2.5),
             Rectangle((0, 3.5), 10, 2.5), Rectangle((10.5, 3.5), 10, 2.5)]
    
    if brain_only == True:
        return boxes[:21]
    else:
        return boxes
#------------------------------------------------------------------------------
#more sort plots
        
def correlate_markers(df, m1, m2, depth=1000, include0=True, L_R=False):
    '''markers: NeuN, TBR1, Olig2, Lhx2, PU1'''
    
    df = df[((df.EXPERIMENT == m1) | (df.EXPERIMENT == m2)) &
            (df.SET_SORT_INCLUDED == True) & (df.DEPTH >= depth)]

    sorter(df)
    
    df['SUM_PRESENCE'] = df.groupby(['POS', 'L_R', 'BRAIN_REGION'])\
                            ['IN_TISSUE'].transform('sum')
    
    if include0 == False:
        df = df[df.SUM_PRESENCE == 2]
    else:
        df = df[df.SUM_PRESENCE != 0]
        df['SUM_COUNT'] = df.groupby(['POS', 'L_R', 'BRAIN_REGION'])\
                            ['IN_TISSUE'].transform('count')
        df= df[df.SUM_COUNT == 2]
    
    
    x = df[df.EXPERIMENT == m1].reset_index().CONSIDERED_MAF_TISSUE
    y = df[df.EXPERIMENT == m2].reset_index().CONSIDERED_MAF_TISSUE
    
    plot_df = df[df.EXPERIMENT == m1].reset_index()
    plot_df['X'] = x**0.5
    plot_df['Y'] = y**0.5

    sns.regplot(x='X', y='Y', data=plot_df,
               scatter_kws={'edgecolor': 'k', 'color': '0.8', 'alpha': 0.5},
               line_kws={'color': 'k'})
    
    if L_R == True:
        sns.regplot(x='X', y='Y', data=plot_df[plot_df.L_R == 'L'],
                    scatter_kws={'edgecolor': 'xkcd:orange',
                                 'color': 'xkcd:orange', 'alpha': 0.2},
                    line_kws={'color': 'xkcd:orange', 'alpha': 0.2,
                              'linestyle': '--'},
                    ci=None)
        sns.regplot(x='X', y='Y', data=plot_df[plot_df.L_R == 'R'],
                    scatter_kws={'edgecolor': 'b', 'color': 'b', 'alpha': 0.2},
                    line_kws={'color': 'b', 'alpha': 0.2, 'linestyle': '--'},
                    ci=None)
    
    plt.xlim(-0.015, 1.)
    plt.ylim(-0.015, 1.)
    
    plt.xticks(ticks=[0., 0.2**0.5, 1.0], labels=['0.0', '0.2', '1.0'])
    plt.yticks(ticks=[0., 0.2**0.5, 1.0], labels=['0.0', '0.2', '1.0'])
    
    name = {'Olig2': 'OLIG2', 'Lhx2': 'LHX2', 'PU1': 'PU.1'}
    
    plt.xlabel('AF (sqrt-t) ' + name.get(m1, m1))
    plt.ylabel('AF (sqrt-t) ' + name.get(m2, m2))
    
    sns.despine(trim=True, offset=5)
    
    ax = plt.gca()
    ax.set(adjustable='box', aspect='equal')
    
    plt.show()
    
    return(sp.stats.spearmanr(plot_df.X, plot_df.Y))

#---
def volcano_left_right(df, depth=1000):
    '''get average af, normalized delta, and p value to plot volcano.'''
    
    df = df[(df.EXPERIMENT.isin(['NeuN', 'Olig2', 'Lhx2'])) &
            (df.SET_SORT_INCLUDED == True) & (df.DEPTH >= depth) &
            (
             ((df.SET_IN_CORTEX) & (df.SET_CORTEX_ONLY == False)) |
             ((df.SET_IN_CORTEX) & (df.SET_LEFT_ONLY == False) &
                                   (df.SET_RIGHT_ONLY == False))
            )]
    df['SUM_SORTS_POS'] = df.groupby('CHR_POS_REF_ALT')['IN_TISSUE']\
                            .transform('sum')
    df = df[df.SUM_SORTS_POS > 0].reset_index()
    
    df['MAX_SORTED_CONSIDERED'] = df.groupby('CHR_POS_REF_ALT')\
                                     ['CONSIDERED_MAF_TISSUE'].transform('max')
    df['LEFT_AVG_CMAF'] = l_r_acmaf(df, l_r='L')
    df['RIGHT_AVG_CMAF'] = l_r_acmaf(df, l_r='R')
    df['NORM_DELTA'] = df.apply(norm_delta, axis=1)
    
    df['NEG_LOG_P'] = neg_log_p(df)
    
    plot_df = df[~(df.CHR_POS_REF_ALT.duplicated())]
    
    #replace values above 10 with 10
    plot_df['NEG_LOG_P_adj'] = plot_df.apply(replace_above_10, axis=1)
    #make markers for values above 10 with triangle
    ecs = ['r' if p > 10 else 'k' for p in plot_df.NEG_LOG_P]
    #make colors for markers abov p=0.05 and below/above -0.5/0.5
    colors = make_colors_volcano(plot_df)
    #make sizes for different af bins
    sizes = make_sizes_volcano(plot_df)
    
    plt.scatter(plot_df.NORM_DELTA, plot_df.NEG_LOG_P_adj, c=colors, s=sizes,
                edgecolors=ecs, alpha=0.7, zorder=100)
    plt.hlines(y=-sp.log10(0.05), xmin=-1, xmax=1, linestyles='--', color='r',
               alpha=0.2)
    plt.vlines(x=[-0.5, 0.5], ymin=0, ymax=10, linestyles='--', color='r',
               alpha=0.2)
    plt.vlines(x=0, ymin=0, ymax=10, linestyles='--', color='k')
    
    plt.xlim(-1.035, 1.035)
    plt.ylim(-.35, 10.35)
    plt.xlabel('Normalized d (R-L)')
    plt.ylabel('-log(P)')
    plt.xticks(ticks=[-1.0, -0.5, 0., 0.5, 1.0])
    
    sns.despine(trim=True, offset=5)
    
    plt.show()
    return plot_df
    

def l_r_acmaf(df, l_r):
    
    df_lat = df[df.L_R == l_r]
    df_grp = df_lat.groupby('CHR_POS_REF_ALT').mean().reset_index()
    
    df_tomerge = df_grp[['CHR_POS_REF_ALT', 'CONSIDERED_MAF_TISSUE']]
    df_tomerge.rename({'CONSIDERED_MAF_TISSUE': 'OUTPUT'}, axis=1,
                      inplace=True)
    
    df_out = pd.merge(df, df_tomerge, how='left')
    
    return df_out.OUTPUT
    
def norm_delta(row):
    
    l = row.LEFT_AVG_CMAF
    r = row.RIGHT_AVG_CMAF
    maxi = max([l, r])
    
    return (r-l)/maxi
    
def neg_log_p(df):
    '''bonferroni corrected!'''
    
    lst = []
    
    for chrposrefalt in df.CHR_POS_REF_ALT.unique():
        
        var_df = df[df.CHR_POS_REF_ALT == chrposrefalt]
        anova = var_df.anova(dv='CONSIDERED_MAF_TISSUE',
                             between=['L_R', 'EXPERIMENT'])
        if 'p-unc' in anova.columns:
            pval = anova.loc[0, 'p-unc']
        else:
            pval = 1
        
        pval = pval*len(df.CHR_POS_REF_ALT.unique())
        
        if pval > 1:
            pval = 1
        
        lst.append((chrposrefalt, -sp.log10(pval)))
    
    pvals = pd.DataFrame(lst, columns=['CHR_POS_REF_ALT', 'OUTPUT'])
    
    df_out = pd.merge(df, pvals, how='left')
    
    return df_out.OUTPUT
    
def replace_above_10(row):
    
    if row.NEG_LOG_P > 10:
        return 10
    else:
        return row.NEG_LOG_P
    
def make_colors_volcano(df):
    
    lst = []
    
    for p, d in zip(df.NEG_LOG_P, df.NORM_DELTA):
        
        if p > -sp.log10(0.05):
            if d < -0.5:
                lst.append('b')
            elif d > 0.5:
                lst.append('xkcd:orange')
            else:
                lst.append('0.8')
        else:
            lst.append('0.8')
    
    return lst
    
def make_sizes_volcano(df):
    
    lst = []
    
    for maf in df.MAX_SORTED_CONSIDERED:
        if maf < 0.02:
            lst.append(25)
        elif maf < 0.05:
            lst.append(50)
        elif maf < 0.1:
            lst.append(100)
        else:
            lst.append(200)
    
    return lst
#---


def loliplot_NeuN_TBR1(df, POS):
    
    '''plot NeuN and TBR1 in a small loliplot'''
    markers = ['NeuN', 'TBR1']
    colors = ['xkcd:brick red', 'xkcd:baby pink']
    
    df = df[(df.EXPERIMENT.isin(markers)) & (df.POS == POS) &
            (((df.L_R == 'L') & (df.BRAIN_REGION == 'T')) |
             ((df.L_R == 'R') & (df.BRAIN_REGION == 'PF')))]
    sorter(df)
    
    df['MAF_'] = (df.CONSIDERED_MAF_TISSUE * df.SET_SORT_INCLUDED)**0.5
    
    maxi, semi, maxi_l, semi_l = get_ticks(df)
    
    for i in range(2):
    
        plt.vlines([0 + i/4, 1 + i/4], ymin=[0, 0],
                   ymax=df[df.EXPERIMENT == markers[i]].MAF_)
        plt.scatter([0 + i/4, 1 + i/4], df[df.EXPERIMENT == markers[i]].MAF_,
                    color=colors[i], s=100, edgecolors='k', zorder=100)
    
    plt.vlines(x=0.625, ymin=0, ymax=maxi**0.5, linestyles='--')
    
    plt.xlim(-0.1, 1.35)
    plt.ylim(-((maxi/100)**0.5), maxi**0.5)
    
    plt.xticks(ticks=[0.125, 1.125], labels=['L-T', 'R-PF'])
    plt.yticks(ticks=[0., semi**0.5, maxi**0.5],
               labels=['0.0', semi_l, maxi_l])
    
    plt.xlabel('')
    plt.ylabel('AF (sqrt-t)')
    plt.title(df.CHR_POS_REF_ALT.unique()[0])
    
    sns.despine(bottom=True, trim=True, offset=5)
    
    plt.show()
    
    
    
    
    
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
    axs[3].set_title(df.CHR_POS_REF_ALT.unique()[0])
                
    plt.show()

#------------------------------------------------------------------------------
#snMPAS

def make_clustermap_sn(df, TLP=True, n_cutoff=21, category='CHR_POS_REF_ALT',
                       metric='sqeuclidean', method='weighted', no01=True,
                       sn_cls = '_05'):
    '''clustermap for single nuclei. has to be sn_cls = '' to work on the old
    annotations.'''
    
    if TLP == True:
        cp = 'NUMBER_CP_TLP' + sn_cls
        present = 'SET_CP_TL_PRESENT' + sn_cls
        number_cell = 'NUMBER_CP_TLP_CELL' + sn_cls
    else:
        cp = 'NUMBER_CP' + sn_cls
        present = 'SET_CELL_POSITIVE' + sn_cls
        number_cell = 'NUMBER_CP_CELL' + sn_cls
        
    df = df[(df[cp] > 0) & (df[cp] < n_cutoff) & (df.TISSUE_CELLS == 'cells')]
    
    if no01 == True:
        df = df[(df[cp] > 1) & (df[number_cell] > 0)]
    
    if category == 'EXPERIMENT':
        a = 'EXPERIMENT'
        b = 'CHR_POS_REF_ALT'
    
    elif category == 'CHR_POS_REF_ALT':
        a = 'CHR_POS_REF_ALT'
        b = 'EXPERIMENT'
    
    df_ = df.groupby([a, b]).mean()[present].unstack()
    
    sns.clustermap(df_.fillna(0), method=method, metric=metric, cmap='Blues')
    
    plt.show()

#---
def make_sorted_sn_plot(df, TLP=True, n_cutoff=21, no01=True, sn_cls='_05',
                        x=40, y=80, cell_type = True):
    '''manually sort variants for plotting based on abundance.'''
    
    if TLP == True:
        cp = 'NUMBER_CP_TLP' + sn_cls
        present = 'SET_CP_TL_PRESENT' + sn_cls
        number_cell = 'NUMBER_CP_TLP_CELL' + sn_cls
    else:
        cp = 'NUMBER_CP' + sn_cls
        present = 'SET_CELL_POSITIVE' + sn_cls
        number_cell = 'NUMBER_CP_CELL' + sn_cls
        
    df = df[(df[cp] > 0) & (df[cp] < n_cutoff) & (df.TISSUE_CELLS == 'cells')]
    
    df = df[['CHR_POS_REF_ALT', 'EXPERIMENT', cp, present, number_cell]]
    
    if no01 == True:
        df = df[(df[cp] > 1) & (df[number_cell] > 0)]
    
    df.sort_values(by=[cp, 'CHR_POS_REF_ALT'], ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.reset_index(inplace=True)
    
    variants = pd.Series(df.CHR_POS_REF_ALT.unique()).reset_index()
    variants.rename({'index':'variant_order', 0: 'CHR_POS_REF_ALT'}, axis=1,
                    inplace=True)
    
    df = pd.merge(df, variants, how='left')
    
    sorter = ['variant_order']
    
    for vo in df.variant_order.unique():
        vostr = str(vo)
        
        df_temp = df[df.variant_order == vo][['EXPERIMENT', present]]
        df_temp[present] = (df_temp[present] == False)
        df_temp.rename({present: vostr}, axis=1, inplace=True)
        df = pd.merge(df, df_temp, how='left')
        
        sorter.append(vostr)
    
    df.sort_values(by=sorter, inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.reset_index(inplace=True)
    
    cells = df[['EXPERIMENT', 'level_0']][~(df.EXPERIMENT.duplicated())]
    cells.rename({'level_0':'cell_order'}, axis=1,
                 inplace=True)
    
    df = pd.merge(df, cells, how='left')
    
    df_ = df[df[present] == True]
    
    if cell_type == False:
        
        colors = 'k'
    
    else:
        
        colors = make_colors_cell_type(df_)
    
    #edgecolors = ['r' if 'NeuN' in exp else 'k' for exp in df.EXPERIMENT]
    
    plt.scatter(df_.variant_order, df_.cell_order, color=colors, marker='s',
                #edgecolors=edgecolors
                )
    
    plt.xlim(-1, x)
    plt.ylim(-1, y)
    
    plt.xlabel('Ranked Variants')
    plt.ylabel('Cells')
    sns.despine(offset=5, trim=True)
    
    plt.show()
    return (df.EXPERIMENT.unique(), df.CHR_POS_REF_ALT.unique())

def make_colors_cell_type(df):
    
    colors = []
    
    for index, row in df.iterrows():
        
        if 'NeuN' in row.EXPERIMENT:
            colors.append((1,0,0,1))
        
        elif 'DAPI' in row.EXPERIMENT:
            colors.append((0,0,0,1))
        
        else:
            colors.append((0,0,1,1))
    
    return colors
#---

def make_bar_sn_plot(df):
    '''note the fake 0.5 are b/c of Illustrator import error for 0 values.'''
    
    df = df[(df.ID.str.contains('T_L_Sml')) & (df.EXPERIMENT == 'bulk') &
            (df.CLADE_FIG != 'ZONK') & (df.SUB_CLADE_FIG == 'Founder')]\
            [['CONSIDERED_MAF_TISSUE', 'VARIANT_ORDER_FIG', 'CLADE_FIG',
              'SUB_CLADE_FIG']].sort_values('VARIANT_ORDER_FIG')
    
    df_sum = pd.DataFrame(data={'CONSIDERED_MAF_TISSUE':
                                [df.CONSIDERED_MAF_TISSUE.sum(), 0.5, 0.5, 0.5,
                                 0.5, 0.5, 0.5, 0.5, 0.5],
                                'VARIANT_ORDER_FIG': [-1 for i in range(9)],
                                'CLADE_FIG': ['Sum', 'I', 'II', 'III', 'IV',
                                              'V', 'VI', 'VII', 'Sum'],
                                'SUB_CLADE_FIG': ['A', 'A', 'A', 'A', 'A', 'A',
                                                  'A', 'A', 'Founder']})
        
    df = df.append(df_sum).sort_values('VARIANT_ORDER_FIG')
    
    df['CONSIDERED_MAF_TISSUE'] = df.CONSIDERED_MAF_TISSUE**0.5
    
    sns.catplot(x='SUB_CLADE_FIG', y='CONSIDERED_MAF_TISSUE', hue='CLADE_FIG',
                data=df, kind='bar',
                palette=sns.color_palette(['0.4', 'b', 'g', 'r', 'y',
                                           'xkcd:purple', 'xkcd:cyan', '0.8']),
                edgecolor='k')
    
    plt.ylim(0, 0.5)
    plt.yticks(ticks=[0, 0.01**0.5, 0.05**0.5, 0.2**0.5, 0.5**0.5],
               labels=['0.00', '0.01', '0.05', '0.20', '0.50'])
    plt.ylabel('L-T-Sml AF (sqrt-t)')
    
    plt.show()


def correlation_variants(df, POS1, POS2, only_25=False, xy_plus=0.05):
    '''given to variants by position, draw correlation across all tissues where
    both are present.'''
    
    df = df[((df.POS == POS1) | (df.POS == POS2)) &
            (df.TISSUE_CELLS == 'tissues') & (df.EXPERIMENT != 'TBR1') &
            ((df.EXPERIMENT.isin(['bulk', 'subsample'])) |
             (df.SET_SORT_INCLUDED == True))]
    
    if only_25 == True:
        df = df[df.EXPERIMENT == 'bulk']
            
    sorter(df)
    
    df['SUM_PRESENCE'] = df.groupby('ID')['IN_TISSUE'].transform('sum')
    
    df = df[df.SUM_PRESENCE == 2]
    
    x = df[df.POS == POS1].reset_index().CONSIDERED_MAF_TISSUE
    y = df[df.POS == POS2].reset_index().CONSIDERED_MAF_TISSUE
    
    plot_df = df[df.POS == POS1].reset_index()
    plot_df['X'] = x
    plot_df['Y'] = y

    sns.regplot(x='X', y='Y', data=plot_df,
               scatter_kws={'edgecolor': 'k', 'color': '0.8', 'alpha': 0.5},
               line_kws={'color': 'k'})
    sns.regplot(x='X', y='Y', data=plot_df[plot_df.EXPERIMENT == 'bulk'],
               scatter_kws={'edgecolor': 'r', 'color': 'r', 'alpha': 0.2},
               line_kws={'color': 'r', 'alpha': 0.2, 'linestyle': '--'},
               ci=None)
    x = plt.gca().get_xlim()[1]
    y = plt.gca().get_ylim()[1]
    
    maxi = max([x, y]) + xy_plus
    
    plt.xlim(0, maxi)
    plt.ylim(0, maxi)
    
    plt.xlabel(str(POS1))
    plt.ylabel(str(POS2))
    
    sns.despine(trim=True, offset=5)
    
    ax = plt.gca()
    ax.set(adjustable='box', aspect='equal')
    
    plt.show()
    return(sp.stats.spearmanr(plot_df.X, plot_df.Y))
    

def correlation_sn(df, method='pearson', only_25=False, cluster=False):
    
    df = df[(df.SET_SN_FIGURE == True) & (df.TISSUE_CELLS == 'tissues') &
            (df.EXPERIMENT != 'TBR1') & (df.CLADE_FIG != 'ZONK') &
            ((df.EXPERIMENT.isin(['bulk', 'subsample'])) |
             (df.SET_SORT_INCLUDED == True))]
    
    if only_25 == True:
        df = df[df.EXPERIMENT == 'bulk']
    
    df['VARIANT'] = df.apply(lambda row: '_'.join([row.CLADE_FIG,
                                                   row.SUB_CLADE_FIG,
                                             str(int(row.VARIANT_ORDER_FIG))]),
                             axis=1)
        
    sorter(df)
    
    df.sort_values(['CLADE_FIG', 'SUB_CLADE_FIG', 'VARIANT_ORDER_FIG'],
                   inplace=True)
    df.reset_index(inplace=True)
    
    df['NUMBERS'] = pd.merge(df,
                             df[~(df.VARIANT.duplicated())]\
                             [['VARIANT', 'VARIANT_ORDER_FIG']].reset_index()\
                             .reset_index()\
                             .drop(['index', 'VARIANT_ORDER_FIG'], axis=1),
                             how='left')['level_0']
    
    df['VARIANT'] = df.apply(lambda row: str(row.NUMBERS + 100) + '_' +
                                         row.VARIANT,
                             axis=1)
    
    piv = pd.pivot_table(df, index=['ID'], values=['CONSIDERED_MAF_TISSUE'],
                             columns=['VARIANT'])
    piv.columns = piv.columns.droplevel()
    
    corr = piv.corr(method=method)
    
    if cluster == False:
        sns.heatmap(corr, vmin=-1, vmax=1, center=0, cmap='RdBu_r',
                    xticklabels=1, yticklabels=1)
    
    else:
        sns.clustermap(corr, vmin=-1, vmax=1, center=0, cmap='RdBu_r',
                       xticklabels=1, yticklabels=1)
    
    ax = plt.gca()
    ax.set(adjustable='box', aspect='equal')
    
    plt.show()
    
    
    


