# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 04:33:04 2019
​
@author: Martin W. Breuss
"""
​
#Import modules
​
import sys, os, gzip
import pandas as pd
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
​
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle, Polygon
​
plt.rcParams['svg.fonttype'] = 'none'
plt.ioff()
​
#Import modules done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Functions
​
#First level - Mosaic Output
#import------------------------------------------------------------------------
​
def import_table(path_to_data):
    '''header of the input table'''
    
    names = ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'ANNO', 'GENE',
             'GNOMAD_FREQ', 'REPEAT_MASKER', 'SEGDUP', 'HOMOPOLYMER',
             'REF_SEQ', 'DINUCLEOTIDE', 'NEAR_INDEL', 'UCSC_RPMSK',
             'REF_COUNT', 'ALT_COUNT', 'MAF', 'LOWER_CI', 'UPPER_CI',
             'CI_IS_GREATER', 'NORMAL_REF_COUNT', 'NORMAL_ALT_COUNT',
             'NORMAL_MAF', 'NORMAL_LOWER_CI', 'NORMAL_UPPER_CI',
             'IN_MOSAICHUNTER', 'IN_MOSAICFORECAST', 'IN_STRELKA2',
             'IN_MUTECT2']
    
    df = pd.read_csv(path_to_data, sep='\t', names=names)
    
    return df
​
​
#annotation--------------------------------------------------------------------
​
def define_mosaic(df):
    '''uses gnomAD AF, repeat info, segdup info, vicinity to homopolymer/di-
    nucleotide repeat/indel, and a variant's presence in both, Strelka2 and
    MuTect2.'''
    
    df['SET_MOSAIC'] = (
                       (df.GNOMAD_FREQ < 0.01) & #redundant, but left in
                       (df.UCSC_RPMSK == 'pass') &
                       (df.REPEAT_MASKER == 0) &
                       (df.SEGDUP == False) &
                       (df.HOMOPOLYMER == False) &
                       (df.DINUCLEOTIDE == False) &
                       (df.NEAR_INDEL == False) &
                       (
                       ((df.IN_STRELKA2 == True) & (df.IN_MUTECT2 == True)) |
                       (df.IN_MOSAICHUNTER == True) | 
                       ((df.IN_MOSAICFORECAST == True) &
                         (df.GNOMAD_FREQ == 0))
                       )
                       )
​
​
def add_indel_info(df):
    '''add info whether a variant is an indel or a long indel (more than 1bp
    difference).'''
    
    df['INDEL'] = (df.REF.str.len() != df.ALT.str.len())
    df['LONG_INDEL'] = ((df.REF.str.len() > 2) | (df.ALT.str.len() > 2))
​
​
#mosaic subset-----------------------------------------------------------------
    
def add_recurrence_info(df):
    '''adds the information whether or not a variant is recurrent.'''
    
    df['DUP'] = df.duplicated('POS', keep=False)
    df['DUP_2nd'] = df.duplicated('POS', keep='first')
​
#Second level - Variant check in all tissues
#import------------------------------------------------------------------------
​
#import can be done with regular read_csv
​
def annotate_order(df):
    '''use a dictionary to annotate all samples in specified order.'''
    
    tissues = ('7614-B-L-PF-2-Lrg', '7614-B-L-PF-2-Sml', '7614-B-L-F-3-Lrg',
               '7614-B-L-F-3-Sml', '7614-B-L-P-4-Lrg', '7614-B-L-P-4-Sml', 
               '7614-B-L-O-5-Lrg', '7614-B-L-O-5-Sml', '7614-B-L-T-1-Lrg',
               '7614-B-L-T-1-Sml', '7614-B-R-PF-2-Lrg', '7614-B-R-PF-2-Sml',
               '7614-B-R-F-3-Lrg', '7614-B-R-F-3-Sml', '7614-B-R-P-4-Lrg',
               '7614-B-R-P-4-Sml', '7614-B-R-O-5-Lrg', '7614-B-R-O-5-Sml',
               '7614-B-R-T-1-Lrg', '7614-B-R-T-1-Sml', '7614-C', '7614-H',
               '7614-L', '7614-K-L', '7614-K-R')
    
    dic = {t:s for t, s in zip(tissues, range(25))}
    
    df['SORTER'] = df.SAMPLE.map(dic)
​
​
def annotate_het_noise_flag(df):
    '''annotate variants that may be noise (flag on tissue and variant basis)
    and those that are het.'''
    
    df['FLAG_HET'] = (df.UPPER > 0.45)
    df['FLAG_NOISE'] = (df.LOWER < 0.001)
    
    df['FLAG_HET_COMBINED'] = df.groupby('CHR_POS_REF_ALT')['FLAG_HET']\
                                .transform('mean')
    
​
def get_mosaics_only(df):
    '''just a reminder for filter used.'''
    
    return df[df.FLAG_HET_COMBINED < 0.5]
​
​
def make_uniqueID_remove_duplicates(df):
    '''only needed if the organs are added not exactly once. should not be
    needed for all versions of output.'''
    
    df['UNIQUE_ID'] = df.apply(lambda row: row['CHR_POS_REF_ALT'] +
                                           row['SAMPLE'], axis=1)
    
    df = df[~df.duplicated('UNIQUE_ID', keep='first')]
    
    return df
​
#Third level - Table returned with annotation for MF etc. and Seq
#Additional information and prepare for analysis
#------------------------------------------------------------------------------
    
def mf_calls_only(df):
    '''only consider variants that were called by mosaicforecast (~80-90%
    validation reported; approximately true for first level validation effort
    as well).'''
    
    return df[df.IN_MOSAICFORECAST == True]
​
def add_seed_info(df):
    '''add a recurrence flag, so analysis on the variant level that do not
    need organ information can be done.'''
    
    df['SEED'] = ~df.duplicated('CHR_POS_REF_ALT', keep='first')
​
​
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#all functions to annotate mutation categories for plotting
def add_mut_cats(df):
    '''extract the mutational category from the sequence snippet.'''
    
    df['TRI_REF'] = df.REF_SEQ.str[7:10]
    
    df['TRI_CAT'] = df.apply(lambda row: categorize_mut(row)[0], axis=1)
    df['CAT'] = df.apply(lambda row: categorize_mut(row)[1], axis=1)
​
​
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
​
    
def rev_comp(seq):
    '''make reverse complement of sequence. simple, assumes that all are upper
    case and regular nucleotides. left without failsaves to make sure that the
    data passed on is as expected.'''
    reverse_comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    rev_seq = ''
    
    for n in seq[::-1]:
        rev_seq += reverse_comp[n]
    
    return rev_seq
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#functions to update the sample info and add sorters
def update_sample(df):
    '''new sample info'''
    
    df['SAMPLE'] = df.apply(get_sample_name, axis=1)
​
def get_sample_name(row):
    '''subfunction for above.'''
    
    sample = '-'.join(['7614', row['ORGAN'], str(row['LEFT_RIGHT']),
                       str(row['BRAIN_REGION']), str(row['BR_NUM']),
                       str(row['LG_SM'])]) #str() needed for nan entries
    
    sample = sample.replace('-nan', '')
    sample = sample.replace('.0', '')
    
    return sample
​
def get_chrom_sorter(df):
    '''make a sortable string from chromosome for sort values.'''
    
    df['CHROM_SORTER'] = df.apply(get_chrom_str, axis=1)
​
def get_chrom_str(row):
    '''use normal string manipulation to get chrN with N bein 01-22 or X.'''
    
    chrom = str(row['CHROM'])
    
    if len(chrom) != 2 and chrom != 'X':
        chrom = '0' + chrom
    
    return ('chr' + chrom)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
def variant_is_present(df):
    '''essentially only the reverse of the noise flag.'''
    
    df['IN_TISSUE'] = df.apply(lambda row: not row['FLAG_NOISE'], axis=1)
    
def annotate_ms_number(df):
    '''add 1-25 based on what will likely be used in ms; i.e. L->R; PF/F/P/O/T;
    C/H/L/K.'''
    
    df['TISSUE_NUMBER_MS'] = df.SORTER + 1
​
def count_positive_tissues(df):
    '''mainly used to remove 0 tissue mosaics, as they are not informative for
    work.'''
    
    df['NUMBER_POSITIVE_TISSUES'] = df.groupby('CHR_POS_REF_ALT')['IN_TISSUE']\
                                        .transform('sum')
​
def remove_0_tissue_vars(df):
    '''s.a.'''
    
    return df[df.NUMBER_POSITIVE_TISSUES != 0]
​
​
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#functions to annotate categories such as lateralized, brain specific etc.
def annotate_expression_pattern(df):
    '''use transform to annotate patterns as bool. note that some of the groups
    require a clean mosaic table, i.e. variants must be determined as positive
    in at least one tissue.'''
    
    grp = df.groupby('CHR_POS_REF_ALT')['IN_TISSUE'] #generates groupby object
    
    df['SET_ONE_TISSUE'] = grp.transform(one_tissue)
    df['SET_ONE_REGION'] = grp.transform(one_region) #cortical regions
    df['SET_CORTEX_ONLY'] = grp.transform(cortex_only)
    df['SET_BRAIN_ONLY'] = grp.transform(brain_only)
    df['SET_KIDNEY_ONLY'] = grp.transform(kidney_only)
    df['SET_LEFT_ONLY'] = grp.transform(left_only)
    df['SET_RIGHT_ONLY'] = grp.transform(right_only)
    df['SET_IN_CORTEX'] = grp.transform(in_cortex)
    df['SET_IN_BRAIN'] = grp.transform(in_brain)
    df['SET_IN_CEREBELLUM'] = grp.transform(in_cerebellum)
    df['SET_IN_HEART'] = grp.transform(in_heart)
    df['SET_IN_LIVER'] = grp.transform(in_liver)
    df['SET_IN_KIDNEY'] = grp.transform(in_kidney)
    df['SET_LEFT_ONLY_BRAIN'] = grp.transform(left_only_br) #can be in organs/cbl
    df['SET_RIGHT_ONLY_BRAIN'] = grp.transform(right_only_br)
​
​
def one_tissue(col):
    
    return sum(col) == 1
​
def one_region(col):
    
    if sum(col[20:]) == 0 and sum(col[:20]) < 3:
        
        region_lst = []
        
        for i in range(0,20,2):
            region_lst.append(bool(sum(col[i:i+2])))
        
        if sum(region_lst) == 1:
            return True
        else:
            return False
    
    else:
        return False
​
def cortex_only(col):
    
    return sum(col[20:]) == 0
​
def brain_only(col):
    
    return sum(col[21:]) == 0
​
def kidney_only(col):
    
    return sum(col[:23]) == 0
​
def left_only(col):
    
    return sum(col[10:23]) == 0 and col.iloc[24] == 0
​
def right_only(col):
    
    return sum(col[:10]) == 0 and sum(col[20:24]) == 0
​
def in_cortex(col):
    
    return sum(col[:20]) > 0
​
def in_brain(col):
    
    return sum(col[:21]) > 0
​
def in_cerebellum(col):
    
    return col.iloc[20] == 1
​
def in_heart(col):
    
    return col.iloc[21] == 1
​
def in_liver(col):
    
    return col.iloc[22] == 1
​
def in_kidney(col):
    
    return sum(col[23:]) > 0
​
def left_only_br(col):
    
    return sum(col[10:20]) == 0 and sum(col[:10]) > 0
​
def right_only_br(col):
    
    return sum(col[:10]) == 0 and sum(col[10:20]) > 0
​
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
​
def get_max_af(df):
    '''get the maximum AF for each variant and the lower/upper coordinates.'''
    
    grp = df.groupby('CHR_POS_REF_ALT') #generates groupby object
    
    df['MAX_MAF'] = grp['MAF'].transform('max')
    df['MAX_MAF_BOOL'] = df.apply(lambda row: row['MAF'] == row['MAX_MAF'],
                                  axis=1)
​
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Functions to define a category label
def define_cat_labels(df):
    '''make several labels that are mutually exclusive. Organ-specific (HLK),
    Organ+brain, Across-brain, lateralized brain, single tissue only.'''
    
    df['CAT_LABEL'] = df.apply(categorize_label, axis=1)
​
def categorize_label(row):
    '''use bools to define label.'''
    
    if row.SET_IN_BRAIN == False: #look at organ only calls
        
        if row.SET_ONE_TISSUE == False:
            
            if row.SET_KIDNEY_ONLY == True:
                return 'both kidneys' #both kidneys still are not one_tissue
            else:
                return 'organs'
        
  ...
