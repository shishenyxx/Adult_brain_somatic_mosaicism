import pandas as pd
import numpy as np
import gzip
import os, sys

HEADER = ["#ID", "CHROM", "POS", "REF", "ALT", "ANNO", "GENE", "GNOMAD_FREQ", \
          "REPEAT_MASKER", "SEGDUP", "HOMOPOLYMER", "REF_SEQ", "DINUCLEOTIDE", \
          "NEAR_INDEL", "UCSC_RPMSK", "REF_COUNT", "ALT_COUNT", "MAF", "LOWER_CI",\
          "UPPER_CI", "CI_IS_GREATER", "NORMAL_REF_COUNT", "NORMAL_ALT_COUNT", \
          "NORMAL_MAF", "NORMAL_LOWER_CI", "NORMAL_UPPER_CI", "NORMAL_CI_IS_GREATER", \
          "TUMOR_IS_BLOOD", "TUMOR_IS_SPERM"]



def summarize_annotations(tumor_id, exon_vf, vf, gnomad_dropped, repeats, \
                          homopolymer, tumor_ci, normal_ci, outfile):
    if ("sperm" in tumor_id) and ("blood" in tumor_id):
       tumor_id = tumor_id.split("_")[0] 
    if "sperm" in tumor_id:
        tumor_is_sperm = 1
        tumor_is_blood = 0
    elif "blood" in tumor_id:
        tumor_is_sperm = 0
        tumor_is_blood = 1
    else:
        tumor_is_sperm = np.nan
        tumor_is_blood = np.nan

    if os.stat(vf).st_size == 0:
        with open(outfile,"w") as f:
            f.write("\t".join(HEADER) + "\n")
        return

    vf_df = pd.read_csv(vf, sep="\t", header=None)
    homopolymer_df = pd.read_csv(homopolymer, sep="\t", header=None)
    tumor_ci_df = pd.read_csv(tumor_ci, sep="\t", header=None)
    normal_ci_df = pd.read_csv(normal_ci, sep="\t", header=None)

    variants = vf_df[[2,3,5,6]].astype(str).apply(lambda x: "_".join(x), axis=1)

    exon_dict = {}
    if os.stat(exon_vf).st_size != 0:
        exon_df = pd.read_csv(exon_vf, sep="\t", header=None)
        exon_variants = exon_df[[3,4,6,7]].astype(str).apply(lambda x: "_".join(x), axis=1)
        exon_dict = dict(zip(exon_variants, zip(exon_df[1], exon_df[2])))

    region = np.copy(vf_df[0])
    gene = np.copy(vf_df[1])
    for i in range(len(variants)):
        variant = variants[i]
        if variant in exon_dict.keys():
            region[i] += ":" + exon_dict[variant][0].replace(" ", "_")
            gene[i] = exon_dict[variant][1]

    gnomad_freqs = np.zeros(len(vf_df))
    gnomad_dict = {}
    if os.stat(gnomad_dropped).st_size != 0:
        dropped_df = pd.read_csv(gnomad_dropped, sep="\t", header=None)
        gnomad_variants = dropped_df[[2,3,5,6]].astype(str).apply(lambda x: "_".join(x), axis=1)
        gnomad_dict = dict(zip(gnomad_variants, dropped_df[1]))

    for i in range(len(gnomad_freqs)):
        if variants[i] in gnomad_dict.keys():
            frequency = gnomad_dict[variants[i]]
            if frequency == ".":
                frequency = -1
            gnomad_freqs[i] = frequency

    repeat = np.zeros(len(vf_df))
    segdup = np.zeros(len(vf_df))
    repeats_dict = {}
    if os.stat(repeats).st_size != 0:
        repeats_df = pd.read_csv(repeats, sep="\t", header=None)
        repeats_df[1] += 1
        repeats_variants = repeats_df[[0,1,3,4]].astype(str).apply(lambda x: "_".join(x), axis=1)
        repeats_dict = dict(zip(repeats_variants, zip(repeats_df[5], repeats_df[6])))
        for i in range(len(vf_df)):
            variant = variants[i]
            repeat[i] = repeats_dict[variant][0]
            segdup[i] = repeats_dict[variant][1]

    chrom = np.copy(vf_df[2])
    pos = np.copy(vf_df[3])
    ref = np.copy(vf_df[5])
    alt = np.copy(vf_df[6])
    ref_seq = np.copy(homopolymer_df[4]).astype(str)
    homopolymer = np.copy(homopolymer_df[5]).astype(int)
    dinucleotide = np.copy(homopolymer_df[6]).astype(int)
    near_indel = np.copy(homopolymer_df[7]).astype(int)
    ucsc_rpmsk = np.copy(homopolymer_df[8]).astype(str)    

    tumor_ref_count = np.copy(tumor_ci_df[4]).astype(str)
    tumor_alt_count = np.copy(tumor_ci_df[5]).astype(str)
    tumor_maf = np.copy(tumor_ci_df[6]).astype(str)
    tumor_lower_ci = np.copy(tumor_ci_df[7]).astype(str)
    tumor_upper_ci = np.copy(tumor_ci_df[8]).astype(str)
    tumor_is_greater = np.copy(tumor_ci_df[9]).astype(str)

    normal_ref_count = np.copy(normal_ci_df[4]).astype(str)
    normal_alt_count = np.copy(normal_ci_df[5]).astype(str)
    normal_maf = np.copy(normal_ci_df[6]).astype(str)
    normal_lower_ci = np.copy(normal_ci_df[7]).astype(str)
    normal_upper_ci = np.copy(normal_ci_df[8]).astype(str)
    normal_is_greater = np.copy(normal_ci_df[9]).astype(str)

    id = np.array([tumor_id] * len(pos))
    is_blood = np.array([tumor_is_blood] * len(pos))
    is_sperm = np.array([tumor_is_sperm] * len(pos))


    final = pd.DataFrame(np.array([
           id, chrom, pos, ref, alt, region, gene, \
           gnomad_freqs, repeat, segdup, \
           homopolymer, ref_seq, dinucleotide, near_indel, ucsc_rpmsk,\
           tumor_ref_count, tumor_alt_count, tumor_maf, tumor_lower_ci, tumor_upper_ci, tumor_is_greater, \
           normal_ref_count, normal_alt_count, normal_maf, normal_lower_ci, normal_upper_ci, normal_is_greater, \
           is_blood, is_sperm]).T, columns=HEADER)

    final.to_csv(outfile, sep="\t", index=False, header=True)


   
def main(argv):
    if len(argv) != 10:
        sys.stderr.write("usage: " + argv[0] + "<tumor_id> <exon_vf> <vf> <gnomad_dropped> <repeats> \
                                                <homopolymer> <tumor_ci> <normal_ci> <output>\n")
        sys.exit(2)
    
    tumor_id = argv[1]
    exon_vf = argv[2]
    vf = argv[3]
    gnomad_dropped = argv[4]
    repeats = argv[5]
    homopolymer = argv[6]
    tumor_ci = argv[7]
    normal_ci = argv[8]
    outfile = argv[9]

    summarize_annotations(tumor_id, exon_vf, vf, gnomad_dropped, repeats,\
                          homopolymer, tumor_ci, normal_ci, outfile)


main(sys.argv)
