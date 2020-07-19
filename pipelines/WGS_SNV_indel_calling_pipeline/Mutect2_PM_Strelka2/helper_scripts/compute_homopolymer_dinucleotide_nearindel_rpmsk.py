import sys
import os
import genome
import genome.db
import re
import util.txtfile
import util.file
import genome.vcf
import gzip

header = ["#CHROM",
          "POS",
          "REF",
          "ALT",
          "REF_SEQ",
          "HOMOPOLYMER",
          "DINUCLEOTIDE",
          "NEAR_INDEL"]


def check_if_in_homopolymer(seq_str):

    pattern = re.compile(r'([ACGT])\1{3,}')
    matches = [m.group() for m in re.finditer(pattern, seq_str)]
    if len(matches) > 0:
        is_homopolymer = 1
    else:
        is_homopolymer = 0

    return is_homopolymer

def check_if_in_dinucleotide_repeat(seq_str):

    pattern = re.compile(r'([ACGT]{2})\1{3,}')
    matches = [m.group() for m in re.finditer(pattern, seq_str)]
    if len(matches) > 0:
        is_dinucleotide = 1
    else:
        is_dinucleotide = 0

    return is_dinucleotide

def make_list_of_germline_indels(path, chrom_dict):
    f = util.file.check_open(path)
    vcf_reader = genome.vcf.VCFReader(f)

    indel_coords = []
    for vcf_row in vcf_reader:
        #print(vcf_row)
        ref_bases = vcf_row.ref_base.split(",")
        alt_bases = vcf_row.alt_base.split(",")

        if sum([len(x) > 1 for x in (ref_bases + alt_bases)]) > 0:
            # this is an indel, get the longest ref allele and make
            # it a coordinate, expanding by 5 bp
            max_size = max([len(x) for x in ref_bases])
            if vcf_row.start <= 5:
                continue
            new_start = vcf_row.start - 5
            new_end = vcf_row.start + max_size + 5

            if not vcf_row.chrom_name.startswith("chr"):
                chrom_name = "chr" + vcf_row.chrom_name
            else:
                chrom_name = vcf_row.chrom_name

            if chrom_name == "chrMT":
                chrom_name = "chrM"
            
            chrom = chrom_dict[str.encode(chrom_name)]
            coord = genome.coord.Coord(chrom, new_start, new_end)
            indel_coords.append(coord)

    f.close()
    genome.coord.sort_coords(indel_coords)
    return indel_coords


def make_repeat_coord_list(path, chrom_dict):
    f = open(path, "r")

    coords = []
    for line in f:
        if line.startswith("#"):
            continue
        words = line.rstrip().split("\t")

        if not words[5].startswith("chr"):
             chrom_name = "chr" + words[5]
        else:
             chrom_name = words[5]

        chrom = chrom_dict[str.encode(chrom_name)]
        start = int(words[6])+1 # ucsc coords are 0-based, convert to 1-based
        end = int(words[7])
        repeat_type = words[11]

        coord = genome.coord.Coord(chrom, start, end)
        coord.repeat_type = repeat_type
        coords.append(coord)

    f.close()

    genome.coord.sort_coords(coords)
    return coords



def main(argv):
    if len(argv) != 5:
        sys.stderr.write("usage: " + argv[0] + " <input_vcf> <germline_vcf_path> <rpmsk_path> <output_file>\n")
        sys.exit(2)
        print("inpur directory is:", bam_dir)

    input_vcf_path = argv[1]
    germline_vcf_path = argv[2]
    rpmsk_path = argv[3]
    output_file = argv[4]


    gdb = genome.db.GenomeDB(assembly="hg19")
    chrom_dict = gdb.get_chromosome_dict()
    seq_track = gdb.open_track("seq")
    print(germline_vcf_path)
    indel_coords = make_list_of_germline_indels(germline_vcf_path, chrom_dict)

    if input_vcf_path.endswith("gz"):
        rfile = gzip.open(input_vcf_path, "rt")
    else:
        rfile = open(input_vcf_path, "r")
   

    repeat_coords = make_repeat_coord_list(rpmsk_path, chrom_dict)
    var_coords = []
    results = []

    for line in rfile:
        if line.startswith("#"):
            continue
        line = line.rstrip().split("\t")
        chrom = "chr" + line[0]
        pos = int(line[1])
        end_pos = pos + len(line[3]) - 1
        ref = line[3]
        alt = line[4]
        b_chrom_name = chrom.encode('UTF-8')
        near_indel = 0

        var_coord = genome.coord.Coord(chrom_dict[b_chrom_name], pos, end_pos)
        var_coords.append(var_coord)

        #homopolymer
        seq_str_9bp = seq_track.get_seq_str(chrom,pos-4,pos+4)
        is_homopolymer = check_if_in_homopolymer(seq_str_9bp)
        #dinucleotide
        seq_str_17bp = seq_track.get_seq_str(chrom, pos-8, pos+8)
        is_dinucleotide = check_if_in_dinucleotide_repeat(seq_str_17bp)
        #near indel
        indel_overlaps = genome.coord.get_coord_overlaps(var_coord, indel_coords)
        if len(indel_overlaps) > 0:
           near_indel = 1
        results.append([chrom, pos, ref, alt, seq_str_17bp, is_homopolymer, is_dinucleotide, near_indel])
    rfile.close()

    genome.coord.sort_coords(var_coords)
    overlaps = genome.coord.get_overlaps(var_coords, repeat_coords)
    repeat_dict = {}
    for i in range(len(overlaps)):
       var_coord = var_coords[i]
       chrom = var_coord.chrom.name.decode()
       pos = var_coord.start
       if len(overlaps[i]) > 0:
           repeat_type = overlaps[i][0].repeat_type
       else:
           repeat_type = "pass"
       repeat_dict[chrom + ":" + str(pos)] = repeat_type


    with open(output_file, "w") as wfile:
        for element in results:
            chrom = element[0]
            pos = element[1]
            repeat_type = repeat_dict[chrom + ":" + str(pos)]
            wfile.write("\t".join(map(str, element + [repeat_type])) + "\n")
    wfile.close()

main(sys.argv)
