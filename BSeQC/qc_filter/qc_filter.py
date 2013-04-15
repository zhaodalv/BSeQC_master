#!/usr/bin/python

# ------------------------------------
#python package
# ------------------------------------
import numpy as np
import logging
import sys
import os

# ------------------------------------
#the read package
# ------------------------------------
from BSeQC.read.read_info import MappingReader as RI
from BSeQC.qc_filter.single_filter import read_trim_single as SF
from BSeQC.qc_filter.paired_filter import read_trim_paired as PF
from BSeQC.qc_filter.duplicate_filter import duplicate_filter as DF


# ------------------------------------
#logging object
# ------------------------------------

logging.basicConfig(level=20,
    format=' %(levelname)-5s @ %(asctime)s: %(message)s ',
    datefmt='%a, %d %b %Y %H:%M:%S',
    stream=sys.stderr,
    filemode="w"
)
info = logging.info
error = logging.error
warning = logging.warning


def filter_sam(sam_inf, strand_t, read_l, single_on, name, s_path, auto, remove_overlap, loc_dict, max_cov,
               not_mapping):
    '''
    Trim the mapping files with the biased positions of every length in every strand,
    which are saved in the variance: strand_t.
    '''
    for sam in sam_inf:
        out_sam = sam[:-4] + '_' + name + '_filter.sam'
        out = open(out_sam, 'w')

        #check the input mapping files
        if sam[-4:].upper() == '.SAM':
            sam_format = 1
            read_inf = open(sam)
        elif sam[-4:].upper() == '.BAM':
            sam_format, read_inf = 1, os.popen('%ssamtools view -h %s' % (s_path, sam))
        else:
            error("The input mapping file is not SAM format or BAM format")

        #scan every read to qc_filter
        for read in read_inf:
            #for sam header
            if read.startswith('@'):
                out.write(read)
                continue

            #Get the read information for trimming
            #If the read isn't unique mapping, we will get a empty list ([]).
            #In: single unique mapping read  Out: [flag,strand,chr,pos,CIGAR,seq,score]
            #In: paired unique mapping read  Out: [flag,strand,chr,pos1,CIGAR,pos2,insert,seq,score]
            read_info = RI(read)
            read_info = read_info.extract_information()
            if len(read_info) == 0:
                if not_mapping:         #keep the not_unique mapping reads
                    out.write(read)
                continue

            if len(loc_dict) > 0: # the --filter_dup has been set True, have to remove duplicate reads
                duplicate, loc_dict = DF(read_info, loc_dict, max_cov, single_on)
            else:
                duplicate = False

            if single_on:
                if auto:
                    if read_l[0] != '':
                        original_length = int(read_l[sam_inf.index(sam)])
                    else:
                        original_length = ''
                    SF(read, strand_t, out, read_info, original_length, duplicate)
                    continue
                else:
                    if not duplicate and len(loc_dict) > 0:
                        out.write(read)                  # not trimming, only output not_duplicate reads
                        continue
            else:
                if auto or remove_overlap:
                    if read_l[0] != '':
                        original_length = [int(i) for i in read_l[sam_inf.index(sam)].split('-')]
                    else:
                        original_length = ''
                    PF(read, strand_t, out, read_info, original_length, auto, remove_overlap, duplicate)
                else:
                    if not duplicate and len(loc_dict) > 0:
                        out.write(read)                  # not trimming, only output not_duplicate reads
        out.close()



	
