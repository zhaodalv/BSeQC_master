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
from BSeQC.read.read_info import read_information as RI


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


def trim_sam(sam_inf, strand_t, read_l, single_on, name, s_path):
    '''
    Trim the mapping files with the biased positions of every length in every strand,
    which are saved in the variance: strand_t.
    '''
    for sam in sam_inf:
        info("Trim the %s..." % sam)
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

        #scan every read to trim
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
                out.write(read)
            else:
                if single_on:
                    if read_l[0] != '':
                        original_length = int(read_l[sam_inf.index(sam)])
                    else:
                        original_length = ''
                    read_trim_single(read, strand_t, out, read_info, original_length)
                else:
                    if read_l[0] != '':
                        original_length = [int(i) for i in read_l[sam_inf.index(sam)].split('-')]
                    else:
                        original_length = ''
                    read_trim_pair(read, strand_t, out, read_info, original_length)
        info("Get the trimmed %s : %s !!" % (sam, out_sam))
        out.close()


def read_trim_single(read, strand_t, out, read_info, original_length):
    '''
    specialise in trimming single-end read
    '''

    strand, pos, CIGAR, seq, score = read_info[1], int(read_info[3]), read_info[4], read_info[5], read_info[6]
    if original_length:
        readlen = original_length
    else:
        readlen = len(seq)
    if strand_t[strand][readlen] == 'NA':
        out.write(read)
        return

    #get the trimming positions
    #note: the length of some mapping reads may be shorter than the original sequence length
    if len(seq) >= (strand_t[strand][readlen][1] + 1):
        exclude_3 = strand_t[strand][readlen][1]
    else:
        exclude_3 = len(seq) - 1
    exclude_5 = strand_t[strand][readlen][0]

    if strand == '-+':
        pos = pos + (len(seq) - exclude_3 - 1)
        seq = seq[::-1][exclude_5:exclude_3 + 1][::-1]
        score = score[::-1][exclude_5:exclude_3 + 1][::-1]
        CIGAR = str(len(seq)) + CIGAR[-1]
    else:
        pos = pos + exclude_5
        seq = seq[exclude_5:exclude_3 + 1]
        score = score[exclude_5:exclude_3 + 1]
        CIGAR = str(len(seq)) + CIGAR[-1]

    #output the new read
    original_read_info = read.rstrip().split('\t')
    out.write('%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s'
              % (original_read_info[0], original_read_info[1], original_read_info[2], pos, original_read_info[4],
                 CIGAR, original_read_info[6], original_read_info[7], original_read_info[8], seq, score))
    for i in range(11, len(original_read_info)):
        out.write('\t%s' % original_read_info[i])
    out.write('\n')
    return


def read_trim_pair(read, strand_t, out, read_info, original_length):
    '''
    specialise in trimming paired-end read
    '''

    strand, pos1, CIGAR, pos2, insert, seq, score =\
    read_info[1], int(read_info[3]), read_info[4], int(read_info[5]), abs(int(read_info[6])), read_info[7], read_info[8]
    if original_length:
        if strand == '++' or strand == '-+':
            readlen = original_length[0]
            p_len = original_length[1]
        else:
            readlen = original_length[1]
            p_len = original_length[0]
    else:
        readlen = len(seq)
        p_len = insert - abs(pos2 - pos1)

    if pos1 == pos2:
        out.write(read)
        return
    '''
    if not strand_t[strand].has_key(readlen):
        out.write(read)
        info(read.rstrip())
        return
    '''
    if strand_t[strand][readlen] == 'NA':
        out.write(read)
        return
    if (pos2 < pos1 and (strand == '++' or strand == '--') or(pos1 < pos2 and (strand == '+-' or strand == '-+'))):
        out.write(read)
        return
    p_dict = {'++': '+-', '+-': '++', '-+': '--', '--': '-+'}
    p_strand = p_dict[strand]
    #get the trimming positions
    #using the other paired read to decide the new pos2 and insert size
    p_exclude_5 = strand_t[p_strand][p_len][0]
    p_exclude_3 = strand_t[p_strand][p_len][1]
    exclude_5 = strand_t[strand][readlen][0]
    exclude_3 = strand_t[strand][readlen][1]

    if strand == '++' or strand == '--':
        pos1_new = pos1 + exclude_5
        pos2_new = pos2 + (p_len - 1 - p_exclude_3)  #the correspond pair read
        pos1_3_new = pos1 + exclude_3
        pos2_3_new = pos2 + (p_len - 1 - p_exclude_5)
        #after trimming if 5' > 3'
        if pos1_new > pos2_new:
            pos2_new = pos1_new
        if pos1_3_new > pos2_3_new:
            pos1_3_new = pos2_3_new
        start = pos1_new - pos1
        end = pos1_3_new - pos1 + 1
        seq = seq[start:end]
        score = score[start:end]
        CIGAR = str(len(seq)) + CIGAR[-1]
        #add the mate read length
        insert = pos2_new - pos1_new + (pos2_3_new - pos2_new + 1)
    else:
        pos1_new = pos1 + (readlen - 1 - exclude_3)
        pos2_new = pos2 + p_exclude_5
        pos1_3_new = pos1 + (readlen - 1 - exclude_5)
        pos2_3_new = pos2 + p_exclude_3
        #after trimming if 5' > 3'
        if pos1_new < pos2_new:
            pos1_new = pos2_new
        if pos1_3_new < pos2_3_new:
            pos2_3_new = pos1_3_new
        start = pos1_new - pos1
        end = pos1_3_new - pos1 + 1
        seq = seq[start:end]
        score = score[start:end]
        CIGAR = str(len(seq)) + CIGAR[-1]
        #add self read length
        insert = -(pos1_new - pos2_new + (pos1_3_new - pos1_new + 1))
    original_read_info = read.rstrip().split('\t')
    out.write('%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%s'
              % (original_read_info[0], original_read_info[1], original_read_info[2], pos1_new, original_read_info[4],
                 CIGAR, original_read_info[6], pos2_new, insert, seq, score))
    for i in range(11, len(original_read_info)):
        out.write('\t%s' % original_read_info[i])
    out.write('\n')
    return
