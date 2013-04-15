#!/usr/bin/python

'''
Use the Poisson distribution to estimate the max coverage based on the sequence depth.

Export the figure of duplicate distribution
'''

# ------------------------------------
#python package
# ------------------------------------
import logging
import numpy as np
from scipy.stats import poisson


# ------------------------------------
#logging object
# ------------------------------------
import sys

logging.basicConfig(level=20,
    format=' %(levelname)-5s @ %(asctime)s: %(message)s ',
    datefmt='%a, %d %b %Y %H:%M:%S',
    stream=sys.stderr,
    filemode="w"
)
info=logging.info
error=logging.error
warning=logging.warning


def duplicate_report(loc_dict,gsize,p_poisson,name):

    total_tag, fre_dict = summary_dup(loc_dict)
    max_dup_cov = cal_max_dup_cov(gsize,total_tag,p_poisson)
    dup_density_plot(fre_dict,max_dup_cov,name)
    return max_dup_cov

def cal_max_dup_cov(gsize,total_tag,p_poisson):
    '''
    Use the Poisson distribution to estimate the max coverage
    The default pvalue is 1e-5
    '''
    poisson_model = poisson(total_tag/gsize)
    max_dup_cov = poisson_model.ppf(1- p_poisson)
    return max_dup_cov

def summary_dup(loc_dict):

    total_tag = 0
    fre_dict = {}
    for chr in loc_dict.keys():
        total_tag += sum(loc_dict[chr].values())
        cal_frequency = np.bincount(loc_dict[chr].values())
        ii = np.nonzero(cal_frequency)[0]
        Frequency = zip(ii,cal_frequency[ii]) # the duplicate read frequency in every chrome: [(coverage: #reads)]
        for f in range(len(Frequency)):
            if fre_dict.has_key(Frequency[f][0]):
                fre_dict[Frequency[f][0]] += Frequency[f][1]
            else:
                fre_dict[Frequency[f][0]] = Frequency[f][1]
    return total_tag, fre_dict

def dup_density_plot(fre_dict, max_dup_cov,name):

    try:
        '''
        The matplotlib, python module, is used for draw the Mbias plot.
        We recommend installing it.
        '''
        import matplotlib
        #matplotlib.use('Agg')
        from matplotlib.lines import Line2D
        import matplotlib.pyplot as plt
        pdf_off = False
    except:
        info('Please install the matplotlib module to draw the Mbias plot./'
             ' If not, you can only get the Mbias table')
        pdf_off = True


    if pdf_off:
        info("Can't import matplotlib package. Skip the Dup_dis.pdf")
    else:
        fre_array = np.array(fre_dict.items())
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(fre_array[:,0], fre_array[:,1], lw=2, c='red')
        ax.set_ylim(0,fre_array[:,1].max())
        #ax.set_xlim(0,fre_array[:,0].max())
        ax.set_xlim(0,25)
        ax.set_ylabel(' # Reads')
        ax.set_xlabel(' read coverage')
        ax.add_line(Line2D([max_dup_cov,max_dup_cov],[0,fre_array[:,1].max()], ls = '--', lw=1))
        ax.annotate('The max coverage: %s' %max_dup_cov, xy = (8,fre_array[:,1].max()/2))
        plt.savefig(name + '_Dup_dis.pdf', format='PDF')
        info('The distribution of duplicate reads has been drawn!')
    return