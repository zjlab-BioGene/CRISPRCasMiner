#!/usr/bin/env python

import argparse
import logging
import os,sys,re
import pkg_resources

from generator import Generator
from prodigal import Prodigal
from hmmer import HMMER
from minced import Minced
from nearcrispr import NearCRISPR
from filters import CRISPRCas
#from plot import Map

########## Arguments ##########
def get_parameters():
    ap = argparse.ArgumentParser()
    # Required
    ap.add_argument('input', help='Input fasta file')
    ap.add_argument('output', help='Prefix for output directory')

    # Optional
    ap.add_argument('--name','-n',type=str,default=None,help='Sample Name')
    ap.add_argument('-t', '--threads', help='Number of parallel processes [%(default)s].', default=4, type=int)
    ap.add_argument('--prodigal', help='Which mode to run prodigal in [%(default)s].', default='meta', type=str, choices=['single','meta'])
    ap.add_argument('--circular', help='Input should be treated as circular.', action='store_true')
    ap.add_argument('--keep_tmp', help='Keep temporary files (prodigal, hmmer, minced).', action='store_true')
    ap.add_argument('--log_lvl', help='Logging level [%(default)s].', default='INFO', type=str, choices=['DEBUG','INFO','WARNING','ERROR'])
    ap.add_argument('--redo_typing', help='Redo the typing. Skip prodigal and HMMER and load the hmmer.tab from the output dir.', action='store_true')
    ap.add_argument('--simplelog', help='No color or progress bar in log.', action='store_true')

    # Data
    apd = ap.add_argument_group('data arguments')
    apd.add_argument('--db', help='Path to database. Only needed if CCMINER_DB environment variable is not set.', default='', type=str)

    # MinCED
    apm = ap.add_argument_group('MinCED arguments')
    apm.add_argument('--seed', help='Seed for sampling when estimating CRISPR stats (see argument above) [%(default)s].', default=42, type=float)
    apm.add_argument('--searchWL', help='MinCED option - Length of search window. Range: 6-9.', default=8, type=int)
    apm.add_argument('--minNR', help='MinCED option - Minimum number of repeats.', default=3, type=int)
    apm.add_argument('--minRL', help='MinCED option - Minimum repeat length.', default=23, type=int)
    apm.add_argument('--maxRL', help='MinCED option - Maximum repeat length.', default=47, type=int)
    apm.add_argument('--minSL', help='MinCED option - Minimum spacer length.', default=26, type=int)
    apm.add_argument('--maxSL', help='MinCED option - Maximum spacer length.', default=50, type=int)
    apm.add_argument('--cctyper_path', help='cctyper path.', type=str) #runze.cai修改,2023.2.24

    # NearCRISPR
    apn = ap.add_argument_group('NearCRISPR arguments')
    #span default改为20 runze.cai 2023.6.15
    apn.add_argument('--span',type=int,default=20,help='number of proteins around CRISPR array')

    # filters
    apf = ap.add_argument_group('Filters arguments')
    apf.add_argument('--database_name','-db_name',default='database',type=str,help='Database name.')
    apf.add_argument('--need_prolen','-prolen',default= 700,type=int,help='need prolen.')
    apf.add_argument('--need_rank','-rank',default= 3,type=int,help='need rank.')
    apf.add_argument('--large_prolen','-large',default= 2000,type=int,help='large prolen.')
    apf.add_argument('--split_prolen_upper','-sup',default= 1500,type=int,help='split prolen upper.')
    apf.add_argument('--split_prolen_lower','-slow',default= 700,type=int,help='split prolen lower.')

    # Plot
    app = ap.add_argument_group('plotting arguments')
    app.add_argument('--expand', help='Expand operons with un-annotated genes. The value determines by how many bp in each end to expand. 0 only fills gaps [%(default)s].', default=0, type=int)
    app.add_argument('--no_plot', help='Do not draw a map of CRISPR-Cas.', action='store_true')
    app.add_argument('--no_grid', help='Do not add grid to plot.', action='store_true')
    
    return ap

def main():
    
    ########## Initialize ##########
    ap = get_parameters()
    master = Generator(ap.parse_args())

    ########## Prodigal ##########
    proteins = Prodigal(master)
    proteins.run_prod()

    ########## Hmmer ##########
    casproteins = HMMER(proteins)
    casproteins.main_hmm()

    ########## CRISPRs ##########
    crispr_obj = Minced(casproteins)
    crispr_obj.run_minced()

    ######### NearCRISPR ########
    nearcrispr = NearCRISPR(crispr_obj)
    nearcrispr.getPotentialCas()

    ######### CRISPR-Cas Filters ########
    candidates = CRISPRCas(nearcrispr)
    candidates.get_candidates()

if __name__ == '__main__':
    
    main()