import os
import logging
import sys
import shutil
import json
import pkg_resources
import subprocess
import pandas as pd
from Bio import SeqIO

class Generator(object):

    def __init__(self, args):
       
        self.fasta = args.input
        self.out = args.output
        self.name = os.path.basename(os.path.abspath(args.output).strip('/')) if args.name == None else args.name
        self.threads = args.threads
        self.prod = args.prodigal
        self.db = args.db
        self.circular = args.circular
        self.keep_tmp = args.keep_tmp
        self.lvl = args.log_lvl
        self.redo = args.redo_typing
        self.noplot = args.no_plot
        self.nogrid = args.no_grid
        self.expand = args.expand
        self.simplelog = args.simplelog
        self.seed = args.seed

        # MinCED
        self.searchWL = args.searchWL
        self.minNR = args.minNR
        self.minRL = args.minRL
        self.maxRL = args.maxRL
        self.minSL = args.minSL
        self.maxSL = args.maxSL
        self.cctyper_path = args.cctyper_path #runze.cai修改 2023.2.27
        
        # NearCRISPR
        self.span = args.span

        # CRISPRCas
        self.database_name = args.database_name
        self.need_prolen = args.need_prolen
        self.need_rank = args.need_rank
        self.large_prolen = args.large_prolen
        self.split_prolen_upper = args.split_prolen_upper
        self.split_prolen_lower = args.split_prolen_lower
        self.out_prolen_max = args.out_prolen_max
        self.out_prolen_min = args.out_prolen_min

        self.any_cas = False
        self.any_operon = False
        self.any_crispr = False

        # Logger
        if self.simplelog:
            logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.lvl)
        else:
            logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.lvl)
        #logging.info('Running CRISPRCasMiner version {}'.format(pkg_resources.require("ccminer")[0].version))

        # Force consistency
        self.out = os.path.join(self.out, '')
        #runze.cai 修改 2023.6.27
        self.prot_path = self.cctyper_path+'/proteins.faa'

        # Check databases
        self.check_db()
        
        # Check input and output
        self.check_out()
        self.check_input()

        # If redo check if any crisprs and operons
        if self.redo:
            if os.path.exists(self.out+'cas_operons.tab') or os.path.exists(self.out+'cas_operons_putative.tab'):
                self.any_operon = True
            if os.path.exists(self.out+'crisprs_all.tab'):
                self.any_crispr = True

        # Write arguments
        da = vars(args)
        f = open(self.out+'arguments.tab', 'w')
        for k, v in da.items():
            f.write('{}:\t{}\n'.format(k, v))
        f.close()

    def check_db(self):
        
        if self.db == '':
            try:
                #self.db = os.environ['CCMINER_DB']
                self.db = '/data/home/runze.cai/data/Analysis/mgnify/ccminer/db/'
            except:
                logging.error('Could not find database directory')
                sys.exit()

        self.pdir = os.path.join(self.db, 'Profiles', '')

        # Look if HMM profile dir exists
        if os.path.isdir(self.pdir):
            for i in os.listdir(self.pdir):
                if not i.lower().endswith('.hmm'):
                    logging.error('There are non-HMM profiles in the HMM profile directory')
                    sys.exit()
        else:
            logging.error('Could not find HMM profile directory')
            sys.exit()

    def check_out(self):

        if not self.redo:
            try:
                os.mkdir(self.out)
            except FileExistsError:
                logging.error('Directory '+self.out+' already exists')
                sys.exit()

    def check_input(self):

        if os.path.isfile(self.fasta):
            self.check_fasta()
        else:
            logging.error('Could not find input file')
            sys.exit()

    def check_fasta(self):
        
        # Get sequence lengths
        with open(self.fasta, 'r') as handle:
            self.len_dict = {}
            self.seq_dict = {}
            for fa in SeqIO.parse(handle, 'fasta'):
                if fa.id in self.len_dict:
                    logging.error('Duplicate fasta headers detected')
                    sys.exit()
                self.len_dict[fa.id] = len(fa.seq)
                self.seq_dict[fa.id] = fa.seq
            
        # Check for numeric headers
        self.num_headers = False
        for i in self.len_dict.keys():
            try:
                dump = float(i)
                self.num_headers = True
            except:
                pass
        
        if self.num_headers:
            logging.warning('Numeric fasta headers detected. A prefix is added to the names')
            new_fasta = open(self.out+'fixed_input.fna', 'w')
            subprocess.run(['sed', 's/^>/>Contig/', self.fasta], stdout = new_fasta)
            new_fasta.close()
            self.fasta = self.out+'fixed_input.fna'
            self.len_dict = {'Contig'+str(key): val for key, val in self.len_dict.items()}
            self.seq_dict = {'Contig'+str(key): val for key, val in self.seq_dict.items()}

    def clean(self):
        if not self.redo:

            if self.num_headers:
                os.remove(self.out+'fixed_input.fna')

            if os.stat(self.out+'hmmer.log').st_size == 0:
                os.remove(self.out+'hmmer.log')
                
            if not self.keep_tmp:
                
                logging.info('Removing temporary files')
                
                shutil.rmtree(self.out+'hmmer')
                
                os.remove(self.out+'minced.out')
                os.remove(self.out+'prodigal.log')
                os.remove(self.out+'proteins.faa')