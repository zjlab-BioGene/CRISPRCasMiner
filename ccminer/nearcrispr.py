import os,sys,re
from time import sleep
import pandas as pd
import subprocess
import logging

class Motifs(object):
    
    def __init__(self,aaseq):
        self.aaseq = aaseq
        
    def searchMotifs(self):
        self.search_motifs()
    
    def search_motifs(self):
        
        def match_pattern(motif,seq):
            ## search for pattern match
            mtf_match = list(re.finditer(motif,seq))
            num = len(mtf_match)
            if num > 0:
                mtfs = []
                poses = []
                for mtf in mtf_match:
                    mtf_str = mtf.group(0)
                    pos = str(mtf.span()).replace(' ','')
                    mtfs.append(mtf_str)
                    poses.append(pos)
            pattern = ';'.join(mtfs) if num > 0 else None
            pos = ';'.join(poses) if num > 0 else None
            return pattern,num,pos
                         
        ## search for Zinfinger: CXXC, CXXC
        self.CXXC, self.CXXCnum, self.CXXCpos = match_pattern(r'C.{2}C',self.aaseq)
        self.RXXXXH, self.RXXXXHnum, self.RXXXXHpos = match_pattern(r'R....H',self.aaseq)
        
class NearCRISPR(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)
    
    def getPotentialCas(self):
        self.fetch_proteins_near_crispr()
        self.get_hmm_annotation()
        self.get_protein_seq()
        self.get_motifs()
        self.write_proteins_near_crispr()
    
    def fetch_proteins_near_crispr(self):
        crisprtab = os.path.join(self.out,'crisprs_calib.tab')
        genetab = os.path.join(self.out,'genes.tab')
        
        ## read CRISPRs
        df1 = pd.read_csv(crisprtab,sep='\t',header=0)
        df1 = df1[['Contig','CRISPR','DRrepeat','N_repeats','Start','End']]
        df1 = df1.set_index('CRISPR').T.to_dict(orient='list')
        crispr = {}
        for key,vlis in df1.items():
            if vlis[0] not in crispr.keys():
                crispr[vlis[0]] = {}
            crispr[vlis[0]][key] = [vlis[3],vlis[4],vlis[1],vlis[2]] # crispr[contig] = [start, end, repeat, n_repeat]
            
        ## read Genes
        df2 = pd.read_csv(genetab,sep='\t',header=0)
        df2['Protein'] = df2['Contig'] + '_' + df2['Pos'].map(str)
        contig, start, end, strand, pos = df2.set_index('Protein').T.to_dict(orient='split')['data']
        proteins = df2.set_index('Protein').T.to_dict(orient='split')['columns']
        starts = {}
        _starts = {}
        ends = {}
        _ends = {}
        strands = {}
        distances = {}
        prior_end = 0
        for i in range(len(contig)):
            if contig[i] not in crispr.keys(): continue
            if contig[i] not in starts.keys():
                starts[contig[i]] = {}
                _starts[contig[i]] = {}
                ends[contig[i]] = {}
                _ends[contig[i]] = {}
                strands[contig[i]] = {}
                distances[contig[i]] = {}
            starts[contig[i]][start[i]] = pos[i]
            _starts[contig[i]][pos[i]] = start[i]
            ends[contig[i]][end[i]] = pos[i]
            _ends[contig[i]][pos[i]] = end[i]
            strands[contig[i]][proteins[i]] = strand[i]
            distances[contig[i]][proteins[i]] = start[i] - prior_end - 1
            prior_end = end[i]
            
        ## Find proteins near CRISPR array: in the same contig, span n proteins, unstream and downstream repectively
        if len(starts.keys()) > 0 :
            header = ['Name','Contig','CRISPR','CRISPRstart','CRISPRend','DRrepeat','N_repeats','Repeat_len',
                    'Protein','PROstart','PROend','Strand','Remark','Rank','PriorDistance','CRISPRdistance']
            nearcris = pd.DataFrame(columns=header)
            for ctg in starts.keys():
                for key,clis in crispr[ctg].items():
                    # clis[0] = cpr_start, clis[1] = cpr_end, clis[2] = repeat
                    down_stream = sorted(starts[ctg].values())[-1]
                    up_stream = sorted(starts[ctg].values())[0]
                    down, up = down_stream, up_stream
                    for s in sorted(starts[ctg].keys()):
                        if s >= clis[1]:
                            if starts[ctg][s] < down:
                                down = starts[ctg][s]
                    for e in sorted(ends[ctg].keys()):
                        if e <= clis[0]:
                            if ends[ctg][e] > up:
                                up = ends[ctg][e]
                    up_count = self.span
                    if down+self.span < down_stream:
                        down_stream = down+self.span
                    if up-self.span+1 > up_stream:
                        up_stream = up-self.span+1
                    else:
                        up_count = up               
                    
                    j = 0 ## upstream rank
                    k = 0 ## downstream rank
                    crisdist = 0
                    for i in range(up_stream,down_stream):
                        ## Protein start & end pos
                        _start = _starts[ctg][i]
                        _end = _ends[ctg][i]
                        ## Upstream or Downstream
                        if _end <= clis[0]:
                            remark = 'up'
                            rank = up_count - j
                            j += 1
                            crisdist = clis[0] - _end -1 
                        elif (_end>clis[0] and _end<=clis[1]) or (_start>=clis[0] and _start<clis[1]) or (_start<=clis[0] and _end>=clis[1]):
                            remark = 'on'
                            rank = 0
                            crisdist = 0
                        else:
                            remark = 'down'
                            k += 1
                            rank = k
                            crisdist = _start - clis[1] - 1
                        nearcris.loc[len(nearcris.index)] = [self.name,ctg,key,clis[0],clis[1],clis[2],clis[3],len(clis[2]),
                                                             ctg+'_'+str(i),_start,_end,strands[ctg][ctg+'_'+str(i)],remark,rank,distances[ctg][ctg+'_'+str(i)],crisdist]
            ## Near CRISPR object
            nearcris = nearcris.drop_duplicates()
            self.nearcrispr = nearcris
            
        else:
            logging.error('No Proteins near CRISPR found! Exit.')
            sys.exit()
        
    def get_hmm_annotation(self):
        # runze.cai 修改于 2023.6.28
        allhmm = os.path.join(self.cctyper_path,'hmmer.tab')
        def annotation(hmm):
            ann = hmm
            if hmm.startswith('Cas12'):
                ann = 'Cas12'
            elif hmm.startswith('Cas13'):
                ann = 'Cas13'
            elif hmm.startswith('Cas9'):
                ann = 'Cas9'
            return ann
                    
        ## Fetch HMM annotation
        hmm = pd.read_csv(allhmm,sep='\t',header=0)
        # runze.cai修改于2023.6.28
        hmm.rename(columns={'ORF':'Protein'},inplace=True)
        hmm['HMM'] = [x.split('_')[0] for x in hmm['Hmm']]
        hmm['Ann'] = [annotation(x) for x in hmm['HMM']]
        hmm['refLen'] = hmm['qlen']
        hmm['hEval'] = hmm['Eval']
        hmm['hScore'] = hmm['score']
        hmm['hCov_seq'] = hmm['Cov_seq']
        hmm['hCov_hmm'] = hmm['Cov_hmm']
        hmm = hmm[['Protein','Ann','HMM','Hmm','refLen','hEval','hScore','hCov_seq','hCov_hmm']]
        ## Choose the best hits
        besthits = pd.DataFrame(columns=['Protein','Ann','HMM','Hmm','refLen','hEval','hScore','hCov_seq','hCov_hmm'])
        for p,df in list(hmm.groupby(by='Protein')):
            df = df.sort_values(by=['hEval','hScore'],ascending=[True,False])
            first_line = list(df.iloc[0])
            first_line[3] = ';'.join(set(df.Hmm))
            besthits.loc[len(besthits.index)] = first_line
        
        ## Merge hmm annotation into nearcrispr
        self.nearcrispr = pd.merge(self.nearcrispr,besthits,left_on=['Protein'],right_on=['Protein'],how='left')
        self.nearcrispr[['Ann','HMM']] =  self.nearcrispr[['Ann','HMM']].fillna('Unknown')
        self.nearcrispr.to_csv(os.path.join(self.out,'nearCRISPR.tab'),sep='\t',index=False)        
    
    def get_protein_seq(self):
        
        ## aaSequence, Length and Completeness of all proteins near CRISPR
        proteinLength = pd.DataFrame(columns=['Protein','Prolen','Complete','aaSeq'])
        allglist = os.path.join(self.out,'all.glist')
        allpro = open(allglist,'w')
        for p in set(self.nearcrispr.Protein):
            allpro.write(p+'\n')
        allpro.close()
        fasta = os.path.join(self.out,'nearCRISPR.faa')
        with open(fasta,'w') as nearcrispr_faa:
            subprocess.run(['seqtk', 'subseq', self.prot_path, allglist], stdout=nearcrispr_faa)
        tab = os.path.join(self.out,'nearCRISPR.tab')
        with open(tab,'w') as nearcrispr_tab:
            subprocess.run(['bioawk', '-cfastx', '{print $name,$seq}', os.path.join(self.out, 'nearCRISPR.faa')], stdout=nearcrispr_tab)
        with open(tab,'r') as nearcrispr_tab:
            for line in nearcrispr_tab:
                p,s = line.split()[0:2]
                if not s.startswith('M'):
                    complete = 'Suspicious'
                elif s.endswith('*'):
                    complete = 'True'
                else:
                    complete = 'False'
                s = s.strip('*')
                length = len(s)
                proteinLength.loc[len(proteinLength.index)] = [p,length,complete,s]
        self.nearcrispr = pd.merge(self.nearcrispr,proteinLength,left_on='Protein',right_on='Protein',how='left')
        ## clean data
        del proteinLength
        for f in [allglist,fasta,tab]:
            os.remove(f)
   
    def get_motifs(self):
        motifs = pd.DataFrame(columns=['Protein','CXXC','CXXCnum','CXXCpos','RXXXXH','RXXXXHnum','RXXXXHpos'])
        
        ## Fetch target effector Cas
        glist = os.path.join(self.out,'target.glist')
        effectors = set(self.nearcrispr[(self.nearcrispr.Ann=='Cas12') | (self.nearcrispr.Ann=='Cas13') | (self.nearcrispr.Ann=='Cas9')].Protein)
        if len(effectors) == 0:
            logging.error('No candidate Cas found near CRISPR!')
        else:
            with open(glist,'w') as target:
                for e in effectors:
                    target.write(e+'\n')
                    
            ## fetch target protein sequence using bioawk
            fasta = os.path.join(self.out, 'putative_cas.faa')
            with open(fasta,'w') as pcas_faa:
                subprocess.run(['seqtk', 'subseq', self.prot_path, glist], stdout=pcas_faa)
            tab = os.path.join(self.out, 'putative_cas.tab')
            with open(tab,'w') as pcas_tab:
                subprocess.run(['bioawk', '-cfastx', '{print $name,$seq}',fasta],stdout=pcas_tab)
                
            ## search Zinfingers motif: CXXC, CXXC
            with open(tab,'r') as pcas_tab:
                for line in pcas_tab:
                    p,s = line.split()[:2]
                    s = s.strip('*')
                    m = Motifs(s)
                    m.searchMotifs()
                    motifs.loc[len(motifs.index)] = [p,m.CXXC,m.CXXCnum,m.CXXCpos,m.RXXXXH,m.RXXXXHnum,m.RXXXXHpos]
            
            ## clean
            os.remove(tab)
            
        self.nearcrispr = pd.merge(self.nearcrispr,motifs,left_on='Protein',right_on='Protein',how='left')
        del motifs
        #os.remove(glist)
                       
    def write_proteins_near_crispr(self):
        ## write near-crispr-proteins table: nearCRISPR.txt
        self.nearcrispr = self.nearcrispr.sort_values(by=['Name','Contig','CRISPRstart','PROstart']).drop_duplicates()
        selcols = ['Name','Contig','CRISPR','DRrepeat','N_repeats','Repeat_len',
                   'PriorDistance','CRISPRdistance','Complete','Protein','Remark','Strand','Rank','Prolen',
                   'HMM','Hmm','CXXCnum','CXXCpos','CXXC','RXXXXHnum','RXXXXHpos','RXXXXH','refLen','hEval','hScore','hCov_seq','hCov_hmm']
        self.nearcrispr[selcols].to_csv(os.path.join(self.out,'nearCRISPR.txt'),sep='\t',index=False)
        
        ## write near-crispr potential unknown proteins:
        