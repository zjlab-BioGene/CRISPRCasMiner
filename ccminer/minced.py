import os,sys,re
import subprocess
import logging
import math
import random
import shutil
import statistics as st
import numpy as np
import pandas as pd

# Define the CRISPR class
class CRISPR(object):

    def __init__(self, cali_out, sequence, crispr):
        self.cali_out = cali_out
        self.sequence = sequence
        self.crispr = crispr
        self.repeats = []
        self.spacers = []

    def setPos(self, start, end, seqence_len):
        self.start = int(start.rstrip())
        self.end = int(end.rstrip())
        self.sequence_len = int(seqence_len.rstrip())
        
    def addRepeat(self, repeat):
        self.repeats.append(repeat.rstrip())
        
    def addSpacer(self, spacer):
        put_spacer = spacer.rstrip()
        if len(put_spacer) > 0:
            self.spacers.append(put_spacer)
            
    def getConsensus(self):
        self.cons = max(set(self.repeats), key = self.repeats.count) 
        
    def check_Repeats(self):
        n = len(self.cons)
        c = len(self.repeats)
        match = [0]*n
        cons_lst = list(self.cons)
        for repeat in self.repeats:
            rlst = list(repeat)
            for i in range(n):
                if i < len(rlst) and rlst[i] == cons_lst[i]:  #runze.cai修改 2023.2.27
                    match[i] += 1
        self.mean_match_rate = round(np.mean(match)/c,2)
               
    def check_Spacers(self):
        self.calibrated = 'No'
        self.calibLen = 0
        fa = os.path.join(self.cali_out,self.crispr+'.fa')
        spacerfa = open(fa,'w')
        for i in range(len(self.spacers)):
            ss = self.spacers[i]
            spacerfa.write('>'+self.crispr+':'+str(i)+'\n'+ss+'\n')
        spacerfa.close()
        ## MSA
        msa = os.path.join(self.cali_out,self.crispr+'.msa')
        os.system("mafft --auto \'%s\' > \'%s\' 2>/dev/null" % (fa,msa))
        msa_seq = os.path.join(self.cali_out,self.crispr+'.txt')
        os.system("bioawk -cfastx \'{print $seq}\' \'%s\' > \'%s\'" % (msa,msa_seq))
        ## Repeat Calibration
        msa_read = open(msa_seq,'r')
        sps = []
        for l in msa_read:
            sps.append(l.strip().upper())
        n = len(sps[0])
        c = len(sps)
        spacer_codes = {}
        
        for i in range(c):
            spacer_cur = sps[i]
            slst = list(spacer_cur)
            for j in range(n):
                if j not in spacer_codes.keys():
                    spacer_codes[j] = []
                #spacer_codes[j][i] = slst[j]  
                spacer_codes[j].append(slst[j])  
  
        spacer_cons = [.0]*n
        codes = []
        for i in sorted(spacer_codes.keys()):
            code = max(set(spacer_codes[i]), key = spacer_codes[i].count)
            codes.append(code)
            rate = (spacer_codes[i].count(code)+1)/(len(spacer_codes[i])+1)
            if code == '-':
                spacer_cons[i] = -1.0
            elif rate >= self.mean_match_rate:
                spacer_cons[i] = rate
            else:
                spacer_cons[i] = 0
        
        ## Forward check
        fwseq = ''
        if spacer_cons[0] > 0:
            for i in range(len(spacer_cons)):
                if spacer_cons[i] > 0:
                    fwseq += codes[i]
                else:
                    break
        ## Backward check
        bkseq = ''
        if spacer_cons[-1] > 0:
            for i in reversed(range(len(spacer_cons))):
                if spacer_cons[i] > 0:
                    bkseq += codes[i]
                else:
                    break
        bkseq = bkseq[::-1]
        ## Update repeats and spacers
        if len(bkseq) > len(fwseq):
            cut_len = len(bkseq)
            self.calibrated = 'Backward'
            self.calibLen = cut_len
            logging.info('%s spacers need calibration backwardly %d bp...' % (self.crispr,cut_len))
            self.cons = bkseq + self.cons
            for i in range(len(self.repeats)):
                self.repeats[i] = bkseq + self.repeats[i]
            for i in range(len(self.spacers)):
                self.spacers[i] = self.spacers[i][0:len(self.spacers[i])-cut_len]
                
        elif len(fwseq) > len(bkseq):
            cut_len = len(fwseq)
            self.calibrated = 'Farward'
            self.calibLen = cut_len
            logging.info('%s spacers need calibration farwardly %d bp...' % (self.crispr,cut_len))
            self.cons += fwseq
            for i in range(len(self.repeats)):
                self.repeats[i] = self.repeats[i] + fwseq
            for i in range(len(self.spacers)):
                self.spacers[i] = self.spacers[i][cut_len:]
                   
    def stats(self):
        if len(self.spacers) > 1:
            self.spacer_len = round(st.mean([len(x) for x in self.spacers]), 1)
            self.spacer_sem = round(st.stdev([len(x) for x in self.spacers])/math.sqrt(len(self.spacers)), 1)
        else:
            self.spacer_identity = 0
            self.spacer_len = len(self.spacers[0])
            self.spacer_sem = 0
        self.repeat_len = round(st.mean([len(x) for x in self.repeats]), 1)

class RNAfold(object):
    
    def __init__(self,out,seq):
        self.seq = str(seq).strip()
        self.out = os.path.abspath(out) # Modified by Wenhui.Li, 2024.01.03
    
    def runRNAfold(self):
        self.run_rnafold()
        self.parce_rnafold()
    
    def run_rnafold(self):
        # runze.cai 2023.6.28备注 self.out只能用绝对路径
        current_seq = os.path.join(self.out,'current_seq')
        self.rnafold_res = os.path.join(self.out,'rnafold_out')
        f = open(current_seq,'w')
        f.write(self.seq)
        f.close()
        ## run RNAfold
        ## runze.cai 修改于 2023.6.15
        tmp_dir = os.getcwd()
        os.chdir(self.out)
        os.system('RNAfold --MEA -p %s > %s' % (current_seq,self.rnafold_res))
        ## remove all .ps files
        for i in os.listdir('./'):
            if i.endswith('ps'):
                os.remove(i)
        os.chdir(tmp_dir)
    
    def parce_rnafold(self):
        
        def gc_content(seq): 
            seq = seq.strip().upper()
            gc = round( (seq.count('G')+seq.count('C'))/len(seq),2 )
            return gc
        
        def string_sub(string,start,end,old,new):
            lst = list(string)
            for i in range(start,end+1):
                if lst[i] == old:
                    lst[i] = new
            string = ''.join(lst)
            return string

        def parce_rnafold_secondary_structure(seq,snd_str):
            seq = seq.strip().upper()
            snd_str = snd_str.strip()
            _left = []
            _right = []
            pair = []
            nLoop, nBubble, nStem, nLink = 0, 0, 0, 0
            loopLen = []
            bubbleLen = []
            stemLen = []
            linkLen = []
        
            for i in re.finditer(r'\(',snd_str):
                pos = int(i.span()[0])
                _left.insert(0,pos)
            
            for i in re.finditer(r'\)',snd_str):
                pos = int(i.span()[0])
                _right.append(pos)
                
            if len(_left) > 0 and len(_right) > 0:
            
                if _left[-1] > 0 or _right[-1] < len(snd_str)-1:
                    nLink += 1
                    linkLen.append(_left[-1] + len(snd_str) -1 - _right[-1])
                
                nStem = 1
                stem_len = 0       
                for k in range(len(_left)):
                    i = _left[k]
                    stem_len += 1
                    left = seq[i:i+1]
                    right = ''
                    for h in re.finditer(r'\)',snd_str):
                        j = h.span()[0]
                        if j > i:
                            right = seq[j:j+1]
                            segment = snd_str[i+1:j]
                            n_dot = segment.count('.')
                            ## Loops and Bubbles
                            if n_dot > 0:
                                if n_dot == len(segment):
                                    nLoop += 1       
                                    loopLen.append(n_dot)
                                    snd_str = string_sub(snd_str,i+1,j-1,'.','-')
                                else:
                                    nBubble += 1
                                    bubbleLen.append(n_dot)
                                    snd_str = string_sub(snd_str,i+1,j-1,'.','-')                    
                            
                            snd_str = string_sub(snd_str,j,j,')',']')                       
                            break     
                    pair.append(''.join(sorted([left,right])))  

                    ## Stems and Links
                    if k+1 < len(_left):
                        i_next = _left[k+1]
                        tmp_seg = snd_str[i_next+1:i]
                        if tmp_seg.count(r')') > 0:
                            nStem += 1
                            stemLen.append(stem_len)
                            stem_len = 0
                            snd_str = string_sub(snd_str,i,len(snd_str)-1,'(','[')
                            
                            for l in reversed(_right):
                                if l < i:
                                    nLink += 1
                                    segment2 = snd_str[l+1:i]
                                    linkLen.append(segment2.count('.'))
                
                stemLen.append(stem_len)                   
            nGC = pair.count('CG')
            nAU = pair.count('AU')
            nGU = pair.count('GU')
            nPair = nGC + nAU + nGU
            stemScore = round((nGC*3 + nAU*2 + nGU*1)/nPair,2) if nPair > 0 else 0.0
            maxLoopLen = max(loopLen) if len(loopLen) > 0 else 0
            sumLoopLen = sum(loopLen) if len(loopLen) > 0 else 0
            maxBubbleLen = max(bubbleLen) if len(bubbleLen) > 0 else 0
            sumBubbleLen = sum(bubbleLen) if len(bubbleLen) > 0 else 0
            maxStemLen = max(stemLen) if len(stemLen) > 0 else 0
            sumStemLen = sum(stemLen) if len(stemLen) > 0 else 0
            maxLinkLen = max(linkLen) if len(linkLen) > 0 else 0
            sumLinkLen = sum(linkLen) if len(linkLen) > 0 else 0        
            return nLoop,maxLoopLen,sumLoopLen,nBubble,maxBubbleLen,sumBubbleLen,nStem,maxStemLen,sumStemLen,nLink,maxLinkLen,sumLinkLen,nGC,nAU,nGU,nPair,stemScore
            
        def read_rnafold(infile):
            if os.path.exists(infile):
                _content = open(infile,'r').readlines()
                ## line 0, sequence
                seq = _content[0].strip()
                length = len(seq)
                ## line 1, Minimum free energy, MFE
                _tmp = _content[1].strip().split()
                mfe_str = _tmp[0]
                mfe = float(re.sub(r'[() ]','',_tmp[-1]))
                del _tmp
                ## line 2, ensemble free energy, EFE
                _tmp = _content[2].strip().split()
                efe_str = _tmp[0]
                efe = float(re.sub(r'[\[\] ]','',_tmp[-1]))
                efe_score = efe_str.count(',')*1 + efe_str.count('|')*2 + (efe_str.count('{')+efe_str.count('}'))*3 + (efe_str.count('(')+efe_str.count(')'))*5
                del _tmp
                ## line 3, centroid structure, CS
                _tmp = _content[3].strip().split()
                cfe = float(re.sub(r'[{]','',_tmp[-2]))
                distance = float(re.sub(r'[d=}]','',_tmp[-1]))
                del _tmp
                ## line 4, Mininum erra accuracy, MEA
                _tmp = _content[4].strip().split()
                mea_str = _tmp[0]
                mea = float(re.sub(r'[{]','',_tmp[-2]))
                accuracy = float(re.sub(r'[MEA=}]','',_tmp[-1]))
                del _tmp
                ## line 5, Frequency of MEF & ensemble diversity
                _tmp1, _tmp2 = _content[5].strip().split(';')
                fre = float(_tmp1.split()[-1])
                diversity = float(_tmp2.split()[-1])
                del _tmp1, _tmp2
                ## gc content, gc content in stem
                gc = gc_content(seq)
                nLoop,maxLoopLen,sumLoopLen,nBubble,maxBubbleLen,sumBubbleLen,nStem,maxStemLen,sumStemLen,nLink,maxLinkLen,sumLinkLen,nGC,nAU,nGU,nPair,stemScore = parce_rnafold_secondary_structure(seq,mfe_str)
                ## Feature part 1: count
                self.GC = gc
                self.nLoop = nLoop
                self.nBubble = nBubble
                self.nStem = nStem
                self.nLink = nLink
                self.nGC = nGC
                self.nAU = nAU
                self.nGU = nGU
                self.nPair = nPair
                self.stemScore = stemScore
                ## Feature part 2: length & ratio
                self.maxLoopLen = maxLoopLen
                self.sumLoopLen = sumLoopLen  
                self.maxBubbleLen = maxBubbleLen
                self.sumBubbleLen = sumBubbleLen
                self.maxStemLen = maxStemLen
                self.sumStemLen = sumStemLen
                self.maxLinkLen = maxLinkLen
                self.sumLinkLen = sumLinkLen
                self.LoopPercent = round(sumLoopLen/length,2)
                self.BubblePercent = round(sumBubbleLen/length,2)
                self.StemPercent = round(sumStemLen/length,2)
                self.LinkPercent = round(sumLinkLen/length,2)
                self.LinkLoopRatio = round(sumLinkLen/sumLoopLen,2) if sumLoopLen > 0 else 0
                self.StemLoopRatio = round(sumStemLen/sumLoopLen,2) if sumLoopLen > 0 else 0
                self.LinkStemRatio = round(sumLinkLen/sumStemLen,2) if sumStemLen > 0 else 0
                ## Feature part 3: free energy
                self.MFE = mfe
                self.MEA = mea
                self.EFE = efe
                self.CFE = cfe
                self.Distance = distance
                self.Diversity = diversity
                self.Accuracy = accuracy
                self.MFEstr = mfe_str
                self.MEAstr = mea_str
                   
            else:
                print("Error: RNAfold failed!")
        
        read_rnafold(self.rnafold_res)   
          
class Minced(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run_minced(self):

        if not self.redo:
            logging.info('Predicting CRISPR arrays with minced')

            # Run minced
            # runze.cai修改，2023.2.27
            
            self.cctyper_Minced_Blast()
            # subprocess.run(['minced',
            #                 '-searchWL', str(self.searchWL),
            #                 '-minNR', str(self.minNR),
            #                 '-minRL', str(self.minRL),
            #                 '-maxRL', str(self.maxRL),
            #                 '-minSL', str(self.minSL),
            #                 '-maxSL', str(self.maxSL),
            #                 self.fasta, 
            #                 self.out+'minced.out'], 
            #                 stdout=subprocess.DEVNULL, 
            #                 stderr=subprocess.DEVNULL)
        
            # Parse
            self.parse_minced()
            self.parse_rnafold()

            # Write results
            self.write_crisprs()
            self.write_spacers()

    #runze.cai修改 2023.2.27 cctyper得到的array传递至ccminer
    def cctyper_Minced_Blast(self):
        os.system('seqkit fx2tab ' + self.fasta + ' -nl > ' + self.out + '/seqLen.tab')
        seq_len = {}
        f_temp = open(os.path.join(self.out, 'seqLen.tab'))
        for t in f_temp.readlines():
            tt = t.strip().split()
            seq_len[tt[0]] = tt[-1]
        f_temp.close()
        try:
            all_data = pd.read_table(os.path.join(self.cctyper_path, 'crisprs_all.tab'), sep='\t')
        except:
            f_out = open(self.out + '/minced.out','w')
            f_out.close()
            return None
        seq_data = {}
        contig = list(all_data.Contig)
        crispr = list(all_data.CRISPR)
        for id in range(len(crispr)):
            seq_data[crispr[id]] = [contig[id] , seq_len[contig[id]] , str(id+1)]
        
        f_out = open(self.out + '/minced.out','w')
        gff_data = open(os.path.join(self.cctyper_path, 'crisprs.gff'))
        outlines = ''
        for lines in gff_data.readlines():
            gfflines = lines.strip().split()
            if gfflines[2] == 'repeat_region':
                if outlines != '':
                    print(outlines , file=f_out)
                    print('Repeats:', file=f_out)
                    print('', file=f_out)
                contig = gfflines[0]
                crispr = gfflines[8].split(';')[0].split('=')[1]
                print('Sequence' + '\t' + '\'' + contig + '\'' + '\t' + '(' + seq_len[contig] + ' bp)', file=f_out)
                print('CRISPR' + '\t' + seq_data[crispr][2] + '\t' + 'Range: ' + gfflines[3] + ' - ' + gfflines[4], file=f_out)
                print('POSITION\tREPEAT\tSPACER', file=f_out)

            if 'REPEAT' in gfflines[8]:
                outlines = ''
                outlines = outlines + gfflines[3]
                seqlines = gfflines[8].split('Note=')[1].split(';')[0]
                outlines = outlines + '\t' + seqlines
                repeat_len = len(seqlines)
            if 'SPACER' in gfflines[8]:
                seqlines = gfflines[8].split('Note=')[1].split(';')[0]
                spacer_len = len(gfflines[0]) 
                outlines = outlines + '\t' + seqlines + '\t' + '[ ' + str(repeat_len) + ', ' + str(spacer_len) + ' ]'
                print(outlines,file=f_out)
                outlines = ''
        # print(outlines + '\t[ ' + str(repeat_len) + ','  + 'NA ]',file=f_out)
        print(outlines , file=f_out)
        print('Repeats: ', file=f_out)
        print('',file=f_out)
        f_out.close()

    def parse_minced(self):
        file = open(os.path.join(self.out,'minced.out'), 'r')#runze.cai修改，2023.2.27
        cali_out = os.path.join(self.out,'MinCEDcalibrate')
        # Make dir
        try:
            os.mkdir(cali_out)
        except:
            shutil.rmtree(cali_out)
            os.mkdir(cali_out)
             
        random.seed(self.seed)
        crisprs = []
        for ll in file:
            # Record sequence accession
            if ll.startswith('Sequence'):
                sequence_current = re.sub("\'","",ll.strip().split()[1])
                sequence_current_len = re.sub('\(','',ll.strip().split()[2])
            # Create instance of CRISPR and add positions
            if ll.startswith('CRISPR'):
                cid = int(ll.split()[1])
                crispr_current = '{}_{}'.format(sequence_current,cid)
                crisp_tmp = CRISPR(cali_out,sequence_current, crispr_current)
                pos = re.sub('.*Range: ', '', ll)
                start = re.sub(' - .*', '', pos)
                end = re.sub('.* - ', '', pos)
                crisp_tmp.setPos(start, end, sequence_current_len)
            # Add Repeats and Spacers to the current instance
            if ll[:1].isdigit():
                lll = ll.split()
                if len(lll) == 7:
                    crisp_tmp.addRepeat(lll[1])
                    crisp_tmp.addSpacer(lll[2])
                if len(lll) == 2:
                    crisp_tmp.addRepeat(lll[1])
            # Save the instance
            if ll.startswith('Repeats'):
                crisp_tmp.getConsensus()
                crisp_tmp.check_Repeats()
                crisp_tmp.check_Spacers()
                crisp_tmp.stats()
                crisprs.append(crisp_tmp)

        file.close()

        ## Clean cali_out
        shutil.rmtree(cali_out)
        
        self.crisprs = crisprs

    def parse_rnafold(self):
        rnafold_out = os.path.join(self.out,'RNAfold')
        # Make dir
        try:
            os.mkdir(rnafold_out)
        except:
            shutil.rmtree(rnafold_out)
            os.mkdir(rnafold_out)
        
        for i in range(len(self.crisprs)):
            crisp = self.crisprs[i]
            rnastr = RNAfold(rnafold_out,crisp.cons)
            rnastr.runRNAfold()
            crisp.name = self.name
            ## Feature part 1: count
            crisp.GC = rnastr.GC
            crisp.nLoop = rnastr.nLoop
            crisp.nBubble = rnastr.nBubble
            crisp.nStem = rnastr.nStem
            crisp.nLink = rnastr.nLink
            crisp.nGC = rnastr.nGC
            crisp.nAU = rnastr.nAU
            crisp.nGU = rnastr.nGU
            crisp.nPair = rnastr.nPair
            crisp.stemScore = rnastr.stemScore
            ## Feature part 2: length & ratio
            crisp.maxLoopLen = rnastr.maxLoopLen
            crisp.sumLoopLen = rnastr.sumLoopLen
            crisp.maxBubbleLen = rnastr.maxBubbleLen
            crisp.sumBubbleLen = rnastr.sumBubbleLen
            crisp.maxStemLen = rnastr.maxStemLen
            crisp.sumStemLen = rnastr.sumStemLen
            crisp.maxLinkLen = rnastr.maxLinkLen
            crisp.sumLinkLen = rnastr.sumLinkLen
            crisp.LoopPercent = rnastr.LoopPercent
            crisp.BubblePercent = rnastr.BubblePercent
            crisp.StemPercent = rnastr.StemPercent
            crisp.LinkPercent = rnastr.LinkPercent
            crisp.LinkLoopRatio = rnastr.LinkLoopRatio
            crisp.StemLoopRatio = rnastr.StemLoopRatio
            crisp.LinkStemRatio = rnastr.LinkStemRatio
            ## Feature part 3: free energy
            crisp.MFE = rnastr.MFE
            crisp.MEA = rnastr.MEA
            crisp.EFE = rnastr.EFE
            crisp.CFE = rnastr.CFE
            crisp.Distance = rnastr.Distance
            crisp.Diversity = rnastr.Diversity
            crisp.Accuracy = rnastr.Accuracy
            crisp.MFEstr = rnastr.MFEstr
            crisp.MEAstr = rnastr.MEAstr
            
            self.crisprs[i] = crisp
        
        ## clean rnafold_out
        shutil.rmtree(rnafold_out)   
        
    def write_crisprs(self):
       
        header = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Name','Contig',
                    'CRISPR','DRrepeat', 'Calibrated','CalibLen',
                    'ContigLen','Start','End',
                    'N_repeats','Repeat_len','Spacer_len_avg','Spacer_len_sem',
                    'Repeat_identity','GC',
                    'nLoop','nBubble','nStem','nLink',
                    'nGC','nAU','nGU','nPair','stemScore',
                    'maxLoopLen','sumLoopLen',
                    'maxBubbleLen','sumBubbleLen',
                    'maxStemLen','sumStemLen',
                    'maxLinkLen','sumLinkLen',
                    'LoopPercent','BubblePercent','StemPercent','LinkPercent',
                    'LinkLoopRatio','StemLoopRatio','LinkStemRatio',
                    'MFE','MEA','EFE','CFE',
                    'Distance','Diversity','Accuracy',
                    'MFEstr','MEAstr')
       
        def write_crispr(handle, cris):
            handle.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(cris.name,cris.sequence,
                    cris.crispr,cris.cons,cris.calibrated,cris.calibLen,
                    cris.sequence_len,cris.start,cris.end,
                    len(cris.repeats),cris.repeat_len,cris.spacer_len,cris.spacer_sem,
                    cris.mean_match_rate,cris.GC,
                    cris.nLoop,cris.nBubble,cris.nStem,cris.nLink,
                    cris.nGC,cris.nAU,cris.nGU,cris.nPair,cris.stemScore,
                    cris.maxLoopLen,cris.sumLoopLen,
                    cris.maxBubbleLen,cris.sumBubbleLen,
                    cris.maxStemLen,cris.sumStemLen,
                    cris.maxLinkLen,cris.sumLinkLen,
                    cris.LoopPercent,cris.BubblePercent,cris.StemPercent,cris.LinkPercent,
                    cris.LinkLoopRatio,cris.StemLoopRatio,cris.LinkStemRatio,
                    cris.MFE,cris.MEA,cris.EFE,cris.CFE,
                    cris.Distance,cris.Diversity,cris.Accuracy,
                    cris.MFEstr,cris.MEAstr))
            
        f = open(self.out+'crisprs_calib.tab', 'w')
        f.write(header)
        for crisp in self.crisprs:
            write_crispr(f, crisp)
        f.close()

    def write_spacers(self):
        
        if len(self.crisprs) > 0:
            os.mkdir(self.out+'spacers')
            for crisp in self.crisprs:
                f = open(self.out+'spacers/{}.fa'.format(crisp.crispr), 'w')
                n = 0
                for sq in crisp.spacers:
                    n += 1
                    f.write('>{}:{}\n'.format(crisp.crispr, n))
                    f.write('{}\n'.format(sq))

                f.close()