import os,sys,re
import pandas as pd
import shutil

class CRISPRCas(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)
    
    def get_candidates(self):

        def read_and_filter_data(file_path, columns):
            with open(file_path, 'r') as file:
                df = pd.read_csv(file, sep='\t')
            return df[columns]

        def parse_tab_certain(sub_all_sample, genome):
            title = ['Contig', 'Operon', 'Prediction', 'CRISPRs']
            co = pd.DataFrame(columns=title)
            data_files = []

            if os.path.exists(os.path.join(sub_all_sample, 'CRISPR_Cas.tab')):
                data_files.append(os.path.join(sub_all_sample, 'CRISPR_Cas.tab'))

            if os.path.exists(os.path.join(sub_all_sample, 'CRISPR_Cas_putative.tab')):
                data_files.append(os.path.join(sub_all_sample, 'CRISPR_Cas_putative.tab'))

            for data_file in data_files:
                data = read_and_filter_data(data_file, title)
                co = pd.concat([co, data])

            genome_dict = {}
            pid_dict = {}
            co['Genome'] = genome

            for index, row in co.iterrows():
                CRISPRid = re.findall(r'\[(.*?)\]', row.CRISPRs)[0].split(',')
                for p in CRISPRid:
                    p = p.replace("'", "")
                    typetag = row.Prediction
                    plevel = 'certain'
                    if 'Unknown' in typetag:
                        plevel = 'putative'
                    if 'Hy' in typetag:##hybrid typetag 20230713
                        plevel = 'hybrid'
                    typelevel = [typetag, plevel]

                    pid_dict[p] = typelevel
                genome_dict[row.Genome] = pid_dict

            return genome_dict

        def ccminer_out(input_data,sample_data,sample,output_path,out_tab):

            try:
                os.mkdir(output_path + '/' + sample)
            except:
                None

            CRISPR_id = list(set(sample_data.CRISPR))
            grouped = sample_data.groupby(sample_data.CRISPR)
            pid_list = open(output_path + '/' + sample + '/ORFid.list','w')
            
            for id in CRISPR_id:
                df = grouped.get_group(id)
                contig = list(df.Contig)[0]
                temp_tab = pd.DataFrame(columns=['DB_name','Sample','Contig','CRISPR','DRrepeat','Remark','Strand','Rank','Prolen','Protein','HMM','typetag','level','pro_seq','gene_seq'])
                ORF_in_cctyper = False
                
                if sample in input_data.keys() and id in input_data[sample].keys():
                    ORF_in_cctyper = True
                
                pid = list(df.Protein)
                rank = list(df.Rank)
                prolen = list(df.Prolen)
                hmm = list(df.HMM)
                strand = list(df.Strand)
                remark = list(df.Remark)
                drrepeat = list(df.DRrepeat)
                
                num_split_temp = []
                ORF_is_typeII = False
                for num in range(len(rank)-1):
                    if re.search(pattern='Cas9|Cas12|Cas13',string=hmm[num]) != None or re.search(pattern='Cas9|Cas12|Cas13',string=hmm[num]) != None:
                        ORF_is_typeII = True
                    if abs(rank[num]-rank[num+1]) == 1 and strand[num] == strand[num+1] and self.split_prolen_lower < prolen[num] < self.split_prolen_upper and self.split_prolen_lower < prolen[num+1] < self.split_prolen_upper:
                        num_split_temp.extend((num,num+1))

                for num in range(len(rank)):
                    if ORF_in_cctyper:
                        temp_tab.loc[len(temp_tab.index)]=[self.database_name,sample,contig,id,drrepeat[num],remark[num],strand[num],rank[num],prolen[num],pid[num],hmm[num],input_data[sample][id][0],input_data[sample][id][1],'','']
                        continue
                    if ORF_is_typeII:
                        temp_tab.loc[len(temp_tab.index)]=[self.database_name,sample,contig,id,drrepeat[num],remark[num],strand[num],rank[num],prolen[num],pid[num],hmm[num],'Unknown','putative','','']
                        continue
                    if num in num_split_temp:
                        temp_tab.loc[len(temp_tab.index)]=[self.database_name,sample,contig,id,drrepeat[num],remark[num],strand[num],rank[num],prolen[num],pid[num],hmm[num],'Unknown','split','','']
                        continue
                    if prolen[num] > self.large_prolen:
                        temp_tab.loc[len(temp_tab.index)]=[self.database_name,sample,contig,id,drrepeat[num],remark[num],strand[num],rank[num],prolen[num],pid[num],hmm[num],'Unknown','large','','']
                        continue
                    if prolen[num] > self.need_prolen:
                        temp_tab.loc[len(temp_tab.index)]=[self.database_name,sample,contig,id,drrepeat[num],remark[num],strand[num],rank[num],prolen[num],pid[num],hmm[num],'Unknown','primary','','']
                        continue
                    temp_tab.loc[len(temp_tab.index)]=[self.database_name,sample,contig,id,drrepeat[num],remark[num],strand[num],rank[num],prolen[num],pid[num],hmm[num],'Unknown','other','','']


                temp_seq = []
                temp_gene = []
                pid = list(df.Protein)
                pid_list = open(output_path + '/' + sample + '/ORFid.list','w')
                for p_id in pid:
                    print(p_id,file=pid_list)
                pid_list.close()
                id_tmp = id.replace('|','_')
                #protein seq
                os.system("seqtk subseq " + self.cctyper_path + '/proteins.faa ' +  output_path + '/' + sample + '/ORFid.list' + ' > ' + output_path + '/' + sample + '/' + id_tmp + '.faa')
                # gene seq
                os.system("seqtk subseq " + self.cctyper_path + '/genes.cds ' +  output_path + '/' + sample + '/ORFid.list' + ' > ' + output_path + '/' + sample + '/' + id_tmp + '.cds')        
                id_array_pro = open(output_path + '/' + sample + '/' + id_tmp + '.faa')
                id_array_gene = open(output_path + '/' + sample + '/' + id_tmp + '.cds')
                for l in id_array_pro.readlines():
                    if l.startswith('>'):
                        continue
                    else:
                        temp_seq.append(l.strip())

                temp_tab['pro_seq'] = temp_seq

                for l in id_array_gene.readlines():
                    if l.startswith('>'):
                        continue
                    else:
                        temp_gene.append(l.strip())
                temp_tab['gene_seq'] = temp_gene
                
                temp_tab = temp_tab[['DB_name','Sample','CRISPR','DRrepeat','Remark','Strand','Rank','Prolen','Protein','HMM','typetag','level','pro_seq','gene_seq']]
                temp_tab = temp_tab[temp_tab['Rank'] <= self.need_rank]
                
                out_tab = pd.concat([out_tab,temp_tab])
            os.system('rm -rf ' + output_path + '/' + sample)
            # shutil.rmtree(output_path + '/' + sample)
            return(out_tab)
            
        if self.nearcrispr.empty:
            outcols = ['DB_name','Sample','CRISPR','DRrepeat','Remark','Strand','Rank','Prolen',
                    'Protein','HMM','typetag','level','pro_seq','gene_seq']
            self.out_tab = pd.DataFrame(columns=outcols)
            self.out_tab.to_csv(os.path.join(self.out,'out.tab'),sep='\t',index=False)
        else:
            outcols = ['DB_name','Sample','CRISPR','DRrepeat','Remark','Strand','Rank','Prolen',
                    'Protein','HMM','typetag','level','pro_seq','gene_seq']
            out_tab = pd.DataFrame(columns=outcols)
            cctyper_dict = parse_tab_certain(self.cctyper_path,self.name)
            self.out_tab = ccminer_out(cctyper_dict,self.nearcrispr,self.name,self.out,out_tab)
        self.out_tab.to_csv(os.path.join(self.out,'out.tab'),sep='\t',index=False)
