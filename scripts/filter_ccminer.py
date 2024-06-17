import os,sys,re
import pandas as pd
import shutil
import argparse

########## Arguments ##########
def get_parameters():
    ap = argparse.ArgumentParser()
    # Required
    ap.add_argument('--cctyper_path', help='Input cctyper path.')
    ap.add_argument('--ccminer_path', help='output ccminer path.')

    # Optional
    ap.add_argument('--gname',type=str,default=None,help='Sample name/Genome name.')
    ap.add_argument('--dbname',type=str,default=None,help='Database name.')
    ap.add_argument('--sele_rank',type=int,default=3,help='Expected maximum distance between effector-proteins and CRISPR array.')
    ap.add_argument('--sele_prolen',type=int,default=700,help='Expected minimum size of effctor-proteins.')
    ap.add_argument('--large_prolen',type=int,default=2000,help='Protein size of (too) large proteins.')
    ap.add_argument('--split_prolen_min',type=int,default=700,help='Minimum protein size of suspected splited protein segments.')
    ap.add_argument('--split_prolen_max',type=int,default=1500,help='Maximum protein size of suspected splited protein segments.')
    
    args = ap.parse_args()
    if (args.cctyper_path == None or args.ccminer_path == None):
        raise Exception('Error: cctyper_path and ccminer_path are required!')
    
    return args

########## Functions ##########
def init_outtab():
    outcols = ['DB_name','Sample','CRISPR','DRrepeat','Remark','Strand','Rank','Prolen',
                    'Protein','HMM','typetag','level','pro_seq','gene_seq']
    out_tab = pd.DataFrame(columns=outcols)
    return out_tab

def read_and_filter_data(file_path, columns):
    with open(file_path, 'r') as file:
        df = pd.read_csv(file, sep='\t')
    return df[columns]

def parse_tab_certain(cctyper_path, gname):
    title = ['Contig', 'Operon', 'Prediction', 'CRISPRs']
    co = pd.DataFrame(columns=title)
    data_files = []

    if os.path.exists(os.path.join(cctyper_path, 'CRISPR_Cas.tab')):
        data_files.append(os.path.join(cctyper_path, 'CRISPR_Cas.tab'))

    if os.path.exists(os.path.join(cctyper_path, 'CRISPR_Cas_putative.tab')):
        data_files.append(os.path.join(cctyper_path, 'CRISPR_Cas_putative.tab'))

    for data_file in data_files:
        data = read_and_filter_data(data_file, title)
        co = pd.concat([co, data])

    genome_dict = {}
    pid_dict = {}
    co['Genome'] = gname

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

def ccminer_out(input_data,nearcrispr,out_tab,args):

    CRISPR_id = list(set(nearcrispr.CRISPR))
    grouped = nearcrispr.groupby(nearcrispr.CRISPR)
    
    for id in CRISPR_id:
        df = grouped.get_group(id)
        contig = list(df.Contig)[0]
        temp_tab = pd.DataFrame(columns=['DB_name','Sample','Contig','CRISPR','DRrepeat','Remark','Strand','Rank','Prolen','Protein','HMM','typetag','level','pro_seq','gene_seq'])
        ORF_in_cctyper = False
        
        if args.gname in input_data.keys() and id in input_data[args.gname].keys():
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
            if abs(rank[num]-rank[num+1]) == 1 and strand[num] == strand[num+1] and args.split_prolen_min < prolen[num] < args.split_prolen_max and args.split_prolen_min < prolen[num+1] < args.split_prolen_max:
                num_split_temp.extend((num,num+1))

        for num in range(len(rank)):
            if ORF_in_cctyper:
                temp_tab.loc[len(temp_tab.index)]=[args.dbname,args.gname,contig,id,drrepeat[num],remark[num],strand[num],rank[num],prolen[num],pid[num],hmm[num],input_data[args.gname][id][0],input_data[args.gname][id][1],'','']
                continue
            if ORF_is_typeII:
                temp_tab.loc[len(temp_tab.index)]=[args.dbname,args.gname,contig,id,drrepeat[num],remark[num],strand[num],rank[num],prolen[num],pid[num],hmm[num],'Unknown','putative','','']
                continue
            if num in num_split_temp:
                temp_tab.loc[len(temp_tab.index)]=[args.dbname,args.gname,contig,id,drrepeat[num],remark[num],strand[num],rank[num],prolen[num],pid[num],hmm[num],'Unknown','split','','']
                continue
            if prolen[num] > args.large_prolen:
                temp_tab.loc[len(temp_tab.index)]=[args.dbname,args.gname,contig,id,drrepeat[num],remark[num],strand[num],rank[num],prolen[num],pid[num],hmm[num],'Unknown','large','','']
                continue
            if prolen[num] > args.sele_prolen:
                temp_tab.loc[len(temp_tab.index)]=[args.dbname,args.gname,contig,id,drrepeat[num],remark[num],strand[num],rank[num],prolen[num],pid[num],hmm[num],'Unknown','primary','','']
                continue
            temp_tab.loc[len(temp_tab.index)]=[args.dbname,args.gname,contig,id,drrepeat[num],remark[num],strand[num],rank[num],prolen[num],pid[num],hmm[num],'Unknown','other','','']


        temp_seq = []
        temp_gene = []
        os.makedirs(os.path.join(args.ccminer_path, args.gname), exist_ok=True)
        pid_list = open(args.ccminer_path + '/' + args.gname + '/ORFid.list','w')
        for p_id in pid:
            print(p_id,file=pid_list)
        pid_list.close()
        id_tmp = id.replace('|','_')
        #protein seq
        os.system("seqtk subseq " + args.cctyper_path + '/proteins.faa ' +  args.ccminer_path + '/' + args.gname + '/ORFid.list' + ' > ' + args.ccminer_path + '/' + args.gname + '/' + id_tmp + '.faa')
        # gene seq
        # os.system("seqtk subseq " + args.cctyper_path + '/genes.cds ' +  args.ccminer_path + '/' + args.gname + '/ORFid.list' + ' > ' + args.ccminer_path + '/' + args.gname + '/' + id_tmp + '.cds')        
        id_array_pro = open(args.ccminer_path + '/' + args.gname + '/' + id_tmp + '.faa')
        # id_array_gene = open(args.ccminer_path + '/' + args.gname + '/' + id_tmp + '.cds')
        for l in id_array_pro.readlines():
            if l.startswith('>'):
                continue
            else:
                temp_seq.append(l.strip())

        temp_tab['pro_seq'] = temp_seq

        # for l in id_array_gene.readlines():
        #     if l.startswith('>'):
        #         continue
        #     else:
        #         temp_gene.append(l.strip())
        # temp_tab['gene_seq'] = temp_gene
        temp_tab['gene_seq'] = 'NoData'
        
        temp_tab = temp_tab[['DB_name','Sample','CRISPR','DRrepeat','Remark','Strand','Rank','Prolen','Protein','HMM','typetag','level','pro_seq','gene_seq']]
        temp_tab = temp_tab[temp_tab['Rank'] <= args.sele_rank]
        
        out_tab = pd.concat([out_tab,temp_tab])

    return(out_tab)

def main():
    args = get_parameters()
    out_tab = init_outtab()
    cctyper_dict = parse_tab_certain(args.cctyper_path, args.gname)
    nearcrispr_file = os.path.join(args.ccminer_path, 'nearCRISPR.txt')
    if os.path.exists(nearcrispr_file):
        nearcrispr = pd.read_csv(nearcrispr_file, sep='\t')
        if nearcrispr.shape[0] > 0: ## is not empty
            out_tab = ccminer_out(cctyper_dict,nearcrispr,out_tab,args)
            out_tab.to_csv(os.path.join(args.ccminer_path, 'out.tab'), sep='\t', index=False)
        else:
            print('NearCRISPR is empty: %s' % args.ccminer_path)
    else:
        print('nearCRISPR.txt file is not exist: %s' % args.ccminer_path)
    
if __name__ == '__main__':
    
    main()
    
"""
## Run example
python filter_ccminer.py \
        --cctyper_path ./CGMR_part128.cctyper \
        --ccminer_path ./CGMR_part128.ccminer \
        --gname CGMR_part128 \
        --dbname CGMR
"""