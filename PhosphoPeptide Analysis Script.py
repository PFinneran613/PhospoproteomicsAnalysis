# -*- coding: utf-8 -*-
"""
Author: Patrick Finneran


Data analysis in:
"Bimodal evolution of Src and Abl kinase substrate specificity revealed using mammalian cell extract as substrate pool"

"""
import pandas as pd
import numpy as np
from classes.phosphoproteomes import ControlProteome, BackgroundProteome,PhosphoProteome,SwissProt,CompareProteomes,DataBaseProteome,Seq_List_Proteome
import seaborn as sns
import matplotlib.ticker as mtick
    
    
    

def MinMax4Heatmap(Proteomes,p=0.05,bias=None):
    #Finds the min and max value so all heatmaps have the same scale
    AminoAcids = ['I', 'Q', 'M', 'D', 'L', 'S', 'T', 'W', 'G','H',
                  'A','E', 'P', 'Y', 'V', 'K','C', 'R', 'N', 'F']
    Positions = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
    a = p/(14*20)
    sig = -np.log10(a/(1-a))
    max_val = ('','','',0)
    min_val = ('','','',100)
    kinases = ['Abl','Src','Anc AS','Anc S1','Anc AST']
    for name in kinases:
        if bias == None:
            PWM = Proteomes[name].PWMs['ln(Enrichment)']
            sigPWM = Proteomes[name].PWMs['LogsOdd']
        if bias == 'length':
            PWM = Proteomes[name].PWMs['ln(Enrichment) LB']
            sigPWM = Proteomes[name].PWMs['LogsOdd LB']
        for pos in Positions:
            for aa in AminoAcids:
                if sigPWM[pos][aa] >= sig:
                    if abs(PWM[pos][aa]) > max_val[3]:
                        max_val = (name,pos,aa,abs(PWM[pos][aa]))
                    elif abs(PWM[pos][aa]) < min_val[3]:
                        min_val = (name,pos,aa,abs(PWM[pos][aa]))
    return max_val,min_val

def Max_Logo(Proteomes,rc):
    #Finds the maximum value so all Logos have the sames scale.
    AminoAcids = ['I', 'Q', 'M', 'D', 'L', 'S', 'T', 'W', 'G','H',
                  'A','E', 'P', 'Y', 'V', 'K','C', 'R', 'N', 'F']
    Positions = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
    cutoff = np.log(rc)
    sums = []
    
    for name in ['Abl','Anc AS','Anc AST','Anc S1','Src']:
        for pos in Positions:
            s2 = 0
            for aa in AminoAcids:
                s = Proteomes[name].PWMs['ln(Enrichment)'][pos][aa]
                if s > cutoff:
                    s2 += s
            sums.append(s2)
    max_val_Logo = max(sums)
    return max_val_Logo


def load_data(p=0.05,rc=1.3):
    # Load Background
    BgData = pd.read_csv('data/BackgroundUnenrichedProteomics.tsv',sep='\t')
    Bg = BackgroundProteome(BgData,'MS Data Background')
    # Load Control
    ControlData = pd.read_csv('data/ControlPhosphoproteomicsData.csv')
    Control = ControlProteome(ControlData)
    # Load Data Sets
    Data = pd.read_csv('data/PhosphoproteomicsData.csv')
    kinases = ['Abl','Anc AS','Anc AST','Anc S1','Src']
    Proteomes = {}
    for name in kinases:
        print('Importing {0}'.format(name))
        Proteomes[name] = PhosphoProteome(name,Data,Bg,Control)
    #Load Short Time Point
    Short_Data = pd.read_csv('data/PhosphoproteomicsSrcShortTimePoint.csv',low_memory=False)
    short_kinases = ['Src 10 Minute']
    for name in short_kinases:
        print('Importing {0}'.format(name))
        Proteomes[name] = PhosphoProteome(name,Short_Data,Bg,Control)
    # Load Database Standards
    sp_Bg = SwissProt()
    DB_data = pd.read_csv('data/Kinase_Substrate_Dataset',sep='\t',skiprows=3)
    db_kinases = ['Abl','Src']
    for name in db_kinases:
        print('Importing {0} PSP'.format(name))
        Proteomes[name+' PSP'] = DataBaseProteome(name,DB_data,sp_Bg)
        
    #Import Manually Curated Data Sets
    manual_data_sets = ['HAKA_MS']
    for name in manual_data_sets:
        seq_fname = 'data/Other Data/' + name + '_Peptides.txt'
        prot_fname = 'data/Other Data/' + name + '_Proteins.txt'
        Proteomes[name] = Seq_List_Proteome(name,seq_fname,prot_fname,sp_Bg)
    # Determine Max and Min values for Logos and Heatmaps
    max_val,min_val = MinMax4Heatmap(Proteomes,p=p)
    max_val_Logo = Max_Logo(Proteomes,rc)
    # Generate figures for each phosphoproteome
    for name in Proteomes.keys():
        Proteomes[name].gen_HeatMap(p=p,max_val=max_val[3],cutoff=min_val[3],ratio_cutoff=rc)
        Proteomes[name].prep_Logo(ratio_cutoff=rc, max_val=max_val_Logo)
        Proteomes[name].gen_BgVenn(venn_type='Seqs')
        Proteomes[name].gen_BgVenn(venn_type='Proteins')
    return Proteomes
        
def compare_protomes(Proteomes,p=0.05,rc=1.3):
    # Initialize Comparisons
    Compare = CompareProteomes(Proteomes,p=p,ratio_cutoff=rc)
    
    # Generate Entropy Figure
    kinases = ['Anc AST','Anc AS','Abl','Anc S1','Src','Background']
    Compare.gen_entropy(kinases,fname='Entropy Figure')
    
    # Venn Diagram Comparison with HAKA_MS Method
    Compare.gen_venn(['Abl','HAKA_MS'],'Seqs')
    Compare.gen_venn(['Abl','HAKA_MS'],'Proteins')
    
    # Determine Max and Min values for Heatmaps
    max_val,min_val = MinMax4Heatmap(Proteomes,p=p)
    
    # Position Dependent Heatmaps
    AA_Pos = [(['E','D','A','P'],-2),(['G','A'],+1),(['E','A'],+2),(['F','I','P'],+3),(['K','G','F'],+4),(['K','A'],+5)]
    for aa,pos in AA_Pos:
        Compare.pos_heatmap(pos, aa,max_val=max_val[3])
    
    
    
    # Graph Certain Residues and Positions
    graph_vals = [('P',3),('P',-2),('A',2),('A',1),('A',-2),('E',-1),('E',2),
          ('F',3),('I',3),('G',4),('V',2),('D',-2),('P',-1),('AV',2),('G',1),
          ('N',-2),('T',-1),('F',3),('I',-1),('L',-1),('V',-1),('IV',-1),
          ('E',-3),('S',-2),('DE',-3),('D',-3),('S',2)]
    for aa,pos in graph_vals:
        Compare.graph_data(aa,pos,PWM_T='ln(Enrichment)') #ln(Enrichment) is actually log base 10
    
    
    
    # Venn Diagrams of Certain Data Sets
    kinase_pairs = [['Abl','Abl PSP'], ['Src', 'Src PSP'],['Src','Abl','Anc AS'],['Src','Anc S1','Anc AS']]
    for kp in kinase_pairs:
        Compare.gen_venn(kp,'Seqs')
        Compare.gen_venn(kp,'Proteins')
    
    # Heatmaps for sequences unique to Src or Abl
    vennpair_Hmaps = ['Src', 'Abl']
    Compare.venn_heatmaps(vennpair_Hmaps,max_val[3])
    return Compare

def main():
    Proteomes = load_data(p=0.05,rc=1.3)
    Compare = compare_protomes(Proteomes)
    return Proteomes, Compare


Proteomes, Compare = main()
