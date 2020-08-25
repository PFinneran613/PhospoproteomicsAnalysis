# -*- coding: utf-8 -*-
"""
Author: Patrick Finneran

PhosphoProteome Classes used for data analysis in:
"Bimodal evolution of Src and Abl kinase substrate specificity revealed using mammalian cell extract as substrate pool"

"""
import pandas as pd
import re
import scipy.stats
import numpy as np
import json
import os
import seaborn as sns
import xlsxwriter

### ControlProteome ###
#      The ControlProteome is to account for phosphorylated tyrosines
#      found phosphorylated before additional kinases are added to the reaction.
class ControlProteome:
    
    def __init__(self,data):
        self.data = data
        self.Seqs = set()
        for i, row in data.iterrows():
            seq = row['Sequence window'].split(';')[0]
            x = int(((len(seq) - 1) / 2) - 7)
            seq = seq[x:-x]
            self.Seqs.add(seq)
        if not os.path.exists('Sequences'):
            os.makedirs('Sequences')
        with open('Sequences/Control Sequences.txt','w') as outfile:
            outfile.write('\n'.join(self.Seqs))


### BackgroundProteome ###
#      The BackgroundProteome is for creating the background distribution of amino acids
#      from a non-phosphorylated data set
class BackgroundProteome:
    
    def __init__(self,data,name):
        self.name = name
        self.data = data
        self.gen_sets()
        self.gen_PWM()
        self.prep_Logo()
        self.gen_heatmap()

    
    def gen_sets(self):
        self.Seqs = set()
        self.Proteins = set()
        self.Peptides = []
        check_ints = [['Unenriched'],[1,2,3]]
        for i, row in self.data.iterrows():
            c = False
            for a in check_ints[0]:
                for b in check_ints[1]:
                    if row['Intensity {0} {1}'.format(a,b)] > 0:
                        c = True
            if c:
                seq = row['Sequence']
                for i in range(len(seq)):
                    if seq[i] =='Y':
                        l = i
                        r = len(seq) - i
                        r = len(row['Sequence']) - i
                        try:
                            if l >= 7 and r > 7:
                                seq_window = seq[i-7:i+8]
                                e = 'Cat 1 Failed'
                            elif l < 7:
                                e = 'Cat 2 Failed'
                                sw = row['N-term cleavage window']
                                seq_window = sw[15+l-7:15+l+8]
                            elif r < 8:
                                sw = row['C-term cleavage window']
                                seq_window = sw[15-r-7:15-r+8]
                                e = 'Cat 3 Failed'
                            if len(seq_window) != 15:
                                print(e)
                                print('l: ',l)
                                print('r: ',r)
                                print(len(seq),seq)
                                print(seq_window)
                                print(sw)
                                return
                            if '_' not in seq_window:
                                self.Seqs.add(seq_window)
                                peptide = row['Sequence']
                                skip = 0
                                for i in range(len(peptide)-1):
                                    if peptide[i] == 'K' or peptide[i] == 'R':
                                        if self.check_site(peptide[i-1:i+2],proteolysis_type=3):
                                            skip = skip + 1
                                temp = {'Peptide':peptide,'Length':len(row['Sequence']),'SeqWindow':seq_window,'Skips':skip}
                                self.Peptides.append(temp)
                                prot = row['Proteins'].split(';')[0]
                                if len(prot)>0:
                                    self.Proteins.add(prot)
                        except:
                            pass
        self.Peptides = pd.DataFrame(self.Peptides)
        if not os.path.exists('Sequences'):
            os.makedirs('Sequences')
        with open('Sequences/Background Sequences.txt','w') as outfile:
            outfile.write('\n'.join(self.Seqs))
    
    def check_site(self,site,proteolysis_type):
        exclusions = []
        if proteolysis_type == 3:
            exclusions = ['CKY', 'DKD', 'CKH', 'CKD', 'KKR', 'RRH', 'RRR', 'CRK', 'DRD', 'RRF', 'KRR']
        if len(site) <= 2:
            return True
        if site[1] in 'KR':
            if site in exclusions:
                return False
            elif site[2] == 'P' and proteolysis_type >= 2:
                return False
            else:
                return True
        else:
            return False
                  
    def gen_PWM(self):
        AminoAcids = ['I', 'Q', 'M', 'D', 'L', 'S', 'T', 'W', 'G','H', 'A',
                      'E', 'P', 'Y', 'V', 'K','C', 'R', 'N', 'F']
        Positions = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
        N = len(self.Seqs)
        self.PWM = {}
        for pos in Positions:
            self.PWM[pos] = {}
            for aa in AminoAcids:
                motif = motif = ['.','.','.','.','.','.','.','Y','.','.','.','.','.','.','.']
                motif[7+pos] = aa
                Motif = re.compile(''.join(motif))
                self.PWM[pos][aa] = len(Motif.findall('\n'.join(self.Seqs))) / N
                
    def gen_heatmap(self):
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set(font_scale=1.5)
        x_vals = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
        y_vals = ['G', 'A', 'V', 'L', 'I', 'W', 'F', 'Y', 'P',  'M', 'C', 'S', 'T', 'N', 'Q', 'D', 'E', 'H', 'K', 'R']
        array_list = []
        i = 0
        for aa in y_vals:
            array_list.append([])
            for pos in x_vals:
                if pos == 0: #and aa == 'Y':
                    array_list[i].append(0)
                else:
                    array_list[i].append(self.PWM[pos][aa]) 
            i = i+1
        Array = np.array(array_list)
        self.PWM_array = Array
        fig, ax = plt.subplots(figsize=(6,6))
        ax = sns.heatmap(Array, xticklabels = x_vals, yticklabels = y_vals, linecolor='black',
                         linewidths = .5, square = True,cmap='Reds',
                         cbar_kws = {'label': 'Residue Frequency'},ax=ax)
        plt.xlabel('Position from P-Tyrosine',size=16)
        plt.ylabel('Amino Acids',size=16)
        plt.yticks(rotation=0,ha='center',size=13)
        plt.xticks(rotation=0,ha='center',size=13)
        ax.set_title('Residue Frequency\nin {0}'.format(self.name))
        fname1 = 'Backgrounds/Residue Frequency in {0}.png'.format(self.name)
        if not os.path.exists('Backgrounds'):
            os.makedirs('Backgrounds')
        plt.savefig(fname1, format='png', dpi=1200, bbox_inches='tight')
        plt.close()
        
    def prep_Logo(self):
        AminoAcids = ['I', 'Q', 'M', 'D', 'L', 'S', 'T', 'W', 'G','H',
                      'A', 'E', 'P', 'Y', 'V', 'K','C', 'R', 'N', 'F']
        Positions = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
        scores = {}
        for pos in Positions:
            scores[pos] = []
            for aa in AminoAcids:
                s = self.PWM[pos][aa]
                scores[pos].append(s)
        Scores = []
        for pos in Positions:
            Scores.append(scores[pos])
        if not os.path.exists('PhosphoProteomes/{0}'.format(self.name)):
            os.makedirs('PhosphoProteomes/{0}'.format(self.name))
        df = pd.DataFrame(Scores)
        df.T.to_csv('PhosphoProteomes/{0}/{0}_Scores.csv'.format(self.name),header=False,index=False)


class PhosphoProteome:
    
    def __init__(self,name,data,Bg,Control):
        self.name = name
        self.data = data
        self.Bg = Bg
        self.Control = Control
        self.gen_sets()
        self.gen_PWMs()
        self.gen_freq_array()
        if not os.path.exists('Sequences'):
            os.makedirs('Sequences')
        with open('Sequences/{} Sequences.txt'.format(self.name),'w') as outfile:
            outfile.write('\n'.join(self.Seqs))
        
        
    def gen_sets(self):
        self.Seqs = set()
        self.Proteins = set()
        self.Peptides = []
        for i, row in self.data.iterrows():
            c = False
            for trial in [1,2,3]:
                if row['Intensity {0} {1}'.format(self.name,trial)] > 0:
                    if row['Localization prob {0} {1}'.format(self.name,trial)] >= 0.7:
                        c = True
            seq = self.check_seq(row['Sequence window'])
            if seq is not None and c == True:
                self.Seqs.add(seq)
                skip = 0
                p = row['Phospho (STY) Probabilities']
                ignore = False
                peptide = ''
                for char in p:
                    if char == ')':
                        ignore = False
                    elif ignore:
                        pass
                    elif char == '(':
                        ignore = True
                    else:
                        peptide = peptide + char
                for i in range(len(peptide)-1):
                    if peptide[i] == 'K' or peptide[i] == 'R':
                        if self.check_site(peptide[i-1:i+2],proteolysis_type=3):
                            skip = skip + 1
                self.Peptides.append({'Peptide':peptide,'Length':len(peptide),'Skips':skip,'SeqWindow':seq})
                prot = row['Proteins'].split(';')[0]
                if len(prot)>0:
                    self.Proteins.add(prot)
        self.Peptides = pd.DataFrame(self.Peptides)
        if not os.path.exists('Sequences'):
            os.makedirs('Sequences')
        with open('Sequences/{0} Sequences.txt'.format(self.name),'w') as outfile:
            outfile.write('\n'.join(self.Seqs))

    def check_seq(self,seq):
        try:
            if len(seq) > 15:
                seq = seq.split(';')[0]
                x = int(((len(seq) - 1) / 2) - 7)
                seq = seq[x:-x]
        except:
            return None
        if seq == '' or seq == 'nan':
            return None
        elif seq[0] == '_' or seq[14] == '_' or seq[7] != 'Y':
            return None
        elif seq in self.Control.Seqs:
            return None
        else:
            return seq
    
    def check_site(self,site,proteolysis_type):
        exclusions = []
        if proteolysis_type == 3:
            exclusions = ['CKY', 'DKD', 'CKH', 'CKD', 'KKR', 'RRH', 'RRR', 'CRK', 'DRD', 'RRF', 'KRR']
        if len(site) <= 2:
            return True
        if site[1] in 'KR':
            if site in exclusions:
                return False
            elif site[2] == 'P' and proteolysis_type >= 2:
                return False
            else:
                return True
        else:
            return False 
    
    def gen_PWMs(self):
        AminoAcids = ['I', 'Q', 'M', 'D', 'L', 'S', 'T', 'W', 'G','H', 'A',
                      'E', 'P', 'Y', 'V', 'K','C', 'R', 'N', 'F']
        Positions = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
        N = len(self.Seqs)
        self.PWMs = {}
        PWM_types = ['Frequency', '-ln(p)','Enrichment Ratio',
                     'Percent Enrichment','ln(Enrichment)','LogsOdd']
        for t in PWM_types:
            self.PWMs[t] = {}
        for pos in Positions:
            for t in PWM_types:
                self.PWMs[t][pos] = {}
            for aa in AminoAcids:
                if pos == 0:
                    self.PWMs['Frequency'][pos][aa] = 0
                    self.PWMs['-ln(p)'][pos][aa] = 0
                    self.PWMs['Enrichment Ratio'][pos][aa] = 0
                    self.PWMs['Percent Enrichment'][pos][aa] = 0
                    self.PWMs['ln(Enrichment)'][pos][aa] = 0
                    self.PWMs['LogsOdd'][pos][aa] = 0
                else:
                    p = self.Bg.PWM[pos][aa]
                    motif = motif = ['.','.','.','.','.','.','.','Y','.','.','.','.','.','.','.']
                    motif[7+pos] = aa
                    Motif = re.compile(''.join(motif))
                    count = len(Motif.findall('\n'.join(self.Seqs)))
                    self.PWMs['Frequency'][pos][aa] = count / N
                    self.PWMs['-ln(p)'][pos][aa] = -np.log(self.PWMs['Frequency'][pos][aa])
                    if p == 0:
                        self.PWMs['Enrichment Ratio'][pos][aa] = 0
                        self.PWMs['Percent Enrichment'][pos][aa] = 0
                    else:
                        self.PWMs['Enrichment Ratio'][pos][aa] = self.PWMs['Frequency'][pos][aa] / p
                        self.PWMs['Percent Enrichment'][pos][aa] = (self.PWMs['Frequency'][pos][aa] - p) / p
                    self.PWMs['ln(Enrichment)'][pos][aa] = np.log(self.PWMs['Enrichment Ratio'][pos][aa])
                    
                    
                    tail1 = scipy.stats.binom_test(count, N, p, alternative='greater')
                    tail2 = scipy.stats.binom_test(count, N, p, alternative='less')
                    self.PWMs['LogsOdd'][pos][aa] = - np.log10(tail1/tail2)
                    


            

    def gen_HeatMap(self,p=0.05,sig=True,max_val=0.45,cutoff=0.0725,bias=None,ratio_cutoff=None):   
        import seaborn as sns
        from matplotlib import pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        import os
        sym = '■'
        rc = cutoff
        
        if sig == True:
            a = p/(14*20)
            sig = -np.log10(a/(1-a))
        elif sig == False:
            sig = 0
        PWM = self.PWMs['ln(Enrichment)']
        sigPWM = self.PWMs['LogsOdd']
        title = 'Residue Preference in {0}\np={1}'.format(self.name,p)
        fname = 'Residue Preference in {0} with {1} cutoff.png'.format(self.name,p)
        if ratio_cutoff != None:
            rc = np.log(ratio_cutoff)
            cutoff = rc
            title = title + '\nCutoff Ratio of {0}'.format(ratio_cutoff)
            fname = fname[:-4] + 'Cutoff Ratio of {0}.png'.format(ratio_cutoff)
        if ratio_cutoff == 0:
            rc = 0
            cutoff=0.0725
            
        Mark = 'Sig'
        Mask = False
        wRange = cutoff/max_val*100
        wRange = int(wRange*10)
        cmap_default = sns.color_palette("RdBu_r",1000)
        cmap_new = []
        blues = cmap_default[100:450]
        reds = cmap_default[550:900]
        b_range = int((1000 - wRange)/2)
        r_range = b_range
        interval = len(blues) / b_range
        for i in range(b_range):
            color = blues[int(interval*i)]
            cmap_new.append(color)
        for i in range(wRange):
            cmap_new.append((1,1,1))
        for i in range(r_range):
            color = reds[int(interval*i)]
            cmap_new.append(color)
        cmap_final = []
        cmap_final = cmap_new
        Custom_cmap = LinearSegmentedColormap.from_list(name='Custom_CMap', colors =cmap_final, N=len(cmap_new))
        sns.set(font_scale=1.5)
        x_vals = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
        y_vals = ['G', 'A', 'V', 'L', 'I', 'W', 'F', 'Y', 'P',  'M', 'C', 'S', 'T', 'N', 'Q', 'D', 'E', 'H', 'K', 'R']
        array_list = []
        annot_list = []
        c_list = []
        i = 0
        for aa in y_vals:
            array_list.append([])
            annot_list.append([])
            c_list.append([])
            for pos in x_vals:
                if pos == 0: #and aa == 'Y':
                    array_list[i].append(0)
                    annot_list[i].append('')
                    c_list[i].append(False)
                else:
                    c_list[i].append(True)
                    if (sigPWM[pos][aa] <sig) and (sigPWM[pos][aa]>-sig):
                        if PWM[pos][aa] > cutoff or PWM[pos][aa] < -cutoff:
                            if Mark == 'Insig':
                                annot_list[i].append(sym)
                            else:
                                annot_list[i].append('')
                        else:
                            annot_list[i].append('')
                        if Mask:
                            array_list[i].append(0)
                        else:
                            array_list[i].append(PWM[pos][aa])
                    else:
                        array_list[i].append(PWM[pos][aa])
                        if Mark == 'Sig':
                            if ratio_cutoff == None:
                                annot_list[i].append(sym)
                            elif abs(PWM[pos][aa]) >= rc:
                                annot_list[i].append(sym)
                            else:
                                annot_list[i].append('')   
                        else:
                            annot_list[i].append('')    
            i = i+1
        Array = np.array(array_list)
        self.PWM_array_lnE = Array
        Annots = np.array(annot_list)
        Center = np.array(c_list)
        fig, ax = plt.subplots(figsize=(6,6))
        ax = sns.heatmap(Array, xticklabels = x_vals, yticklabels = y_vals , vmax=max_val, vmin=-max_val, linecolor='black',
                         linewidths = .5, cmap=Custom_cmap, square = True, annot=Annots,
                         annot_kws={'ha':'center','va':'center',"size": 10,'fontname':'Arial','color':'black','fontdict':{'color':'black'}},#'color':'#929591'
                         fmt='', cbar_kws = {'label': 'Disfavored  -  Enrichment  -  Favored'},ax=ax)
        sns.heatmap(Array,cmap=['#929591','#929591'],mask=Center,cbar=False,annot=False,linewidth=1, linecolor='black', xticklabels = x_vals, yticklabels = y_vals)
        plt.xlabel('Position from P-Tyrosine',size=16)
        plt.ylabel('Amino Acids',size=16)
        plt.yticks(rotation=0,ha='center',size=13)
        plt.xticks(rotation=0,ha='center',size=13)
        ax.set_title(title,size=20)
        if not os.path.exists('PhosphoProteomes/{0}'.format(self.name)):
            os.makedirs('PhosphoProteomes/{0}'.format(self.name))
        fname1 = 'PhosphoProteomes/{0}/{1}'.format(self.name,fname)
        plt.savefig(fname1, format='png', dpi=600, bbox_inches='tight')
        fname1b = 'PhosphoProteomes/{0}/{1}'.format(self.name,fname)[:-4]+'.eps'
        plt.savefig(fname1b, format='eps', dpi=1200, bbox_inches='tight')
        plt.close()
        
        fig, ax = plt.subplots(figsize=(6,6))
        ax = sns.heatmap(Array, xticklabels = x_vals, yticklabels = y_vals , vmax=max_val, vmin=-max_val, linecolor='black',
                         linewidths = .5, cmap=Custom_cmap, square = True, annot=Annots,
                         annot_kws={'ha':'center','va':'center',"size": 10,'fontname':'Arial','color':'black','fontdict':{'color':'black'}},#'color':'#929591'
                         fmt='', cbar=False,ax=ax)
        sns.heatmap(Array,cmap=['#929591','#929591'],mask=Center,cbar=False,annot=False,linewidth=1, linecolor='black', xticklabels = x_vals, yticklabels = y_vals)
        plt.xlabel('Position from P-Tyrosine',size=16)
        plt.ylabel('Amino Acids',size=16)
        plt.yticks(rotation=0,ha='center',size=13)
        plt.xticks(rotation=0,ha='center',size=13)
        fname2 = 'PhosphoProteomes/{0}/{1}'.format(self.name,fname)[:-4] + 'For Combo.eps'
        plt.savefig(fname2, format='eps', dpi=1200, bbox_inches='tight')
        plt.close()
        
        
        
        
    def gen_BgVenn(self,venn_type='Seqs'):
        from matplotlib_venn import venn2,venn2_circles
        from matplotlib import pyplot as plt
        if venn_type == 'Seqs':
            label1 = '{0}\nSubstrates\n{1}'.format(self.name,len(self.Seqs))
            label2 = 'Background\nSubstrates\n{0}'.format(len(self.Bg.Seqs))
            labels = [label1,label2]
            sets = [self.Seqs, self.Bg.Seqs]
            title = '{0} Substrates\ncompared to Background'.format(self.name)
            fname = 'PhosphoProteomes/{0}/{0} Sequences Compared to Background Venn.png'.format(self.name)
        if venn_type == 'Proteins':
            label1 = '{0}\nProteins\n{1}'.format(self.name,len(self.Proteins))
            label2 = 'Background\nProteins\n{0}'.format(len(self.Bg.Proteins))
            labels = [label1,label2]
            sets = [self.Proteins, self.Bg.Proteins]
            title = '{0} Proteins\ncompared to Background'.format(self.name)
            fname = 'PhosphoProteomes/{0}/{0} Proteins Compared to Background Venn.png'.format(self.name)
        plt.figure()
        ax = plt.gca()
        v = venn2(sets, set_labels=labels)
        try:
            v.get_patch_by_id('10').set_color('red')
            v.get_patch_by_id('01').set_color('blue')
            v.get_patch_by_id('11').set_color('magenta')
        except:
            pass
        c = venn2_circles(sets, linestyle='solid')
        c[0].set_edgecolor('red')
        c[1].set_edgecolor('blue')
        plt.title(title)
        h, l = [],[]
        Sets = ['10','11','01']
        for i in Sets:
            # remove label by setting them to empty string:
            try:
                count = v.get_label_by_id(i).get_text()
                v.get_label_by_id(i).set_text("")
                # append patch to handles list
                h.append(v.get_patch_by_id(i))
                # append count to labels list
                l.append(count)
            except:
                pass
        ax.legend(handles=h, labels=l, title="Counts",loc="upper left", bbox_to_anchor=(1,1))
        if not os.path.exists('PhosphoProteomes/{0}'.format(self.name)):
            os.makedirs('PhosphoProteomes/{0}'.format(self.name))
        plt.savefig(fname, format='png', dpi=600, bbox_inches='tight')
        plt.close()
        
    def prep_Logo(self,ratio_cutoff,max_val=None):
        AminoAcids = ['I', 'Q', 'M', 'D', 'L', 'S', 'T', 'W', 'G','H',
                      'A', 'E', 'P', 'Y', 'V', 'K','C', 'R', 'N', 'F']
        Positions = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
        cutoff = np.log(ratio_cutoff)
        if max_val == None:
            scores = {}
            sums = {}
            for pos in Positions:
                sums[pos] = 0
                scores[pos] = []
                for aa in AminoAcids:
                    s = self.PWMs['ln(Enrichment)'][pos][aa]
                    if s > cutoff:
                        sums[pos] = sums[pos] + s
                        scores[pos].append(s)
                    else:
                        scores[pos].append(0)
            scores[0][13] = max(sums.values())
            Scores = []
            for pos in Positions:
                Scores.append(scores[pos])
        else:
            Scores = []
            for pos in Positions:
                scores = []
                for aa in AminoAcids:
                    s = self.PWMs['ln(Enrichment)'][pos][aa]
                    if pos == 0 and aa == 'Y':
                        scores.append(max_val)
                    elif s > cutoff:
                        scores.append(s)
                    else:
                        scores.append(0)
                Scores.append(scores)
        if not os.path.exists('PhosphoProteomes/{0}'.format(self.name)):
            os.makedirs('PhosphoProteomes/{0}'.format(self.name))
        df = pd.DataFrame(Scores)
        df.T.to_csv('PhosphoProteomes/{0}/{0}_Scores.csv'.format(self.name),header=False,index=False)
        
    def gen_freq_array(self):
        x_vals = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
        y_vals = ['G', 'A', 'V', 'L', 'I', 'W', 'F', 'Y', 'P',  'M', 'C', 'S', 'T', 'N', 'Q', 'D', 'E', 'H', 'K', 'R']
        array_list = []
        i = 0
        for aa in y_vals:
            array_list.append([])
            for pos in x_vals:
                if pos == 0: #and aa == 'Y':
                    array_list[i].append(0)
                else:
                    array_list[i].append(self.PWMs['Frequency'][pos][aa]) 
            i = i+1
        Array = np.array(array_list)
        self.PWM_array = Array
        pass

class DataBaseProteome(PhosphoProteome):
    
    def __init__(self,name,data,Bg):
        self.name = name
        self.data = data
        self.Bg = Bg
        self.gen_sets()
        self.gen_PWMs()
        self.name = name + ' PSP'
        
    def gen_sets(self):
        self.Seqs = set()
        self.Proteins = set()        
        for i,row in self.data.iterrows():
            if row['KINASE'] == self.name:
                if row['SUB_ORGANISM'] != 'human':
                    pass
                elif row['SITE_+/-7_AA'][0] == '_' or row['SITE_+/-7_AA'][-1] == '_':
                    pass
                else:
                    self.Seqs.add(row['SITE_+/-7_AA'].upper())
                    prot = row['SUB_ACC_ID'].split(';')[0]
                    if len(prot)>0:
                        self.Proteins.add(prot)
        if not os.path.exists('Sequences'):
            os.makedirs('Sequences')
        with open('Sequences/{0} PSP Sequences.txt'.format(self.name),'w') as outfile:
            outfile.write('\n'.join(self.Seqs))
        
class Seq_List_Proteome(PhosphoProteome):
    
    def __init__(self,name,seq_fname,prot_fname,Bg):
        self.name = name
        with open(seq_fname) as infile:
            self.Seqs = set([seq.rstrip() for seq in infile.readlines() if len(seq.rstrip()) == 15])
        with open(prot_fname) as infile:
            self.Proteins = set([prot.rstrip() for prot in infile.readlines() if len(prot.rstrip()) > 0])
        self.Bg = Bg
        self.gen_PWMs()      


class SwissProt(BackgroundProteome):
    
    def __init__(self):
        self.name = 'Swiss Prot'
        self.gen_sets()
        self.gen_PWM()
        if not os.path.exists('Sequences'):
            os.makedirs('Sequences')
        with open('Sequences/Swiss Prot Sequences.txt','w') as outfile:
            outfile.write('\n'.join(self.Seqs))
        self.prep_Logo()
        self.gen_heatmap()
            
            
    
    def gen_sets(self):
        try:
            with open('data/SwissProt/Proteome.json','r') as infile:
                self.Proteome = json.load(infile)
            with open('data/SwissProt/Proteins.json','r') as infile:
                self.Proteins = set(json.load(infile))
            with open('data/SwissProt/Sequences.json','r') as infile:
                self.Seqs = set(json.load(infile))
        except:
            self.Seqs = set()
            self.Proteins = set()
            self.Proteome = {}
            fname = 'data/uniprot_sprot.fasta'
            organism = 'OS=Homo sapiens'
            use = False
            with open(fname, 'r') as file:
                lines = file.readlines()
            for i in lines:
                if i.startswith('>'):
                    use = False
                    if organism in i: 
                        line = i.split('|')
                        try: 
                            name = line[1]
                        except: 
                            name = line[0]
                        self.Proteome[name] = {'sequence': ''}
                        use = True
                elif use:
                    self.Proteome[name]['sequence'] += i.strip()
            for entry in self.Proteome.keys():
                seq = self.Proteome[entry]['sequence']
                index=-1
                length=len(seq)
                for aa in seq:
                    index=index+1
                    if index>7 and index<(length-7):
                        if aa == 'Y':
                            start=index-7
                            end=index+8
                            peptide=seq[start:end]
                            self.Seqs.add(peptide)
                            self.Proteins.add(entry)
            if not os.path.exists('data/SwissProt'):
                os.makedirs('data/SwissProt')
            with open('data/SwissProt/Proteome.json','w') as outfile:
                json.dump(self.Proteome,outfile)
            with open('data/SwissProt/Proteins.json','w') as outfile:
                json.dump(list(self.Proteins),outfile)
            with open('data/SwissProt/Sequences.json','w') as outfile:
                json.dump(list(self.Seqs),outfile)
            if not os.path.exists('Sequences'):
                os.makedirs('Sequences')
            with open('Sequences/SwissProt Sequences.txt','w') as outfile:
                outfile.write('\n'.join(self.Seqs))
                
                
        


class CompareProteomes:
    
    def __init__(self, Proteomes,p=0.05,ratio_cutoff=1.25):
        self.Proteomes = Proteomes
        self.p = p
        self.ratio_cutoff = ratio_cutoff
        if not os.path.exists('PhosphoProteomes/Comparisons'):
            os.makedirs('PhosphoProteomes/Comparisons')


    def gen_venn(self,kinases,venn_type):
        labels = []
        sets = []
        title = venn_type + ' of\n'
        fname = 'PhosphoProteomes/Comparisons/' + venn_type + ' of '
        for name in kinases:
            if venn_type == 'Seqs':
                sets.append(set([seq.split('-')[0] for seq in self.Proteomes[name].Seqs]))
                label = '{0}\nSubstrates\n{1}'.format(self.Proteomes[name].name,len(self.Proteomes[name].Seqs))
            if venn_type == 'Proteins':
                sets.append(set([prot.split('-')[0] for prot in self.Proteomes[name].Proteins]))
                label = '{0}\nSubstrates\n{1}'.format(self.Proteomes[name].name,len(self.Proteomes[name].Proteins))
            labels.append(label)
            title = title + name + ' and '
            fname = fname + name + ' and '
        title = title[0:-4] +'Venn Diagram'
        fname = fname[0:-4] +'Venn Diagram.svg'
        
        
        if len(sets) < 2:
            pass
        elif len(sets) == 2:
            from matplotlib_venn import venn2,venn2_circles
            from matplotlib import pyplot as plt
            plt.figure()
            ax = plt.gca()
            v = venn2(sets, set_labels=labels)
            try: v.get_patch_by_id('10').set_color('red')
            except: pass
            try: v.get_patch_by_id('01').set_color('blue')
            except: pass
            try: v.get_patch_by_id('11').set_color('magenta')
            except: pass
            c = venn2_circles(sets, linestyle='solid')
            c[0].set_edgecolor('red')
            c[1].set_edgecolor('blue')
            plt.title(title)
            h, l = [],[]
            Sets = ['10','11','01']
            for i in Sets:
                # remove label by setting them to empty string:
                try:
                    count = v.get_label_by_id(i).get_text()
                    v.get_label_by_id(i).set_text("")
                    # append patch to handles list
                    h.append(v.get_patch_by_id(i))
                    # append count to labels list
                    l.append(count)
                except:
                    pass
            ax.legend(handles=h, labels=l, title="Counts",loc="upper left", bbox_to_anchor=(1,1))
            plt.savefig(fname,bbox_inches='tight', format = 'svg', dpi = 1200)
            plt.close()
            
        elif len(sets) == 3:
            from matplotlib_venn import venn3, venn3_circles
            from matplotlib import pyplot as plt
            plt.figure()
            ax = plt.gca()
            v = venn3(sets,set_labels = labels)
            try: v.get_patch_by_id('100').set_color('red')
            except: pass
            try: v.get_patch_by_id('001').set_color('green')
            except: pass
            try: v.get_patch_by_id('010').set_color('blue')
            except: pass
            try: v.get_patch_by_id('101').set_color('yellow')
            except: pass
            try: v.get_patch_by_id('011').set_color('turquoise')
            except: pass
            try: v.get_patch_by_id('110').set_color('purple')
            except: pass
            try: v.get_patch_by_id('111').set_color('grey')
            except: pass
            circles = venn3_circles(sets, linestyle='solid')
            circles[0].set_edgecolor('red')
            circles[1].set_edgecolor('blue')
            circles[2].set_edgecolor('green')
            plt.title(title)
            h, l = [],[]
            Sets = ['100','010','001','110','011','101','111']
            for i in Sets:
                # remove label by setting them to empty string:
                try:
                    count = v.get_label_by_id(i).get_text()
                    v.get_label_by_id(i).set_text("")
                    # append patch to handles list
                    h.append(v.get_patch_by_id(i))
                    # append count to labels list
                    l.append(count)
                except:
                    pass
            ax.legend(handles=h, labels=l, title="Counts",loc="upper left", bbox_to_anchor=(1,1))
            plt.savefig(fname,bbox_inches='tight', format = 'eps', dpi = 1200)
            plt.close()
        
    def gen_entropy(self,kinases,fname='Entropy Figure',IC=False):
        import matplotlib.pyplot as plt
        E_Vals = {}
        Annots_Full = []
        Masks_Full = []
        for name in kinases:
            E_Vals[name] = {}
            annots = []
            mask = []
            if name == 'Background':
                PWM = self.Proteomes[kinases[0]].Bg.PWM
            else:
                PWM = self.Proteomes[name].PWMs['Frequency']
            for pos in PWM.keys():
                probs = []
                for aa in PWM[pos].keys():
                    probs.append(PWM[pos][aa])
                e = scipy.stats.entropy(probs)
                if IC:
                    e = np.log(20) - e
                
                if pos != 0:
                    annots.append('{:.2f}'.format(e))
                    mask.append(True)
                    E_Vals[name][pos] = e
                else:
                    annots.append('pY')
                    mask.append(False)
                    E_Vals[name][pos] = 0
            Annots_Full.append(annots)
            Masks_Full.append(mask)
        
        Mask = np.array(Masks_Full)
        Annots = np.array(Annots_Full)
        df = pd.DataFrame(E_Vals)
        df = df.T
        df = df.reindex(kinases)
        min_val = 6
        max_val = 0
        for i in df.columns:
            if i == 0:
                pass
            else:
                for val in df[i]:
                    if val == 'pY':
                        pass
                    else:
                        if val < min_val:
                            min_val = val
                        if val > max_val:
                            max_val = val
        fig, ax = plt.subplots(figsize=(5,2)) 
        sns.heatmap(df,vmin=min_val, vmax=max_val,cmap='Reds_r',linewidth=1,ax=ax,
                    annot=False,fmt='', square=True,annot_kws={"size": 12},
                    cbar_kws={'shrink':.7})
        sns.heatmap(df,cmap=['#929591','#929591'],mask=Mask,cbar=False,annot=Annots,fmt='',annot_kws={"size": 12},linewidth=1, linecolor='#929591',xticklabels=1, yticklabels=1)
        plt.yticks(rotation=0,fontsize=12,fontweight='bold')
        plt.xticks(fontsize=12,fontweight='bold')
        plt.title('Sequence Entropy of Substrates\nat Each Position',fontweight="bold",fontsize=14)
        cbar = ax.collections[0].colorbar
        # here set the labelsize by 20
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label('Entropy',size=13,fontweight='bold')
        if not os.path.exists('Output'):
            os.makedirs('Output')
        plt.savefig('Output/{0}.png'.format(fname), format='png', dpi=600, bbox_inches='tight')
        plt.savefig('Output/{0}.eps'.format(fname), format='eps', dpi=600, bbox_inches='tight')
        plt.show()
        plt.close()
        return df
        
    def pos_heatmap(self,pos,AAs='All',max_val=.45,sig=True):
        from matplotlib import pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        sym = '■'
        if sig == True:
            p = 0.05
            a = p/(14*20)
            sig = -np.log10(a/(1-a))
        elif sig == False:
            sig = 0
        if self.ratio_cutoff == 0:
            cutoff = 0
        else:
            cutoff = np.log(self.ratio_cutoff)
        
        wRange = cutoff/max_val*100
        wRange = int(wRange*10)
        cmap_default = sns.color_palette("RdBu_r",1000)
        cmap_new = []
        blues = cmap_default[100:450]
        reds = cmap_default[550:900]
        b_range = int((1000 - wRange)/2)
        r_range = b_range
        interval = len(blues) / b_range
        for i in range(b_range):
            color = blues[int(interval*i)]
            cmap_new.append(color)
        for i in range(wRange):
            cmap_new.append((1,1,1))
        for i in range(r_range):
            color = reds[int(interval*i)]
            cmap_new.append(color)
        cmap_final = []
        cmap_final = cmap_new
        Custom_cmap = LinearSegmentedColormap.from_list(name='Custom_CMap', colors =cmap_final, N=len(cmap_new))
        x_vals = ['Abl', 'Anc AS', 'Anc AST', 'Anc S1', 'Src']
        y_vals = AAs
        if AAs == 'All':
            y_vals = ['G', 'A', 'V', 'L', 'I', 'W', 'F', 'Y', 'P',  'M', 'C', 'S',
                      'T', 'N', 'Q', 'D', 'E', 'H', 'K', 'R']
        array_list = []
        annot_list = []
        i = 0
        for aa in y_vals:
            array_list.append([])
            annot_list.append([])
            for dname in x_vals:
                p_val = self.Proteomes[dname].PWMs['LogsOdd'][pos][aa]
                e_val = self.Proteomes[dname].PWMs['ln(Enrichment)'][pos][aa]
                annot_list[i].append('{:.2f}'.format(e_val))
                array_list[i].append(self.Proteomes[dname].PWMs['ln(Enrichment)'][pos][aa]) 
            i = i+1
        Array = np.array(array_list)
        Annots = np.array(annot_list)
        fig, ax = plt.subplots(figsize=(2.1,2))
        x_vals = ['Abl', 'Anc AS', 'Anc AST', 'Anc S1', 'Src']
        ax = sns.heatmap(Array, xticklabels = x_vals, yticklabels = y_vals , vmax=max_val, vmin=-max_val, linecolor='black',
                         linewidths = 1, cmap=Custom_cmap, square = True, annot=Annots, fmt='',
                         cbar=False,ax=ax,annot_kws={'size':9})
                        # annot_kws={'ha':'center','va':'center',"size": 10,'fontname':'Arial','color':'#929591'}
        plt.ylabel('Amino Acids',size=14,fontweight='bold')
        plt.yticks(rotation=0,ha='center',size=12,fontweight='bold')
        plt.xticks(rotation=90,ha='center',size=12,fontweight='bold') 
        ax.set_title('Preference at {:+0}'.format(pos),size=14,fontweight='bold')
        d = 'Heatmap_Pos'
        if not os.path.exists('Output/{0}'.format(d)):
            os.makedirs('Output/{0}'.format(d))
        plt.savefig('Output/Heatmap at {:+0}.eps'.format(pos), format='eps', dpi=1200, bbox_inches='tight')
        plt.close()
        
        
    def graph_data(self,aa,pos,PWM_T='ln(Enrichment)',plot_sig=False):
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mtick
        plt.style.use('seaborn-whitegrid')
        plt.style.use('seaborn-paper')
        c1 = (44/255,169/255,225/255)#'#92cad1'
        c2 = (195/255,44/255,59/255)#'#e9724d'
        marker_size = 45
        line_size = 3
        f_size = 12
        t_size = 14
        fig = plt.figure(figsize=(2.5,2))
        fig.tight_layout()
        
        ax = plt.subplot(111)
        ax.set_xmargin(.075)
        ax.set_ymargin(.15)
        if len(aa) == 1:
            y1_vals = np.array([self.Proteomes['Anc AST'].PWMs[PWM_T][pos][aa],
                                self.Proteomes['Anc AS'].PWMs[PWM_T][pos][aa],
                                self.Proteomes['Anc S1'].PWMs[PWM_T][pos][aa],
                                self.Proteomes['Src'].PWMs[PWM_T][pos][aa]]).astype(np.double) 
            y2_vals = np.array([self.Proteomes['Anc AST'].PWMs[PWM_T][pos][aa],
                                self.Proteomes['Anc AS'].PWMs[PWM_T][pos][aa],
                                None,
                                self.Proteomes['Abl'].PWMs[PWM_T][pos][aa]]).astype(np.double) 
        else:
            y1_vals, y2_vals = self.calc_graph_vals_motif(aa,pos,PWM_T)
        
    
        
        y1mask = np.isfinite(y1_vals)
        y2mask = np.isfinite(y2_vals)
        x_vals = np.array([1,2,3,4])
        x_ticks = ['Anc AST', 'Anc AS', 'Anc S1 / None', 'Src / Abl']
        plt.plot(x_vals[y1mask], y1_vals[y1mask], color=c1, linestyle='-',linewidth=line_size,marker='o',ms=3)
        plt.scatter(x_vals[y1mask], y1_vals[y1mask], color=c1, label='Src Evolution', s=marker_size)
        plt.scatter(x_vals[y2mask], y2_vals[y2mask], color=c2, label='Abl Evolution', s=marker_size)
        plt.plot(x_vals[y2mask], y2_vals[y2mask], color=c2, linestyle='-',linewidth= line_size)
        
        
        
        if PWM_T == 'Percent Enrichment':
            plt.ylabel('Percent Enriched', fontsize=f_size)
        if PWM_T == 'ln(Enrichment)':
            plt.ylabel('Normalized\nLog Probability', fontsize=f_size,fontweight="bold")
        else:
            plt.ylabel('Enrichment Ratio', fontsize=f_size)
        plt.xlabel('Time', fontsize=f_size,fontweight="bold")
        plt.title('Preference for {} at {:+0}'.format(aa,pos), fontsize=t_size,fontweight="bold")
        
        plt.plot(x_vals[y1mask][0:2], y1_vals[y1mask][0:2], color=c2, linestyle='-',linewidth=line_size)
        plt.plot(x_vals[y1mask][0:2], y1_vals[y1mask][0:2], color=c1, linestyle='--',dashes=(3, 3),linewidth=line_size,marker='o',ms=7,mew=2,mec=c2)

        ### Annotate Points with Kinase Labels
        o = 9
        if y1_vals[0] > y1_vals[1]:
            oAST = o
            oAS = -o
        else:
            oAST = -o
            if y1_vals[1]>y1_vals[2]:
                oAS = o
            elif y1_vals[1]>y2_vals[3]:
                oAS = o
            else:
                oAS = -o
        plt.annotate('AncAST', (1, y1_vals[0]),xytext=(13, oAST),textcoords='offset points',horizontalalignment='center',verticalalignment='center',fontweight="bold",fontsize=f_size)
        plt.annotate('AncAS', (2, y1_vals[1]),xytext=(0, oAS),textcoords='offset points',horizontalalignment='center',verticalalignment='center',fontweight="bold",fontsize=f_size)
        plt.annotate('AncS1', (3, y1_vals[2]),xytext=(0, -o),textcoords='offset points',horizontalalignment='center',verticalalignment='center',fontweight="bold",fontsize=f_size)
        plt.annotate('Src', (4, y1_vals[3]),xytext=(0, -o),textcoords='offset points',horizontalalignment='center',verticalalignment='center',fontweight="bold",fontsize=f_size)
        plt.annotate('Abl', (4, y2_vals[3]),xytext=(0, -o),textcoords='offset points',horizontalalignment='center',verticalalignment='center',fontweight="bold",fontsize=f_size)
        
        
        plt.xticks([], [])
        plt.tick_params(labelsize=f_size)
        plt.grid(False)
        if PWM_T == 'Percent Enrichment':
            fmt = '{x:,.0%}'
            yticks = mtick.StrMethodFormatter(fmt)
            ax=plt.gca()
            ax.yaxis.set_major_formatter(yticks)
        d = 'Graphs'
        if not os.path.exists('Output/{0}'.format(d)):
            os.makedirs('Output/{0}'.format(d))
        plt.savefig('Output/Graphs/{0} at {1}.eps'.format(aa,pos), format='eps', dpi=1200, bbox_inches='tight')
        plt.close()
    
    def calc_graph_vals_motif(self,aa,pos,PWM_T):
        motif = ['.', '.', '.', '.', '.', '.', '.', 'Y', '.', '.', '.', '.', '.', '.', '.']
        motif[pos + 7] = '[' + aa + ']'
        motif = re.compile(''.join(motif))
        bg_p = len(motif.findall('\n'.join(self.Proteomes['Abl'].Bg.Seqs))) / len(self.Proteomes['Abl'].Bg.Seqs)
        bg_lp = np.log(bg_p)
        Scores = {}
        dnames = ['Anc AST','Anc AS','Anc S1','Src','Abl']
        for dname in dnames:
            p = len(motif.findall('\n'.join(self.Proteomes[dname].Seqs))) / len(self.Proteomes[dname].Seqs)
            if PWM_T == 'ln(Enrichment)':
                Scores[dname] = np.log(p) - bg_lp
            if PWM_T == 'Enrichment':
                Scores[dname] = p / bg_p
            if PWM_T == 'Percent Enrichment':
                Scores[dname] = (p - bg_p) / bg_p 
        y1_vals = np.array([Scores['Anc AST'],Scores['Anc AS'],Scores['Anc S1'],Scores['Src']]).astype(np.double) 
        y2_vals = np.array([Scores['Anc AST'],Scores['Anc AS'], None, Scores['Abl']]).astype(np.double) 
        return y1_vals, y2_vals
        
    
    def venn_heatmaps(self,kinases,max_val):
        AminoAcids = ['I', 'Q', 'M', 'D', 'L', 'S', 'T', 'W', 'G','H',
                      'A', 'E', 'P', 'Y', 'V', 'K','C', 'R', 'N', 'F']
        Positions = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
        
        for i in [(0,1),(1,0),'Both']:
            if i == 'Both':
                Seqs = self.Proteomes[kinases[0]].Seqs | self.Proteomes[kinases[1]].Seqs
                title = 'Residue Preference For Substrates\nIn Both {0} and {1}'.format(kinases[0],kinases[1])
                fname = 'PhosphoProteomes/Comparisons/Heatmap {0} and {1}.eps'.format(kinases[0],kinases[1])
            else:
                Seqs = self.Proteomes[kinases[0]].Seqs - self.Proteomes[kinases[1]].Seqs
                title = 'Residue Preference For {0} Substrates\nNot In {1}'.format(kinases[i[0]],kinases[i[1]])
                fname = 'PhosphoProteomes/Comparisons/Heatmap {0} not {1}.eps'.format(kinases[i[0]],kinases[i[1]])
            PWM = {}
            sigPWM = {}
            for pos in Positions:
                PWM[pos] = {}
                sigPWM[pos] = {}
                for aa in AminoAcids:
                    p = self.Proteomes[kinases[0]].Bg.PWM[pos][aa]
                    motif = motif = ['.','.','.','.','.','.','.','Y','.','.','.','.','.','.','.']
                    motif[7+pos] = aa
                    Motif = re.compile(''.join(motif))
                    count = len(Motif.findall('\n'.join(Seqs)))
                    N = len(Seqs)
                    freq = count / N
                    try:
                        PWM[pos][aa] = np.log(freq/p)
                    except:
                        PWM[pos][aa] = 0
                    tail1 = scipy.stats.binom_test(count, N, p, alternative='greater')
                    tail2 = scipy.stats.binom_test(count, N, p, alternative='less')
                    sigPWM[pos][aa] = - np.log10(tail1/tail2)
            self.gen_HeatMap(PWM,sigPWM,title,fname,p=p,max_val=max_val)
        
        
        
    
    def gen_HeatMap(self,PWM,sigPWM,title,fname,p=0.05,max_val=0.45):   
        import seaborn as sns
        from matplotlib import pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        sym = '■'
        a = self.p/(14*20)
        sig = -np.log10(a/(1-a))
        rc = np.log(self.ratio_cutoff)
        cutoff = rc
        Mark = 'Sig'
        Mask = False
        wRange = cutoff/max_val*100
        wRange = int(wRange*10)
        cmap_default = sns.color_palette("RdBu_r",1000)
        cmap_new = []
        blues = cmap_default[100:450]
        reds = cmap_default[550:900]
        b_range = int((1000 - wRange)/2)
        r_range = b_range
        interval = len(blues) / b_range
        for i in range(b_range):
            color = blues[int(interval*i)]
            cmap_new.append(color)
        for i in range(wRange):
            cmap_new.append((1,1,1))
        for i in range(r_range):
            color = reds[int(interval*i)]
            cmap_new.append(color)
        cmap_final = []
        cmap_final = cmap_new
        Custom_cmap = LinearSegmentedColormap.from_list(name='Custom_CMap', colors =cmap_final, N=len(cmap_new))
        sns.set(font_scale=1.5)
        x_vals = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
        y_vals = ['G', 'A', 'V', 'L', 'I', 'W', 'F', 'Y', 'P',  'M', 'C', 'S', 'T', 'N', 'Q', 'D', 'E', 'H', 'K', 'R']
        array_list = []
        annot_list = []
        c_list = []
        i = 0
        for aa in y_vals:
            array_list.append([])
            annot_list.append([])
            c_list.append([])
            for pos in x_vals:
                if pos == 0: #and aa == 'Y':
                    array_list[i].append(0)
                    annot_list[i].append('')
                    c_list[i].append(False)
                else:
                    c_list[i].append(True)
                    if (sigPWM[pos][aa] <sig) and (sigPWM[pos][aa]>-sig):
                        if PWM[pos][aa] > cutoff or PWM[pos][aa] < -cutoff:
                            if Mark == 'Insig':
                                annot_list[i].append(sym)
                            else:
                                annot_list[i].append('')
                        else:
                            annot_list[i].append('')
                        if Mask:
                            array_list[i].append(0)
                        else:
                            array_list[i].append(PWM[pos][aa])
                    else:
                        array_list[i].append(PWM[pos][aa])
                        if Mark == 'Sig':
                            if abs(PWM[pos][aa]) >= rc:
                                annot_list[i].append(sym)
                            else:
                                annot_list[i].append('')   
                        else:
                            annot_list[i].append('')    
            i = i+1
        Array = np.array(array_list)
        Annots = np.array(annot_list)
        Center = np.array(c_list)
        
        fig, ax = plt.subplots(figsize=(6,6))
        ax = sns.heatmap(Array, xticklabels = x_vals, yticklabels = y_vals , vmax=max_val, vmin=-max_val, linecolor='black',
                         linewidths = .5, cmap=Custom_cmap, square = True, annot=Annots,
                         annot_kws={'ha':'center','va':'center',"size": 10,'fontname':'Arial','color':'black','fontdict':{'color':'black'}},#'color':'#929591'
                         fmt='', cbar=False,ax=ax)
        sns.heatmap(Array,cmap=['#929591','#929591'],mask=Center,cbar=False,annot=False,linewidth=1, linecolor='black', xticklabels = x_vals, yticklabels = y_vals)
        plt.xlabel('Position from P-Tyrosine',size=16)
        plt.ylabel('Amino Acids',size=16)
        plt.yticks(rotation=0,ha='center',size=13)
        plt.xticks(rotation=0,ha='center',size=13)
        plt.savefig(fname, format='eps', dpi=1200, bbox_inches='tight')
        plt.close()
        
        
        
