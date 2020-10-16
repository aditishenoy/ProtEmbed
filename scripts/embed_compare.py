import numpy as np
import pandas as pd
import argparse

import sys
import os
from os import path
from tqdm import tqdm
import h5py
import torch

from scipy import stats

import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

matplotlib.use('agg')

def main():

    parser = argparse.ArgumentParser(description="Metrics to compare embeddings")
    parser.add_argument("--pairs", default=None, type=str, help="")
    parser.add_argument("--emd_dirr1", default=None, type=str, help="")
    parser.add_argument("--emd_dirr2", default=None, type=str, help="")
    parser.add_argument("--comment", type=str, help="")

    args = parser.parse_args()
    path = "/home/a/aditi/pfs/current_projects/ProtEmbed/data/"

    pairs1 = open(path+args.pairs, 'r')
    lines_pairs1 = pairs1.readlines()
 
    mod = swap(lines_pairs1)
    plst=[]
    clst= []  
    ''' 
    
    for i in lines_pairs1:
       
       i = i.split()
       P1 = i[0]
       P2 = i[1]

       epath1 = os.path.join(path+args.emd_dirr1 + P1 + '.npz')
       emd1 = np.load(epath1)

       epath2 = os.path.join(path+args.emd_dirr2 + P2 + '.npz')
       emd2 = np.load(epath2)

       f1 = emd1.files
       for item in f1:
         emd1 = (emd1[item])
       f2 = emd2.files
       for item in f2:
         emd2 = (emd2[item])'''



    '''
    df = pd.DataFrame()
    df['p'] = plst
    df['c']= clst
    sns.jointplot(df)
    plt.xlabel('plst')
    plt.ylabel('clst')
    #name = args.type + '_pcor_hp_' + args.comment
    plt.savefig('1.png', bbox_inches='tight')
    plt.close()
    sns.displot(a, kde=True)
    plt.xlabel('Correlation values between homology pairs (Spearman)')
    plt.ylabel('Count')
    name = 'Homology_cval_hp_'+ args.comment
    plt.savefig('{}.png'.format(name), bbox_inches='tight')
    plt.close()
    sns.distplot(clst)
    plt.xlabel('Correlation values between homology pairs (Correlate Numpy)')
    plt.ylabel('Count')
    #name = args.type +'_pval_hp_'+ args.comment
    plt.savefig('clst_np.png', bbox_inches='tight')
    plt.close()'''


def stack_embed(emd1, emd2):
    print (emd1.shape, emd2.shape)   

    plt.plot(emd1)
    plt.plot(emd2)
    plt.savefig('stacked.png', bbox_inches='tight')


def acc_len(emd1, emd2, P1, P2, clst):
 
   df = pd.read_csv('/home/a/aditi/pfs/current_projects/ProtEmbed/data/EXP-2-MultipleProtein/Protset1_len', sep = ",", names=['ID', 'Seq_len'])
   length  = (df.loc[df['ID']==P1]['Seq_len'].values)[0]

   if int(length) < 200:
       clst.append(emd1)

   return clst
   
   
def moving_average(emd1, emd2):
    emd1 = pd.Series(emd1)
    a = emd1.rolling(window=10).mean()
    print (a.shape)
    emd2 = pd.Series(emd2)
    b = emd2.rolling(window=10).mean()
    print (b.shape)

    plt.plot(a)
    plt.plot(b)
    plt.savefig('ma4.png', bbox_inches='tight')


def homology(emd1, emd2, plst, clst):

    p_val = stats.spearmanr(emd1, emd2)
    print (p_val)

    corr = np.correlate(emd1, emd2)
    print (corr)
       
    plst.append(p_val[0])
    clst.append(corr[0])
    print (clst)

    return plst, clst

    """

    plt.xlabel('hp 1 emd values')
    plt.ylabel('hp 2 emd values')
    name = args.type + '_pcor_hp_' + args.comment
    plt.savefig('{}.png'.format(name), bbox_inches='tight')
    plt.close()

    sns.distplot(diff_lst)
    plt.xlabel('Correlation values between homology pairs')
    plt.ylabel('Count')
    name = args.type +'_pval_hp_'+ args.comment
    plt.savefig('{}.png'.format(name), bbox_inches='tight')
    plt.close()"""

def swap(lines):
      
    for i in lines:
        i = i.split()
        P1 = i[0]
        print (P1)

        epath1 = os.path.join('/home/a/aditi/pfs/current_projects/ProtEmbed/data/EXP-2-MultipleProtein/ProtSet1/' + P1 + '.npz')
        emd1 = np.load(epath1)

        epath2 = os.path.join('/home/a/aditi/pfs/current_projects/ProtEmbed/data/EXP-2-MultipleProtein/mod_5/' + P1 + '_5' + '.npz')
        emd2 = np.load(epath2)

        f1 = emd1.files
        for item in f1:
           emd1 = (emd1[item])
        f2 = emd2.files
        for item in f2:
           emd2 = (emd2[item])
       
        p_val = stats.spearmanr(emd1, emd2)
        print (p_val)

        plst.append(p_val[0])

    sns.displot(plst, kde=True)
    plt.xlabel('Correlation values between original and modified pairs')
    plt.ylabel('Count')
    name = args.type +'_pval_mod_'+ args.comment
    plt.savefig('{}.png'.format(name), bbox_inches='tight')
    plt.close()

    """
    if args.type == "Swap-NC":

        print ("Swapped NC residues")

        p_val = stats.spearmanr(emd1, emd2)
        print (p_val)
        sns.distplot(p_val)
        plt.xlabel('Correlation values between original and swapped NC terminal residue')
        plt.ylabel('Count')
        name = args.type +'_pval_swapnc_'+ args.comment
        plt.savefig('{}.png'.format(name), bbox_inches='tight')
        plt.close()

        sns.scatterplot(emd1,emd2)
        plt.xlabel('Original emd values')
        plt.ylabel('Swapped NC emd values')
        name = args.type + '_pcor_swapnc_' + args.comment
        plt.savefig('{}.png'.format(name), bbox_inches='tight')
        plt.close()

    if args.type == "Split":

        print ("Shuffled residues")

        p_val = stats.spearmanr(emd1, emd2)
        print (p_val)
        
        sns.scatterplot(emd1,emd2)
        plt.xlabel('Original emd values')
        plt.ylabel('split emd values')
        name = args.type + '_pcor_split_' + args.comment
        plt.savefig('{}.png'.format(name), bbox_inches='tight')
        plt.close()

    if args.type == "RamdomShuffle":

        print ("Ramdomly shuffled residues")

        p_val = stats.spearmanr(emd1, emd2)
        print (p_val)
        
        sns.scatterplot(emd1,emd2)
        plt.xlabel('Original emd values')
        plt.ylabel('Randomly shuffled emd values')
        name = args.type + '_pcor_Rshuf_' + args.comment
        plt.savefig('{}.png'.format(name), bbox_inches='tight')
        plt.close()
 
        std1 = stats.tstd(emd1, axis=0)
         std2 = stats.tstd(emd2, axis=0)
         diff_std = abs((std1)-(std2))
         print (diff_std)

 
        for E1, E2 in zip(os.listdir(args.emd_dirr1), os.listdir(args.emd_dirr2)):
        if E1.endswith('.npz') and E2.endswith('.npz'):
            emb1 = np.load(args.emd_dirr1+E1)
            emb2 = np.load(args.emd_dirr2+E2)
            diff = emd1 - emd2
            trel = stats.ttest_rel(emd1, emd2, axis=1)
            print (t_rel)
            diff_lst.append(diff_std)
            """

def read_fasta_sequence(afile, query_id=''):
    seq_dict = {}
    header = ''
    seq = ''

    for aline in afile:
        aline = aline.strip()

        if aline.startswith('>'):
            if header != '' and seq != '':
                if header in seq_dict:
                    seq_dict[header].append(seq)
                else:
                    seq_dict[header] = [seq]
            seq = ''
            if aline.startswith('>%s' % query_id) and query_id !='':
                header = query_id
            else:
                header = aline[1:]
        else:
            #aline_seq = aline.translate(None, '.-').upper()
            seq += aline
    return seq


if __name__ == '__main__':
    main()

