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
    parser.add_argument("--type", type=str, help="What type of pairs are being compared? (homology-based, NC-terminal, split-proteins)")
    parser.add_argument("--pairs", default=None, type=str, help="")
    parser.add_argument("--emd_dirr1", default=None, type=str, help="")
    parser.add_argument("--emd_dirr2", default=None, type=str, help="")
    parser.add_argument("--comment", type=str, help="")

    args = parser.parse_args()
    diff_lst=[]

    if args.type == "Homology":
      print ("Homologous protein pairs")
      
      pairs1 = open(args.pairs, 'r')
      lines_pairs1 = pairs1.readlines()
      
      for i in lines_pairs1:
         i = i.split()
         P1 = i[0]
         P2 = i[1]
         #print (P1, P2)

         epath1 = os.path.join(args.emd_dirr1 + P1 + '.npz')
         emd1 = np.load(epath1)

         epath2 = os.path.join(args.emd_dirr2 + P2 + '.npz')
         emd2 = np.load(epath2)

         f1 = emd1.files
         for item in f1:
            emd1 = (emd1[item])
         f2 = emd2.files
         for item in f2:
            emd2 = (emd2[item])
         #print (emd1.shape)
         #print (emd2.shape)

         p_val = stats.spearmanr(emd1, emd2)
         #print (p_val)
     
         '''std1 = stats.tstd(emd1, axis=0)
         std2 = stats.tstd(emd2, axis=0)
         diff_std = abs((std1)-(std2))
         print (diff_std)'''

         diff_lst.append(p_val)
 
         '''for E1, E2 in zip(os.listdir(args.emd_dirr1), os.listdir(args.emd_dirr2)):
         if E1.endswith('.npz') and E2.endswith('.npz'):
            emb1 = np.load(args.emd_dirr1+E1)
            emb2 = np.load(args.emd_dirr2+E2)
            diff = emd1 - emd2
            trel = stats.ttest_rel(emd1, emd2, axis=1)
            print (t_rel)
            diff_lst.append(diff_std)'''
 
         sns.scatterplot(emd1,emd2)

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
      plt.close()

    if args.type == "Non-Homology":
      print ("Non-homologous protein pairs")
      
      pairs1 = open(args.pairs, 'r')
      lines_pairs1 = pairs1.readlines()
      
      for i in lines_pairs1:
         i = i.split()
         P1 = i[0]
         P2 = i[1]

         epath1 = os.path.join(args.emd_dirr1 + P1 + '.npz')
         emd1 = np.load(epath1)

         epath2 = os.path.join(args.emd_dirr2 + P2 + '.npz')
         emd2 = np.load(epath2)

         f1 = emd1.files
         for item in f1:
            emd1 = (emd1[item])
         f2 = emd2.files
         for item in f2:
            emd2 = (emd2[item])
         
         p_val = stats.spearmanr(emd1, emd2)
     
         diff_lst.append(p_val)

         sns.scatterplot(emd1,emd2)
      plt.xlabel('nhp 1 emd values')
      plt.ylabel('nhp 2 emd values')
      name = args.type + '_pcor_nhp_' + args.comment
      plt.savefig('{}.png'.format(name), bbox_inches='tight')
      plt.close()
 
      sns.distplot(diff_lst)
      plt.xlabel('Correlation values between non-homologous pairs')
      plt.ylabel('Count')
      name = args.type +'_pval_nhp_'+ args.comment
      plt.savefig('{}.png'.format(name), bbox_inches='tight')
      plt.close()

    if args.type == "Swap-motifs":
      
      pairs1 = open(args.pairs, 'r')
      lines_pairs1 = pairs1.readlines()
      
      for i in lines_pairs1:
         i = i.split()
         P1 = i[0]

         print (P1)

         epath1 = os.path.join(args.emd_dirr1 + P1 + '.npz')
         emd1 = np.load(epath1)
         try:
            epath2 = os.path.join(args.emd_dirr2 + P1 + '_10' + '.npz')
            emd2 = np.load(epath2)
         except FileNotFoundError:
            continue

         f1 = emd1.files
         for item in f1:
            emd1 = (emd1[item])
         f2 = emd2.files
         for item in f2:
            emd2 = (emd2[item])
         
         p_val = stats.spearmanr(emd1, emd2)
     
         diff_lst.append(p_val)

         sns.scatterplot(emd1,emd2)

      plt.xlabel('original emd values')
      plt.ylabel('modified emd values')
      name = args.type + '_pcor_mod_' + args.comment
      plt.savefig('{}.png'.format(name), bbox_inches='tight')
      plt.close()
 
      sns.distplot(diff_lst)
      plt.xlabel('Correlation values between original and modified pairs')
      plt.ylabel('Count')
      name = args.type +'_pval_mod_'+ args.comment
      plt.savefig('{}.png'.format(name), bbox_inches='tight')
      plt.close()


      '''

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
        plt.close()'''


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

