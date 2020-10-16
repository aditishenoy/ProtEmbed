import argparse
import numpy as np
import pandas as pd
import sys
import os
import json
import random 

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

    return header, seq
     
def get_seqlen(fastafolder):
   df = pd.DataFrame()
   seqlist = []   
   flist = []
   for filename in os.listdir(fastafolder):
        if filename.endswith('.fa'):
            f = open(fastafolder+filename, 'r')
            flist.append(filename[:-3])
            seq = (read_fasta_sequence(f))
            seqlen = len(seq)
            seqlist.append(seqlen)

   a = pd.Series(flist)
   b = pd.Series(seqlist)
   df = pd.DataFrame({'ID' : a, 'Seq_len' : b})
   print (df)

   return df

def get_multiples(seqlen, pointer):
     ind_list = []
     count = 0
     ind_list.append(0)
     for i in range(seqlen):
         if count < pointer:
             count=count+1
         else:
             ind_list.append(i)
             count = 0
     return ind_list

def remove_at(seq, start_pos, end_pos):
    return seq[:start_pos] + seq[end_pos:]

def actual_swap(seq, ran1, ran2, frame):
       
        tmp1 = seq[ran1:ran1+frame+1] 
        tmp2 = seq[ran2:ran2+frame+1] 
        new_seq_1 = remove_at(seq, ran1, ran1+len(tmp1))
        new_seq_2 = remove_at(new_seq_1, ran2-len(tmp1), ran2)
        new_seq3 = new_seq_2[:ran1] + tmp2 + new_seq_2[ran1:]
        new_seq4 = new_seq3[:ran2] + tmp1 + new_seq3[ran2:]
        return new_seq4


def shuffle_seqs(fastafolder, frame, outputfolder):

    for filename in os.listdir(fastafolder):

        if filename.endswith('.fa'):
            f = open(fastafolder+filename, 'r')
 
            header, seq = (read_fasta_sequence(f))

            seqlen = len(seq)
            if seqlen < 50:
                continue
            else:
                ind_list = get_multiples(seqlen, frame)

                for i in range(10):
                    ran1 = random.choice(ind_list)
                    ran2 = random.choice(ind_list)
                    if ran2 > ran1:
                        seq = actual_swap(seq, ran1, ran2, frame)

                with open(outputfolder+'{}'.format(filename[:-3]+'_'+str(frame+1)+'.fa'), 'w') as f:
                       print ('>'+header, file=f)
                       print(seq, file=f)

       
def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta_folder', type=str, help="")
    parser.add_argument('--frame', type=int, help="")
    parser.add_argument('--output', type=str, help="")
    args = parser.parse_args()
    print (args.fasta_folder)

    #df = get_seqlen(args.fasta_folder)
    #df.to_csv('Protset2_len', index = False)
     
    shuf5 = shuffle_seqs(args.fasta_folder, args.frame, args.output)
    

if __name__ == "__main__":
    main()

