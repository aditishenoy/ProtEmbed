import argparse
import numpy as np
import pandas as pd
import sys
import os
import json

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_dataset(xml_file, length):
    ids = []
    seq_len = []
    scl = []
    sct = []
    seq = []

    df = pd.DataFrame()
    print (type(length))
    
    for index, record in enumerate(SeqIO.parse(xml_file, "uniprot-xml")):
         ids.append(record.id)
         seq_len.append(record.annotations['sequence_length'])
         seq.append(record.seq)
         if (len(ids)) >= length:
             break

    a = pd.Series(ids)
    b = pd.Series(seq_len)
    c = pd.Series(seq)
    df = pd.DataFrame({'ID' : a, 'Seq_len' : b, 'Sequence' : c})
    print (df)
    df.to_csv('uniprot_2000.csv') 
    return df


def get_ids(xml_file):

   for index, record in enumerate(SeqIO.parse(xml_file, "uniprot-xml")):
       with open('swissprot_ids', 'a') as f:
           print(record.id, file=f)

def main():

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('--xml_file', type=str, help="")
    parser.add_argument('--l', type=int, help="")
    args = parser.parse_args()

    print (args.xml_file)
    print (args.l)
    #df = get_dataset_from_xml(args.xml_file, args.l)
    df = get_ids(args.xml_file)

if __name__ == "__main__":
    main()

