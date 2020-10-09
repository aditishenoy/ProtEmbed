import argparse
import numpy as np
import pandas as pd
import sys
import os
import json

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_ids(xml_file, output):
   for index, record in enumerate(SeqIO.parse(xml_file, "uniprot-xml")):
       with open('{}'.format(output), 'a') as f:
           print(record.id, file=f)

def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--xml_file', type=str, help="")
    parser.add_argument('--output', type=str, help="")
    args = parser.parse_args()

    df = get_ids(args.xml_file, args.output)

if __name__ == "__main__":
    main()

