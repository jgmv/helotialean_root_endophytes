#!/usr/bin/python

'''
    File name: filter_fasta_seqs.py
    Author: Jose G. Macia-Vicente
    Date created: 2018-07-25
    Date last modified: 2018-07-25
    Python Version: 2.7
'''


from Bio import SeqIO
import re
import sys
import argparse


parser = argparse.ArgumentParser(
    description='Filter fasta sequences by size and rename.')

parser.add_argument(
    "i",
    help = 'input fasta file')

parser.add_argument(
    "-o", help = 'output file',
    type=str,
    default='output.fasta')

parser.add_argument(
    "-l", help = 'minimum sequence length',
    type=int,
    default=150)

parser.add_argument(
    "-e", help = 'regex expression (default: "\_(.*?)\_")',
    type=str,
    default='\_(.*?)\_')

parser.add_argument(
    "-p", help = 'regex pattern position (default: 0)',
    type=int,
    default=0)


args = parser.parse_args()


if1 = args.i
of1 = args.o
desc1 = args.e
desc2 = args.p
minLength = args.l

count = 0

handle = open(if1)
seqOut = open(of1, "w")

for seq_record in SeqIO.parse(handle, "fasta"):
	if len(seq_record.seq) > minLength:
		strain = re.findall(desc1, seq_record.id)[desc2]
		seqOut.write(">"+strain+"\n")
		seqOut.write(str(seq_record.seq.lower())+"\n")		

		count +=1

print str(count)+" sequences extracted"

handle.close()
seqOut.close()
