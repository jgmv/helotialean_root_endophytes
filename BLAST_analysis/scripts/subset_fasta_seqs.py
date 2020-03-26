#!/usr/bin/python

'''
    File name: subset_fasta_seqs.py
    Author: Jose G. Macia-Vicente
    Date created: 2019-08-04
    Python Version: 2.7
'''


from Bio import SeqIO
import re
import sys
import argparse


parser = argparse.ArgumentParser(
    description='Subset fasta sequences from a file with sequence names.')

parser.add_argument(
    "i",
    help = 'input fasta file')

parser.add_argument(
    "f",
    help = 'file with sequence names')

parser.add_argument(
    "-o", help = 'output file',
    type=str,
    default='output.fasta')

args = parser.parse_args()

if1 = args.i
if2 = args.f
of1 = args.o

count = 0

handle = open(if1)
names = open(if2).read().splitlines()
seqOut = open(of1, "w")

processed = []
for seq_record in SeqIO.parse(handle, "fasta"):
	if seq_record.id in names:
		seq = seq_record.id
		seqOut.write(">"+seq+"\n")
		seqOut.write(str(seq_record.seq.lower())+"\n")
		processed.append(seq)
		count +=1

for i in names:
    if i not in processed:
		print i+" not in fasta file."

print str(count)+" sequences extracted"

handle.close()
seqOut.close()
