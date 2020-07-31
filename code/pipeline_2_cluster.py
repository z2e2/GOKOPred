#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 17:32:06 2020

@author: taha
"""

import pandas as pd
import numpy as np
import os
import pickle
import sys

InputFile = sys.argv[1]
OUTPUT = sys.argv[2]
NCPU = int(sys.argv[3])


# InputFile = '/home/taha/Desktop/GO/Test/TEST.pkl'
# OUTPUT = '/home/taha/Desktop/GO/Train_Test_Split'
# NCPU = 4

def Acc2Ind(ACC):
    IND = []
    for acc in ACC:
        IND.append(DATA.index[DATA['accession']==acc].tolist()[0])
    return  sorted(IND)


print('CONVERTING INPUT TO FASTA...')
DATA = pickle.load(open(InputFile,'rb'))


f = open(OUTPUT+'/KOGO.faa','w+')

for counter,ii in enumerate(DATA.index):
    
    f.write('>'+DATA.loc[ii,'accession']+'\n')
    f.write(DATA.loc[ii,'sequence']+'\n')
    
f.close()

FASTA = OUTPUT+'/KOGO.faa'
print('DONE!\n')



# Read sequence names
print('READING SEQUENCE NAMES...')
f = open(FASTA)
HEADERS = []
while 1:
    line = f.readline()
    
    if line.startswith('>'):
        HEADERS.append(line[1:-1])
        
    if not line:
        break

print(str(len(HEADERS))+' sequences were found')
print('DONE!\n')

# Build database
print('BUILDING BLAST DATABASE...')
commandDB = 'makeblastdb -in '+FASTA+' -dbtype prot -parse_seqids'
os.system(commandDB)
print('DONE!\n')

# BLAST everything against everything
print('BLASTING ALL SEQUENCES AGAINST ALL...')
commandBLAST = 'blastp -db '+FASTA+' -query '+FASTA+' -out '+OUTPUT+'/BlastRes.tsv -outfmt 7 -num_threads '+str(NCPU)
os.system(commandBLAST)
print('DONE!\n')

# Reading results
print('READING THE BLAST RESULTS...')
COL = ['query acc.ver', 'subject acc.ver', 'identity', 'alignment length',\
       'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end',\
           'evalue', 'bit scor']

RES = pd.read_csv(OUTPUT+'/BlastRes.tsv', sep='\t',header=None,names=COL)

NoNeighbor = list(HEADERS)
WithNeighbor = []

for ii in range(RES.shape[0]):
    
    if RES.loc[ii,'query acc.ver'].startswith('#'):
        continue
    
    # NEW HIT
    
    # ignore self-match:
    if RES.loc[ii,'query acc.ver'] == RES.loc[ii,'subject acc.ver']:
        continue
    
    # ignore less than 50% identitiy
    if RES.loc[ii,'identity'] < 50:
        continue
    
    if RES.loc[ii,'query acc.ver'] in NoNeighbor:
        NoNeighbor.remove(RES.loc[ii,'query acc.ver'])
        
    if not RES.loc[ii,'query acc.ver'] in WithNeighbor:
        WithNeighbor.append(RES.loc[ii,'query acc.ver'])
    
    #print(RES.loc[ii,'query acc.ver'],RES.loc[ii,'subject acc.ver'],RES.loc[ii,'identity'])


# Write results to disk
with open(OUTPUT+'/WithNeighbor.txt', 'w') as filehandle:
    for listitem in WithNeighbor:
        filehandle.write('%s\n' % listitem)

with open(OUTPUT+'/NoNeighbor.txt', 'w') as filehandle:
    for listitem in NoNeighbor:
        filehandle.write('%s\n' % listitem)

DATA.loc[Acc2Ind(WithNeighbor),:].to_pickle(OUTPUT+'/TRAIN.pkl')

DATA.loc[Acc2Ind(NoNeighbor),:].to_pickle(OUTPUT+'/TEST.pkl')

print('DONE!\n')



# report:
print(str(len(HEADERS))+' sequences were found')
print(str(len(WithNeighbor))+' sequences have neighbors')
print(str(len(NoNeighbor))+' sequences no not have neighbors')
print('Results were save in ' + OUTPUT)
