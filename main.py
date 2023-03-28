#!/usr/bin/env python3
# coding=utf-8

import pandas as pd
import pysam
import os
import subprocess

# Check whether the specified path exists or not
subprocess.run([
        "rm", "-rf", './motifRes'
    ])
subprocess.run([
        "rm", "-rf", './motifResMeme'
    ])
if not os.path.exists('motifRes'):

   # Create a new directory because it does not exist
   os.makedirs('motifRes')

# Check whether the specified path exists or not
if not os.path.exists('motifResMeme'):

   # Create a new directory because it does not exist
   os.makedirs('motifResMeme')


def getControl(path):
    df = pd.read_csv(path,sep='\t')
    df2 = df.loc[df['Diff events'] == 'No']
    return df2

def getChanged(path):
    df = pd.read_csv(path,sep='\t')
    df2 = df.loc[df['Diff events'] == 'Yes']
    return df2

def getExonPosition(row, exon):
    ret = {'name':'','chr':'','start':'','end':'','strand':''}
    ret['chr'] = row[exon].split(':')[0][3:]
    ret['start'] = int(row[exon].split(':')[1].split('-')[0])
    ret['end'] = int(row[exon].split(':')[1].split('-')[1])
    ret['name'] = row[exon]    
    return ret

def getRNA(seq, strand):
    if strand == '-':
        seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
        seq = seq.upper()
        seq = seq[::-1]
    else:
        seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
        seq = seq.upper()
    seq=seq.replace("T", "U")
    #print(seq)
    return seq

pathList = ['results/PSI_FIVE.txt','results/PSI_MXE.txt','results/PSI_RI.txt','results/PSI_SE.txt','results/PSI_THREE.txt',]

resList = {
    'se_exon_up200_upRegulated':{'ctr':'','sig':''},
    'se_exon_down200_upRegulated':{'ctr':'','sig':''},
    'five_exonA_down200_upRegulated':{'ctr':'','sig':''},
    'five_exonB_down200_upRegulated':{'ctr':'','sig':''},
    'three_exonA_up200_upRegulated':{'ctr':'','sig':''},
    'three_exonB_up200_upRegulated':{'ctr':'','sig':''},
    'mxe_exonA_up200_upRegulated':{'ctr':'','sig':''},
    'mxe_exonA_down200_upRegulated':{'ctr':'','sig':''},
    'mxe_exonB_up200_upRegulated':{'ctr':'','sig':''},
    'mxe_exonB_down200_upRegulated':{'ctr':'','sig':''},
    'ri_exonA_down200_upRegulated':{'ctr':'','sig':''},
    'ri_exonB_up200_upRegulated':{'ctr':'','sig':''},
    'se_exon_up200_dnRegulated':{'ctr':'','sig':''},
    'se_exon_down200_dnRegulated':{'ctr':'','sig':''},
    'five_exonA_down200_dnRegulated':{'ctr':'','sig':''},
    'five_exonB_down200_dnRegulated':{'ctr':'','sig':''},
    'three_exonA_up200_dnRegulated':{'ctr':'','sig':''},
    'three_exonB_up200_dnRegulated':{'ctr':'','sig':''},
    'mxe_exonA_up200_dnRegulated':{'ctr':'','sig':''},
    'mxe_exonA_down200_dnRegulated':{'ctr':'','sig':''},
    'mxe_exonB_up200_dnRegulated':{'ctr':'','sig':''},
    'mxe_exonB_down200_dnRegulated':{'ctr':'','sig':''},
    'ri_exonA_down200_dnRegulated':{'ctr':'','sig':''},
    'ri_exonB_up200_dnRegulated':{'ctr':'','sig':''},
}
genome = pysam.Fastafile('Mus_musculus.GRCm38.dna.primary_assembly.fa')

#PSI_SE.txt upstream 200
tmpControl = getControl('./results/PSI_SE.txt')
tmpSig = getChanged('./results/PSI_SE.txt')
for index, row in tmpControl.iterrows():
    ret = getExonPosition(row,'exon')
    try:
        seq = genome.fetch(ret['chr'], ret['start']-200, ret['start'])
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    resList['se_exon_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    resList['se_exon_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'

for index, row in tmpSig.iterrows():
    ret = getExonPosition(row,'exon')
    try:
        seq = genome.fetch(ret['chr'], ret['start']-200, ret['start'])
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    if row['dPSI']>0:
        resList['se_exon_up200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
    else:
        resList['se_exon_up200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
    
    


#PSI_SE.txt downstream 200
tmpControl = getControl('./results/PSI_SE.txt')
tmpSig = getChanged('./results/PSI_SE.txt')
for index, row in tmpControl.iterrows():
    ret = getExonPosition(row,'exon')
    try:
        seq = genome.fetch(ret['chr'], ret['end'], ret['end']+200)
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    resList['se_exon_down200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    resList['se_exon_down200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'

for index, row in tmpSig.iterrows():
    ret = getExonPosition(row,'exon')
    try:
        seq = genome.fetch(ret['chr'], ret['end'], ret['end']+200)
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    if row['dPSI']>0:
        resList['se_exon_down200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
    else:
        resList['se_exon_down200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'

#five exona down200
tmpControl = getControl('./results/PSI_FIVE.txt')
tmpSig = getChanged('./results/PSI_FIVE.txt')
for index, row in tmpControl.iterrows():
    ret = getExonPosition(row,'exon_a')
    try:
        seq = genome.fetch(ret['chr'], ret['end'], ret['end']+200)
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    resList['five_exonA_down200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    resList['five_exonA_down200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'

for index, row in tmpSig.iterrows():
    ret = getExonPosition(row,'exon_a')
    try:
        seq = genome.fetch(ret['chr'], ret['end'], ret['end']+200)
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    if row['dPSI']>0:
        resList['five_exonA_down200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
    else:
        resList['five_exonA_down200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'

#five exonb down200
tmpControl = getControl('./results/PSI_FIVE.txt')
tmpSig = getChanged('./results/PSI_FIVE.txt')
for index, row in tmpControl.iterrows():
    ret = getExonPosition(row,'exon_b')
    try:
        seq = genome.fetch(ret['chr'], ret['end'], ret['end']+200)
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    resList['five_exonB_down200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    resList['five_exonB_down200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    
for index, row in tmpSig.iterrows():
    ret = getExonPosition(row,'exon_b')
    try:
        seq = genome.fetch(ret['chr'], ret['end'], ret['end']+200)
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    if row['dPSI']>0:
        resList['five_exonB_down200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
    else:
        resList['five_exonB_down200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'

#three exon a up 200
tmpControl = getControl('./results/PSI_THREE.txt')
tmpSig = getChanged('./results/PSI_THREE.txt')
for index, row in tmpControl.iterrows():
    ret = getExonPosition(row,'exon_a')
    try:
        seq = genome.fetch(ret['chr'], ret['start']-200, ret['start'])
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    resList['three_exonA_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    resList['three_exonA_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'

for index, row in tmpSig.iterrows():
    ret = getExonPosition(row,'exon_a')
    try:
        seq = genome.fetch(ret['chr'], ret['start']-200, ret['start'])
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    if row['dPSI']>0:
        resList['three_exonA_up200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
    else:
        resList['three_exonA_up200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
        

#three exon b up 200
tmpControl = getControl('./results/PSI_THREE.txt')
tmpSig = getChanged('./results/PSI_THREE.txt')
for index, row in tmpControl.iterrows():
    ret = getExonPosition(row,'exon_b')
    try:
        seq = genome.fetch(ret['chr'], ret['start']-200, ret['start'])
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    resList['three_exonB_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    resList['three_exonB_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'

for index, row in tmpSig.iterrows():
    ret = getExonPosition(row,'exon_b')
    try:
        seq = genome.fetch(ret['chr'], ret['start']-200, ret['start'])
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    if row['dPSI']>0:
        resList['three_exonB_up200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
    else:
        resList['three_exonB_up200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'

#mxe exon a up 200
tmpControl = getControl('./results/PSI_MXE.txt')
tmpSig = getChanged('./results/PSI_MXE.txt')
for index, row in tmpControl.iterrows():
    ret = getExonPosition(row,'exon_a')
    try:
        seq = genome.fetch(ret['chr'], ret['start']-200, ret['start'])
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    resList['mxe_exonA_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    resList['mxe_exonA_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'

for index, row in tmpSig.iterrows():
    ret = getExonPosition(row,'exon_a')
    try:
        seq = genome.fetch(ret['chr'], ret['start']-200, ret['start'])
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    if row['dPSI']>0:
        resList['mxe_exonA_up200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
    else:
        resList['mxe_exonA_up200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'

#mxe exon a down 200
tmpControl = getControl('./results/PSI_MXE.txt')
tmpSig = getChanged('./results/PSI_MXE.txt')
for index, row in tmpControl.iterrows():
    ret = getExonPosition(row,'exon_a')
    try:
        seq = genome.fetch(ret['chr'], ret['end'], ret['end']+200)
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    resList['mxe_exonA_down200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    resList['mxe_exonA_down200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'

for index, row in tmpSig.iterrows():
    ret = getExonPosition(row,'exon_a')
    try:
        seq = genome.fetch(ret['chr'], ret['end'], ret['end']+200)
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    if row['dPSI']>0:
        resList['mxe_exonA_down200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
    else:
        resList['mxe_exonA_down200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'

#mxe exon b up 200
tmpControl = getControl('./results/PSI_MXE.txt')
tmpSig = getChanged('./results/PSI_MXE.txt')
for index, row in tmpControl.iterrows():
    ret = getExonPosition(row,'exon_b')
    try:
        seq = genome.fetch(ret['chr'], ret['start']-200, ret['start'])
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    resList['mxe_exonB_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    resList['mxe_exonB_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'

for index, row in tmpSig.iterrows():
    ret = getExonPosition(row,'exon_a')
    try:
        seq = genome.fetch(ret['chr'],ret['start']-200, ret['start'])
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    if row['dPSI']>0:
        resList['mxe_exonB_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    else:
        resList['mxe_exonB_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'

#mxe exon b down 200
tmpControl = getControl('./results/PSI_MXE.txt')
tmpSig = getChanged('./results/PSI_MXE.txt')
for index, row in tmpControl.iterrows():
    ret = getExonPosition(row,'exon_b')
    try:
        seq = genome.fetch(ret['chr'], ret['end'], ret['end']+200)
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    resList['mxe_exonB_down200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    resList['mxe_exonB_down200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'

for index, row in tmpSig.iterrows():
    ret = getExonPosition(row,'exon_b')
    try:
        seq = genome.fetch(ret['chr'], ret['end'], ret['end']+200)
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    if row['dPSI']>0:
        resList['mxe_exonB_down200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
    else:
        resList['mxe_exonB_down200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'

#ri exon a down 200
tmpControl = getControl('./results/PSI_RI.txt')
tmpSig = getChanged('./results/PSI_RI.txt')
for index, row in tmpControl.iterrows():
    ret = getExonPosition(row,'exon_a')
    try:
        seq = genome.fetch(ret['chr'], ret['end'], ret['end']+200)
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    resList['ri_exonA_down200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    resList['ri_exonA_down200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'

for index, row in tmpSig.iterrows():
    ret = getExonPosition(row,'exon_a')
    try:
        seq = genome.fetch(ret['chr'], ret['end'], ret['end']+200)
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    if row['dPSI']>0:
        resList['ri_exonA_down200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
    else:
        resList['ri_exonA_down200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'

#ri exon b up 200
tmpControl = getControl('./results/PSI_RI.txt')
tmpSig = getChanged('./results/PSI_RI.txt')
for index, row in tmpControl.iterrows():
    ret = getExonPosition(row,'exon_b')
    try:
        seq = genome.fetch(ret['chr'], ret['start']-200, ret['start'])
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    resList['ri_exonB_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
    resList['ri_exonB_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'

for index, row in tmpSig.iterrows():
    ret = getExonPosition(row,'exon_b')
    try:
        seq = genome.fetch(ret['chr'], ret['end']-200, ret['end'])
    except:
        print('error fetching '+ret['chr'])
    res = getRNA(seq, ret['strand'])
    if row['dPSI']>0:
        resList['ri_exonB_up200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
    else:
        resList['ri_exonB_up200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'

for item in resList:
    with open('motifRes/'+item+'_ctr.txt', 'w') as f:
        f.write(resList[item]['ctr'])
    with open('motifRes/'+item+'_sig.txt', 'w') as f:
        f.write(resList[item]['sig'])
        
with open('motifRes/total_ctr.txt', 'w') as f:
    for item in resList:
        f.write(resList[item]['ctr'])
with open('motifRes/total_sig.txt', 'w') as f:
    for item in resList:
        f.write(resList[item]['sig'])

    
for item in resList:
    print('now running '+item)
    print('xstreme --p motifRes/'+item+'_sig.txt --n motifRes/'+item+'_ctr.txt --oc motifResMeme/'+item+' --rna --m /opt/dbs/meme/motif_databases/motifs.meme --meme-p 32 --minw 4 --minw 8')
    os.system('xstreme --p motifRes/'+item+'_sig.txt --n motifRes/'+item+'_ctr.txt --oc motifResMeme/'+item+' --rna --m /opt/dbs/meme/motif_databases/motifs.meme --meme-p 32 --minw 4 --minw 8  > /dev/null  &')


os.system('xstreme --p motifRes/total_sig.txt --n motifRes/total_ctr.txt --oc motifResMeme/total --rna --m /opt/dbs/meme/motif_databases/motifs.meme --meme-p 32 --minw 4 --minw 8 > /dev/null &')