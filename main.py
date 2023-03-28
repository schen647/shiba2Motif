import argparse
import pandas as pd
import pysam
import os
import subprocess

class mainEventLoop:
    def __init__(self, args):
        self.args = args
        self.tmp = args.tmp
        self.out = args.out
        self.input_five = args.input_five
        self.input_mxe = args.input_mxe
        self.input_ri = args.input_ri
        self.input_se = args.input_se
        self.input_three = args.input_three
        self.ref = args.ref
        self.xstremePath = args.xstremePath
    
    def getControl(self,path):
        df = pd.read_csv(path,sep='\t')
        df2 = df.loc[df['Diff events'] == 'No']
        return df2

    def getChanged(self,path):
        df = pd.read_csv(path,sep='\t')
        df2 = df.loc[df['Diff events'] == 'Yes']
        return df2

    def getExonPosition(self,row, exon):
        ret = {'name':'','chr':'','start':'','end':'','strand':''}
        ret['chr'] = row[exon].split(':')[0][3:]
        ret['start'] = int(row[exon].split(':')[1].split('-')[0])
        ret['end'] = int(row[exon].split(':')[1].split('-')[1])
        ret['name'] = row[exon]    
        ret['strand'] = row['strand']
        return ret

    def getRNA(self,seq, strand):
        #print(seq,strand)
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

    def getSeq(self,chrN, start, end,shift, strand):
        genome = pysam.Fastafile(self.ref)
            #print('checking')
            #print(chrN, start, end, shift, strand)
        if strand == '+' and shift > 0:
            startPos = end
            endPos = end+shift
        if strand == '-' and shift > 0:
            startPos = start-shift
            endPos = start
        if strand == '+' and shift < 0:
            startPos = start+shift
            endPos = start
        if strand == '-' and shift < 0:
            startPos = end
            endPos = end-shift
        return genome.fetch(chrN, startPos, endPos)
    
    def launch(self):

        subprocess.run([
                "rm", "-rf", self.tmp +'/motifRes'
            ])
        subprocess.run([
                "rm", "-rf", self.out+'./motifResMeme'
            ])
        if not os.path.exists(self.tmp +'/motifRes'):

        # Create a new directory because it does not exist
            os.makedirs(self.tmp +'/motifRes')

        # Check whether the specified path exists or not
        if not os.path.exists(self.out+'./motifResMeme'):

        # Create a new directory because it does not exist
            os.makedirs(self.out+'./motifResMeme')
        # parse the input list into 5 files
        #pathList = ['results/PSI_FIVE.txt','results/PSI_MXE.txt','results/PSI_RI.txt','results/PSI_SE.txt','results/PSI_THREE.txt',]
        pathList= [
            self.input_five,
            self.input_mxe,
            self.input_ri,
            self.input_se,
            self.input_three,
        ]

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

        #PSI_SE.txt upstream 200
        tmpControl = self.getControl(self.input_se)
        tmpSig = self.getChanged(self.input_se)
        for index, row in tmpControl.iterrows():
            ret = self.getExonPosition(row,'exon')
            try:
            #seq = genome.fetch(ret['chr'], ret['start']-200, ret['start'])
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],-200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
            #print(res)
                resList['se_exon_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                resList['se_exon_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            #print('added')
            except:
                print('error fetching '+ret['chr'])


        for index, row in tmpSig.iterrows():
            ret = self.getExonPosition(row,'exon')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],-200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                if row['dPSI']>0:
                    resList['se_exon_up200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
                else:
                    resList['se_exon_up200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])

            
            

        #PSI_SE.txt downstream 200
        tmpControl = self.getControl(self.input_se)
        tmpSig = self.getChanged(self.input_se)
        for index, row in tmpControl.iterrows():
            ret = self.getExonPosition(row,'exon')
            try:
                #print('checking')
                #print(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                resList['se_exon_down200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                resList['se_exon_down200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        for index, row in tmpSig.iterrows():
            ret = self.getExonPosition(row,'exon')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                if row['dPSI']>0:
                    resList['se_exon_down200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
                else:
                    resList['se_exon_down200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        #five exona down200
        tmpControl = self.getControl(self.input_five)
        tmpSig = self.getChanged(self.input_five)
        for index, row in tmpControl.iterrows():
            ret = self.getExonPosition(row,'exon_a')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                resList['five_exonA_down200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                resList['five_exonA_down200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        for index, row in tmpSig.iterrows():
            ret = self.getExonPosition(row,'exon_a')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                if row['dPSI']>0:
                    resList['five_exonA_down200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
                else:
                    resList['five_exonA_down200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        #five exonb down200
        tmpControl = self.getControl(self.input_five)
        tmpSig = self.getChanged(self.input_five)
        for index, row in tmpControl.iterrows():
            ret = self.getExonPosition(row,'exon_b')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                resList['five_exonB_down200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                resList['five_exonB_down200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])

            
        for index, row in tmpSig.iterrows():
            ret = self.getExonPosition(row,'exon_b')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                if row['dPSI']>0:
                    resList['five_exonB_down200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
                else:
                    resList['five_exonB_down200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        #three exon a up 200
        tmpControl = self.getControl(self.input_three)
        tmpSig = self.getChanged(self.input_three)
        for index, row in tmpControl.iterrows():
            ret = self.getExonPosition(row,'exon_a')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],-200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                resList['three_exonA_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                resList['three_exonA_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        for index, row in tmpSig.iterrows():
            ret = self.getExonPosition(row,'exon_a')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],-200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                if row['dPSI']>0:
                    resList['three_exonA_up200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
                else:
                    resList['three_exonA_up200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])

                

        #three exon b up 200
        tmpControl = self.getControl(self.input_three)
        tmpSig = self.getChanged(self.input_three)
        for index, row in tmpControl.iterrows():
            ret = self.getExonPosition(row,'exon_b')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],-200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                resList['three_exonB_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                resList['three_exonB_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        for index, row in tmpSig.iterrows():
            ret = self.getExonPosition(row,'exon_b')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],-200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                if row['dPSI']>0:
                    resList['three_exonB_up200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
                else:
                    resList['three_exonB_up200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        #mxe exon a up 200
        tmpControl = self.getControl(self.input_mxe)
        tmpSig = self.getChanged(self.input_mxe)
        for index, row in tmpControl.iterrows():
            ret = self.getExonPosition(row,'exon_a')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],-200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                resList['mxe_exonA_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                resList['mxe_exonA_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        for index, row in tmpSig.iterrows():
            ret = self.getExonPosition(row,'exon_a')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],-200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                if row['dPSI']>0:
                    resList['mxe_exonA_up200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
                else:
                    resList['mxe_exonA_up200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        #mxe exon a down 200
        tmpControl = self.getControl(self.input_mxe)
        tmpSig = self.getChanged(self.input_mxe)
        for index, row in tmpControl.iterrows():
            ret = self.getExonPosition(row,'exon_a')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                resList['mxe_exonA_down200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                resList['mxe_exonA_down200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        for index, row in tmpSig.iterrows():
            ret = self.getExonPosition(row,'exon_a')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                if row['dPSI']>0:
                    resList['mxe_exonA_down200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
                else:
                    resList['mxe_exonA_down200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        #mxe exon b up 200
        tmpControl = self.getControl(self.input_mxe)
        tmpSig = self.getChanged(self.input_mxe)
        for index, row in tmpControl.iterrows():
            ret = self.getExonPosition(row,'exon_b')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],-200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                resList['mxe_exonB_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                resList['mxe_exonB_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        for index, row in tmpSig.iterrows():
            ret = self.getExonPosition(row,'exon_a')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],-200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                if row['dPSI']>0:
                    resList['mxe_exonB_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                else:
                    resList['mxe_exonB_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        #mxe exon b down 200
        tmpControl = self.getControl(self.input_mxe)
        tmpSig = self.getChanged(self.input_mxe)
        for index, row in tmpControl.iterrows():
            ret = self.getExonPosition(row,'exon_b')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                resList['mxe_exonB_down200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                resList['mxe_exonB_down200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        for index, row in tmpSig.iterrows():
            ret = self.getExonPosition(row,'exon_b')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                if row['dPSI']>0:
                    resList['mxe_exonB_down200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
                else:
                    resList['mxe_exonB_down200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        #ri exon a down 200
        tmpControl = self.getControl(self.input_ri)
        tmpSig = self.getChanged(self.input_ri)
        for index, row in tmpControl.iterrows():
            ret = self.getExonPosition(row,'exon_a')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                resList['ri_exonA_down200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                resList['ri_exonA_down200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        for index, row in tmpSig.iterrows():
            ret = self.getExonPosition(row,'exon_a')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                if row['dPSI']>0:
                    resList['ri_exonA_down200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
                else:
                    resList['ri_exonA_down200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        #ri exon b up 200
        tmpControl = self.getControl(self.input_ri)
        tmpSig = self.getChanged(self.input_ri)
        for index, row in tmpControl.iterrows():
            ret = self.getExonPosition(row,'exon_b')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],-200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                resList['ri_exonB_up200_upRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
                resList['ri_exonB_up200_dnRegulated']['ctr']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        for index, row in tmpSig.iterrows():
            ret = self.getExonPosition(row,'exon_b')
            try:
                seq = self.getSeq(ret['chr'], ret['start'], ret['end'],-200,ret['strand'])
                res = self.getRNA(seq, ret['strand'])
                if row['dPSI']>0:
                    resList['ri_exonB_up200_upRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
                else:
                    resList['ri_exonB_up200_dnRegulated']['sig']+='>'+ret['name']+'\n'+res+'\n'
            except:
                print('error fetching '+ret['chr'])


        for item in resList:
            with open(self.tmp+'/motifRes/'+item+'_ctr.txt', 'w') as f:
                f.write(resList[item]['ctr'])
            with open(self.tmp+'/motifRes/'+item+'_sig.txt', 'w') as f:
                f.write(resList[item]['sig'])
                
        with open(self.tmp+'/motifRes/total_ctr.txt', 'w') as f:
            for item in resList:
                f.write(resList[item]['ctr'])
        with open(self.tmp+'/motifRes/total_sig.txt', 'w') as f:
            for item in resList:
                f.write(resList[item]['sig'])

            
        for item in resList:
            print('now running '+item)
            print(self.xstremePath+' --p '+self.tmp+'/motifRes/'+item+'_sig.txt --n '+self.tmp+'/motifRes/'+item+'_ctr.txt --oc '+self.out+'/motifResMeme/'+item+' --rna --m '+self.ref+' --meme-p 32 --minw 4 --minw 8')
            os.system(self.xstremePath+' --p '+self.tmp+'/motifRes/'+item+'_sig.txt --n '+self.tmp+'/motifRes/'+item+'_ctr.txt --oc '+self.out+'/motifResMeme/'+item+' --rna --m '+self.ref+' --meme-p 32 --minw 4 --minw 8 > /dev/null  &')


        os.system(self.xstremePath+' --p  '+self.tmp+'/motifRes/total_sig.txt --n '+self.tmp+'/motifRes/total_ctr.txt --oc '+self.out+'/motifResMeme/'+item+' --rna --m '+self.ref+' --meme-p 32 --minw 4 --minw 8 > /dev/null &')


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--tmp",
        type=str,
        required=False,
        default="/dev/shm",
        help="Path to the tmp dir for intermediate files",
    )


    parser.add_argument(
        "--out",
        type=str,
        required=True,
        #default="/dev/shm"
        help="Where the output will be saved",
    )

    parser.add_argument(
        "--input_five",
        type=str,
        required=True,
        #default="/dev/shm"
        help="results file with pattern \"FIVE\" produced by shiba",
    )

    parser.add_argument(
        "--input_mxe",
        type=str,
        required=True,
        #default="/dev/shm"
        help="results file with pattern \"MXE\" produced by shiba",
    )

    parser.add_argument(
        "--input_ri",
        type=str,
        required=True,
        #default="/dev/shm"
        help="results file with pattern \"RI\" produced by shiba",
    )

    parser.add_argument(
        "--input_se",
        type=str,
        required=True,
        #default="/dev/shm"
        help="results file with pattern \"SE\" produced by shiba",
    )

    parser.add_argument(
        "--input_three",
        type=str,
        required=True,
        #default="/dev/shm"
        help="results file with pattern \"THREE\" produced by shiba",
    )

    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        #default="/dev/shm"
        help="path to the reference genome",
    )

    parser.add_argument(
        "--xstremePath",
        type=str,
        required=False,
        default="xstreme",
        help="path to the meme binary",
    )

    
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    #a = os.system(args.xstremePath)
    #print('===')
    #print(a)
    #print('===')
    assert os.system(args.xstremePath+' > /dev/null 2>&1')==512, "Unable to reach meme, please specify a path to the meme binary"

    # check if reference genome can be loaded
    assert pysam.FastaFile(args.ref) is not None, "Unable to load reference genome, please provide a file that can be parsed as a fasta reference genome"

    # check if all the input files exist
    assert os.path.exists(args.input_five) and os.path.exists(args.input_mxe) and os.path.exists(args.input_ri) and os.path.exists(args.input_se) and os.path.exists(args.input_three), "Unable to find input file"

    # create output folder if it doesn't exist
    if not os.path.exists(args.out):
        os.makedirs(args.out)
    
    # check if the tmp folder exists, if not, use /tmp, if still not, create /tmp
    if not os.path.exists(args.tmp):
        args.tmp = "/tmp"
        if not os.path.exists("/tmp"):
            os.makedirs("/tmp")          

    program = mainEventLoop(args)
    program.launch()  
            
