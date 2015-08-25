#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python ReadsCount.py --input Coldstress.table
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
Experiment	NB	NB	NB	HEG4	HEG4	HEG4	EG4	EG4	EG4
C_day1_hr1	1	2	3	4	5	6	7	8	9
C_day1_hr3	10	11	12	13	14	15	16	17	18
C_day2_hr1	19	20	21	22	23	24	25	26	27
C_day2_hr3	28	29	30	31	32	33	34	35	36
C_day3_hr1	37	38	NA	40	41	42	43	44	45
C_day3_hr3	46	47	NA	49	50	51	52	53	54
T_day1_hr1	55	56	57	58	59	60	61	62	63
T_day1_hr3	64	65	66	67	68	69	70	71	72
T_day2_hr1	73	74	75	76	77	78	79	80	81
T_day2_hr3	82	83	84	85	86	87	88	89	90
T_day3_hr1	91	92	93	94	95	96	97	98	99
T_day3_hr3	100	101	102	103	104	105	106	107	108
'''
def readtable(infile):
    bamdir = '/rhome/cjinfeng/Rice/RNAseq/Danforth/rice_cold/bam'
    data = defaultdict(lambda : defaultdict(lambda : defaultdict(str)))
    header = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith('Experiment'): 
                header = re.split(r'\t',line)
            elif(len(line) > 2):
                unit = re.split(r'\t',line)
                counts = []
                samples = []
                for i in range(1,len(unit)):
                    data[unit[0]][header[i]] = unit[i]
                    if unit[i] != 'NA':
                        prefix = '%s_%s' %(header[i], unit[i])
                        bam = '/rhome/cjinfeng/Rice/RNAseq/Danforth/rice_cold/bam/%s_%s.bam' %(header[i], unit[i])
                        count = '%s_%s.count' %(header[i], unit[i]) 
                        counts.append(count)
                        samples.append(header[i])
                        hts_count(bam, prefix)
                merge_count(counts, samples, unit[0])
                DESeq(unit[0], samples)
                #print ','.join(counts)
                #print ','.join(samples)
    return data

def hts_count(bamfile, prefix):
    cmd = '/usr/local/bin/samtools view %s | /opt/Python/2.7.3/bin/htseq-count --stranded=no - /rhome/cjinfeng/BigData/02.Transcription/Transcriptome/bin/DiffExpression_Stress/MSU7.gene.gtf > %s.count 2> %s.log' %(bamfile, prefix, prefix) 
    countfile = '%s.count' %(prefix)
    if not os.path.isfile(countfile):
        os.system(cmd)

def merge_count(counts, samples, test):
    nums = range(1,len(counts)*2+1)
    cols = [x for x in nums if x % 2 == 0]
    cmd1 = 'echo %s | sed \'s/ /\t/g\' > %s_NB_HEG4_EG4.header' %('\t'.join(samples), test)
    cmd2 = 'paste %s | grep "^LOC" | cut -f1,%s | cat %s_NB_HEG4_EG4.header - > %s_NB_HEG4_EG4.count' %(' '.join(counts), ','.join(map(str,cols)), test, test)
    #print cmd1
    #print cmd2
    os.system(cmd1)
    os.system(cmd2)

def DESeq(test, samples):
    line = ['"' + str(x) + '"' for x in samples]
    sampleline = ','.join(line)
    Rcmd = '''
library(DESeq) # version 1.9.11
library(VennDiagram)
pdf("''' + test + '''.DESeq.pdf") 
# Read in data ------------------------------------------------------------
 
## Read in the data making the row names the first column
counttable <- read.table("''' + test + '''_NB_HEG4_EG4.count", header=T, row.names=1)
head(counttable)
 
## Make metadata data.frame
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c(''' + sampleline + '''))

meta$condition <- relevel(meta$condition, ref="NB")
 
# DESeq -------------------------------------------------------------------
 
## Make a new countDataSet
d <- newCountDataSet(counttable, meta$condition)
 
## Estimate library size and dispersion
d <- estimateSizeFactors(d)
d <- estimateDispersions(d)
plotDispEsts(d, main="DESeq: Per-gene dispersion estimates")
 
## Principal components biplot on variance stabilized data, color-coded by condition-librarytype
print(plotPCA(varianceStabilizingTransformation(d), intgroup=c("condition")))
 
## Calling differential expression by N normal test
dtable1= nbinomTest(d, "NB", "HEG") 
#dtable <- dtable[order(dtable$padj), ]
write.table(dtable1,file="''' + test  + '''.DESeq.NB_HEG4.test",sep="\\t") 

dtable2= nbinomTest(d, "NB", "EG4") 
write.table(dtable2,file="''' + test  + '''.DESeq.NB_EG4.test",sep="\\t")

dtable2= nbinomTest(d, "HEG", "EG4") 
write.table(dtable2,file="''' + test  + '''.DESeq.HEG4_EG4.test",sep="\\t") 

dev.off() 
'''
    outfile = '%s.DESeq.R' %(test)
    ofile = open(outfile, 'w')
    print >> ofile, Rcmd
    ofile.close() 
    cmd = 'cat %s | R --slave' %(outfile)
    os.system(cmd)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    readtable(args.input)

if __name__ == '__main__':
    main()

