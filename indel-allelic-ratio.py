#! /usr/bin/env python3
import sys, math, os, time
from optparse import OptionParser
from collections import defaultdict

#Want the 200kb leading up to leseion, the first 200kb, the whole lesion, the last 200kb in lesion, and the following 200kb
#gen the in gene for all five of these, if borcder has bases, use it, if at either end, use periods as placeholders
#use percentage for all 5.

"""
File takes as input:
1. A file of intervals to consider and samples to consider, supplied to this program with the -l option
2. A file of allele-specific read depth at predetermined sites across the genome for some
   population. This file is the output of CallAllelesAB.py and is supplied to this program
   with the -a option.
3. The name of the desired output file, supplied to this program with the -o option.

For each sample x interval specified in the interval file, calculate the coverage of parent
"A" allele at all SNP loci contained within that interval. Also compare parent "A" coverage
at the same interval for all other samples.

Example usage:
python <path>/<to>/indel-allelic-ratio.py -a alleles_2019_0104_LOP.txt -l 2019_0129_LOP_uniq_lesions.tsv -o alleles_2019_0129_LOP_uniq_lesions.tsv
"""

usage = ""
parser = OptionParser(usage=usage)
parser.add_option("-a", "--allelecalls", dest="a", help="allele calls file")
parser.add_option("-l", "--lesions", dest="l", help="lesion file")
parser.add_option("-o", "--out", dest="o", help="output file")
#parser.add_option("-k", "--key", dest="filekey", default = "IFG", help="\"IFG\" or \"GWR\" ")


(opt, args) = parser.parse_args()


#split to length
def split(l, n):
   return(list(splitter(l, n)))

def splitter(l, n):
    i = 0
    chunk = l[:n]
    while chunk:
        yield chunk
        i += n
        chunk = l[i:i+n]



genes = defaultdict(list)

#li = os.listdir(os.getcwd())

#get ifg and gwr files, but keep separate
#fileset = filter(lambda x: opt.filekey in x and x.endswith('.txt') and 'AlleleCalls' in x, li)
#fileset.sort()
#print fileset
#oldname -> newname
#names = defaultdict(str)

#load in correct names
#f = open(opt.n)
#headname = f.readline()
#for l in f:
#   x = l.split('\t')
#   origname = x[0]
#   newname = x[5].replace('\n','')
#   if origname in names.keys() and names[origname] != newname:
#      parser.error("ERROR, same base name, different new name")
#   if newname == '':
#      x[424343]
#   names[origname] = newname
#   
#f.close()

#data array, sample -> data
piledat = defaultdict(lambda: defaultdict(list))


#first read in all data to get average 
chromlist = []
slist = set() # sample list?

f = open(opt.a)
header = f.readline()
head = header[:-1].split('\t')
temp = head[5::6] # starting at index column 5, keep every 6th column
# subnames = map(lambda x: x.replace("Snptype-", ''), temp) # this gets a list of sample names as an iterable map object rather than a list
subnames = list(map(lambda x: x.replace("Snptype-", ''), temp)) # returns the above as a list
#t = split(head[5:],6)
ct = 0
for l in f:
   ct+=1
   if ct % 10000 == 8:
#       print ct, p1 # some python2 stoofs
        print(ct, p1) # WELCOME TO THE FUTURE
#      if ct > 50000:
#         break
   x = l[:-1].split('\t') # strips the trailing '\n' then returns the rest of the data as list
   slist.add(len(x)) # should be going through allelecalls, gets length of each line
   p1 = x[1:5] # list of pos ref A B
   rest = split(x[5:],6) # list of lists, per sample, list of Snptype, SNP1, SNP2, TotalCov, CovA.
   chrom = x[0]
   if chrom not in chromlist:
      chromlist.append(chrom)
#    for snum in range(len(rest)):
   for index, name in enumerate(subnames):
      vals = rest[index]
      # oname = subnames[snum] # this isn't working, map object is not subscriptable. May be a python3 issue?
      # vals = rest[snum]
      if vals == ['.', '.', '.', '.', '.', '.']:
         continue
      #newname = names[oname]
      # newname = oname.replace("upop-",'') # probably don't need this
      if vals[4] == '.':
         vals[4] = 0
      vals = [int(vals[3]), float(vals[4])]
      # piledat[newname][chrom].append(p1+vals) # need to modify since upop not applicable
#       piledat[oname][chrom].append(p1+vals)
      piledat[name][chrom].append(p1+vals)

f.close()
# done parsing callAllelesAB.py output file, can now parse piledat for intervals of interest

databysample = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
# load in lesions by chromosome
f = open(opt.l)
leshead = f.readline()

# this loop will not parse in the lesion file that Meric shared with me.
# For my analyses, I used a modified lesion file. Name: Path:
for l in f:
   x = l[:-1].split('\t')
   # sample = x[0]
   sample = x[4]
   #####
   # for ck in piledat.keys(): # keys are sample names
   for ck in piledat.keys():
#       if x[0] in ck: # has to change if I keep my lesion file format
      if x[4] in ck: # can I change this to sample?
         sample = ck[:]
         chrom = x[0]
         start = x[1]
         end = x[2]
#       if x[0].replace('_','-')+'-' in ck: # can set up lesion file to get rid of this. However it's done, needs to match sample names from parsed CallAlleles table
#          sample = ck[:]
         #print x[:4], ck
#          chrom = x[1]
#          start = x[2]
#          end = x[3]
         for z in x:
            if z == '':
               x[x.index(z)] = '.'
#          databysample[sample][chrom][start+'-'+end] = [[sample]+x[1:],[]]
         databysample[sample][chrom][start+'-'+end] = [x[:],[]] # chrom was missing from output
         #print sample, chrom, start, end, [sample]+x[1:]

f.close()

# left off here, going to lunch.
# what does this part do?
notinlesion = defaultdict(list)
for sample in piledat.keys(): # i.e., most of the samples do not have a lesion. Hence, the loads of empty dictionaries.
   for chrom in piledat[sample].keys():
      temp = piledat[sample][chrom]
      chromlesions = databysample[sample][chrom].keys()
      for sub in temp:
         lesionhit = list(filter(lambda z: int(sub[0]) >= int(z.split('-')[0]) and int(sub[0]) <= int(z.split('-')[1]), chromlesions)) # this is the source of the crash, filter object is not subscriptable. This is a Python versioning issue.
         if lesionhit == []:
            notinlesion[sample].append(sub) # sample x interval is appended to growing list of those not affected by lesion
         else:
            databysample[sample][chrom][lesionhit[0]][1].append(sub) # this line crashes the loop. See line 154
            #print databysample[sample][chrom][lesionhit[0]][0],sample,sub
         if len(lesionhit) > 1:
            x[23232]

# doesn't crash in Python3
notinvals = defaultdict(list)
for sample in notinlesion.keys():
   temp = notinlesion[sample]
   sumSNPCovoverperD = sum(map(lambda z: z[4]*float(z[5]), temp)) # foreach SNP in those covering lesions, count how many match A allele
   allsnpcov = sum(map(lambda z: z[4], temp)) # count all reads supporting SNP in lesion
   perdnotlesion = sumSNPCovoverperD/float(allsnpcov) # calculate percent A
   numsnps = len(temp) # calculate number of SNPs in region.
   notinvals[sample] = [numsnps, round(perdnotlesion,3)]

o = open(opt.o, 'w')
# Crashes when running through loop. Probably a Python2 vs Python3 issue
o.write(leshead[:-1]+'\t' + '\t'.join(['NumSNPLes', 'PerALes', 'NumSNPsNotLes', 'PerANotLes'])+'\n') # prints 79 to console and doesn't write anything to file. Have to close before stuff shows up.
ind1 = sorted(databysample.keys()) # does wrapping around a sorted() function work? Yes it does.
# ind1.sort() # doesn't work, is a dict_keys object instead of a sortable list. No longer necessary with edits to line above.

# Line 16 of loop crashes. TypeError can only concatenate list (not "map") to list.
for sample in ind1:
   for chrom in databysample[sample].keys():
      for lesion in databysample[sample][chrom].keys():
         temp = databysample[sample][chrom][lesion]
         if temp[-1] != []:
            sumSNPCovoverperD = sum(map(lambda z: z[4]*float(z[5]), temp[1]))
            allsnpcov = sum(map(lambda z: z[4], temp[1]))
            perdlesion = sumSNPCovoverperD/float(allsnpcov)
            numsnps = len(temp[1])
            oset = [numsnps, round(perdlesion,3)]
         else:
            oset = ['0', '0.0']
         nset = notinvals[sample]
         if nset == []:
            nset = ['0', '.']
         oline = temp[0] + list(map(lambda k: str(k), oset+nset)) # wrap map in list call to correct.
         o.write('\t'.join(oline)+'\n')

o.close()