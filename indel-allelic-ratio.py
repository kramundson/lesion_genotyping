#! /usr/bin/env python2.6
import sys, math, os, time
from optparse import OptionParser
from collections import defaultdict

#Want the 200kb leading up to leseion, the first 200kb, the whole lesion, the last 200kb in lesion, and the following 200kb
#gen the in gene for all five of these, if borcder has bases, use it, if at either end, use periods as placeholders
#use percentage for all 5.



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
slist = set()

f = open(opt.a)
header = f.readline()
head = header[:-1].split('\t')
temp = head[5::6]
subnames = map(lambda x: x.replace("Snptype-", ''), temp)
#t = split(head[5:],6)
ct = 0
for l in f:
   ct+=1
   if ct % 10000 == 8:
      print ct, p1
#      if ct > 50000:
#         break
   x = l[:-1].split('\t')
   slist.add(len(x))
   p1 = x[1:5]
   rest = split(x[5:],6)
   chrom = x[0]
   if chrom not in chromlist:
      chromlist.append(chrom)
   for snum in range(len(rest)):
      oname = subnames[snum]
      vals = rest[snum]
      if vals == ['.', '.', '.', '.', '.', '.']:
         continue
      #newname = names[oname]
      newname = oname.replace("upop-",'')
      if vals[4] == '.':
         vals[4] = 0
      vals = [int(vals[3]), float(vals[4])]
      piledat[newname][chrom].append(p1+vals)

f.close()

databysample = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
# load in lesions by chromosome
f = open(opt.l)
leshead = f.readline()
for l in f:
   x = l[:-1].split('\t')
   sample = x[0]
   #####
   for ck in piledat.keys():
      if x[0].replace('_','-')+'-' in ck:
         sample = ck[:]
         #print x[:4], ck
         chrom = x[1]
         start = x[2]
         end = x[3]
         for z in x:
            if z == '':
               x[x.index(z)] = '.'
         databysample[sample][chrom][start+'-'+end] = [[sample]+x[1:],[]]
         #print sample, chrom, start, end, [sample]+x[1:]

f.close()


notinlesion = defaultdict(list)

for sample in piledat.keys():
   for chrom in piledat[sample].keys():
      temp = piledat[sample][chrom]
      chromlesions = databysample[sample][chrom].keys()
      for sub in temp:
         lesionhit = filter(lambda z: int(sub[0]) >= int(z.split('-')[0]) and int(sub[0]) <= int(z.split('-')[1]), chromlesions)
         if lesionhit == []:
            notinlesion[sample].append(sub)
         else:
            databysample[sample][chrom][lesionhit[0]][1].append(sub)
            #print databysample[sample][chrom][lesionhit[0]][0],sample,sub
         if len(lesionhit) > 1:
            x[23232]

notinvals = defaultdict(list)
for sample in notinlesion.keys():
   temp = notinlesion[sample]
   sumSNPCovoverperD = sum(map(lambda z: z[4]*float(z[5]), temp))
   allsnpcov = sum(map(lambda z: z[4], temp))
   perdnotlesion = sumSNPCovoverperD/float(allsnpcov)
   numsnps = len(temp)
   notinvals[sample] = [numsnps, round(perdnotlesion,3)]

o = open(opt.o, 'w')
o.write(leshead[:-1]+'\t'.join(['NumSNPLes', 'PerALes', 'NumSNPsNotLes', 'PerANotLes'])+'\n')
ind1 = databysample.keys()
ind1.sort()
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
         oline = temp[0] + map(lambda k: str(k), oset+nset)
         o.write('\t'.join(oline)+'\n')


o.close()



         
         
   

































