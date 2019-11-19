import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --sync input.sync > output.af "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'''
''')

#########################################################   CODE   #########################################################################

parser.add_option("--sync", dest="IN", help="Input file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

def sync2freqh(x,ALL):
    ''' convert string in SYNC format to dictionary of freqencies where x is a string in sync format'''
    from collections import defaultdict as d
    if x==".:.:.:.:.:." or x=="0:0:0:0:0:0":
        return "na","na"
    nuc=["A","T","C","G"]
    counts=[int(X) for X in x.split(":")[:4]]
    if sum(counts)==0:
        return ({"A":0.0,"T":0.0,"C":0.0,"G":0.0},0)
    CO=dict(zip(*[nuc,counts]))
    for k,v in CO.items():
        if k not in ALL:
            del CO[k]
    h=d(float)
    if sum(CO.values())==0:
        return "na","na"
    for k,v in CO.items():
        h[k]=v/float(sum(CO.values()))

    return h,sum(CO.values())

def all_alleles(v):
    ''' returns most common strings'''
    from collections import Counter
    nuc=v.replace("N","")
    if len(nuc)==0 or len(set(nuc))<2:
        return "NA"
    return zip(*Counter(nuc).most_common())[0]

def sync2string(x):
    ''' convert sync format to string of nucleotides  where x is a string in sync format '''
    string=""
    alleles=["A","T","C","G"]
    if x==".:.:.:.:.:." or x=="0:0:0:0:0:0":
        return "na"
    ah=dict(zip(alleles,[int(X) for X in x.split(":")[:4]]))
    for k,v in ah.items():
        string+=v*k
    return string

def load_data(x):
  ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
  import gzip
  if x=="-":
      y=sys.stdin
  elif x.endswith(".gz"):
      y=gzip.open(x,"rt", encoding="latin-1")
  else:
      y=open(x,"r", encoding="latin-1")
  return y

for l in load_data(options.IN):

    # split in columns
    a=l.rstrip().split()
    ID=a[0]+":"+a[1]

    ## make long string of all summed-up alleles across all populations
    AlleleStrings=""
    for pop in a[3:]:
        Counts=sync2string(pop)
        if Counts=="na":
            continue
        AlleleStrings+=Counts

    ## identify two most common alleles
    MA=all_alleles(AlleleStrings)[:2]

    ## skip if position invariate
    if MA=="NA":
        continue
    fl=[]

    ## calculate and print frequencies of major alleles
    for pop in a[3:]:
        FH=sync2freqh(pop,MA)
        if "na" in FH:
            fl.append("NA")
            continue
        fl.append(FH[0][MA[0]])
    print("\t".join(a[:2])+"\t"+"/".join(MA[:2])+"\t"+"\t".join([str(x) for x in fl]))
