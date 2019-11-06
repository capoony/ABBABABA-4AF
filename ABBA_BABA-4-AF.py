import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python3 %prog --input file --output file "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'''
          /\
         /  \
        /\   \
       /  \   \
      /\   \   \
     /  \   \   \
    H1  H2  H3  H4
''')

#########################################################   CODE   #########################################################################

parser.add_option("--AlleleFrequencies", dest="AF", help="Input file; First three columns: Chrom, Pos, MajorAllele/MinorAllele")
parser.add_option("--output", dest="OUT", help="Output prefix")
parser.add_option("--order", dest="ORDER", help="column positions of populations H1,H2,H3,H4 in input file")
parser.add_option("--SNPs", dest="SNP", help="number of SNPS (overrules window and stepsize)",default="NA")

(options, args) = parser.parse_args()
parser.add_option_group(group)

def D1(x,XL,YL,fXL,fYL,N):
    '''calculate nominator X and demoninator Y a la Soraggi et al. 2018, assumes X=BABA-ABBA '''
    A,B,C,D=x
    X=(A-B)*(C-D)
    Y=(A+B-2*A*B)*(C+D-2*C*D)
    XL.append(X)
    fXL[N].append(X)
    YL.append(Y)
    fYL[N].append(Y)

def D2(x,XL,YL,fXL,fYL,N):
    '''calculate nominator X and demoninator Y a la Durand et al 2011, assumes X=ABBA-BABA'''
    A,B,C,D=x
    X=(1-A)*B*C*(1-D)-A*(1-B)*C*(1-D)
    Y=(1-A)*B*C*(1-D)+A*(1-B)*C*(1-D)
    XL.append(X)
    fXL[N].append(X)
    YL.append(Y)
    fYL[N].append(Y)



def blockJackEven(fx,fy,D):
    ''' see Busing et al. 1999 '''

    ## at first retain only windows that don't result in window-wise D="NA"
    fxnew=[]
    fynew=[]
    for i in range(len(D)):
        if D[i]== "NA":
            continue
        fxnew.append(fx[i])
        fynew.append(fy[i])

    ## then calculate overall expected D
    Dest=sum([sum(x) for x in fxnew])/float(sum([sum(y) for y in fynew]))

    # calculate D when removing one group at a time:
    Dmean=[]
    for i in range(len(fxnew)):
        fxprime=fxnew[:i] + fxnew[i+1 :]
        fyprime=fynew[:i] + fynew[i+1 :]
        Dmean.append(sum([sum(x) for x in fxprime])/float(sum([sum(y) for y in fyprime])))

    ## g is the number of independent groups
    g=float(len(fxnew))

    ## caluclate the mean of all group-wise estimates of D
    meanJack=sum(Dmean)/float(len(Dmean))

    ## calculate the jackknife estimate and variance of D and z-score
    DjackEst=g*Dest-(g-1)*meanJack
    DjackVar=((g-1)/g)*sum([(Di-Dest)**2 for Di in Dmean])
    Z=Dest/(DjackVar**0.5)

    ## return
    return Dest,DjackEst,DjackVar,Z

def test0(x,y,z):
    ''' calculate the ratio of x and y if bin contains enough SNPs and if y is not 0, then append to list z '''

    if sum(y)==0 or len(x)!=SNPs:
        xy="NA"
    else:
        xy=str(sum(x)/float(sum(y)))
        if float(xy)>1 or float(xy)<-1:
            xy="NA"
    z.append(xy)
    return xy

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


SNPs=int(options.SNP)
ORDER=[int(x) for x in options.ORDER.split(",")]

Chr=""

# Store the position of every SNP
SNPwindows=[]

## to calculate a global value
fSNPX1,fSNPX2=d(list),d(list)
fSNPY1,fSNPY2=d(list),d(list)

## to calculate window-wise values
SNPX1,SNPX2=[],[]
SNPY1,SNPY2=[],[]

## make list to store window-wise values of D
D1L,D2L=[],[]

## open Outputfile for window-wise estimates
outWindow=open(options.OUT+"_"+str(SNPs)+"SNPs.D","w")

## Counter of SNP groups
N=0

outWindow.write("Chromosome\tStartPos\tAveragePos\tEndPos\tWindowSize\tNoOfSNPs\tSoraggi\tDurand\n")

for l in load_data(options.AF):
    a=l.rstrip().split()

    ## only keep allele frequencies of populations used for ABBA BABA and retain order
    AFs=[a[x] for x in ORDER]
    if Chr=="":
        Chr=a[0]

    ## continue if one AL missing
    if "NA" in AFs or "na" in AFs:
        continue

    AFs=[float(x) for x in AFs]

    ## continue if site is not informative according to Sorragi et al., i.e. if H1==H2 or H3==H4
    if AFs[0]==AFs[1] or AFs[2]==AFs[3]:
        continue

    # condition AFs for the ancestral allele which is more common in the outgroup population
    if AFs[-1]<0.5:
        AFs=[1-x for x in AFs]

    # if the end of a window is reached or a new contig starts, calculate the window-specific D, print it to the outputfile and move on
    if Chr!=a[0] or len(SNPX1)==SNPs:

        #print(Chr,a[0],len(SNPX1),SNPwindows)

        # calculate D
        D1s=test0(SNPX1,SNPY1,D1L)
        D2s=test0(SNPX2,SNPY2,D2L)

        if len(SNPX1)!=0:
            START=min(SNPwindows)
            END=max(SNPwindows)
            MIDDLE=sum(SNPwindows)/float(len(SNPwindows))
            LENGTH=END-START

            # only print the values of D if the number of SNPs is equivalent to the number of SNPs defined as the binsize
            outWindow.write(Chr+"\t"+"\t".join([str(x) for x in [START,MIDDLE,END,LENGTH]])+"\t"+str(len(SNPX1))+"\t"+D1s+"\t"+D2s+"\n")

        # reset the lists
        SNPwindows=[]
        SNPX1,SNPX2=[],[]
        SNPY1,SNPY2=[],[]
        Chr=a[0]

        SNPwindows.append(int(a[1]))

        # store the numerator and denominator of D
        D1([float(x) for x in AFs],SNPX1,SNPY1,fSNPX1,fSNPY1,N)
        D2([float(x) for x in AFs],SNPX2,SNPY2,fSNPX2,fSNPY2,N)

        ## increase counter by one to start a new group
        N+=1

        continue

    SNPwindows.append(int(a[1]))
    Chr=a[0]

    D1([float(x) for x in AFs],SNPX1,SNPY1,fSNPX1,fSNPY1,N)
    D2([float(x) for x in AFs],SNPX2,SNPY2,fSNPX2,fSNPY2,N)

D1s=test0(SNPX1,SNPY1,D1L)
D2s=test0(SNPX2,SNPY2,D2L)

if len(SNPX1)!=0:

    START=min(SNPwindows)
    END=max(SNPwindows)
    MIDDLE=sum(SNPwindows)/float(len(SNPwindows))
    LENGTH=END-START

    outWindow.write(Chr+"\t"+"\t".join([str(x) for x in [START,MIDDLE,END,LENGTH]])+"\t"+str(len(SNPX1))+"\t"+D1s+"\t"+D2s+"\n")


## open Outputfile for genome-wide estimates
outGW=open(options.OUT+"_GenomeWide.D","w")

### calculate Jacknife stat
outGW.write("Method\tDest\tDjackEst\tDjackVar\tZ\n")

Z1=blockJackEven(fSNPX1,fSNPY1,D1L)
outGW.write("SoraggiEtAl\t"+"\t".join([str(x) for x in Z1])+"\n")

Z2=blockJackEven(fSNPX2,fSNPY2,D2L)
outGW.write("DurandEtAl\t"+"\t".join([str(x) for x in Z2])+"\n")
