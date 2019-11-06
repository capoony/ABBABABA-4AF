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

parser.add_option("--AlleleFrequencies", dest="AF", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--order", dest="ORDER", help="column positions of populations H1,H2,H3,H4")
parser.add_option("--window", dest="W", help="window-size")
parser.add_option("--step", dest="S", help="step-size")
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


def D3(x,XL,YL,fXL,fYL,N):
    '''calculate nominator X and demoninator Y a la Martin et al 2014, assumes X=ABBA-BABA'''
    A,B,C,D=x
    if C>B:
        A,C,B,D=x
    X=(1-A)*B*C*(1-D)-A*(1-B)*C*(1-D)
    Y=(1-A)*B*B*(1-D)-A*(1-B)*B*(1-D)
    XL.append(X)
    fXL[N].append(X)
    YL.append(Y)
    fYL[N].append(Y)

def blockJackEven(fx,fy,D):
    ''' see Busing et al. 1999 '''

    ## at first retain windows that don't result in window-wise D="NA"
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

    ## calculate the jackknife estimate and variance of D
    DjackEst=g*Dest-(g-1)*meanJack
    DjackVar=((g-1)/g)*sum([(Di-Dest)**2 for Di in Dmean])

    ## return
    return Dest,DjackEst,DjackVar,Dest/(DjackVar**0.5)

def test0(x,y,z):
    ''' calculate the ratio of x and y if enough SNPs and if y is not 0, then append to list z '''

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
ORDER=[*map(int,options.ORDER.split(","))]

Chr=""

# Store the position of every SNP
SNPwindows=[]

## to calculate a global value
fSNPX1,fSNPX2,fSNPX3=d(list),d(list),d(list)
fSNPY1,fSNPY2,fSNPY3=d(list),d(list),d(list)

## to calculate window-wise values
SNPX1,SNPX2,SNPX3=[],[],[]
SNPY1,SNPY2,SNPY3=[],[],[]

## make list to store window-wise values of D
D1L,D2L,D3L=[],[],[]

out=open(options.OUT,"w")
N=0
for l in load_data(options.AF):
    a=l.rstrip().split()

    ## only keep allele frequencies of populations used for ABBA BABA and retain order
    AFs=[a[x] for x in ORDER]
    if Chr=="":
        Chr=a[0]

    ## continue if one AL missing
    if "NA" in AFs or "na" in AFs:
        continue

    AFs=[*map(float,AFs)]

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
        D3s=test0(SNPX3,SNPY3,D3L)

        if len(SNPX1)!=0:
            START=min(SNPwindows)
            END=max(SNPwindows)
            MIDDLE=sum(SNPwindows)/float(len(SNPwindows))
            LENGTH=END-START

            # only print the values of D if the number of SNPs is equivalent to the number of SNPs defined as the binsize
            out.write(Chr+"\t"+"\t".join([*map(str,[START,MIDDLE,END,LENGTH])])+"\t"+str(len(SNPX1))+"\t"+D1s+"\t"+D2s+"\t"+D3s+"\n")

        # reset the lists
        SNPwindows=[]
        SNPX1,SNPX2,SNPX3=[],[],[]
        SNPY1,SNPY2,SNPY3=[],[],[]
        Chr=a[0]

        SNPwindows.append(int(a[1]))

        # store the numerator and denominator of D
        D1([*map(float,AFs)],SNPX1,SNPY1,fSNPX1,fSNPY1,N)
        D2([*map(float,AFs)],SNPX2,SNPY2,fSNPX2,fSNPY2,N)
        D3([*map(float,AFs)],SNPX3,SNPY3,fSNPX3,fSNPY3,N)

        ## increase counter by one to start a new group
        N+=1

        continue

    SNPwindows.append(int(a[1]))
    Chr=a[0]

    D1([*map(float,AFs)],SNPX1,SNPY1,fSNPX1,fSNPY1,N)
    D2([*map(float,AFs)],SNPX2,SNPY2,fSNPX2,fSNPY2,N)
    D3([*map(float,AFs)],SNPX3,SNPY3,fSNPX3,fSNPY3,N)

D1s=test0(SNPX1,SNPY1,D1L)
D2s=test0(SNPX2,SNPY2,D2L)
D3s=test0(SNPX3,SNPY3,D3L)

if len(SNPX1)!=0:

    START=min(SNPwindows)
    END=max(SNPwindows)
    MIDDLE=sum(SNPwindows)/float(len(SNPwindows))
    LENGTH=END-START

    out.write(Chr+"\t"+"\t".join([*map(str,[START,MIDDLE,END,LENGTH])])+"\t"+str(len(SNPX1))+"\t"+D1s+"\t"+D2s+"\t"+D3s+"\n")

### calculate Jacknife stat
print("Method\tDest\tDjackEst\tDjackVar\tZ")

Z1=blockJackEven(fSNPX1,fSNPY1,D1L)
print("SoraggiEtAl\t"+"\t".join([*map(str,Z1)]))

Z2=blockJackEven(fSNPX2,fSNPY2,D2L)
print("DurandEtAl\t"+"\t".join([*map(str,Z2)]))

Z3=blockJackEven(fSNPX3,fSNPY3,D3L)
print("MartinEtAl\t"+"\t".join([*map(str,Z3)]))
