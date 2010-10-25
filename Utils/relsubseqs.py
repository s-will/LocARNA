#!/usr/bin/python

## enumerate the relevant sub sequences for a crossing RNA (dp format)

import re

from optparse import OptionParser


## ============================================================
## function definitions

def parseCmdLine():
    "parse the command line arguments"

    usage = "usage: %prog [options] input_file"
    parser=OptionParser(usage=usage)

    parser.add_option("-p", "--minprob", dest="minprob",
                      default=1e-6,
                      help="minimal accepted probability, used for filtering")

    parser.add_option("-s", "--strategy", dest="strategy",
                      default="right",
                      help="strategy for reduction (left,right,leftright)")

    parser.add_option("-i", "--informat", dest="informat",
                      default="dp",
                      help="input format")

    (options,args)=parser.parse_args()
    
    if len(args) != 1:
        parser.error("incorrect number of arguments")
        
    return (options,args[0])


def readdp(filename,minprob):
    """
    read a dp file (output of RNAfold -p)
    and return tuple of sequence and list of basepairs with probabilities
    """

    input=open(filename,"r")
    Lines=input.readlines()

    state=''
    sequence=''
    pairs=[]
    
    uboxpat=re.compile('(\d+) (\d+) (\d\.\d+) ubox')
    for line in Lines:
        #print line
        if state=='sequence':
            if re.match('\) } def', line):
                state=''
            else:
                line=line[:-2]
                sequence = sequence + line
        
        elif re.match('/sequence.*',line):
            state='sequence'
        else:
            ubox=uboxpat.match(line)
            if ubox:
                uboxgr=ubox.groups()
                p=float(uboxgr[2])*float(uboxgr[2])
                if (p>=float(minprob)):
                    pairs.append((int(uboxgr[0]),int(uboxgr[1]),p))
        
    return sequence,pairs

def readsimple(filename,minprob):
    """
    read a simple format for describing RNAs with crossing secondary structure
    and return tuple of sequence and list of basepairs with probabilities
    """

    input=open(filename,"r")
    Lines=input.readlines()

    state='sequence'
    sequence=''
    pairs=[]
    
    pat=re.compile('(\d+) (\d+)')
    for line in Lines:
        #print line
        if state=='sequence':
            if re.match('^$',line):
                state=''
            else:
                line=line[:-1]
                sequence = sequence + line
        else:
            pair=pat.match(line)
            if pair:
                pairgr=pair.groups()
                pairs.append((int(pairgr[0]),int(pairgr[1])))

    return sequence,pairs



def incident_pairs_from_right(start,end,pairs):
    """
    returns list of basepairs that have right end 'end'
    and left end >= 'start'
    """
    return filter((lambda p: p[0]>=start and p[1]==end),pairs)         

def incident_pairs_from_left(start,end,pairs):
    """
    symmetric version of incident_pairs_from_right
    """
    return filter((lambda p: p[0]==start and p[1]<=end),pairs)         



def decideChoiceKlein(start,end,l_incpairs,r_incpairs,pairs):
    "make the choice for reducing left or right --- much in the way of Klein"
    def max_dist(pairs):
        sizes=[p[1]-p[0] for p in pairs]
        sizes.append(0)
        size=reduce((lambda x,y: max(x,y)), sizes)
        return size

    lsize=max_dist(l_incpairs)
    rsize=max_dist(r_incpairs)
    
    if (lsize<rsize): return 'left'
    else: return 'right'


def count_enclosing_pairs(pos,pairs):
    "count the number of pairs that enclose a position >pos<"
    return len(filter((lambda p: p[0]<=pos and p[1]>=pos),pairs))

def decideChoice(start,end,pairs):
    "make the choice for reducing left or right --- new experimental way"
    
    lsize=len(filter((lambda p: p[0]<=start and p[1]>start),pairs)) #count_enclosing_pairs(start,pairs)
    rsize=len(filter((lambda p: p[0]<end and p[1]>=end),pairs)) #count_enclosing_pairs(end,pairs)
    
    if (lsize<rsize): return 'left'
    else: return 'right'


    
def reduce_subseq(start,end,pairs,seen,strategy):
    "recursively reduce to subsequences from left or right or using sort of Kleins strategy"

    if end<start: return
    if seen.has_key((start,end)):
        # print 'seen',start,end
        return

    seen[(start,end)]=True

    if (start-1,end+1) in pairs:
        print " "*(start-1+2)+"("+"."*(end-start+1)+")"

    print " "*(start-1+3-len(str(start)))+str(start)+"-"*(end-start+1)+str(end)


    ## make a choice, whether to reduce from left or right

    r_incpairs=incident_pairs_from_right(start,end,pairs)
    l_incpairs=incident_pairs_from_left(start,end,pairs)

    if strategy=='right' or strategy=='left': theChoice=strategy
    else:
        # theChoice=decideChoiceKlein(start,end,l_incpairs,r_incpairs,pairs)
        theChoice=decideChoice(start,end,pairs)
    
    if theChoice=='right' :
        for p in r_incpairs:
            reduce_subseq(start,p[0]-1,pairs,seen,strategy)
            reduce_subseq(p[0]+1,p[1]-1,pairs,seen,strategy)
            
        reduce_subseq(start,end-1,pairs,seen,strategy)

    elif theChoice=='left':
        for p in l_incpairs:
            reduce_subseq(p[0]+1,p[1]-1,pairs,seen,strategy)
            reduce_subseq(p[1]+1,end,pairs,seen,strategy)
            
        reduce_subseq(start+1,end,pairs,seen,strategy)

## ============================================================
## main program

(options,filename)=parseCmdLine()


if options.informat=='dp':
    (seq,pairs)=readdp(filename,options.minprob)
else:
    (seq,pairs)=readsimple(filename,options.minprob)
    
seqlen=len(seq)
print "seqlen =", seqlen, ", #pairs =",len(pairs)
print
print "   "+seq


seen={}
reduce_subseq(1,seqlen,pairs,seen,options.strategy)
