import pysam
import argparse
import gzip
import time
from itertools import repeat
import multiprocessing as mp

def parseArgs():
    '''
    function to parse args in main
    '''
    parser=argparse.ArgumentParser(description='find where the insert occurs')
    parser.add_argument('-1', '--r1', type=str, required=True,
                        help='paired r1 reads')
    parser.add_argument('-2', '--r2', type=str, required=True,
                        help='paired r2 reads')
    parser.add_argument('-3', '--r3', type=str, required=False,
                        help='unpaired r1 reads')
    parser.add_argument('-4', '--r4', type=str, required=False,
                        help='unpaired r2 reads')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output csv file')
    parser.add_argument('-r', '--ref', type=str, required=True,
                        help='reference fasta with >cre and >insert seqs only')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose outputs in csv')
    parser.add_argument('-t', '--threads', type=int, required=True,
                        help='number of threads to use')
    args=parser.parse_args()
    return args


def fasta_dict(reffile):
    '''
    read fasta into dictionary
    '''
    fa=pysam.FastaFile(reffile)
    tigs=fa.references
    fastadict={x:fa.fetch(x) for x in tigs}
    return fastadict


def read_fastq(fastqfile):
    '''
    read fastq
    return dict like {readname:seq}
    '''
    if fastqfile.endswith('.gz'):
        f=gzip.open(fastqfile, 'rt')
        contents=f.read().splitlines()
    else:
        f=open(fastqfile, 'rt')
        contents=f.read().splitlines()
    names=contents[0::4]
    seqs=contents[1::4]
    reads=list(zip(names, seqs))
    return reads


def revcomp(seq):
    '''
    reverse complement a string
    https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases=list(seq)
    newbases=[complement[base] for base in bases]
    revbases=''.join(newbases)[::-1]
    return revbases
    
    
def make_scar_flanks(fastadict, size):
    '''
    figure out what scar sequences to look for at each position
    size is how many bases per side - so size 6 results in a 12mer
    '''
    scarpos={}
    for i in range(0, len(fastadict['cre'])-size-size+1):
        frontcre=fastadict['cre'][i:i+size]
        backcre=fastadict['cre'][i-5+size:i-5+size+size]
        pos=str(i+size)
        scarpos[pos]=[frontcre, backcre]
    return scarpos


def get_scar_seq(readseq, scarstart, scarend):
    '''
    given read seqeunce, and the scar starts and ends
    return the full scar sequence
    '''
    edge1=readseq.index(scarstart)
    edge2=readseq.index(scarend)
    if edge1+len(scarstart)<edge2:
        scarseq=readseq[edge1+len(scarstart):edge2]
        fullscarseq=readseq[edge1:edge2+len(scarend)]
        scarinfo=[scarseq, fullscarseq]
    else:
        scarinfo=[]
    return scarinfo
                        

def scan_for_scars(reads, scarpos, fastadict, pair, verbose):
    '''
    scan through a read for evidence of scar
    takes list [readname, seq]
    takes dict {pos:[frontfwd, frontrev, ...]}
    returns list of lists [name, pos, insert_end, orient, readpair]
    '''
    if verbose:
        print('there are '+ str(len(reads))+ ' reads to process')
    insertsites=[]
    dups=0
    count=0
    start_time=time.time()
    frontins=fastadict['insert'][0:15]
    backins=fastadict['insert'][-15:]
    insseqs=[frontins, backins, revcomp(frontins), revcomp(backins)]
    for i in reads:
        if verbose:
            count+=1
            if count % 25000 == 0:
                elapsed=time.time()-start_time
                print('done '+str(count)+' reads in ' + str(elapsed) + ' seconds')
        name=i[0].split(' ')[0][1:]
        readseq=i[1]
        potentialsites=[]
        for pos in scarpos:
            seqs=scarpos[pos]+[revcomp(x) for x in scarpos[pos]]
            if any(x in readseq for x in seqs): ##hit on one of the potential scar sequences
                if not any(x in readseq for x in insseqs): ##no insert hit
                    if seqs[0] in readseq and seqs[1] in readseq:
                        scarinfo=get_scar_seq(readseq, seqs[0], seqs[1])
                        if len(scarinfo)>1:
                            potentialsites.append([name, pos, scarinfo[0], scarinfo[1], pair])
                    elif seqs[2] in readseq and seqs[3] in readseq:
                        scarinfo=get_scar_seq(readseq, seqs[2], seqs[3])
                        if len(scarinfo)>1:
                            potentialsites.append([name, pos, scarinfo[0], scarinfo[1], pair])
        added=False
        if len(potentialsites)==1:
            insertsites.append(potentialsites[0])
            added=True
        elif len(potentialsites)>1:
            for site in potentialsites:
                if site[2].startswith('TG') and site[2].endswith('CA'):
                    insertsites.append(site)
                    added=True
        if len(potentialsites)>1 and not added:
            if verbose:
                for site in potentialsites:
                    site[4]='multi'
                    insertsites.append(site)
                    dups+=1
                    print(potentialsites)
    if verbose:
        print(str(dups)+' reads supported more than one insert site')
    return insertsites
                
    
def scan_fastq(fqinfo, scarpos, fastadict, L, verbose):
    '''
    run through process of scanning a fastq
    adds insertsites to manager list
    done this way to parallel
    '''
    fastqfile=fqinfo[0]
    pair=fqinfo[1]
    reads=read_fastq(fastqfile)    
    insertsites=scan_for_scars(reads, scarpos, fastadict, pair, verbose)
    L+=insertsites


##for testing
##reffile='/mithril/Data/NGS/projects/dunlop_insert/refs/construct1.fa'
##r1file='/mithril/Data/NGS/projects/dunlop_insert/run3/trimmed/NT278_fwd_paired.fq.gz'
##r2file='/mithril/Data/NGS/projects/dunlop_insert/run3/trimmed/NT278_rev_paired.fq.gz'

def main(reffile, pairedr1, pairedr2, unpairedr1, unpairedr2, outfile, threads, verbose):
    fastadict=fasta_dict(reffile)    
    scarpos=make_scar_flanks(fastadict, 12)

    fqinfo=[[pairedr1, 'r1'], [pairedr2, 'r2']]
    if unpairedr1 is not None:
        fqinfo.append([unpairedr1, 'r1'])
    if unpairedr2 is not None:
        fqinfo.append([unpairedr2, 'r2'])
    
    manager=mp.Manager()
    L=manager.list()
    pool=mp.Pool(threads)
    pool.starmap(scan_fastq, zip(fqinfo, repeat(scarpos), repeat(fastadict),  repeat(L), repeat(verbose)))

    ##ensure unique
    uniquesites={}
    for i in L:
        if i[0] not in uniquesites:
            uniquesites[i[0]]=i[1:]
    
    with open(outfile, 'w') as f:
        for i in uniquesites:
            uniqueinfo=[i]+uniquesites[i]
            f.write(','.join(uniqueinfo)+'\n')
    pool.close()
    pool.join()

    
if __name__ == "__main__":
    args=parseArgs()
    main(args.ref, args.r1, args.r2, args.r3, args.r4, args.out, args.threads, args.verbose)
