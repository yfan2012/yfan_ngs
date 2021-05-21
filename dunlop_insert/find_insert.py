import pysam
import argparse
import numpy

def parseArgs():
    '''
    function to parse args in main
    '''
    parser=argparse.ArgumentParser(description='find where the insert occurs')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='aligned reads')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output csv file')
    parser.add_argument('-r', '--ref', type=str, required=True,
                        help='reference fasta with >cre and >insert seqs only')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose outputs in csv')
    args=parser.parse_args()
    return args



def count_multiples(bamfile, chrom):
    '''
    prelim count of aligments per read
    not used in main, just as sanity check
    '''
    counts=[0,0,0,0]
    bam=pysam.AlignmentFile(bamfile, 'rb')
    baminfo=pysam.IndexedReads(bam)
    baminfo.build()
    read_ids_all=[]
    for readrecord in bam.fetch():
        read_ids_all.append(readrecord.query_name)
    read_ids=set(read_ids_all)
    for name in read_ids:
        readrecord=baminfo.find(name)
        aligns=[0,0]
        for read in readrecord:
            if read.is_read1 and read.reference_name==chrom:
                aligns[0]+=1
            if read.is_read2 and read.reference_name==chrom:
                aligns[1]+=1
        for i in aligns:
            if i <= 2:
                counts[i]+=1
            else:
                counts[3]+=1
    return counts

def real_query_start(ctups):
    '''
    get real query start given hard clipping
    take read.cigartuples and read.query_alignment_start
    '''
    realstart=0
    for i in ctups:
        if i[0]!=0:
            realstart+=i[1]
        else:
            break
    return realstart


class insertsinfo:
    '''
    store alignment info
    for each alignment [refstart, refend, reverse, readstart, readend]
    '''
    def __init__(self, name, readrecord):
        self.readname=name
        self.r1_ins=[]
        self.r1_cre=[]
        self.r2_ins=[]
        self.r2_cre=[]
        self.mapped=self.fillinfo(readrecord)
    def fillinfo(self, readrecord):
        for read in readrecord:
            if not read.is_unmapped:
                realstart=real_query_start(read.cigartuples)
                realend=realstart+read.query_alignment_length
                if read.is_read1:
                    if read.reference_name == 'insert':
                        self.r1_ins.append([read.reference_start, read.reference_end, read.is_reverse, realstart, realend])
                    elif read.reference_name == 'cre':
                        self.r1_cre.append([read.reference_start, read.reference_end, read.is_reverse, realstart, realend])
                elif read.is_read2:
                    if read.reference_name == 'insert':
                        self.r2_ins.append([read.reference_start, read.reference_end, read.is_reverse, realstart, realend])
                    elif read.reference_name == 'cre':
                        self.r2_cre.append([read.reference_start, read.reference_end, read.is_reverse, realstart, realend])
        return('done')
            
                

def find_inserts(bamfile):
    '''
    take alignment info and get info on the interfaces
    '''
    bam=pysam.AlignmentFile(bamfile, 'rb')
    baminfo=pysam.IndexedReads(bam)
    baminfo.build()
    read_ids_all=[]
    for readrecord in bam.fetch():
        read_ids_all.append(readrecord.query_name)
    read_ids=set(read_ids_all)
    interfaces=[]
    for name in read_ids:
        readrecord=baminfo.find(name)    
        ins=insertsinfo(name, readrecord)
        if len(ins.r1_ins)==1 and len(ins.r1_cre)==1:
            interfaces.append(ins)
        if len(ins.r2_ins)==1 and len(ins.r2_cre)==1:
            interfaces.append(ins)
    return interfaces


def filter_direction(interfaces):
    '''
    make sure that cre and ins are facing the same direction
    '''
    cre_positions=[]
    for i in interfaces:
        if len(i.r1_ins)>0 and len(i.r1_cre)>0:
            if i.r1_ins[0][2]==i.r1_cre[0][2]:
                cre_positions.append([i.readname, 'r1']+i.r1_cre[0]+i.r1_ins[0])
        if len(i.r2_ins)>0 and len(i.r2_cre)>0:
            if i.r2_ins[0][2]==i.r2_cre[0][2]:
                cre_positions.append([i.readname, 'r2']+i.r2_cre[0]+i.r2_ins[0])
    return cre_positions


def fasta_dict(reffile):
    '''
    read fasta into dictionary
    '''
    fa=pysam.FastaFile(reffile)
    tigs=fa.references
    fastadict={x:fa.fetch(x) for x in tigs}
    return fastadict


def filter_cre(cre_positions, fastadict):
    '''
    filter down cre_positions
    report only one of the read pairs if they give the same info
    '''
    inslen=len(fastadict['insert'])
    filt_cre={}
    for i in cre_positions:
        if i[0] not in filt_cre:
            if i[7]==0 or i[8]==inslen:
                filt_cre[i[0]]=i[1:]
    return filt_cre


def correct_starts(filt_cre):
    '''
    adjust cre positions if the alignment is ambiguous using reference
    check matches of arbitrary 10bp extended. Should be enough for uniqueness
    assumes that alignment will extend as far as possible - so overlaps will happen if possible
    '''
    newfiltcre={}
    for i in filt_cre:
        if filt_cre[i][6]==0:
            adjust=filt_cre[i][5]-filt_cre[i][9]
            if adjust>0:
                filt_cre[i][2]-=adjust
                print(i, filt_cre[i][0], filt_cre[i][6], adjust)
            newfiltcre[i]=filt_cre[i]
        elif filt_cre[i][7]==985:
            adjust=filt_cre[i][10]-filt_cre[i][4]
            if adjust>0:
                filt_cre[i][1]+=adjust
                print(i, filt_cre[i][0], filt_cre[i][6], adjust)
            newfiltcre[i]=filt_cre[i]
    return newfiltcre



def main():
    args=parseArgs()
    fastadict=fasta_dict(args.ref)
    interfaces=find_inserts(args.bam)
    cre_positions=filter_direction(interfaces)
    filt_cre=filter_cre(cre_positions, fastadict)
    filt_cre_corrected=correct_starts(filt_cre)
    with open(args.out, 'w') as f:
        for i in filt_cre_corrected:
            info=filt_cre_corrected[i]
            if args.verbose:
                towrite=[i]+info
            else:
                ##if left it was the left insert
                if filt_cre_corrected[i][6]==0:
                    towrite=[i, info[0], info[3], info[2], info[6]]
                else:
                    towrite=[i, info[0], info[3], info[1], info[6]]
            f.write(','.join([str(j) for j in towrite])+'\n')
    f.close()

    
if __name__ == "__main__":
        main()
