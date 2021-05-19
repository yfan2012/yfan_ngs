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


class insertsinfo:
    '''
    store alignment info
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
                if read.is_read1:
                    if read.reference_name == 'insert':
                        self.r1_ins.append([read.reference_start, read.reference_end, read.is_reverse, read.query_alignment_start, read.query_alignment_end])
                    elif read.reference_name == 'cre':
                        self.r1_cre.append([read.reference_start, read.reference_end, read.is_reverse, read.query_alignment_start, read.query_alignment_end])
                elif read.is_read2:
                    if read.reference_name == 'insert':
                        self.r2_ins.append([read.reference_start, read.reference_end, read.is_reverse, read.query_alignment_start, read.query_alignment_end])
                    elif read.reference_name == 'cre':
                        self.r2_cre.append([read.reference_start, read.reference_end, read.is_reverse, read.query_alignment_start, read.query_alignment_end])
        return('done')
            
                

def find_inserts(bamfile):
    '''
    find the cre insert location for each read
    [readname, readpair, creposition, orientation_match, called_by_both_pairs]
    doesn't bother with any read that doesn't have exactly 1 cre alignment and 1 insert alignment
    gets insert position by looking at shared cre/insert alignment bounds
    checks for alignment orientation of both alignments
    '''
    cre_positions=[]
    
    ##read in bam
    bam=pysam.AlignmentFile(bamfile, 'rb')
    baminfo=pysam.IndexedReads(bam)
    baminfo.build()
    
    ##get read names
    read_ids_all=[]
    for readrecord in bam.fetch():
        read_ids_all.append(readrecord.query_name)
    read_ids=set(read_ids_all)
    
    ##get info on each read
    interfaces=[]
    for name in read_ids:
        readrecord=baminfo.find(name)    
        ins=insertsinfo(name, readrecord)
        if len(ins.r1_ins)>0 and len(ins.r1_cre)>0:
            interfaces.append(ins)
        if len(ins.r2_ins)>0 and len(ins.r2_cre)>0:
            interfaces.append(ins)
            
    ##get positions, check that cre and ins on a single reads goes in the same direction
    cre_positions=[]
    for i in interfaces:
        if len(i.r1_ins)==1 and len(i.r1_cre)==1:
            if i.r1_ins[0][2]==i.r1_cre[0][2]:
                cre_positions.append([i.readname, 'r1']+i.r1_cre[0]+i.r1_ins[0])
        if len(i.r2_ins)==1 and len(i.r2_cre)==1:
            if i.r2_ins[0][2]==i.r2_cre[0][2]:
                cre_positions.append([i.readname, 'r2']+i.r2_cre[0]+i.r2_ins[0])
                
    return cre_positions



def filter_cre(cre_positions):
    '''
    filter down cre_positions
    '''
    filt_cre={}
    for i in cre_positions:
        if i[0] not in filt_cre:
            filt_cre[i[0]]=i[1:12]
    return filt_cre
                
    
def main():
    args=parseArgs()
    cre_positions=find_inserts(args.bam)
    filt_cre=filter_cre(cre_positions)
    with open(args.out, 'w') as f:
        for i in filt_cre:
            info=filt_cre[i]
            if args.verbose==reu:
                towrite=[i]+info
            else:
                ##if left it was the left insert
                if filt_cre[i][6]<10:
                    towrite=[i, info[0], info[2], info[6]]
                else:
                    towrite=[i, info[0], info[1], info[7]]
            f.write(','.join([str(j) for j in towrite])+'\n')
    f.close()

    
if __name__ == "__main__":
        main()
