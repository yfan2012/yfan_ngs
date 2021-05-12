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

    ##analyze each read
    for name in read_ids:
        readrecord=baminfo.find(name)    
        r1={'cre':[], 'insert':[]}
        r2={'cre':[], 'insert':[]}
        for read in readrecord:
            if not read.is_unmapped:
                if read.is_read1:
                    if read.reference_name == 'insert':
                        r1['insert'].append([read.reference_start, read.reference_end, read.is_reverse, read.query_alignment_start, read.query_alignment_end])
                    elif read.reference_name == 'cre':
                        r1['cre'].append([read.reference_start, read.reference_end, read.is_reverse, read.query_alignment_start, read.query_alignment_end])
                elif read.is_read2:
                    if read.reference_name == 'insert':
                        r2['insert'].append([read.reference_start, read.reference_end, read.is_reverse, read.query_alignment_start, read.query_alignment_end])
                    elif read.reference_name == 'cre':
                        r2['cre'].append([read.reference_start, read.reference_end, read.is_reverse, read.query_alignment_start, read.query_alignment_end])
        ##calculate cre position if the read had one alignment to cre and one to insert
        ##if cre/insert boundary doesn't line up, just give it a None position
        creposr1=None
        creposr2=None
        orientr1=None
        orientr2=None
        if len(r1['insert'])==1 and len(r1['cre'])==1:
            if r1['insert'][0][3] == r1['cre'][0][4]:
                creposr1=r1['cre'][0][1]
            elif r1['insert'][0][4] == r1['cre'][0][3]:
                creposr1=r1['cre'][0][0]
                
            orientr1=r1['cre'][0][2]==r1['insert'][0][2]
            
        if len(r2['insert'])==1 and len(r2['cre'])==1:
            if r2['insert'][0][3] == r2['cre'][0][4]:
                creposr2=r2['cre'][0][1]
            elif r2['insert'][0][4] == r2['cre'][0][3]:
                creposr2=r2['cre'][0][0]
                
            orientr2=r2['cre'][0][2]==r2['insert'][0][2]
            
        ##only write out one read if the cre call was the same
        if creposr1==creposr2 and creposr1!=None:
            cre_positions.append([name, 'read1', creposr1, orientr1, True])
        else:
            cre_positions.append([name, 'read1', creposr1, orientr1, False])
            cre_positions.append([name, 'read2', creposr2, orientr2, False])
            
    return cre_positions

            

    
def main():
    args=parseArgs()
    cre_positions=find_inserts(args.bam)
    with open(args.out, 'w') as f:
        for i in cre_positions:
            f.write(','.join([str(j) for j in i])+'\n')
    f.close()

    
if __name__ == "__main__":
        main()
