import pysam
import argparse
import numpy

def parseArgs():
    '''
    function to parse args in main
    '''
    parser=argparse.ArgumentParser(description='find where the insert occurs')
    parser.add_argument('-r', '--ref', type=str, required=True,
                        help='reference aligned to')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='aligned reads')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output csv file')
    args=parser.parse_args()
    return args


def fasta_dict(reffile):
    '''
    from fastafile, get a dict of seqname:seq
    '''
    fa=pysam.FastaFile(reffile)
    tigs=fa.references
    fastadict={ x:fa.fetch(x) for x in tigs }
    return fastadict

def find_inserts(bamfile):
    '''
    check for edges of cre alignments. if it lines up with any insert alignments, note the cre position
    readname, read1, read_idx_ins, cre_idx_ins, 
    '''
insert_positions=[]

bam=pysam.AlignmentFile(bamfile, 'rb')
baminfo=pysam.IndexedReads(bam)
baminfo.build()

read_ids_all=[]
for readrecord in bam.fetch():
    read_ids_all.append(readrecord.query_name)
read_ids=set(read_ids_all)
    
for name in read_ids:
    readrecord=baminfo.find(name)
    ''
    
r1=np.array([])
r2=np.array([])
for read in readrecord:
    if not read.is_unmapped:
        if read.is_read1:
                np.append(r1, [read.reference_name, read.reference_start, read.reference_end, read.query_alignment_start, read.query_alignment_end])
                r1['cre'].extend([read.query_alignment_start, read.query_alignment_end])

        elif read.is_read2:
            if read.reference_name == 'cre':
                r2['cre'].extend([read.query_alignment_start, read.query_alignment_end])
            elif read.reference_name == 'insert':
                r2['insert'].extend([read.query_alignment_start, read.query_alignment_end])
    
                

            
        
##print(' '.join([read.reference_name, read.query_name, str(read.is_read1),  str(read.is_reverse), str(read.reference_start), str(read.reference_end), str(read.query_alignment_start), str(read.query_alignment_end), read.cigarstring]))


    
        
        
    
def main():
    args=parseArgs()

    pysam
    ##set up parallel jobs for each read
    manager=mp.Manager()
    q=manager.Queue()
    pool=mp.Pool(args.threads)
    watcher=pool.apply_async(listener, (q, args.out))


    ##submit jobs that will accumulate in the watcher
    jobs=[]
    for i in fast5s:
        job=pool.apply_async(read_mods, (i, args.bam, motifpos, motifs, args.pos, q))
        jobs.append(job)

        
    for job in jobs:
        job.get()
    q.put('Done now, ty 4 ur service')
    pool.close()
    pool.join()
