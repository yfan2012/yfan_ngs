import pysam
import argparse
from find_insert import fasta_dict

##for testing
##bamfile='/mithril/Data/NGS/projects/dunlop_insert/run3/check/NT296_insplas.sorted.bam'

def parseArgs():
    '''
    function to parse args in main
    '''
    parser=argparse.ArgumentParser(description='count alignment states')
    parser.add_argument('-b', '--bamfile', type=str, required=True,
                        help='aligned reads')
    args=parser.parse_args()
    return args

def get_readinfo(read, crerange):
    if read.is_unmapped:
        return 'unmapped'
    if read.reference_name=='PTKEI-DEST-CRE-LOXP-SFGFP':
        if (read.reference_start > crerange[0] and read.reference_start < crerange[1]) or (read.reference_end > crerange[0] and read.reference_end < crerange[1]):
            return 'cre'
        else:
            return 'plas'
    if read.reference_name=='INSERT':
        return 'ins'
    
def classify_read(r1info):
    uniquenames=list(set(r1info))
    if len(uniquenames)!=len(r1info):
        ##recombination
        if len(uniquenames)==1 and uniquenames[0]=='plas':
            return 8
        elif len(uniquenames)==1 and uniquenames[0]=='cre':
            return 9
        elif len(uniquenames)==1 and uniquenames[0]=='ins' :
            return 10
        else:
            return 11
    elif len(uniquenames)==1:
        if uniquenames[0]=='unmapped':
            return 0
        if uniquenames[0]=='plas':
            return 1
        if uniquenames[0]=='ins':
            return 2
        if uniquenames[0]=='cre':
            return 3
    elif len(uniquenames)==2:
        if 'ins' in uniquenames and 'cre' in uniquenames:
            return 4
        elif 'ins' in uniquenames and 'plas' in uniquenames:
            return 5
        elif 'cre' in uniquenames and 'plas' in uniquenames:
            return 6
    else:
        return 7 

class statesinfo:
    '''
    store info on alignment states
    '''
    def __init__(self, read_ids, baminfo, crerange):
        '''
        states info:
        unmapped - totally unmapped
        plas only - read had one alignment, and it was plasmid, not covering cre
        ins only - read had one alignment, and it was ins
        cre only - read had one alignment, and it touched the cre section of plasmid
        ins/cre - read had two alignments, insert and cre section of plasmid
        ins/plas - read had two alignments, insert and non-cre section of plasmid
        cre/plas - read had two separate alignments, one touched cre and the other was plasmid but did not touch cre
        multiple - read aligned to three different chrs
        plas_recomb - read had multiple separate alignments to plas only
        cre_recomb - read had multiple separate alignments to cre only
        ins_recomb - read had multiple separate alignments to insert only
        other_recomb - read has multiple separate alignments to something, and aligned to another chr
        ##[unmapped, plas only, ins only, cre only, ins/cre, ins/plas, cre/plas, multiple, recombined]
        '''
        self.r1_states=[0,0,0,0,0,0,0,0,0,0,0,0]
        self.r2_states=[0,0,0,0,0,0,0,0,0,0,0,0]
        self.info=self.getstates(read_ids, baminfo, crerange)
    def getstates(self, read_ids, baminfo, crerange):
        for name in read_ids:
            readrecord=baminfo.find(name)
            r1info=[]
            r2info=[]
            for read in readrecord:
                if read.is_read1:
                    r1info.append(get_readinfo(read, crerange))
                if read.is_read2:
                    r2info.append(get_readinfo(read, crerange))
            try:
                state=classify_read(r1info)
                self.r1_states[classify_read(r1info)]+=1
            except TypeError:
                print('weirdness for '+ name)
            try:
                state=classify_read(r2info)
                self.r2_states[state]+=1
            except TypeError:
                print('weirdness for '+ name)
        return('done')
            

def count_states(bamfile):
    '''
    take in bamfile
    get info how many reads are in each alignment state
    '''
    bam=pysam.AlignmentFile(bamfile, 'rb')
    baminfo=pysam.IndexedReads(bam)
    baminfo.build()
    
    all_read_ids=[]
    for i in bam.fetch():
        all_read_ids.append(i.query_name)
    read_ids=set(all_read_ids)

    crerange=check_plas()
    
    states=statesinfo(read_ids, baminfo, crerange)
    return states


def check_plas():
    '''
    check if cre is in the plasmid sequence
    not used in main, just for manually finding cre coordinates in the plasmid
    '''
    reffile1='/mithril/Data/NGS/projects/dunlop_insert/refs/construct1_plas.fa'
    reffile2='/mithril/Data/NGS/projects/dunlop_insert/refs/construct2_plas.fa'
    ref1=fasta_dict(reffile1)
    ref2=fasta_dict(reffile2)
    start1=ref1['PTKEI-DEST-CRE-LOXP-SFGFP'].find(ref1['CRE'])
    start2=ref2['PTKEI-DEST-CRE-LOXP-SFGFP'].find(ref2['CRE'])
    len1=len(ref1['CRE'])
    len2=len(ref2['CRE'])
    crerange=[start1,start1+len1]
    return crerange
    
def main():
    args=parseArgs()
    crerange=check_plas()
    states=count_states(args.bamfile)
    ##print('#unmapped,plas_only,ins_only,cre_only,ins_cre,ins_plas,cre_plas,multiple,recombined')
    print(','.join([str(x) for x in states.r1_states]))
    print(','.join([str(x) for x in states.r2_states]))

if __name__ == "__main__":
    main()
