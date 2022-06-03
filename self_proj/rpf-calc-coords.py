#!/usr/bin/env python
import pysam
import pandas as pd

MINIMUM_CDS_COVERAGE = 0.4

def load_annotations(filename):
    tranno = pd.read_csv(filename, sep='\t',
                         names=['seqname', 'start', 'end', 'name', 'score', 'strand'], index_col=0)
    return tranno.apply(lambda row: (row['start'], row['end']), axis=1).to_dict()

def process(bam, anno):
    anno = load_annotations(anno)

    for aln in pysam.AlignmentFile(bam):
        trid = aln.reference_name
        if trid not in anno:
            continue

        cdsstart, cdsend = anno[trid]

        if cdsstart == cdsend:
            region = 'ncRNA'
            reltostart = aln.reference_start
            reltostop = reltostart - aln.reference_length
        else:
            reltostart = aln.reference_start - cdsstart
            reltostop = aln.reference_start - cdsend

            startborder = aln.reference_length * (MINIMUM_CDS_COVERAGE - 1)
            stopborder = -aln.reference_length * MINIMUM_CDS_COVERAGE

            if reltostart < startborder:
                region = 'UTR5'
            elif reltostop > stopborder:
                region = 'UTR3'
            else:
                region = 'CDS'

        print(aln.query_name, aln.reference_name, aln.reference_start, aln.reference_length,
              region, reltostart, reltostop, sep='\t')

if __name__ == '__main__':
    import sys

    bam = sys.argv[1]
    anno = sys.argv[2]

    process(bam, anno)
