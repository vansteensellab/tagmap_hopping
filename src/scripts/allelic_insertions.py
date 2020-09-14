#!/usr/bin/env python

import argparse
import pandas
import subprocess
import io
from Bio import pairwise2
from Bio import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Do an in silico intergration of DNA elements")
parser.add_argument('--alleles', '-a', nargs=2,
                    help=('VCF file with allelic information'))
parser.add_argument('--vcf', '-v', help=('VCF output file for bam consensus'))
parser.add_argument('--bam', '-b', nargs=2,
                    help=('bam files with tagmentation'))
parser.add_argument('--ref', '-r', help='reference genome (fasta)')
parser.add_argument('--insert', '-i', help='insertion sites')
parser.add_argument('--label', '-l', nargs=2,
                    help=('labels for the alleles'))
parser.add_argument('--out', '-o',
                    help='output file for the allele calls')
args = parser.parse_args()

insertion_pd = pandas.read_csv(args.insert, sep='\t')


##TODO: make it optional to use this setting.
## this is wheter to change from UCSC annotation to Ensembl annotation
## (from chr1:2000-3000 to 1:2000-3000)
is_grc_vcf = True

command_pileup = ['samtools', 'mpileup', '-v', '-f', args.ref] + args.bam

pileup_sub = subprocess.Popen(command_pileup, stdout=subprocess.PIPE)

command_vcf =  ['bcftools', 'call', '--ploidy', '1', '-mv', '-Oz', '-o', args.vcf, '/dev/stdin']
vcf_sub = subprocess.call(command_vcf, stdin=pileup_sub.stdout)

index_sub = subprocess.call(['bcftools', 'index', args.vcf])


def run_bcftools(vcf, ref, loc, is_grc_vcf):
    command_sam = ['samtools', 'faidx', ref]
    command_bcf = ['bcftools', 'consensus']
    samtools_sub = subprocess.Popen(command_sam + [loc], stdout=subprocess.PIPE)
    if is_grc_vcf:
        sed_sub = subprocess.Popen(["sed", "s/>chr/>/"],
                                   stdin=samtools_sub.stdout,
                                   stdout=subprocess.PIPE)
        in_pipe = sed_sub
    else:
        in_pipe = samtools_sub
    bcf_sub = subprocess.Popen(command_bcf + [vcf], stdin=in_pipe.stdout,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = bcf_sub.communicate()
    for record in SeqIO.parse(io.StringIO(out.decode("utf-8")), "fasta"):
        yield(record)

with open(args.out, 'w') as f_out:
    print('\t'.join(['seqnames', 'start', 'end'] + args.label + ['call']),
          file=f_out)
    for i, row in insertion_pd.iterrows():
        loc_tuple = (row.seqnames, row.start, row.end)
        loc = '%s:%i-%i' % loc_tuple

        new_fa = list(run_bcftools(args.vcf, args.ref, loc, False))[0]

        score_vec = []
        for vcf in args.alleles:
            allelic_fa = list(run_bcftools(vcf, args.ref, loc, is_grc_vcf))[0]
            aln = pairwise2.align.globalms(str(new_fa.seq).upper(),
                                           str(allelic_fa.seq).upper(),
                                           2, -1, -0.5, -0.1,
                                           score_only=True)

            score_vec.append(aln)

        if score_vec[0] > score_vec[1]:
            call = args.label[0]
        elif score_vec[1] > score_vec[0]:
            call = args.label[1]
        else:
            call = 'Ambiguous'
        out_list = list(loc_tuple) + score_vec + [call]
        print('\t'.join(str(o) for o in out_list), file=f_out)
