#!/usr/bin/env python

import argparse
import pandas
from Bio import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Do an in silico intergration of DNA elements")
parser.add_argument('--sites', '-s',
                    help=('A file with sites in which to integrate DNA ',
                          'elements. Must contain headers with the provided',
                          'column names'))
parser.add_argument('--names', '-n', nargs=2,
                    help=('column names for start and end position with site'
                          'of intergration in the center'))
parser.add_argument('--ref', '-r', help='reference genome (fasta)')
parser.add_argument('--insert', '-i', help='sequence of the insert (fasta)')
parser.add_argument('--chain-out', '-c',
                    help=('output chain file linking locations in in silico'
                          'genome back to the reference'))
parser.add_argument('--genome-out', '-g',
                    help='output file for the in silico edited genome (fasta)')
parser.add_argument('--gff-out', '-f',
                   help='output file for a gff with integrations on insilico genome')
parser.add_argument('--overhang', '-o',
                    help='sequence of recognition site (e.g TTAA, TA)')
parser.add_argument('--ori', choices=['fwd', 'rev'],
                     help='orientation ')

args = parser.parse_args()

with open(args.insert, "rU") as ref_in:
    home_list = [seq for seq in SeqIO.parse(args.insert, "fasta")]

home_seq = home_list[0].seq


overhang = args.overhang.lower() if args.overhang not in ['-', '0'] else ''
if args.ori == 'fwd':
    fwd_seq = home_seq.lower()
else:
    fwd_seq = str(Seq.reverse_complement(home_seq)).lower()

if not fwd_seq.startswith(overhang):
    fwd_seq = ''.join((overhang, fwd_seq))
o_revcomp = str(Seq.reverse_complement(overhang))
if not fwd_seq.endswith(o_revcomp):
    fwd_seq = ''.join((fwd_seq, o_revcomp))

rev_seq = str(Seq.reverse_complement(fwd_seq))


insertion_pd = pandas.read_csv(args.sites, sep='\t')



chain_pattern = ('chain {ref_len} {seqname} {ref_len} + 0 {ref_len} {seqname} '
                 '{new_len} + 0 {new_len} {id}')
chain_id = 1
gff_pattern = ('{seqname}\ttagmap_hopping\ttransposon\t{start}\t{end}\t'
               '{score}\t{strand}\t.\tID={id}')

ref_in = open(args.ref, "rU")
genome_out = open(args.genome_out, "w")
chain_out = open(args.chain_out, "w")
gff_out = open(args.gff_out, "w")

for refseq in SeqIO.parse(ref_in, "fasta"):
    if refseq.id in insertion_pd['seqnames'].values:
        seq_insertions = insertion_pd[insertion_pd['seqnames']==refseq.id]

        start_list = seq_insertions[args.names[1]].tolist()
        if overhang=='':
            start_list = [start - 1 for start in start_list] # because 0 length bed regions are not accepted
        
        start_list.insert(0,0)
        end_list = seq_insertions[args.names[0]].tolist()
        end_list.append(len(refseq))
        sub_seq_list = [refseq.seq[start:end] for start, end
                        in zip(start_list, end_list)]


        insert_list = [fwd_seq if s =='+' else rev_seq
                       for s in seq_insertions['strand'].tolist()]
        for i in range(0,len(insert_list)):
            sub_seq_list.insert(i*2+1, insert_list[i])
        new_seq=''.join([str(s) for s in sub_seq_list])
        new_record = SeqIO.SeqRecord(id=refseq.id, description=refseq.id,
                                     seq=Seq.Seq(new_seq))
        SeqIO.write(new_record, genome_out, "fasta")
        chain_str = chain_pattern.format(ref_len = len(refseq),
                                         seqname = refseq.id,
                                         new_len = len(new_seq),
                                         id = chain_id)
        print(chain_str, file=chain_out)
        for start, end in zip(start_list, end_list):
            if (end!=len(refseq)):
                print('%i %i %i' % (end-start, len(overhang), len(fwd_seq)),
                      file=chain_out)
            else:
                print(len(refseq) - start, file=chain_out)
        i = 0
        for j, row in seq_insertions.iterrows():
            new_start = row[args.names[0]] + i * (len(fwd_seq) - len(overhang)) + 1
            i += 1
            new_end = new_start + len(fwd_seq)
            p = 0 if 'p_adj' not in row.index else row['p_adj']
            print(gff_pattern.format(seqname = row['seqnames'], start=new_start,
                                     end=new_end, score=p,
                                     strand=row['strand'], id=j), file=gff_out)
    else:
        SeqIO.write(refseq, genome_out, "fasta")
        chain_str = chain_pattern.format(ref_len = len(refseq),
                                         seqname = refseq.id,
                                         new_len = len(refseq),
                                         id = chain_id)
        print(chain_str, file=chain_out)
        print(len(refseq), file=chain_out)
    chain_id += 1

ref_in.close()
genome_out.close()
chain_out.close()
gff_out.close()
