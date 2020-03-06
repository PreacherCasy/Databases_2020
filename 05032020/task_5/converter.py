#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import click

def convert(gb_file:str, fa_file:str):
    with open(gb_file, 'r') as gb_handle:
        with open(fa_file, 'w') as fa_handle:
            file = SeqIO.read(gb_handle, 'genbank')
            save_list = []
            for feature in file.features:
                if feature.type == 'CDS':
                    save_list.append(SeqRecord(id=feature.qualifiers['product'][0],\
seq=Seq(feature.qualifiers['translation'][0])))
            SeqIO.write(save_list, fa_handle, 'fasta')

@click.command()
@click.option('-i')
@click.option('-o')

def conv(i, o):
    convert(i, o)

if __name__ == '__main__':
    conv()

