#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
from Bio import SeqRecord
import click
import os
import subprocess

Entrez.email, Entrez.tool = 'yu.malovichko@arriam.ru', 'My Custom Script'
devnull = open(os.devnull, 'wb')

def download_ref(out:str):

    """Find the newest reference assembly for coronavira and download it"""

    handle = Entrez.esearch(db='nucleotide', term=f'coronavirus[ALL] AND 2020[PDAT] AND srcdb_refseq[PROP] ', \
retmode = 'xml', sort='relevance')
    uids = Entrez.read(handle)['IdList']
    uid = uids[0]
    fetch = Entrez.efetch(db='nucleotide', id=uid, rettype='fasta', retmode='txt')
    assembly = SeqIO.read(fetch, 'fasta')

    with open(out, 'w') as handle:
        SeqIO.write(assembly, out, 'fasta')

    return assembly

def fetch_relatives(entry:SeqRecord.SeqRecord, od:str='.'):

    """Find matching genomes by BLASTing 2020 nCoV reference against the RefSeq database""""

    seq = entry.seq
    result_handle = NCBIWWW.qblast("blastn", "refseq_genomes", seq)
    records = list(NCBIXML.read(result_handle).alignments)
    for record in records:
        query = record.title.split('|')[3]
        rel = SeqIO.read(Entrez.efetch(db='nucleotide', id=query, retmode='txt'), 'genbank')
        SeqIO.write(rel, os.path(os.getcwd(), f'{rel.id}.gb'), 'genbank')

def extract_CDS(i:str):

    """extract coding sequnces and write them in FASTA format"""

    from converter import convert
    for file in os.listdir(i):
        converter(i, i.replace('.gb', '.fa'))

@click.command()
@click.option('-i', help='A reference file', default=None)
@click.option('-o', help='Path to reference', default=None)
@click.option('-od', help='A folder for downloaded sequences')


def prepare_refs():
    if i is not None:
        ref = i
    else:
        ref = download_ref(o)
    fetch_relatives(od)
    extract_CDS(od)

if __name__ == '__main__':
    prepare_refs()
