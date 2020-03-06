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

def fetch_relatives(entry:SeqRecord.SeqRecord):
    seq = entry.seq
    result_handle = NCBIWWW.qblast("blastn", "refseq_genomes", seq)
    records = list(NCBIXML.read(result_handle).alignments)
    for record in records:
        print(record.title)

def extract_CDS(i:str):
    from converter import convert
    for file in os.listdir(i):
        converter(i, i.replace('.gb', '.fa'))


def brew_db():
    pass

def blast_n_sort():
    pass

#@click.command
#print(return_XML())

ref = download_ref('ref.fasta')
print('Written to ref.fasta')
fetch_relatives(ref)
#mine_references()
