#!/usr/bin/env python3


from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Entrez
from hmmer_API import upload_seq
import click

Entrez.email, Entrez.tool = 'yu.malovichko@arriam.ru', 'Assignment 2'

def blaster(query:str, entrez:str='Escherichia coli[ORGN]'):
    seq = upload_seq(query)
    result_handle = NCBIWWW.qblast("blastp", "nr", seq, entrez_query=entrez,\
 url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi')
    records = list(NCBIXML.read(result_handle).alignments)
    print('Alignment complete')
    print(len(records))
    for al in records:
        print(al)

blaster('dummy.fa')
