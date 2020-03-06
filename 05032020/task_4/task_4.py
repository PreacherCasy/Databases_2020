#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Entrez
from Bio import SeqIO
import click
import pandas as pd

Entrez.email, Entrez.tool = 'yu.malovichko@arriam.ru', 'My Custom Script'

def return_XML(query:str, species:str='Homo sapiens'):
    handle = Entrez.esearch(db='nucleotide', term=f'{query}[GENE] AND {species}[ORGANISM]', retmode = 'xml')

    return handle


def return_pivot(query:str, species:str='Homo sapiens'):
    handle = return_XML(query, species)
    uids, accs, lengths = [], [], []
    for uid in Entrez.read(handle)['IdList']:
        uids.append(uid)
        entry = Entrez.read(Entrez.esummary(db='nucleotide', id=uid))[0]
        acc, length = entry['Caption'], entry['Length']
        accs.append(acc)
        lengths.append(length)
    output = pd.DataFrame({'NCBI ID': uids, 'NCBI Accession': accs, 'Sequence Length': lengths})

    return output


def return_seqs(query:str, species:str='Homo sapiens')->list:
    handle = return_XML(query, species)
    id_list = Entrez.read(handle)['IdList']
    ID, save_list = ','.join(id_list), []
    for record in SeqIO.parse(Entrez.efetch(db='nucleotide', id=ID, rettype='fasta', retmode='txt'), 'fasta'):
        save_list.append(record)

    return save_list


def fetch_seqs_by_PMID(query:str)->list:
    handle = Entrez.elink(dbfrom='pubmed', db='nucleotide', id=query)
    uids = Entrez.read(handle)[0]['LinkSetDb'][0]['Link']
    ID, save_list = ','.join(list(map(lambda x: list(x.values())[0], uids))), []
    for entry in SeqIO.parse(Entrez.efetch(db='nucleotide', id=ID, rettype='fasta', retmode='txt'), 'fasta'):
        save_list.append(entry)
    return save_list

@click.command()
@click.option('-q', prompt='Provide a query name', help='Query name (trivial name, NCBI accession or PMID accession (for fetch_seqs mode))')
@click.option('-sp', default='Homo sapiens', help='Query species name (default: Homo sapiens)')
@click.option('-m', type=click.Choice(['d', 'p', 'r', 'f']), help=('Choose mode to launch the script in. Possible values are: ' +
'd (default mode; returns XML search output for given query and species), ' +
'p (pivot table mode; returns a table containing NCBI ID, accession number and sequence length for each match), ' + 
'r (retrieval mode; returns a FASTA file with sequences which matched the provided query), ' +
'f (fetching mode; fetches and retrieves all NCBI Nucleotide links from the paper accessed by the provided PMID'), default='d')
@click.option('-o', help='Path to output file', default=None)

def execute(q:str, sp:str, m:str, o:str):
    query, species = q, sp
    if m == 'd':
        out = Entrez.read(return_XML(query, species))
    elif m == 'p':
        out = return_pivot(query, species)
        if o is not None:
            out.to_csv(o, sep='\t', index=False)
    elif m == 'r':
        out = return_seqs(query, species)
        if o is not None:
            with open(o, 'w') as handle:
                SeqIO.write(out, handle, 'fasta')
    elif m == 'f':
        out = fetch_seqs_by_PMID(query)
        if o is not None:
            with open(o, 'w') as handle:
                SeqIO.write(out, handle, 'fasta')
    if o is None:
        click.echo(out)

    return out


if __name__ == '__main__':
    execute()

