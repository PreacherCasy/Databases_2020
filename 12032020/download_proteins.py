#!/usr/bin/env python3

from Bio import Entrez
from Bio import ExPASy
from Bio import PDB
from Bio import SeqIO
from Bio.ExPASy import ScanProsite
from collections import defaultdict
from urllib import request
from xml.etree import ElementTree
import click
import subprocess

Entrez.email, Entrez.tool = 'yu.malovichko@arriam.ru', 'Assignment2'

def get_NCBI_ID(query:str='(PCNA[title] OR "proliferating cell nuclear antigen"[title]) AND "Homo Sapiens"[Organism] AND srcdb_refseq[PROP]')->str:
    """Get an NCBI accession numer by a valid Entrez query"""
    handle = Entrez.esearch(db='protein', term=query, retmode='xml', retmax=100000, sort='relevance')
    res = Entrez.read(handle)['IdList'][0]
    
    return res

def fetch_NCBI_ref(query:str='33239451'):
    """Fetch a RefSeq entry by the accession number"""
    entry = Entrez.efetch(db='protein', id=query, rettype='fasta', retmode='txt')
    
    return SeqIO.read(entry, 'fasta')


def get_UniProt_ID(query:str='33239451')->str:
    """Get a UniProt (SwissProt if available) accession from a previously fetched RefSeq entry"""
    handle = Entrez.elink(dbfrom='protein', db='protein', id=query)
    links = Entrez.read(handle)
    for link in links[0]['LinkSetDb']:
        if link['LinkName'].find('uniprot') > -1:
            res = link['Link'][0]['Id']
    res_handle = Entrez.esummary(db='protein', id=res, retmode='xml', report='full')
    res_sum = Entrez.read(res_handle)
    uni_id = res_sum[0]['Caption']

    return uni_id

def fetch_UniProt(query:str='33239451', flag:str='full')->str:
    """Fetch a UniProt entry (SwissProt if available)"""
    with request.urlopen(f'http://www.uniprot.org/uniprot/{query}.xml') as handle:
        record = SeqIO.read(handle, 'uniprot-xml')

    if flag == 'subunit':
        return '\n'.join(record.annotations['comment_subunit'][0].split('. '))
    elif flag == 'description':
        return record.annotations['comment_function'][0]	
    elif flag == 'full':
        return record

def download_ProSite_motifs(query:str='P12004')->defaultdict:
    """Performs ExPASy ProSite search for molecular signatures of the query protein. Accepts UniProt ID as a sole argument"""
    in_handle = ScanProsite.scan(seq=query)
    reader = ScanProsite.read(in_handle)
    storage_dict = defaultdict(dict)
    for motif in reader:
        storage_dict[motif['signature_id']]['start'] = motif['start']
        storage_dict[motif['signature_id']]['stop'] = motif['stop']

    return storage_dict

def get_PDB_ID(query:str='33239451', num:int=0)->str:
    """Get a PDB model ID for a provided NCBI accession number"""
    handle = Entrez.elink(dbfrom='protein', db='protein', id=query)
    links = Entrez.read(handle)
    for link in links[0]['LinkSetDb']:
        if link['LinkName'].find('uniprot') > -1:
            res = link['Link'][0]['Id']
    res_handle = Entrez.efetch(db='protein', id=res, retmode='txt', rettype='gb')
    res_sum = SeqIO.read(res_handle, 'gb')
    attrs = res_sum.annotations['db_source'].split(', ')
    pdbs = [x for x in attrs if x.find('PDBsum') == 0]
    pdbs = list(map(lambda x: x.replace('PDBsum:', ''), pdbs))

    if num >= len(pdbs):
        num = len(pdbs) - 1

    return pdbs[num]


def download_PDB_struct(query:str='1AXC'):
    """Obtain PDB structural models"""
    pdbl = PDB.PDBList()
    pdbl.retrieve_pdb_file(query, pdir='.')


@click.command()
@click.option('-eq', help='Entrez query to search a protein of interest by', default='(PCNA[title] OR "proliferating cell nuclear antigen"[title]) AND "Homo Sapiens"[Organism] AND srcdb_refseq[PROP]')
@click.option('-up', type=click.Choice(['full', 'description', 'subunit']), help='Information inferred from the matching UniProt entry. Allowed options are: '+\
'full (fetch the full UniProt entry; does not stack with pu option, use stream redirection instead); description (extract only data related to protein function); subunit (extract only data related to monomer protein structure and interactions', default=None)
@click.option('-p', help='Path to protein sequence file', default=None)
@click.option('-pu', help='Path to UniProt entry file', default=None)
@click.option('-pro', help='Fetch ProSite motifs (results are directed to stdout', is_flag=True)
@click.option('-pdb', help='Download the PDB structural model', is_flag=True)
@click.option('-num', help='Number of the model to download', default=0)


def fetch_data(eq:str, up:str, p:str, pu:str, pro:bool, pdb:bool, num:int):
    res = get_NCBI_ID(eq)
    entry = fetch_NCBI_ref(res)
    if p is not None:
        with open(p, 'a') as handle:
            SeqIO.write(entry, handle, 'fasta')
    if up is not None:
        unip = fetch_UniProt(get_UniProt_ID(res), up)
        if up != 'full':
            if pu is not None:
                with open(pu, 'a') as handle:
                    handle.write(unip)
        else:
            print(unip)
    if pro:
        print(download_ProSite_motifs(get_UniProt_ID(res)))
    if pdb:
        struct = get_PDB_ID(res, num)
        download_PDB_struct(struct)

if __name__ == '__main__':
    fetch_data()



