#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Entrez
from Bio import SeqIO
from collections import defaultdict
import click
import pandas as pd
import re

Entrez.email, Entrez.tool = 'yu.malovichko@arriam.ru', 'My Custom Script'

def sort_assemblies(query:str, first:int=None)->pd.core.frame.DataFrame:
    handle = Entrez.esearch(db='assembly', term=f'{query}[Organism] OR {query}[All Fields]', \
 retmode='xml', retmax=20000)
    uids, naive_dict = Entrez.read(handle)['IdList'], defaultdict(int)
    if first is not None:
        uids = uids[:first]
    for uid in uids:
        sum_handle = Entrez.esummary(db='assembly', id=uid, retmode='txt', report='full')
        sum_report = Entrez.read(sum_handle)
        environment = sum_report['DocumentSummarySet']['DocumentSummary'][0]['SpeciesName']
        if environment.find('metagenome') > -1:
            environment = environment.replace(' metagenome', '')
            length = int(re.findall('"total_length" sequence_tag="all">[0-9]*</Stat>', sum_report['DocumentSummarySet']['DocumentSummary'][0]['Meta'])\
[0].lstrip('"total_length" sequence_tag="all">').rstrip('</Stat>'))
            naive_dict[environment] += length
    output = pd.DataFrame({'Environment': list(naive_dict.keys()), 'Size, Gb': list(map(lambda x: x / 1000000, list(naive_dict.values())))})
    output = output.sort_values(by='Size, Gb')
    return output

@click.command()
@click.option('-q', prompt='Enter a query organism name', help='Organism name to be searched for')
@click.option('-mx', help='Integer defining first n target entries to be processed', default=None)
@click.option('-o', help='Path to tsv output', default=None)

def pivot(q, mx, o):
    out = sort_assemblies(q, mx)
    if o:
        out.to_csv(o, sep='\t', index=False)
    else:
        click.echo(out)

if __name__ == '__main__':
    pivot()

