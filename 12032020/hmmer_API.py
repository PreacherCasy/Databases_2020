#!/usr/bin/env python3

from Bio import SeqIO
import click
import urllib


class SmartRedirectHandler(urllib.request.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return headers


def upload_seq(file:str)->str:
    with open(file, 'r') as handle:
        entry = SeqIO.read(handle, 'fasta')
    return str(entry.seq)


def submit(seq:str, db:str='pfam'):

    opener = urllib.request.build_opener(SmartRedirectHandler())
    urllib.request.install_opener(opener)

    parameters = {
    'hmmdb': f'{db}',
    'seq':f'>Seq\n{seq}'
}
    enc_params = urllib.parse.urlencode(parameters).encode('utf-8')

    submission = urllib.request.Request('https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan', enc_params)

    results_url = urllib.request.urlopen(submission)
    res_params = {'format': 'tsv'}

    enc_res_params = urllib.parse.urlencode(res_params).encode('utf-8')
    modified_res_url = results_url._headers[6][1].replace('results', 'download') + '?' + enc_res_params.decode()

    results_request = urllib.request.Request(modified_res_url)
    data = urllib.request.urlopen(results_request)


    return data.read().decode('utf-8')


@click.command()
@click.option('-i', prompt='Provide an input FASTA file', type=click.Path(exists=True))
@click.option('-db', default='pfam')
@click.option('-o', type=click.Path(exists=False), default=None)

def submit_hmmscan_job(i:str, db:str, o:str):
    seq = upload_seq(i)
    res = submit(seq, db)
    if o is not None:
        with open(o, 'a') as handle:
            handle.write(res)
    else:
        click.echo(res)

if __name__ == '__main__':
   submit_hmmscan_job()

