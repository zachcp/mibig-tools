import json
import click
#import glob

#import pandas as pd
from Bio import SeqIO, SeqUtils
#from toolz import assoc, dissoc, thread_first,partial
#from multiprocessing import pool


@click.command()
@click.option("--jsonfile",type=click.File('r'), prompt=True, help="Mibig json file")
@click.option("--gbkfile", type=click.File('r'), prompt=True, help="Mibig GBK file")
def general_cluster_data(jsonfile, gbkfile):
    """
    return a list of dictionaries with clusterinformation
    """
    clusterdata = json.load(open(jsonfile))
    gbks = SeqIO.parse(open(gbkfile),'genbank')

    results = []
    for gbk in gbks:
        results.append(
            {"biosyn_class": ';'.join(clusterdata['general_params']['biosyn_class']),
             "mibig_accession": gbk.annotations['accessions'][0],
             "version": gbk.annotations['sequence_version'],
             "organism": gbk.annotations['source'],
             "organism_taxa": ';'.join(gbk.annotations['taxonomy']),

             "cluster_length":  len(gbk.seq),
             "cluster_GC": SeqUtils.GC(gbk.seq),
             "accession": gbk.annotations['comment'].split()[-1]})
    return results

# jsonfiles = glob.glob("mibigdata/*.json")
#
# results = []
# for jsonfile in jsonfiles:
#     gbkfile = jsonfile[:-4] + "1.final.gbk"
#     results.append(general_cluster_data(jsonfile,gbkfile))
#
#
# df = pd.DataFrame([f for r in results for f in r])
# df.to_csv("Clusters.csv",index=False)
#general_cluster_data('mibig.secondarymetabolites.org/repository/BGC0000001/BGC0000001.json',
#                     'mibig.secondarymetabolites.org/repository/BGC0000001/BGC0000001.1.final.gbk')