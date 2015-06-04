import json

from Bio import SeqIO, SeqUtils

def general_cluster_data(data):
    """
    return a list of dictionaries with clusterinformation
    """
    jsonfile, gbkfile = data
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
