import json
import click
#import glob

#import pandas as pd
from Bio import SeqIO, SeqUtils
#from toolz import assoc, dissoc, thread_first,partial
#from multiprocessing import pool


# def get_domains(gbk):
#     proteins = [p for p in gbk.features if p.type =="CDS"]
#     get_nucl = lambda x: str(x.extract(gbk).seq)

#     for p in proteins:
#         print p
#         datadict = {"cluster": gbk.annotations['accessions'][0],
#                     "proteinid": p.qualifiers['protein_id'],
#                     "location": p.location,
#                     "nucleotide": get_nucl(p),
#                     "protein": p.qualifiers['translation'][0],
#                     "product": p.qualifiers['product'][0]
#                    }


# t1 = "NRPS/PKS Domain: PKS_Docking_Nterm (4-31). E-value: 8.6e-08. Score: 22.9;"
# t2 = "NRPS/PKS Domain: PKS_KS (35-460). E-value: 3.1e-185. Score: 607.4;"
# t3 = "NRPS/PKS Domain: PKS_AT (570-867). E-value: 4.5e-105. Score: 342.5;"
# t4 = "NRPS/PKS Domain: PKS_KR (1132-1312). E-value: 1.7e-55. Score: 179.0;"
# t5 = "NRPS/PKS Domain: ACP (1425-1497). E-value: 9.5e-25. Score: 78.0;"
# t6 = "NRPS/PKS Domain: PKS_KS (1520-1946). E-value: 1.2e-180. Score: 592.2;"
# t7 = "NRPS/PKS Domain: PKS_AT (2056-2356). E-value: 4.3e-112. Score: 365.6;"


def process_secmet_domain(st, dna, accession, proteinid):
    """
    Convert something that looks like this:
        "NRPS/PKS Domain: PKS_Docking_Nterm (4-31). E-value: 8.6e-08. Score: 22.9;"

       'NRPS/PKS Domain: AMP-binding (30-457). E-value: 1.8e-119. Score: 390.2;'

    Into:
        {'Domain': 'PKS_Docking_Nterm',
        'E-value': 8.6e-08,
        'Nucleotide_End': 93,
        'Nucleotide_Start': 10,
        'Protein_End': 31,
        'Protein_start': 4,
        'Score': 22.9}

    """
    fulldomain = st.split("Domain:")[1].split(".")[0]
    domain = fulldomain.split("(")[0].strip()
    start  = int(fulldomain.split(")")[0].split("(")[1].split("-")[0])
    end    = int(fulldomain.split(")")[0].split("(")[1].split("-")[1])
    evalue = float(st.split('E-value:')[1].split(". Score:")[0].strip())
    score  = float(st.split("Score:")[1].strip()[:-1])

    return {"Accession": accession,
            "Protein_ID": proteinid,
            "Domain": domain,
            "Protein_start" : start,
            "Nucleotide_Start": (start*3)-2,
            "Protein_End": end,
            "Nucleotide_End": end*3,
            "E-value": evalue,
            "Score": score,
            "DNA-Sequence": dna[(start*3)-2 : end*3]}

#@click.command()
#@click.option("--gbk", type=click.File('r'), help="gbkfile")
def process_secmet(gbk):
    "process the protein qualifiers key for secondary metabolism"
    proteins = [p for p in gbk.features if p.type =="CDS" and "sec_met" in p.qualifiers.keys()]

    secmetdata = []
    for p in proteins:
        #set the default values
        clustersubtype = ""
        clustertype= ""
        clusterkind=""
        if 'protein_id' in p.qualifiers.keys():
            proteinid = p.qualifiers['protein_id'][0]
        else:
            proteinid = ""

        # pass through the list of secmet and
        # get the data that applies to all domains
        secmet = p.qualifiers['sec_met']
        for x in secmet:
            if x.startswith("Type:"):
                clustertype =  x.split(':')[1].strip()
            if x.startswith("Kind:"):
                clusterkind =  x.split(':')[1].strip()
            if x.startswith("NRPS/PKS subtype:"):
                clustersubtype =  x.split(':')[1].strip()

        #gather all the per-domain data in one spot
        domaindata = [d for d in secmet if "Domain:" in d]
        domainfunc = partial(process_secmet_domain,
                             #clustype=clustertype,
                             #cluskind=clusterkind,
                             #clussubtype=clustersubtype,
                             dna = str(p.extract(gbk).seq),
                             accession = gbk.annotations['accessions'][0],
                             proteinid = proteinid)

        for d in domaindata:
            secmetdata.append(domainfunc(d))

    return pd.DataFrame(secmetdata)

# gbks = glob.glob("mibigdata/*.final.gbk")
#
# results = []
# for genbank in gbks:
#     #print genbank
#     gbkrecords = open(genbank,'r')
#     for record in SeqIO.parse(gbkrecords,'genbank'):
#         results.append( process_secmet(record))
#
# df = pd.concat(results)
# df.to_csv("Domains.csv", index=False)