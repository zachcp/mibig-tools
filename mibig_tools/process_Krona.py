import os
import glob
import json
import pandas as pd
import click

@click.command()
@click.option("--mibigfolder", type=click.Path(exists=True), prompt=True, help="Mibig data folder")
@click.option("--outfile", type=click.File("w"), prompt=True, help="Mibig data folder")
def process_krona(mibigfolder, outfile):
    """
    example; process_krona("mibigdata", "mibig_class.txt")

    """
    jsonfiles = glob.glob("{}/*.json".format(mibigfolder))

    records = []

    for j in jsonfiles:
        data = json.load(open(j,'r'))
        count = 1
        accession = data['general_params']['mibig_accession']
        biosyn_class = data['general_params']['biosyn_class']
        biosyn = []

        if 'NRP' in data['general_params'].keys():
            try:
                biosyn_subclass =  data['general_params']['NRP']['subclass']
                biosyn.append(biosyn_subclass)
            except:
               pass


        if 'Polyketide' in data['general_params'].keys():
            try:
                biosyn_subclass =  data['general_params']['Polyketide']['pk_subclass']
                biosyn.append(biosyn_subclass)
            except:
                pass

        if 'RiPP' in data['general_params'].keys():
            try:
                biosyn_subclass =  data['general_params']['Polyketide']['ripp_subclass']
                biosyn.append(biosyn_subclass)
            except:
                pass

        if 'Terpene' in data['general_params'].keys():
            try:
                biosyn_subclass =  data['general_params']['Polyketide']['terpene_subclass']
                biosyn.append(biosyn_subclass)
            except:
                pass


        if 'Sacharide' in data['general_params'].keys():
            try:
                biosyn_subclass =  data['general_params']['Polyketide']['saccharide_subclass']
                biosyn.append(biosyn_subclass)
            except:
                pass

        if 'Alkaloid' in data['general_params'].keys():
            try:
                biosyn_subclass =  data['general_params']['Polyketide']['alkaloid_subclass']
                biosyn.append(biosyn_subclass)
            except:
                pass

        if 'Other' in data['general_params'].keys():
            try:
                biosyn_subclass =  data['general_params']['Polyketide']['other_subclass']
                biosyn.append(biosyn_subclass)
            except:
                pass

        records.append({"count": count,
                        "accession": accession,
                        "biosyn_class": "-".join(biosyn_class),
                        "biosyn_subclass": "-".join(biosyn).replace(" ","_") })

    df = pd.DataFrame(records)
    df = df.groupby(["biosyn_class","biosyn_subclass"]).sum().reset_index()
    df[['count','biosyn_class','biosyn_subclass']].to_csv(outfile, sep="\t", index=False)

    return df

