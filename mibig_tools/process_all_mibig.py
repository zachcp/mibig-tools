
import os
import glob
import pandas as pd
import click
from Bio import SeqIO
from .process_cluster import general_cluster_data
from .process_domains import process_secmet
from .process_nrps    import process_NRP
from .process_pks     import process_Polyketide

@click.command()
@click.option("--gbk_folder", type=click.Path(), prompt=True, help="Mibig GBK data folder")
@click.option("--json_folder", type=click.Path(), prompt=True, help="Mibig JSON data folder")
@click.option("--clusterdata", type=click.STRING, default="clusters.csv", help="Cluster Outputfile (csv)")
@click.option("--domaindata", type=click.STRING, default="domains.csv", help="Domainwh Outputfile (csv)")
@click.option("--nrpsdata", type=click.STRING, default="NRPS.csv", help="NRPS Outputfile (csv)")
@click.option("--pksdata", type=click.STRING, default="PKS.csv", help="PKS Outputfile (csv)")
@click.option("--fnafolder", type=click.STRING, default="fnaout", help="output directory for FNA files")
def process_mibig_cluster_folder(gbk_folder, json_folder, clusterdata, domaindata, nrpsdata, pksdata, fnafolder):
    """
    Assumes that mibig cluster data (and only mibig cluster data) is in the
    given folder. Mibig data consists of two files:

     - json: e.g. BGC00001.json
     - gbk:  e.g. BGC00001.1.final.gbk

    This script will process the entire folder an output the followign files:
        clusterdata: clusterinfo
        domaindata: domain info (from the GBK file)
        nrpsdata: NRPS-specific Info
        pksdata: PKS-specific Info
    """
    jsonfiles = glob.glob("{}/*.json".format(json_folder))
    gbkfiles = glob.glob("{}/*.gbk".format(gbk_folder))

    jsons = set([os.path.basename(jf).replace(".json", "") for jf in jsonfiles])
    gbks  = set([os.path.basename(gb).replace(".gbk", "") for gb in gbkfiles])

    targets = jsons.intersection(gbks)

    datafiles =  [ ("{}/{}.json".format(json_folder, target), "{}/{}.gbk".format(gbk_folder, target)) for target in targets]


    def makeclusters():
        ## Process Clusters
        ####################################################################################
        ####################################################################################

        # results = []
        # for jsonfile in jsonfiles:
        #     gbkfile = jsonfile[:-4] + "1.final.gbk"
        #     results.append(general_cluster_data(jsonfile,gbkfile))
        #
        #
        # df = pd.DataFrame([f for r in results for f in r])
        # df.to_csv("Clusters.csv",index=False)
        #general_cluster_data('mibig.secondarymetabolites.org/repository/BGC0000001/BGC0000001.json',
        #

        clusterinfo = [general_cluster_data(data) for data in  datafiles]
        clusterdf   = pd.DataFrame([f for r in clusterinfo for f in r]) #general_cluster_data returns a list
        clusterdf.to_csv(clusterdata, index=False)


    # def makedomains():
    #     ## Process Domains
    #     ####################################################################################
    #     ####################################################################################
    #     #there can be more than one gbk per cluster. This will pull out ALL the cluster info
    #     # should we use only the first cluster?
    #
    #     allgbks = [gbk for (json,gbk) in datafiles]
    #     domaininfo = []
    #     for genbank in allgbks:
    #         #print(genbank)
    #         gbkrecords = open(genbank,'r')
    #         for record in SeqIO.parse(gbkrecords,'genbank'):
    #             #print(process_secmet(record))
    #             domaininfo.append(process_secmet(record))
    #
    #     domains = pd.concat(domaininfo)
    #     #add a UniqueID columns
    #     #print(domains.head())
    #     domains['UniqueID'] =  domains.Protein_ID + "." + domains.Nucleotide_Start.astype(str)
    #     #eliminate dubplicates based on UniqueID (groupby and take first)
    #     domains = domains.groupby('UniqueID').head(1)
    #     domains.to_csv(domaindata,index=False)

    ## Process NRPs
    ####################################################################################
    ####################################################################################
    def makeNRPS():
        NRPs = [process_NRP(js) for (js, gbk) in datafiles]
        NRP_dfs = [pd.DataFrame(nrp) for nrp in NRPs if nrp]
        print("There are {} jsonfiles but only {} with NRP sections".format(len(jsonfiles), len(NRP_dfs)))

        nrpsdf = pd.concat(NRP_dfs)
        nrpsdf.to_csv(nrpsdata,index=False)
        # print df.shape

    ## Process PKs
    ####################################################################################
    ####################################################################################
    def makePKS():
        PKSs = [process_Polyketide(js) for (js, gbk) in datafiles]
        PKS_dfs = [pd.DataFrame(pks) for pks in PKSs if pks]
        print("There are {} jsonfiles but only {} with PKS sections".format(len(jsonfiles), len(PKS_dfs)))
        pksdf = pd.concat(PKS_dfs)
        pksdf.to_csv(pksdata,index=False)
        # print df.shape

    # ## Write FNAs
    # ####################################################################################
    # ####################################################################################
    # def write_fnafile(domain_df, outfile):
    #     with open(outfile,'w') as f:
    #         for idx, row in domain_df.iterrows():
    #             f.write(">{}\n{}\n".format(row.UniqueID,row.DNA_Sequence))
    #
    # def writeFNAs():
    #     if not os.path.exists(fnafolder): os.mkdir(fnafolder)
    #     df = pd.read_csv(domaindata)
    #     domains = list(pd.unique(df.Domain))
    #     for domain in domains:
    #         domain_df = df[df.Domain == domain]
    #         write_fnafile(domain_df, outfile="{}/{}.fna".format(fnafolder,domain))
    #

    # processing scripts as functions to use local scope/avoid memory issues
    makeclusters()
    #makedomains()
    makeNRPS()
    makePKS()
    #writeFNAs()
