
import os
import glob
import pandas as pd
from process_cluster import general_cluster_data
from process_domains import process_secmet
from process_nrps import process_NRP
from process_pks import process_Polyketide

def process_mibig_cluster_folder(mibigfolder,
                                 clusterdata = "Clusters.csv",
                                 domaindata  = "Domains.csv",
                                 nrpsdata ="NRPs.csv",
                                 pksdata ="PKSs.csv"):
    """
    Assumes that mibig cluster data (and only mibig cluster data) is in the
    given folder. Mibig data consists of two files:

     - json: e.g. BGC00001.json
     - gbk:  e.g. BGC00001.1.final.gbk




    :return:
    """
    jsonfiles = glob.glob("{}/*.json".format(mibigfolder))
    gbkfiles = [jsonfile[:-4] + "1.final.gbk" for jsonfile in jsonfiles]

    for gbkfile in gbkfiles:
        assert os.path.exists(gbkfile)

    datafiles = zip(jsonfiles,gbkfiles)

    def makeclusters():
        ## Process Clusters
        ####################################################################################
        ####################################################################################

        clusterinfo =  [general_cluster_data(data) for data in  datafiles]
        clusterdf = pd.DataFrame([f for r in clusterinfo for f in r]) #general_cluster_data returns a list
        clusterdf.to_csv(clusterdata, index=False)

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

    def makedomains():
        ## Process Domains
        ####################################################################################
        ####################################################################################
        #there can be more than one gbk per cluster. This will pull out ALL the cluster info
        # should we use only the first cluster?
        allgbks = glob.glob("{}/*.final.gbk".format(mibigfolder))
        domaininfo = []
        for genbank in allgbks:
            gbkrecords = open(genbank,'r')
            for record in SeqIO.parse(gbkrecords,'genbank'):
                domaininfo.append( process_secmet(record))

        domains = pd.concat(domaininfo)
        domains.to_csv(domaindata,index=False)


    ## Process NRPs
    ####################################################################################
    ####################################################################################
    def makeNRPS():
        NRPs = [process_NRP(t) for t in jsonfiles]
        NRP_dfs = [pd.DataFrame(nrp) for nrp in NRPs if nrp]
        print "There are {} jsonfiles but only {} with NRP sections".format(len(jsonfiles), len(NRP_dfs))

        nrpsdf = pd.concat(NRP_dfs)
        nrpsdf.to_csv(nrpsdata,index=False)
        # print df.shape

    ## Process PKs
    ####################################################################################
    ####################################################################################
    def makePKS():
        PKSs = [process_Polyketide(t) for t in jsonfiles]
        PKS_dfs = [pd.DataFrame(pks) for pks in PKSs if pks]
        print "There are {} jsonfiles but only {} with PKS sections".format(len(jsonfiles), len(PKS_dfs))
        pksdf = pd.concat(PKS_dfs)
        pksdf.to_csv(pksdata,index=False)
        # print df.shape

    # processing scripts as funcitons to use local scope
    makeclusters()
    makedomains()
    makeNRPS()
    makePKS()
