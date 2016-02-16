import os
import json
import pprint

from toolz import assoc, thread_first


#Process NRPs
pp = pprint.PrettyPrinter(indent=2)

def condv(key, dic, default="None"):
    "Helper Function to handle missing keys"
    try:
        val = dic[key]
        if type(val) == list:
            if not val: #emptylist
                return ""
            if len(val) == 1:
                return val[0]
            else:
                return ";".join(val)
        else:
            return dic[key]
    except:
        return default

def process_NRP(jsonfile):
    """process the NRP field of a MIBIG JSON Document

    iterate to the innermost level (nrps->nrps_genes->nrps_module->a_subst_spec)
    and genetate a flat disctionary for each module/adomain.

    Use the condv funcitonto handle missing key-value pairs.

    Return a list of these flat files.

    {u'nrps_release_type': u'Hydrolysis',
     u'nrps_thioesterase': [],
     u'lin_cycl_nrp': u'Cyclic',
     u'nrps_te_type': u'Unknown',
     u'subclass': u'Other
     u'nrps_genes': [{u'nrps_gene': u'ACZ55945.1'
                     u'nrps_module': [
                        {u'cdom_subtype': u'N/A',
                        u'module_nr': u'0',
                        u'nrps_mod_doms': u'None',
                        u'a_substr_spec': {u'evidence_a_spec': u'None',
                                           u'epimerized': False,
                                           u'nonprot_adom_spec': u'DVEDIGSVAK',
                                           u'aa_type': u'Nonproteinogenic'},
                        u'nrps_mod_skip_iter': u'Neither'},
                        .....]
                        }]
    }
    """
    j =json.load(open(jsonfile,'r'))

    defaultvalues = {"aa_subcluster": "NA",
                        "aa_type": "Unknown",
                        "cdom_subtype": "Unknown",
                        "epimerized": "NA",
                        "evidence_a_spec": "NA",
                        "lin_cycl_nrp": "Linear",
                        "module_nr": "NA",
                        "nonprot_adom_spec": "NA",
                        "nrps_gene": "NA",
                        "nrps_mod_doms": "NA",
                        "nrps_mod_skip_iter": "NA",
                        "nrps_release_type": "Unknown",
                        "nrps_te_type": "Unknown",
                        "nrps_thioesterase": "NA",
                        "other_spec": "NA",
                        "prot_adom_spec": "NA",
                        "subclass": "Other"}



    #check for NRPS
    try:
        nrps = j['general_params']['NRP']
    except:
        #don't do anything with non-NRP JSON
        #print "No NRP in {}".format(jsonfile)
        return

    #check for NRPS genes:
    try:
        genes = nrps['nrps_genes']
    except:
        #in this case there is only simple data i.e. cyclic vs. linear
        print("jsonfile {} has only a simple record".format(jsonfile))
        data = defaultvalues
        for key, value in nrps.items():
            data = assoc(data, key, value)
        return [data]

    #gene information but no module information:
    try:
        #check that genes have modules
        modules = [gene['nrps_module'] for gene in genes]
    except:
        print("handle the case where there are genes but no modules")
        data = defaultvalues
        for key, value in nrps.items():
            if key != "nrps_genes":
                data = assoc(data, key, value)

        alldata = []
        for gene in genes:
            tempdata = data
            for k,v in gene.items():
                tempdata = assoc(tempdata,k,v)
            alldata.append(tempdata)
        return alldata

    #process all of the full records
    allvalues = []
    for gene in genes:
        for module in gene['nrps_module']:
            try:
                assert('a_substr_spec' in module.keys())
                adomain = module['a_substr_spec']
                data = thread_first(defaultvalues,
                    #adomain level info
                    (assoc, 'prot_adom_spec', condv('prot_adom_spec', adomain)),
                    (assoc, 'nonprot_adom_spec', condv('nonprot_adom_spec', adomain)),
                    (assoc, 'evidence_a_spec', condv('evidence_a_spec, adomain',adomain)),
                    (assoc, 'aa_type', condv('aa_type', adomain,  default="Unknown")),
                    (assoc, 'epimerized', condv('epimerized', adomain)),
                    (assoc, 'other_spec', condv('other_spec', adomain)),
                    (assoc, 'a_multiple_spec', condv('a_multiple_spec', adomain)),
                    (assoc, 'aa_subcluster', condv('aa_subcluster', adomain)),

                    #modulelevel info
                    (assoc, "cdom_subtype", condv('cdom_subtype',module, default="Unknown")),
                    (assoc, "module_nr", condv('module_nr',module)),
                    (assoc, "nrps_mod_doms", condv('nrps_mod_doms', module)),
                    (assoc, "nrps_other_mod_dom", condv('nrps_other_mod_dom', module)),
                    (assoc, "nrps_mod_skip_iter", condv('nrps_mod_skip_iter', module)),
                    (assoc, "nrps_evidence_skip_iter", condv( "nrps_evidence_skip_iter", module)),

                    #genelevel info
                    (assoc, "nrps_gene", condv('nrps_gene',gene)),
                    #clusterlevel info
                    (assoc, "lin_cycl_nrp", condv('lin_cycl_nrp',nrps, default="Unknown")),
                    (assoc, "nrps_release_type", condv('nrps_release_type', nrps, default="Unknown")),
                    (assoc, "nrps_te_type", condv('nrps_te_type', nrps, default="Unknown")),
                    (assoc, "nrps_thioesterase", condv('nrps_thioesterase', nrps)),
                    (assoc, "subclass", condv('subclass', nrps, default="Other")),
                    (assoc, "mibigaccession", j['general_params']['mibig_accession']))


                allvalues.append(data)
                print("Data: \n\n\n\n")
                print(data)
                print(adomain)
                print(module)
                print(nrps)
                print("end \n\n\n")
            except:
                print("error in jsonfile:{}\n \
                       gene:{}\n \
                       module:{}\n \
                       alldata:{}\n".format(jsonfile, gene['nrps_gene'], module,nrps))
    return allvalues
