import os
import json
import pprint

from toolz import assoc, thread_first


#Process KS
pp = pprint.PrettyPrinter(indent=2)

def condv(key, dic, default="None"):
    "Helper Funciton to handle missing keys"
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
    
def process_Polyketide(jsonfile):
    """process the Polyketide field of a MIBIG JSON Document
    
    iterate to the innermost level (nrps->nrps_genes->nrps_module->a_subst_spec)
    and genetate a flat disctionary for each module/adomain. 
    
    Use the condv funcitonto handle missing key-value pairs.
    
    Return a list of these flat files.
    
    { u'lin_cycl_pk': u'Cyclic',
      u'mod_pks_genes': [ { u'mod_pks_gene': u'AEK75502.1',
                             u'pks_module': [ { u'at_substr_spec': u'Unknown',
                                                u'module_nr': u'4',
                                                u'pks_domains': [ u'Ketosynthase',
                                                                  u'Acyltransferase',
                                                                  u'Ketoreductase',
                                                                  u'Dehydratase',
                                                                  u'Thiolation (ACP/PCP)'],
                                                u'pks_mod_doms': u'None',
                                                u'pks_mod_skip_iter': u'Neither'},
                                            ....]
                            },
                            ......
                            ],
      u'pk_subclass': u'Other',
      u'pks_release_type': u'Macrolactonization',
      u'pks_subclass': [u'Modular type I'],
      u'pks_te_type': u'None',
      u'starter_unit': u'Acetyl-CoA'}
    """
    j =json.load(open(jsonfile,'r'))
    
    defaultvalues = {k:"NA" for k in  
                       [u'lin_cycl_pk',
                        u'mod_pks_genes',
                        u'mod_pks_gene',
                        u'pks_module',
                        u'at_substr_spec',
                        u'at_multiple_spec',
                        u'evidence_at_spec',
                        u'kr_stereochem',
                        u'module_nr',
                        u'pk_subclass',
                        u'pks_domains',
                        u'pks_mod_doms',
                        u"pks_other_mod_dom"
                        u'pks_mod_skip_iter',
                        u'pks_evidence_skip_iter',
                        u'pks_release_type',
                        u'pks_subclass',
                        u'pks_te_type',
                        u'starter_unit',
                        u"nr_iterations",
                        u"iterative_subtype",
                        u"iter_cycl_type",
                        u"trans_at",
                        u"ketide_length",
                        u"cyclases",
                        u"pks_te_type",
                        u"pks_thioesterase"]}    
    
    #check for PKS
    try:
        polyketide = j['general_params']['Polyketide']
        #pp.pprint(polyketide)
    except:
        #don't do anythign with non-NRP JSON
        #print "No NRP in {}".format(jsonfile)
        return
    
    #check for PKS genes:
    try:
        genes = polyketide['mod_pks_genes']
    except:
        #in this case there is only simple data i.e. cyclic vs. linear
        print("jsonfile {} has only a simple record".format(jsonfile))
        data = defaultvalues
        for key, value in polyketide.items():
            data = assoc(data, key, condv(value, polyketide))
        data = thread_first(data,
                            (assoc, "mibigaccession", j['general_params']['mibig_accession']),
                            (assoc, "mibigaccession_file", os.path.basename(jsonfile)[:-5]))
        return [data]
        
    #gene information but no module information:
    try:
        #check that genes have modules
        modules = [gene['pks_module'] for gene in genes]
    except:
        print("handle the case where there are genes but no modules")
        data = defaultvalues
        for key, value in polyketide.items():
            if key != "mod_pks_genes":
                data = assoc(data, key, condv(value,polyketide))

        alldata = []
        for gene in genes:
            tempdata = data
            for k,v in gene.items():
                tempdata = assoc(tempdata,k,condv(v,gene))
            alldata.append(tempdata)    
        return alldata
            
    #process all of the full records
    allvalues = [] 
    for gene in genes:
        for module in gene['pks_module']:
            try:
                data = thread_first(defaultvalues,     
                    #modulelevel info
                    (assoc, 'at_substr_spec', condv('at_substr_spec',module)),
                    (assoc, 'at_multiple_spec', condv('at_multiple_spec',module)), 
                    (assoc, 'evidence_at_spec', condv('evidence_at_spec',module)), 
                    (assoc, 'kr_stereochem', condv('kr_stereochem',module)), 
                    (assoc, 'module_nr', condv('module_nr',module)),
                    (assoc, 'pks_mod_doms', condv('pks_mod_doms',module)),
                    (assoc, 'pks_mod_skip_iter', condv('pks_mod_skip_iter',module)), 
                    (assoc, 'pks_domains', condv('pks_domains',module)),
                    (assoc, 'pks_other_mod_dom', condv('pks_other_mod_dom',module)),
                    (assoc, 'pks_evidence_skip_iter', condv('pks_evidence_skip_iter',module)),
                    #genelevel info
                    (assoc, "mod_pks_gene", condv('nrps_gene',gene)),
                    #clusterlevel info
                    (assoc, "lin_cycl_pk", condv('lin_cycl_pk',polyketide)),
                    (assoc, 'pk_subclass', condv('pk_subclass', polyketide)),
                    (assoc, 'pks_release_type', condv('pks_release_type', polyketide)),
                    (assoc, 'pks_subclass', condv('pks_subclass', polyketide)),
                    (assoc, 'pks_te_type', condv('pks_te_type', polyketide)),
                    (assoc, 'starter_unit', condv('starter_unit', polyketide)),
                                    
                    (assoc, 'pks_thioesterase', condv('pks_thioesterase', polyketide)),
                    (assoc, 'cyclases', condv('cyclases', polyketide)),
                    (assoc, 'ketide_length', condv('ketide_length', polyketide)),
                    (assoc, 'starter_unit', condv('starter_unit', polyketide)),
                    (assoc, 'trans_at', condv('trans_at', polyketide)),
                    (assoc, 'iter_cycl_type', condv('iter_cycl_type', polyketide)),
                    (assoc, 'iterative_subtype', condv('iterative_subtype', polyketide)),
                    (assoc, 'nr_iterations', condv('nr_iterations', polyketide)),                                    
                    (assoc, "mibigaccession", j['general_params']['mibig_accession']),
                    (assoc, "mibigaccession_file", os.path.basename(jsonfile)[:-5]))  
                
                allvalues.append(data)
            except:
                print("error in jsonfile:{}\n \
                       gene:{}\n \
                       module:{}\n \
                       alldata:{}\n".format(jsonfile, gene['mod_pks_gene'], module,polyketide))
    return allvalues
