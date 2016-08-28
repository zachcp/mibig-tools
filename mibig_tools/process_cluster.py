import json
from Bio import SeqIO, SeqUtils
from toolz import merge

def general_cluster_data(data):
    """
    return a list of dictionaries with cluster information
    """
    jsonfile, gbkfile = data
    clusterdata = json.load(open(jsonfile))
    gbks = SeqIO.parse(open(gbkfile),'genbank')

    compounds = get_compounds(clusterdata)

    results = []
    for gbk in gbks:
        for compound in compounds:
            clusterdict = {
                "biosyn_class": ';'.join(clusterdata['general_params']['biosyn_class']),
                 "mibig_accession": gbk.annotations['accessions'][0],
                 "version": gbk.annotations['sequence_version'],
                 "organism": gbk.annotations['source'],
                 "organism_taxa": ';'.join(gbk.annotations['taxonomy']),
                 "cluster_length":  len(gbk.seq),
                 "cluster_GC": SeqUtils.GC(gbk.seq),
                 "accession": gbk.annotations['comment'].split()[-1],
                 "molecule": compound
            }

        taxonomydata = process_taxonomy(';'.join(gbk.annotations['taxonomy']))

        resultdata = merge(clusterdict, taxonomydata)
        results.append(resultdata)
    return results


def get_compounds(BGCjson):
    "pull out compound information form teh Mibig JSON"
    try:
        compounddicts = BGCjson['general_params']['compounds']
        compounds     = [cd["compound"] for cd in compounddicts]
        return compounds
    except:
        return list()


def process_taxonomy(tax):
    """
    Needs to Handle the Taxonomy from NCBI

    This parser was handmade for MIBIG and woks by verifying the taxonomy of each MiBiG taxonomy field.
    """

    try:
        taxonomy = tax.split(';')
    except:
        taxonomy = "Unknown" * 7

    # Handle By Taxonomy Length
    ######################################################################
    ######################################################################
    if len(taxonomy) == 2:
        tax1 =  ['Bacteria', 'environmental samples']
        tax2 =  ['Bacteria', 'Unknown']
        if taxonomy == tax1 or tax2:
            # extended taxonomy
            k,p = taxonomy
            c,o,f,g = [p for i in range(4)]

        else:
            print(len(taxonomy))
            print(taxonomy)

    elif len(taxonomy) == 4:
        tax1 = ['Bacteria', 'Cyanobacteria', 'Oscillatoriales', 'Trichodesmium']
        tax2 = ['Bacteria', 'Proteobacteria', 'Deltaproteobacteria', 'Candidatus Entotheonella']
        tax3 = ['Bacteria', 'Cyanobacteria', 'Stigonematales', 'Hapalosiphon']
        if taxonomy == tax1 or tax3:
            # extended taxonomy
            k,p,o,g = taxonomy
            c,f = ["Unknown","Unknown"]
        elif taxonomy == tax2:
            k,p,c,g = taxonomy
            o,f= ["Unknown","Unknown"]
        else:
            print(len(taxonomy))
            print(taxonomy)

    elif len(taxonomy) == 5:
        tax1 = ['Bacteria', 'Cyanobacteria', 'Oscillatoriophycideae', 'Oscillatoriales', 'Kamptonema']
        tax2 = ['Bacteria', 'Cyanobacteria', 'Oscillatoriophycideae', 'Oscillatoriales', 'Moorea']
        if taxonomy == tax1 or tax2:
            # extended taxonomy
            k,p,c,o,f = taxonomy
            g ="Unknown"

        else:
            print(len(taxonomy))
            print(taxonomy)

    elif len(taxonomy) == 6:
        tax1 = ['Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Pseudomonadales', 'Pseudomonadaceae', 'Pseudomonas']
        tax2 = ['Bacteria', 'Firmicutes', 'Bacilli', 'Bacillales', 'Paenibacillaceae', 'Brevibacillus']
        if taxonomy == tax1 or tax2:
            # extended taxonomy
            k,p,c,o,f,g = taxonomy

        else:
            print(len(taxonomy))
            print(taxonomy)

    elif len(taxonomy) == 7:
        #classic taxonomy
        tax1 = ['Bacteria', 'Actinobacteria', 'Actinobacteridae', 'Actinomycetales', 'Micromonosporineae', 'Micromonosporaceae', 'Verrucosispora']
        tax2 = ['Bacteria', 'Actinobacteria', 'Actinobacteridae', 'Actinomycetales', 'Pseudonocardineae', 'Pseudonocardiaceae', 'Kutzneria']
        if taxonomy == tax1 or tax2:
            k,p,subclass,o,suborder,f,g =  taxonomy
            c="Actinobacteria"

        else:
            print(len(taxonomy))
            print(taxonomy)

    elif len(taxonomy) == 8:
        tax1 = ['Bacteria', 'Actinobacteria', 'Actinobacteridae', 'Actinomycetales', 'Streptomycineae', 'Streptomycetaceae', 'Streptomyces', 'Streptomyces albidoflavus group']
        tax2 = ['Bacteria', 'Actinobacteria', 'Actinobacteridae', 'Actinomycetales', 'Corynebacterineae', 'Mycobacteriaceae', 'Mycobacterium', 'Mycobacterium avium complex (MAC)']
        if taxonomy == tax1 or tax2:
            # extended taxonomy
            k,p,c,o,suborder,f,g,subgenus = taxonomy

        else:
            print(len(taxonomy))
            print(taxonomy)

    elif len(taxonomy) == 9:
        tax1 = ['Eukaryota', 'Fungi', 'Dikarya', 'Ascomycota', 'Pezizomycotina', 'Leotiomycetes', 'Helotiales', 'Sclerotiniaceae', 'Botrytis']
        tax2 = ['Eukaryota', 'Fungi', 'Dikarya', 'Ascomycota', 'Pezizomycotina', 'Leotiomycetes', 'Helotiales', 'Helotiaceae', 'Glarea']
        if taxonomy == tax1 or tax2:
            # extended taxonomy
            superkingdom,k,subkindgdom,p,subphylum,c,o,f,g = taxonomy

        else:
            print(len(taxonomy))
            print(taxonomy)

    elif len(taxonomy) == 10:
        tax1 = ['Eukaryota', 'Fungi', 'Dikarya', 'Ascomycota', 'Pezizomycotina', 'Eurotiomycetes', 'Eurotiomycetidae', 'Eurotiales', 'Aspergillaceae', 'Aspergillus']
        tax2 = ['Eukaryota', 'Fungi', 'Dikarya', 'Ascomycota', 'Pezizomycotina', 'Eurotiomycetes', 'Eurotiomycetidae', 'Eurotiales', 'Aspergillaceae', 'Monascus']
        if taxonomy == tax1 or tax2:
            # extended taxonomy
            superkingdom,k,subkingdom,p,subphlum,c,subclass,o,f,g = taxonomy
        else:
            print(len(taxonomy))
            print(taxonomy)

    elif len(taxonomy) == 11:
        tax1 = ['Eukaryota', 'Fungi', 'Dikarya', 'Ascomycota', 'Pezizomycotina', 'Dothideomycetes', 'Pleosporomycetidae', 'Pleosporales', 'Pleosporineae', 'Pleosporaceae', 'Alternaria']
        tax2 = ['Eukaryota', 'Fungi', 'Dikarya', 'Ascomycota', 'Pezizomycotina', 'Sordariomycetes', 'Hypocreomycetidae', 'Hypocreales', 'Nectriaceae', 'Fusarium', 'Fusarium fujikuroi species complex']
        if taxonomy == tax1 or tax2:
            # extended taxonomy doth=??
            superkingdom,k,subkingdom,p,subphylum,c,subclass,o,suborder,f,g = taxonomy
        else:
            print(len(taxonomy))
            print(taxonomy)

    elif len(taxonomy) == 12:
        tax1 = ['Eukaryota', 'Fungi', 'Dikarya', 'Ascomycota', 'Pezizomycotina', 'Dothideomycetes', 'Pleosporomycetidae', 'Pleosporales', 'Pleosporineae', 'Pleosporaceae', 'Alternaria', 'Alternaria alternata group']
        tax2 = ['Eukaryota', 'Fungi', 'Dikarya', 'Ascomycota', 'Pezizomycotina', 'Dothideomycetes', 'Pleosporomycetidae', 'Pleosporales', 'Pleosporineae', 'Didymellaceae', 'mitosporic Didymellaceae', 'Phoma']
        tax3 = ['Eukaryota', 'Fungi', 'Dikarya', 'Ascomycota', 'Pezizomycotina', 'Dothideomycetes', 'Pleosporomycetidae', 'Pleosporales', 'Pleosporineae', 'Leptosphaeriaceae', 'Leptosphaeria', 'Leptosphaeria maculans complex']
        if taxonomy == tax1 or tax3:
            superkingdom,k,subkingdom,p,subphylum,c,subclass,o,suborder,f,g,subgenus = taxonomy
        elif taxonomy == tax2:
            superkingdom,k,subkingdom,p,subphylum,c,subclass,o,suborder,f,_,g = taxonomy

        else:
            print(len(taxonomy))
            print(taxonomy)

    elif len(taxonomy) == 14:
        tax1 = ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Poales', 'Poaceae', 'BEP clade', 'Ehrhartoideae', 'Oryzeae', 'Oryza']
        tax2 = ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Poales', 'Poaceae', 'PACMAD clade', 'Panicoideae', 'Andropogoneae', 'Sorghum']
        if taxonomy == tax1 or tax2:
            # extended taxonomy doth=??
            superkingdom,k,p,__,___,____,____,c,o,f,x,subfamily,tribe,g = taxonomy
        else:
            print(len(taxonomy))
            print(taxonomy)

    elif len(taxonomy) == 16:
        tax1 = ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'eudicotyledons', 'Gunneridae', 'Pentapetalae', 'rosids', 'malvids', 'Brassicales', 'Brassicaceae', 'Camelineae', 'Arabidopsis']
        if taxonomy == tax1:
            # extended taxonomy doth=??
            superkingdom,k,p,__,___,____,____,x,xx,xxx,subclass,xxxx,o,f,tribe,g = taxonomy
            c="Unknown"

        else:
            print(len(taxonomy))
            print(taxonomy)

    else:
        print(len(taxonomy))
        print(taxonomy)
        return {"kingdom": "Unknown",
                "phylum": "Unknown",
                "class": "Unknown",
                "order": "Unknown",
                "family": "Unknown",
                "genus": "Unknown"}

    # Return
    ######################################################################
    ######################################################################

    return {"kingdom": k,
            "phylum": p,
            "class": c,
            "order": o,
            "family": f,
            "genus": g}
