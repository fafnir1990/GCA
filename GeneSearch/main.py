from Bio.KEGG import REST
from bioservices import Reactome
import csv
from input import inp

#gene_list = ['POLD1', 'POLE3', 'ABO', 'TP53']
gene_list = inp
specie = "hsa"


human_pathways = REST.kegg_list("pathway", specie).read()
human_pathways_dict = {}
repair_pathways = []
repair_pathways_dict = {}
for line in human_pathways.rstrip().split("\n"):
    entry, description = line.split("\t")
    human_pathways_dict[entry] = description
    if "DNA" in description:
        repair_pathways.append(entry)
        repair_pathways_dict[entry] = description


rejected = []
gene_dict = dict((gene,[]) for gene in gene_list)

i = 0
len_ = len(human_pathways_dict.keys())
for pathway in human_pathways_dict.keys():
    i += 1
    print str(i) + ' // ' + str(len_)
    #print pathway
    pathway_file = REST.kegg_get(pathway).read()
    current_section = None
    for line in pathway_file.rstrip().split("\n"):
        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section
        if current_section == "GENE":
            try:
                gene_identifiers, gene_description = line[12:].split("; ")
                gene_id, gene_symbol = gene_identifiers.split()
                if gene_symbol in gene_dict.keys():
                    gene_dict[gene_symbol].append(pathway)
            except:
                rejected.append(pathway)

print "!!     KEGG     !!"
with open('KEGG_stats.csv', 'wb') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow(['GENE', 'HITS'])
    for gene in gene_dict.keys():
        nop = len(gene_dict[gene])
        print 'Gene code: ' + gene
        print 'Number of pathways: ' + str(nop)
        spamwriter.writerow([str(gene), str(nop)])
        for pathway in gene_dict[gene]:
            print human_pathways_dict[pathway]
        print '\n'


print "!!     REACTOME     !!"
react = Reactome()
#"Currently only human pathways will be returned from this method"

with open('REACTOME_stats.csv', 'wb') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow(['GENE', 'HITS'])
    for gene in gene_list:
        pathway_list = react.query_hit_pathways(gene)
        nop = len(pathway_list)
        spamwriter.writerow([str(gene), str(nop)])
        for p in pathway_list:
            print gene +': '+ p["displayName"] +' , '+ p["stableIdentifier"]["displayName"]








dicta= {'POLD1': ['path:hsa00240', 'path:hsa03410', 'path:hsa03420', 'path:hsa03430', 'path:hsa03030', 'path:hsa00230', 'path:hsa03440', 'path:hsa05166'], 'POLE3': ['path:hsa00240', 'path:hsa03410', 'path:hsa03420', 'path:hsa03030', 'path:hsa00230', 'path:hsa05166'], 'ABO': ['path:hsa00601']}


