import json
from tqdm import tqdm
import sys
import os

if len(sys.argv) < 5:
    exit("run: python analyse.py <species_1_og2genes> <species_2_og2genes> <idx_pairwise> <idx_namesMap>")

species_1_file_name = sys.argv[1] # "data/species_1__odb10v0_OG2genes.tab"
species_2_file_name = sys.argv[2] #"data/species_2__odb10v0_OG2genes.tab"

pairwise_file = sys.argv[3] #"/home/mabuelanin/Desktop/kprocessor/refseq_orthodb/skipmers_effect/idx_KMERS_SKIPMERS_species_1_species_2_.tsv"
names_map_file = sys.argv[4] #"./KMERS_SKIPMERS/species_1_species_2_/idx_KMERS_SKIPMERS_species_1_species_2_.namesMap"


#######################################################################
#                   Parsing kProcessor NamesMap File                  #
#######################################################################

print("Parsing namesMap")
names_map = {}
with open(names_map_file) as namesMap:
    for line in namesMap:
        line = line.strip().split("\t")
        gene_name = line[0]
        gene_id = line[1]
        names_map[gene_id] = gene_name

#######################################################################
#                   Parsing kCluster pairwise File                    #
#######################################################################

print("Parsing pairwise")
pairwise = {}
with open(pairwise_file) as pw:
    next(pw)  # Skip headers
    for line in pw:
        line = line.strip().split("\t")
        gene1_id = line[0]
        gene1 = names_map[gene1_id].split("|")[0]
        gene2_id = line[1]
        gene2 = names_map[gene2_id].split("|")[0]
        similarity = float(line[-1])
        
        if similarity not in pairwise:
            pairwise[similarity] = [[gene1, gene2]]
        else:
            pairwise[similarity].append([gene1, gene2])



#######################################################################
#                   Parsing species Ortho Group to Genes              #
#                   ex: 0at100953	9606_0:000408                     #
#                   gene2OGs[9606_0:000408] = [0at100953]             #
#######################################################################

gene2OGs = {}  # dict{key: odb_geneID, value: [Og_groups]}

for species_file in [species_1_file_name, species_2_file_name]:
    print("Parsing %s" % os.path.basename(species_file))
    with open(species_file, 'r') as species:
        for line in species:
            line = line.strip().split("\t")
            odb_gene_group = line[0]
            odb_gene_id = line[1]
            if odb_gene_id in gene2OGs:
                gene2OGs[odb_gene_id].append(odb_gene_group)
            else:
                gene2OGs[odb_gene_id] = [odb_gene_group]


relation_type = ""
total_count = {}

#######################################################################
#                           Preparing Levels                          #
#######################################################################

levels1 = ['2759', '33208', '7742', '32523', '40674', '9347', '314146', '9443', '314295', '9604', '9606']
levels2 = ['2759', '33208', '7742', '32523', '40674', '9347', '314146', '9443', '314294', '9544']

levels = []
for i in range(min(len(levels1), len(levels2))):
    if levels1[i] == levels2[i]:
        levels.append(levels1[i])
    else:
        break

levels_col = "\t".join(["L-" + str(i) for i in range(len(levels))])

new_pairwise_file = open("sk21_new_pairwise.tsv", 'w')
new_pairwise_header = "\t".join(["sim%", "gene1", "gene2", "Ortho", levels_col]) + "\n"
new_pairwise_file.write(new_pairwise_header)


def get_tax(ogs):
    tax = []
    for og in ogs:
        tax.append(og.split("at")[1])
    
    return tax


for sim, records in tqdm(pairwise.items()):
    for record in records:
        gene1 = record[0]
        gene2 = record[1]
        tax1 = gene1.split("_")[0]
        tax2 = gene2.split("_")[0]

        if tax1 == tax2:
            relation_type = "para"
        else:
            relation_type = "ortho"

        if relation_type == "para":
            continue

        gene1_groups = set(gene2OGs[gene1])
        gene2_groups = set(gene2OGs[gene2])
        intersection = get_tax(list(gene1_groups.intersection(gene2_groups)))
        
        
        matched = len(intersection) > 0

        levels_alignment = [0] * len(levels)
        for i in intersection:
            if i in levels:
                levels_alignment[levels.index(i)] = 1
        
        levels_alignment = "\t".join(list(map(str, levels_alignment)))

        # print("Gene1 : %s" % gene1)
        # print("Gene2 : %s" % gene2)
        # print("Matched: %s" % str(matched))
        # print("relation_type : %s" % relation_type)
        # print("Gene1 Groups: %s" % (",".join(gene1_groups)))
        # print("Gene2 Groups: %s" % (",".join(gene2_groups)))
        # print("Intersection: %s" % (",".join(list(intersection))))

        new_pairwise_line = "\t".join([str(sim), gene1, gene2, str(matched), levels_alignment]) + "\n"
        new_pairwise_file.write(new_pairwise_line)
