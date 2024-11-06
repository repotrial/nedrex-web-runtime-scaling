import os.path
from collections import defaultdict

from random import randrange

disorders = {}
disorder_names ={}
disease_gene_edges=defaultdict(set)

num_per_step = 5
steps = [10,25,50,100,250,500,1000]
allowed_divergence_per_step_as_fraction=0.15
all_steps_matching=defaultdict(set)

def save_disorders(node):
    disorders[node["primaryDomainId"]]=node
    disorder_names[node["primaryDomainId"]]=normalize_disorder_name(node["displayName"])

def normalize_disorder_name(name):
    name = name.lower()
    return name.replace(" ", "_").replace("/","&").replace(",",".")

def save_gene_disease_edges(edge):
    gene = edge["sourceDomainId"]
    disorder = edge["targetDomainId"]
    disease_gene_edges[disorder].add(gene)

import requests

# nedrex_api = 'https://api.web.nedrex.net/api/'
nedrex_api = 'http://localhost:8090/api/'

db_version = requests.get(nedrex_api + "getMetadata").json()["repotrial"]["version"]

userId = "03649658-2e71-46db-8cb5-8517a616ebdb"
# userId = "c35a4be1-6a93-4074-ad74-b82c3bd020e2"
disease_gene_network_build_payload= {
    "nodes": {
        "gene": {
            "filters": []
        },
        "disorder": {
            "filters": []
        }
    },
    "edges": {
        "GeneAssociatedWithDisorder": {
            "filters": []
        }
    },
    "connectedOnly": True,
    "interactions": {
        "ProteinInteractsWithProtein": False,
        "GeneInteractsWithGene": False
    },
    "options": {
        "nodes": {
            "codingGenesOnly": False,
            "approvedDrugsOnly": False,
            "filterElementDrugs": True
        },
        "edges": {
            "experimentalInteraction": True,
            "disorderParents": False,
            "drugTargetsWithAction": False,
            "disorderAssociationCutoff": 0,
            "extendGGI": False,
            "extendPPI": False
        }
    },
    "uid": userId
}


graph_id = requests.post(nedrex_api + "getGraphInfo", json=disease_gene_network_build_payload).json()["id"]


if not os.path.exists("./graph"):
    os.mkdir("./graph")
os.system(f"curl -o ./graph/graph.graphml {nedrex_api}downloadGraph?gid={graph_id}")

import networkx as nx
from collections import defaultdict

def translate_genes(entrez_id):
    payload = {
        "name": "gene",
        "query": entrez_id.split(".")[1] if "." in entrez_id else entrez_id,
        "typeCount": "gene"
    }
    sid = requests.post(nedrex_api + "getSuggestions", json=payload).json()["suggestions"][0]["sid"]
    payload = {
    "sourceType": "gene",
    "targetType": "gene",
    "sugId": sid,
    "noloop":True
}
    return int(requests.post(nedrex_api + "getConnectedNodes", json=payload).json()[0]["id"])


def analyze_graphml(file_path):
    # Read the GraphML file
    G = nx.read_graphml(file_path)

    # Find all nodes that start with "mondo."
    mondo_nodes = [node for node in G.nodes() if node.startswith('mondo.')]

    # Create a dictionary to store mondo nodes and their neighbors
    mondo_neighbors = defaultdict(set)

    # For each mondo node, get its immediate neighbors
    for mondo_node in mondo_nodes:
        print(G.nodes[mondo_node].get('displayName', ''))
        disorders[mondo_node] = G.nodes[mondo_node]
        disorder_names[mondo_node] = normalize_disorder_name(G.nodes[mondo_node].get('displayName', ''))
        # Get all neighbors
        neighbors = list(G.neighbors(mondo_node))
        disease_gene_edges[mondo_node].update(neighbors)

analyze_graphml("./graph/graph.graphml")

for disease,genes in disease_gene_edges.items():
    num_genes = len(genes)
    for step in steps:
        num_genes_frac = num_genes/step
        if 1-allowed_divergence_per_step_as_fraction <= num_genes_frac <=1+allowed_divergence_per_step_as_fraction:
            all_steps_matching[step].add(disease)
            break

for step, diseases in all_steps_matching.items():
    print(f"Step: {step} [{len(diseases)}]")
    for disease in diseases:
        print(f"\t{disease}: {disorder_names[disease]} [{len(disease_gene_edges[disease])}]")

index_file = os.path.join("./seed_files","index.txt")
index_fw = open(index_file,"w")

def write_seed_file(file_name, seeds):
    index_fw.write(file_name+"\n")
    with open(file_name, "w") as f:
        f.write("#EntrezID\n")
        seeds = {translate_genes(entrez) for entrez in seeds}
        for seed in seeds:
            f.write(f"{seed}\n")

for step in steps:
    print(f"Step {step}")
    wd = os.path.join("./seed_files",str(step))
    os.makedirs(wd, exist_ok=True)
    diseases = list(all_steps_matching[step])
    print(diseases)
    target = num_per_step
    while target > 0 and len(diseases) > 0:
        target -=1
        pos = randrange(0, len(diseases))
        disease = diseases[pos]
        diseases.remove(disease)
        seeds = disease_gene_edges[disease]
        write_seed_file(os.path.join(wd,f"{disorder_names[disease]}-{len(seeds)}seeds.tsv"),seeds)

index_fw.close()