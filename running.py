import os

import time
import pandas as pd
import random
import requests

nedrex_api = 'https://api.web.nedrex.net/api/'
# nedrex_api = 'http://localhost:8090/api/'

db_version = requests.get(nedrex_api + "getMetadata").json()["repotrial"]["version"]

userId = "0ec00f73-3715-49e1-b271-079b73cf5a34"
# params = {"elements": True, "type": "gene", "addInteractions": True, "nodesOnly": False}
# tissue = None


# test_file = "./seed_files/10/transient_ischemic_attack-10seeds.tsv"


def read_seed_files(file):
    seeds = []
    with open(file) as f:
        for line in f:
            if line[0] == "#":
                continue
            seeds.append(int(line.strip()))
    return seeds


# payload["nodes"] = read_seed_files(test_file)

def run_request(payload):
    try:
        print(payload)
        start_time = time.time()
        job = requests.post(nedrex_api + "submitJob", json=payload).json()
        jid = job["jid"]
        print(job)

        while (job["state"] not in ["DONE", "ERROR"]):
            time.sleep(0.5)
            job = requests.get(nedrex_api + "getJob?id=" + jid).json()
            print(job)
        end_time = time.time()
        if job["state"] == "DONE":
            return end_time - start_time
        else:
            return 0
    except Exception as e:
        print(e)
        return -1


#
#
arguments = {
    "module_identification": {"diamond": {
        "experimentalOnly": True,
        "addInteractions": True,
        "nodesOnly": False,
        "n": 200,
        "alpha": 1,
        "pcutoff": -3,
        "type": "gene"
    }, "must": {
        "experimentalOnly": True,
        "addInteractions": True,
        "nodesOnly": False,
        "penalty": 0,
        "multiple": True,
        "trees": 10,
        "maxit": 10,
        "type": "gene"
    }, "kpm": {
        "experimentalOnly": True,
        "addInteractions": True,
        "nodesOnly": False,
        "k": 1,
        "type": "gene"
    }, "domino": {
        "experimentalOnly": True,
        "addInteractions": True,
        "nodesOnly": False,
        "type": "gene"
    }, "robust": {
        "experimentalOnly": True,
        "addInteractions": True,
        "nodesOnly": False,
        "trees": 30,
        "initFract": 0.25,
        "threshold": 0.1,
        "reductionFactor": 0.9,
        "type": "gene"
    }},
    "drug_prioritization": {"trustrank": {
        "experimentalOnly": True,
        "interactions": True,
        "addInteractions": True,
        "nodesOnly": False,
        "direct": False,
        "approved": False,
        "damping": 0.85,
        "type": "gene",
        "topX": 100,
        "elements": False
    }, "centrality": {
        "experimentalOnly": True,
        "interactions": True,
        "addInteractions": True,
        "nodesOnly": False,
        "direct": True,
        "approved": True,
        "type": "gene",
        "topX": 100,
        "elements": False
    }
    }
}
file = "./seed_files/10/hyperthyroidism-9seeds.tsv"
index_file = "./seed_files/index.txt"
files = []

with open(index_file) as f:
    for line in f:
        files.append(line.strip())


results = {}


def add_results(goal, algorithm, file, exec_time):
    step = file.split("/")[2]
    case = file.split("/")[3]
    if goal not in results:
        results[goal] = {}
    if algorithm not in results[goal]:
        results[goal][algorithm] = {}
    if step not in results[goal][algorithm]:
        results[goal][algorithm][step] = {}
    if case not in results[goal][algorithm][step]:
        results[goal][algorithm][step][case] = set()
    results[goal][algorithm][step][case].add(exec_time)


stats = list()
for file in files:
    for goal, algorithms in arguments.items():
        for algorithm, params in algorithms.items():
            def_payload = {
                "userId": userId,
                "dbVersion": db_version,
                "algorithm": algorithm,
                "goal": goal,
                "params": {
                    "experimentalOnly": True,
                    "addInteractions": True,
                    "nodesOnly": False,
                    "type": "gene"
                },
                "tissue": "all",
                "selection": True,
                "experimentalOnly": True,
            }
            def_payload["params"].update(params)
            genes = read_seed_files(file)
            def_payload["nodes"] = genes
            exec_time = run_request(def_payload)
            print(f"Runtime for {algorithm} is {exec_time}")
            stats.append({'algorithm': algorithm, 'category': goal, 'runtime': exec_time, 'seed_size_step': file.split("/")[2], 'seeds': genes})

df = pd.DataFrame.from_records(stats)
df.to_csv('./nedrex_web_runtime.tsv', index=False, sep='\t')

