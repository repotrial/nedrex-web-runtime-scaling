import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = pd.read_csv("./nedrex_web_runtime.tsv", sep="\t")
data = data[data['runtime'] != -1]
sns.set(style="whitegrid")

plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 12
})

targets = ['module_identification', 'drug_prioritization']
custom_ticks = [10,25,50,100,250]

titles = {
    'module_identification': "Runtime Scaling for Module Identification",
    'drug_prioritization': "Runtime Scaling for Drug Prioritization"
}

algorithm_label_map = {
    "diamond": "DIAMOnD",
    "must": "MuST",
    "kpm": "KeyPathwayMiner",
    "domino": "DOMINO",
    "robust": "ROBUST",
    "centrality": "Closeness Centrality",
    "trustrank": "TrustRank"
}

data['algorithm'] = data['algorithm'].map(algorithm_label_map)

for target in targets:
    target_data = data[data['category'] == target]

    plt.figure(figsize=(10, 6))
    ax = plt.subplot(1, 1, 1)

    sns.lineplot(data=target_data, x='seed_size_step', y='runtime', hue='algorithm', style='algorithm',
                 markers=True, dashes=False, ci='sd', err_style='band', ax=ax)

    ax.set_ylabel('Runtime (seconds)', fontsize=14)
    ax.set_xlabel('Number of Seed Genes', fontsize=14)
    ax.set_title(titles[target], fontsize=16)
    ax.set_yscale('log')
    ax.set_xticks(custom_ticks)
    ax.set_xticklabels(custom_ticks)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend(title='Algorithm', title_fontsize='13', fontsize='11', loc='upper left')

    plt.tight_layout()
    plt.savefig(f'./nedrex_runtime_scaling_{target}.png', dpi=600)

    # plt.show()

# saved_plots = [f'./drugstone_runtime_scaling_{target}.png' for target in targets]