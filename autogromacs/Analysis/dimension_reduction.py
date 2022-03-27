import os
from umap import umap_ as umap

import MDAnalysis as MDA
from MDAnalysis.analysis.dihedrals import Ramachandran

from sklearn.manifold import TSNE

import matplotlib.pyplot as plt
import seaborn as sns


RETRIEVE_LAST = False


def plot_2D(values, labels, output_filepath):

    plt.clf()
    p1 = sns.scatterplot(
        x=values[:, 0],
        y=values[:, 1],
    )

    for i, label in enumerate(labels):
        V = 0.5
        p1.text(
            *values[i],
            label,
            horizontalalignment='left',
            size='small',
            color='black',
            weight='semibold',
            rotation=18
        )

    plt.title('UMAP projection of a MD simulation', fontsize=12)
    plt.tight_layout()
    plt.legend(labels)

    plt.savefig(f"{output_filepath}.png")

    plt.clf()


def umap_reduce(data):
    reducer = umap.UMAP()
    return reducer.fit_transform(data)

def tsne_reduce(data):
    tsne = TSNE(n_components=2, random_state=0)
    return tsne.fit_transform(data)

Reduce = tsne_reduce

