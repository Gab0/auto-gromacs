
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

from . import autoGromacs

from multiprocessing import Process, Pipe

ForceFields = [
    "amber03",
    "amber12sb",
    "amber14sb",
    "amber99",
    "amber96",
    "amber94",
    "charmm27",
    "gromos53a5",
    "gromos54a7",
]

ForceFields = ForceFields[:3]#["charmm36"]

class AutoGromacsArgs():
    itp = None
    ligand = None
    verbose = False
    quiet = True
    resume = False
    ligand_topology = ""

    def __init__(self, p, ff, runmd=False, w="work/"):
        self.protein = p
        self.FF = ff
        self.runmd = runmd
        self.working_dir = w


def read_directory():
    for F in os.listdir():
        if F.endswith(".pdb"):
            yield F


def execute(arguments, Files):
    results = np.zeros(shape=(len(Files), len(ForceFields)))
    for i, File in enumerate(Files):
        for f, FF in enumerate(ForceFields):
            arguments = AutoGromacsArgs(
                File,
                FF,
                runmd=arguments.runmd,
                w=f"work_{i + 1}_{f + 1}/"
            )

            w = Process(target=autoGromacs.run_pipeline, args=(arguments,))
            w.start()
            w.join()
            results[i, f] = w.exitcode

    return results


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--runmd', action="store_true", default=False)

    return parser.parse_args()


def main():
    arguments = parse_arguments()

    Files = list(read_directory())

    print("Files loaded:")
    print("\n".join(Files))
    print()
    print("Force Fields:")
    print("\n".join(ForceFields))
    print()
    input("OK?")

    results = execute(arguments, Files)

    show_results(results, Files, ForceFields)


def show_results(results, Files, ForceFields):
    print(results)

    fig, ax = plt.subplots()
    im = ax.imshow(results)

    # We want to show all ticks...
    ax.set_yticks(np.arange(len(Files)))
    ax.set_xticks(np.arange(len(ForceFields)))
    # ... and label them with the respective list entries
    ax.set_yticklabels(Files)
    ax.set_xticklabels(ForceFields)

    plt.tight_layout()
    plt.savefig("results.png")


if __name__ == "__main__":
    main()
