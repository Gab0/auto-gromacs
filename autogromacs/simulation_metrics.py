
import re
import os
import argparse

import pandas as pd
from Bio.SeqUtils import ProtParam
from antigen_protocol.Mutation import StructureMutator


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", dest="directory")
    parser.add_argument("-r", dest="reference_structure")
    parser.add_argument("-o", dest="output_file")

    return parser.parse_args()


def get_sequence_properties(seq):
    analyzer = ProtParam.ProteinAnalysis(seq)
    return {
        "IP": analyzer.isoelectric_point(),
        "Gravy": analyzer.gravy()
    }


def process_categorized_simulation_name(name):
    replacements = {
        "DUMMY": "Artificial",
        "mutation": "Natural"
    }

    for word, repl in replacements.items():
        name = re.sub(word, repl + " ", name)

    return name


def get_seq(w):
    return "".join([q[1] for q in w])


def sort_simulation_rule(identifier):

    numbering = re.findall(r"[\d-]+", identifier)
    if not numbering:
        numbering = ""
    else:
        numbering = numbering[0]
    prefix = identifier.replace(numbering + r".*", "")
    prefix = prefix.strip()

    number = float(numbering.replace("-", "."))

    return prefix, number


def analyze_directory(options):
    relevant_files = {
        "structure": "protein.pdb",
    }

    reference_sequence = StructureMutator.loadStructureSequence(options.reference_structure)
    print(reference_sequence)
    data = []

    subdirectories = sorted(
        os.listdir(options.directory),
        key=sort_simulation_rule
    )

    for subdirectory in subdirectories:
        structure_path = os.path.join(
            options.directory,
            subdirectory,
            relevant_files["structure"]
        )

        entry_sequence = get_seq(StructureMutator.loadStructureSequence(structure_path))
        print(entry_sequence)

        len_a = len(entry_sequence)
        len_b = len(reference_sequence)
        if len_a != len_b:
            continue

        try:
            structure_mutations = list(StructureMutator.extract_mutations(
                reference_sequence,
                entry_sequence
            ))

            # FIXME: Solve StructureMutator so this is not needed!
            for mut in structure_mutations:
                mut.fromAA = mut.fromAA[-1]

        except IndexError:
            if len_a == len_b:
                raise

            print("ERROR: Incompatible protein lengths." +
                  f"{len_a} vs {len_b}"
                  )

        w = {
            "Nome": process_categorized_simulation_name(subdirectory),
            "Variações": "; ".join([
                m.show()
                for m in structure_mutations
            ])
        }

        additional_properties = get_sequence_properties(entry_sequence)
        w.update(additional_properties)

        data.append(w)

    print("---")
    df = pd.DataFrame(data)
    for col in df.columns:
        if df[col].dtype == "float":
            df[col] = df[col].round(3)

    df.to_csv(options.output_file, index=False)


def main():
    options = parse_arguments()
    analyze_directory(options)
