
from typing import List, Dict, Any
import os

import MDAnalysis as mda

from antigen_protocol.Mutation import structure_name


def access_dict_key_fallback(dct: Dict[str, Any], keys: List[str]) -> Any:
    for target_key in keys:
        try:
            return dct[target_key]
        except KeyError:
            pass
        raise KeyError(f"Unable to get any of\n{dct}\n{keys}.")


def extract_mutation_labels(recipe_filepath: str, identifiers: List[str]) -> List[str]:

    variant_map = read_mutation_recipes(recipe_filepath)

    def get_variant_map(label):
        access_dict_key_fallback(
            variant_map, [label, structure_name.process_simulation_name(label)]
        )

    return [get_variant_map(label) for label in identifiers]


def read_mutation_recipes(filepath: str) -> Dict[str, str]:

    def read_recipe(recipe_file: str) -> str:
        with open(os.path.join(filepath, recipe_file), encoding="utf-8") as f:
            return " ".join(list(f.readlines()))

    return {
        recipe_file: read_recipe(recipe_file)
        for recipe_file in os.listdir(filepath)
    }


def get_label(u: mda.Universe) -> str:
    """Extract a label to identify a single universe."""

    identifier = clean_universe_prefix(u.filename)

    if False:
        t = u.trajectory.totaltime
        T = "t=%.2fns" % (t / 1000)
        return identifier + " " + T

    else:
        return structure_name.process_simulation_name(identifier)


def clean_universe_prefix(p: str) -> str:
    """Transforms the universe filepath into an indentifier."""
    try:
        p = p.split("/")[-2]
    except IndexError:
        pass

    p = p.replace("work", "")
    p = p.strip("_")

    return p
