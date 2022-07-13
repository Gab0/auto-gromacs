
from typing import List
import MDAnalysis as mda

from pymol import cmd


def generate_filename(method: str, identifier: List[str], label: str) -> str:
    """ Generate a filename for an output pdb file."""
    segments = [method, *identifier, label]
    return "_".join(segments) + ".pdb"


def build_snapshot(us: List[mda.Universe], identifier: List[str], labels: List[str]):
    resolution_divisor = 5

    pdb_filenames = []

    for universe, label in zip(us, labels):
        fname = generate_filename("movie", identifier, label)
        atom_group = universe.select_atoms(sel="protein")

        with mda.Writer(fname, multiframe=True) as output_pdb:
            for idx, ts in enumerate(universe.trajectory):
                if not idx % resolution_divisor:
                    output_pdb.write(atom_group)

        #fname = generate_filename("snapshot", identifier, label)
        #group.write(fname)
        pdb_filenames.append(fname)

    for pdb in pdb_filenames:
        cmd.load(pdb)

    cmd.extra_fit(
        selection="(all)",
        #"name CA",
        #reference="1ake",
        method="super",
        object="aln_super"
    )

    fname = generate_filename("snapshot", identifier, "all")
    cmd.save(fname)


#def build_movie(us: List[mda.Universe], identifier: List[str], labels: List[str]):
