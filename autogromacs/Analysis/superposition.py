
from typing import List
import MDAnalysis as mda

from pymol import cmd


def build_snapshot(us: List[mda.Universe], identifier: List[str], labels: List[str]):
    ags = map(lambda u: u.select_atoms(sel="protein"), us)

    pdb_filenames = []
    for group, label in zip(ags, labels):
        fname = f"snapshot_{'_'.join(identifier)}_{label}.pdb"
        group.write(fname)
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

    cmd.save("output.pdb")
