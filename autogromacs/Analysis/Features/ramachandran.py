
import MDAnalysis as mda


from MDAnalysis.analysis.dihedrals import Ramachandran

from autogromacs.Analysis import mdplots, constants


def analyze(u: mda.Universe):

    sel = u.select_atoms(constants.STANDARD_SELECTION)
    return Ramachandran(sel).run()


def plot(rama, label: str, output_filepath: str):
    mdplots.plot_ramachandran(
        rama,
        label,
        output_filepath
    )
