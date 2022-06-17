
import argparse

from . import mdp_control


def parse_arguments():
    """ Parse CLI arguments. """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-l',
        '--ligand',
        help='Input a ligand file [*.gro]'
    )

    parser.add_argument(
        '-i',
        '--itp',
        dest="ligand_topology",
        help='Input a ligand topology file [*.itp]'
    )

    parser.add_argument(
        '-p',
        '--protein',
        default="protein.pdb",
        help='Input a protein file'
    )

    parser.add_argument(
        '-w',
        '--wdir',
        dest="working_dir",
        help='Working Directory of project (default:work)',
        default='work'
    )

    parser.add_argument(
        '-v',
        '--verbose',
        help='Loud and Noisy[default]',
        action="store_true"
    )

    parser.add_argument(
        '-q',
        '--quiet',
        help='Be very quiet.',
        action="store_true"
    )

    parser.add_argument(
        '--force-field',
        dest="FF",
        default="amber03",
        help="Gromacs force field to use."
    )

    parser.add_argument(
        '--water',
        dest="solvent",
        default="spce",
        help="Select gromacs water model."
    )

    parser.add_argument(
        "--dummy",
        action="store_true",
        help="Do not run simulations (for debugging)."
    )

    parser.add_argument(
        '--norunmd',
        dest="runmd",
        action="store_false",
        default=True
    )

    parser.add_argument(
        '-R',
        '--resume',
        action="store_true"
    )

    parser.add_argument(
        '--gpu',
        action="store_true",
        help="Use GPU on MD steps. Requres GROMACS compiled with GPU support."
    )

    parser.add_argument(
        '-S',
        action="store_true",
        dest="refresh_mdp",
        help="Replace the mdp in the working directory with another."
    )

    parser.add_argument(
        "--box-size",
        type=float,
        default=1.0,
        help="Solvation box size."
    )

    parser.add_argument(
        "-P",
        "--postprocess-only",
        action="store_true",
        help="Execute postprocessing steps only."
    )

    parser.add_argument(
        "-A",
        "--analysis-only",
        action="store_true",
        help="Execute analysis steps only."
    )

    parser.add_argument(
        "--remove-dir",
        dest="RemoveDirectory",
        action="store_true",
        help="Force removal of existing working directory."
    )

    parser.add_argument(
        "--hpc",
        action="store_true",
        help="This flag tunes the run to HPC environments, " +
        "by disabling critical GROMACS flags."
    )

    parser.add_argument(
        "--gpu-offload",
        type=str,
        default="",
        help="Custom calculation methods to offset to gpu, as a comma separated list." +
        "Will override all other GPU-related options. Example: 'pme,update'"
    )

    parser.add_argument(
        "--ntomp",
        type=int,
        default=0,
        help="Force a specific ntomp value for the MD runs."
    )

    parser.add_argument(
        "--ntmpi",
        type=int,
        default=0,
        help="Force a specific ntmpi value for the MD runs."
    )

    parser.add_argument(
        "--nname",
        default="CL",
        help="The identifier of the negative ion."
    )

    parser.add_argument(
        "--pname",
        default="NA",
        help="The identifier of the positive ion."
    )

    parser.add_argument(
        "--salt-concentration",
        type=float,
        default=0.0,
        help="Additional solute concentration in M."

    )

    mdp_control.add_option_override(parser, "MD", "dt")
    mdp_control.add_option_override(parser, "MD", "nsteps")
    mdp_control.add_option_override(parser, "MD", "nstlist")

    mdp_control.add_option_override(parser, "NVT", "nsteps")

    mdp_control.add_option_override(parser, "NPT", "nsteps")

    mdp_control.add_option_override(parser, "IONS", "nsteps")

    return parser.parse_args()
