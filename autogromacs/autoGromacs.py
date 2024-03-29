#!/bin/python
from typing import Tuple, List, Union, Optional

import enum
import shutil
import sys
import os
import re
import subprocess
import datetime
import pathlib
import itertools

from .core.messages import welcome_message
from .core import settings

from . import cli_arguments, mdp_control

IONS_MDP = "ions.mdp"
EM_MDP = "em.mdp"
NVT_MDP = "nvt.mdp"
NPT_MDP = "npt.mdp"

STEP_NO = 0


def welcome():
    """
    Prints out a welcome message, license info and the version.
    """
    print(welcome_message)


def build_step_title(step_no: int, description: str):
    """ Build a message to identify a single pipeline step. """
    now = datetime.datetime.now().strftime("%H:%M:%S")
    return f"\t[{now}]> STEP {step_no}: {description}"


def pipeline_step(func):
    """ Decorator to help organize the pipeline steps. """

    global STEP_NO
    STEP_NO += 1

    def inner(pipeline, *args):
        step_no = STEP_NO
        print(build_step_title(step_no, func.__doc__))
        func(pipeline, step_no, func.__name__, *args[1:])

    return inner


def handle_error(error_code, step_no, log_file=None):
    """ Handles error codes. """
    if error_code == 0:
        print(f'STEP {step_no} Completed!')
        print("")
        return
    else:
        if log_file is not None:
            with open(log_file, encoding="utf8") as f:
                print(f.read())

        if error_code < 0:
            print(f"HEADS UP: Killed by signal {-error_code} ;{{")
            sys.exit(error_code)
        else:
            print("\n")
            print(f"HEADS UP: Command failed for Step {step_no}: return code {error_code}.")
            sys.exit(error_code)


def calculate_multithread_parameters(arguments):
    """ Calculate multithreading parameters. """
    ncores = os.cpu_count() / 2

    max_omp = 4
    if arguments.ntomp:
        ntomp = arguments.ntomp
    else:
        ntomp = os.getenv("OMP_NUM_THREADS")

    if not ntomp:
        ntomp = min(max_omp, ncores)

    if arguments.ntmpi:
        ntmpi = arguments.ntmpi
    else:
        ntmpi = max(1, round(ncores / max_omp))

    return [
        "-ntmpi", str(ntmpi),
        "-ntomp", str(ntomp)
    ]


class GromacsExecutables():
    """ List GROMACS' 'executables'. """
    gmx_commands = [
        "pdb2gmx",
        "solvate",
        "mdrun",
        "covar",
        "grompp",
        "genion",
        "do_dssp",
        "nmeig",
        "anaeig",
        "trjconv"
    ]

    def __init__(self):
        for gmx_command in self.gmx_commands:
            setattr(self, gmx_command, settings.g_prefix + gmx_command)


class GromacsSimulation():
    """Main controller of the pipeline."""
    bashlog = None

    def __init__(self, arguments):
        self.protein_file_path = arguments.protein

        self.ligand_file_path = arguments.ligand
        if not self.ligand_file_path:
            pass

        self.ligand_topology_file_path = arguments.ligand_topology
        if not self.ligand_topology_file_path:
            self.ligand_topology_file_path = 'posre.itp'

        self.working_dir = arguments.working_dir

        self.module_dir = os.path.dirname(__file__)

        self.verbose = arguments.verbose
        self.quiet = arguments.quiet
        self.dummy = arguments.dummy

        self.maxwarn = str(2)
        # A user cant use both the verbose and the quiet flag together
        if self.verbose is True and self.quiet is True:
            print('Can\'t use both the verbose and quiet flags together')
            sys.exit()

        self.gromacs = GromacsExecutables()

        # FIXME: Organize these variables:
        self.downsample_prefix = "mdf"

        self.STEP_NO = 0


    def to_wd(self, f, subdir: List[str] = []) -> str:
        basedir = os.path.join(self.working_dir, *subdir)
        pathlib.Path(basedir).mkdir(parents=True, exist_ok=True)

        return os.path.join(basedir, f)

    def path_state_file(self):
        return self.to_wd("md.cpt")

    def fix_includes(self, fpath):
        with open(fpath, encoding="utf8") as f:
            contents = f.read()

        BASE = '#include "'
        QUERY = BASE + os.path.join(self.working_dir, "*")

        print(f"Replace {QUERY} with {BASE}\n\t@{fpath}?")
        if re.findall(QUERY, contents):
            output = re.sub(
                QUERY,
                BASE,
                contents
            )

            print("Ok.")

            with open(fpath, 'w', encoding="utf8") as f:
                f.write(output)
        else:
            print("Query not found!")

    def path_log(self, code, extra: str = ""):
        """Builds the filename for log files."""
        msg = ""
        if extra:
            msg = "_" + extra
        return self.to_wd(f"step_{code}{msg}.log")

    def run_process(
            self,
            step_no: Union[str, int],
            step_name: str,
            command: Union[List[str], str],
            log_file=None,
            stdin_input: Optional[str] = None):
        """Execute an external process."""

        if not isinstance(step_no, str):
            step_no = str(step_no)

        print("INFO: Attempting to execute " + step_name +
              " [STEP:" + step_no + "]")

        if isinstance(command, list):
            command = " ".join(command)

        # Assure logging is enabled;
        if ">" not in command and log_file is not None:
            command += f" >> {log_file} 2>&1"

        if self.bashlog is not None:
            self.bashlog.write(f"{command}\n")

        if self.dummy:
            return

        with subprocess.Popen(command,
                              stdin=subprocess.PIPE, shell=True) as ret:

            if stdin_input is not None:
                proc_input = stdin_input.encode("utf-8")
                ret.communicate(proc_input)

            ret.wait()

            handle_error(ret.returncode, step_no, log_file)

    @pipeline_step
    def gather_files(self, step_no, step_name, arguments):
        if self.ligand_file_path and not os.path.isfile(
                self.ligand_file_path):
            print('Ligand file not found at ', self.ligand_file_path)

        elif any((
                    not self.ligand_topology_file_path,
                    not os.path.isfile(self.ligand_topology_file_path)
                )):
            print('Ligand Topology file not found at ',
                  self.ligand_topology_file_path)

        elif not os.path.isfile(self.protein_file_path):
            print('Protein file not found at ', self.protein_file_path)
            sys.exit()

        else:
            print('All data files found.')

        pathlib.Path(self.working_dir).mkdir(parents=True, exist_ok=True)

        self.bashlog = open(self.to_wd('bashlog'), 'w')

        print("CHEERS: Working Directory " + self.working_dir +
              " created Successfully")
        print("Moving the files to the Working Directory.")

        shutil.copy2(self.protein_file_path, self.to_wd('protein.pdb'))

        if self.ligand_file_path:
            shutil.copy2(
                self.ligand_file_path,
                self.to_wd('ligand.pdb')
            )
        # shutil.copy2(self.ligand_topology_file_path,
        #             self.working_dir + 'ligand.itp')

    @pipeline_step
    def pdb2gmx_coord(self, step_no, step_name, arguments):
        if self.pdb2gmx_proc(step_no, step_name, arguments, "protein"):
            return 1
        if self.ligand_file_path:
            if self.pdb2gmx_proc(step_no, step_name, arguments, "ligand"):
                return 1
            return self.prepare_system()

        return 0

    def pdb2gmx_proc(self, step_no, step_name, arguments, TARGET):
        assert (TARGET in ["protein", "ligand"])

        print("-> STEP 1: Initiating Procedure to generate topology for %s."
              % TARGET)

        step_no = "1"
        step_name = "Topology Generation"
        log_file = self.path_log(step_no, TARGET)

        POSRE_PATH = self.to_wd("posre.itp")
        TOPOL_PATH = self.to_wd("topol.top")
        command = [
            self.gromacs.pdb2gmx,
            "-f", self.to_wd(TARGET + ".pdb"),
            "-o", self.to_wd(TARGET + ".gro"),
            "-ignh",
            "-i", POSRE_PATH,
            "-p", TOPOL_PATH,
            "-ff", arguments.FF,
            "-water", arguments.solvent
        ]

        command = " ".join(command)
        self.run_process(step_no, step_name, command, log_file)

        # FIX TOPOLOGY INCLUDE PATHS
        # Gromacs won't consider we're creating files in another folder
        # When processing 'includes' for topology files.
        self.fix_includes(TOPOL_PATH)

    @pipeline_step
    def prepare_system(self):
        """ Merge two molecules (PDB files). """
        sys.exit()
        start_from_line = 3  # or whatever line I need to jump to

        # TODO: WHAT IS THIS?
        protein = self.to_wd("protein.gro")
        system = self.to_wd("system.gro")
        ligand = self.to_wd("ligand.gro")

        protein_file = open(protein, "r")
        ligand_file = open(ligand, "r")
        system_file = open(system, 'a')

        # get the last line of protein
        # get the count of Protein and Ligand files
        protien_lines_count = len(protein_file.readlines())
        ligand_lines_count = len(ligand_file.readlines())

        # print protien_lines_count
        # print ligand_lines_count
        # count of the system
        # TODO: Better name
        system_count = protien_lines_count + ligand_lines_count - 6
        protein_file.close()
        ligand_file.close()

        # open files for reading
        protein_file = open(protein, "r")
        ligand_file = open(ligand, "r")

        system_file.write(
            "System.gro Designed by autogromacs\n")
        system_file.write(str(system_count) + "\n")

        line_counter = 1
        for line in protein_file:
            if line_counter in range(start_from_line,
                                     protien_lines_count):  # start_from_line :
                # print line
                system_file.write(line)
            line_counter += 1
        protein_file.close()

        line_counter = 1
        for line in ligand_file:
            if line_counter in range(start_from_line, ligand_lines_count):
                # print line
                system_file.write(line)
            line_counter += 1

        # get the last line of protein [the coordinates of the center]
        protein_file = open(protein, "r")
        last_line = protein_file.readlines()[-1]
        # print last_line
        system_file.write(last_line)
        print("CHEERS: system.gro WAS GENERATED SUCCESSFULLY")

        f1 = open(self.to_wd('topol.top', 'r'), encoding="utf8")
        f2 = open(self.to_wd('topol_temp.top', 'w'), encoding="utf8")

        f1.close()
        f2.close()

        # swaping the files to get the original file
        f1 = open(self.to_wd('topol.top', 'w'))
        f2 = open(self.to_wd('topol_temp.top', 'r'))
        for line in f2:
            f1.write(line)

        # f1.write("UNK        1\n")
        f1.close()
        f2.close()
        os.unlink(self.to_wd('topol_temp.top'))
        print("INFO: Topology File Updated with Ligand topology info ")
        print("CHEERS: STEP[2] SUCCESSFULLY COMPLETED :)\n\n\n")

    @pipeline_step
    def solvate_complex(self, step_no, step_name, arguments):
        """Define the solvation box for the complex."""
        print(">STEP3 : Initiating Procedure to Solvate Complex")
        editconf = settings.g_prefix + "editconf"
        step_no = "3"
        step_name = "Defining the Solvation Box"

        log_file = self.path_log(step_no)

        # -f system.gro
        command = [
            editconf,
            "-f", self.to_wd("protein.gro"),
            "-o", self.to_wd("newbox.gro"),
            "-bt", "cubic",
            "-d", str(arguments.box_size),
            "-c"
        ]

        self.run_process(step_no, step_name, " ".join(command), log_file)

    @pipeline_step
    def solvate_box(self, step_no, step_name, arguments):
        """Solvate the complex."""
        print(">STEP4 : Initiating Procedure to Solvate Complex")
        genbox = settings.g_prefix + "genbox"
        step_no = "4"
        step_name = "Solvating the Box"

        solvent_file = "spc216.gro"

        command = [
            genbox,
            "-cp", self.to_wd("newbox.gro"),
            "-p", self.to_wd("topol.top"),
            "-o", self.to_wd("solv.gro")
        ]

        if self.run_process(step_no, step_name, " ".join(command)):
            return 1

        command = [
            self.gromacs.solvate,
            "-cp", self.to_wd("newbox.gro"),
            "-p", self.to_wd("topol.top"),
            "-cs", solvent_file,
            "-o", self.to_wd("solv.gro")
        ]

        return self.run_process(step_no, step_name, " ".join(command))

    @pipeline_step
    def add_ions(self, step_no, step_name, arguments):
        """Add ions to neutralize the complex."""

        print(">STEP5 : Initiating Procedure to Add Ions & Neutralise the "
              "Complex")

        step_no = "5"
        step_name = "Check Ions "

        log_file = self.path_log(step_no)

        mdp_control.load_mdp(self, arguments, IONS_MDP)

        command = [
            self.gromacs.grompp,
            "-f", self.to_wd(IONS_MDP),
            "-c", self.to_wd("solv.gro"),
            "-p", self.to_wd("topol.top"),
            "-o", self.to_wd("ions.tpr"),
            "-po", self.to_wd("mdout.mdp"),
            "-maxwarn", self.maxwarn
        ]
        command = " ".join(command)
        self.run_process(step_no, step_name, command, log_file)

        step_no = "6"
        step_name = "Neutralizing the complex"
        command = [
            self.gromacs.genion,
            "-s", self.to_wd("ions.tpr"),
            "-o", self.to_wd("solv_ions.gro"),
            "-p", self.to_wd("topol.top"),
            "-neutral",
            "-conc", str(arguments.salt_concentration),
            "-nname", arguments.nname,
            "-pname", arguments.pname
        ]

        self.run_process(step_no, step_name, " ".join(command), stdin_input="13")

        print("DOUBLE CHEERS: SUCCESSFULLY PREPARED SYSTEM FOR SIMULATION")

    @pipeline_step
    def minimize_prepare(self, step_no, step_name, arguments):
        """Prepare the files for the Energy Minimization steps."""

        command = [
            self.gromacs.grompp,
            "-f", self.to_wd(IONS_MDP),
            "-c", self.to_wd("solv_ions.gro"),
            "-p", self.to_wd("topol.top"),
            "-o", self.to_wd("em.tpr"),
            "-po", self.to_wd("mdout.mdp"),
            "-maxwarn 3"
        ]

        self.run_process(step_no, step_name, command)

    @pipeline_step
    def minimize(self, step_no, step_name, arguments):
        """Execute the Energy Minimization simulation."""

        command = self.base_mdrun(arguments, allow_gpu=False, file_prefix="em")

        self.run_process(step_no, step_name, command)

    @pipeline_step
    def nvt(self, step_no, step_name, arguments):
        """Equilibration for volume and temperature."""

        if not os.path.isfile(self.to_wd(NVT_MDP)):
            mdp_control.load_mdp(self, arguments, NVT_MDP)

        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = [
            self.gromacs.grompp,
            "-f", self.to_wd("nvt.mdp"),
            "-c", self.to_wd("em.gro"),
            "-r", self.to_wd("em.gro"),
            "-p", self.to_wd("topol.top"),
            "-o", self.to_wd("nvt.tpr"),
            "-po", self.to_wd("mdout.mdp"),
            "-maxwarn 3"
            ]

        self.run_process(step_no, step_name, command)

        step_no = "10"
        step_name = "NVT Equilibration"
        command = self.base_mdrun(arguments, file_prefix="nvt")
        command += ["-deffnm", self.to_wd("nvt")]

        self.run_process(step_no, step_name, command)

    @pipeline_step
    def npt(self, step_no, step_name, arguments):
        """Equilibration for pressure and temperature."""
        grompp = settings.g_prefix + "grompp"
        step_no = "11"
        step_name = "Preparing files for NPT Equilibration"

        if not os.path.isfile(self.to_wd(NPT_MDP)):
            mdp_control.load_mdp(self, arguments, NPT_MDP)

        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = [
            grompp,
            "-f", self.to_wd(NPT_MDP),
            "-c", self.to_wd("nvt.gro"),
            "-r", self.to_wd("nvt.gro"),
            "-p", self.to_wd("topol.top"),
            "-o", self.to_wd("npt.tpr"),
            "-po", self.to_wd("mdout.mdp"),
            "-maxwarn", "3"
        ]

        if not arguments.dummy:
            self.run_process(step_no, step_name, command)

        step_no = "12"
        step_name = "NPT Equilibration"
        command = self.base_mdrun(arguments, file_prefix="npt")

        if not arguments.dummy:
            self.run_process(step_no, step_name, command)

    @pipeline_step
    def initmd(self, step_no, step_name, arguments):
        """ Prepare the production run. """
        print(">STEP13 : Initiating the Production Run")
        grompp = settings.g_prefix + "grompp"
        step_no = "13"
        step_name = "Preparing files for NPT Equilibration"

        if not (os.path.isfile(self.to_wd("md.mdp"))):
            mdp_control.load_mdp(self, arguments, "md.mdp")

        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = [
            grompp,
            "-f", self.to_wd("md.mdp"),
            "-c", self.to_wd("npt.gro"),
            "-p", self.to_wd("topol.top"),
            "-o", self.to_wd("md.tpr"),
            "-po", self.to_wd("mdout.mdp"),
            "-maxwarn", "3"
        ]

        command = " ".join(command)

        self.run_process(step_no, step_name, command, self.path_log(step_no))

    def base_mdrun(self, arguments, allow_gpu=True, file_prefix="md"):
        """ Generate the base arguments for a MD session. """

        mdrun_arguments_extensions = {
            "-s": ".tpr",
            "-c": ".gro",
            "-o": ".trr",
            "-e": ".edr",
            "-x": ".xtc",
            "-g": ".log",
            "-mtx": ".mtx"
        }

        command = [self.gromacs.mdrun, "-v"]

        for arg, ext in mdrun_arguments_extensions.items():
            command += [arg, self.to_wd(file_prefix + ext)]

        if allow_gpu:
            command += get_gpu_arguments(arguments.gpu,
                                         arguments.hpc,
                                         arguments.gpu_offload)

        command += calculate_multithread_parameters(arguments)

        command = fix_pme_ranks(command)

        command += [
            "-cpo", self.path_state_file()
        ]

        return command

    @pipeline_step
    def main_md(self, step_no, step_name, arguments):
        """ Execute the production MD simulation. """
        step_no = "14"
        step_name = f"Running producion MD for {self.protein_file_path}."
        print(build_step_title(step_no, step_name))

        command = self.base_mdrun(arguments)

        command = " ".join(command)
        log_file = self.path_log(step_no)
        if not arguments.dummy:
            self.run_process(step_no, step_name, command, log_file)

    @pipeline_step
    def continue_mdrun(self, step_no, step_name, arguments):
        """ Resumes a MD simulation. """
        critical_file = self.path_state_file()
        simulation_log = self.path_log("14")

        step_no = "14"
        log_file = self.path_log(step_no)

        if arguments.refresh_mdp:
            mdp_control.load_mdp(self, arguments, "md.mdp")
            command = [
                self.gromacs.grompp,
                "-f", self.to_wd("md.mdp"),
                "-c", self.to_wd("md.tpr"),
                "-p", self.to_wd("topol.top"),
                "-t", self.path_state_file()
            ]

            step_no = "14.1"
            self.run_process(
                step_no,
                "Update simulation",
                " ".join(command),
                self.path_log(step_no)
            )

        if not os.path.isfile(critical_file):
            print(f"RESUME {critical_file} FILE NOT FOUND.")
            return

        command = self.base_mdrun(arguments) + [
            "-cpi", critical_file,
            "-append"
        ]

        with open(simulation_log, encoding="utf8") as f:
            content = f.read()
            cli_command = re.findall(r"gmx [ \-\d\w/\.]+", content)

        chosen = sorted(cli_command, key=len, reverse=True)[0]

        cli_command = chosen.split(" ") + command
        command = " ".join(command)

        print(command)
        self.run_process(
            step_no,
            "Resume simulation",
            command,
            log_file
        )

    @pipeline_step
    def skip_frames(self, step_no, step_name, arguments):
        """
        Extract a shorter trajectory from the original
        by skipping frames.
        """
        step_no = "SKIP"

        command = [
            self.gromacs.trjconv,
            "-f", self.to_wd("md.trr"),
            "-o", self.to_wd("md5.trr"),
            "-skip", str(5)
        ]

        self.run_process(
            step_no,
            "Skip frames",
            command,
            self.path_log(step_no)
        )

    @pipeline_step
    def solve_periodic_boundaries(self, step_no, step_name, arguments):
        """ Remove the periodic boundary conditions. """
        command = [
            self.gromacs.trjconv,
            "-pbc", "mol",
            "-ur", "compact",
            "-f", self.to_wd("md5.trr"),
            "-s", self.to_wd("md.tpr"),
            "-o", self.to_wd(self.downsample_prefix + ".trr")
        ]

        step_no = "REMOVE_PBC"
        self.run_process(
            step_no,
            "Resolve Periodic Boundary Conditions",
            command,
            self.path_log(step_no),
            stdin_input="0"
        )

    @pipeline_step
    def secondary_structure_analysis(self, step_no, step_name, arguments):
        """ Executes DSSP analysis. """
        step_no = "ANALYSIS"
        command = [
            self.gromacs.do_dssp,
            "-f", self.to_wd(self.downsample_prefix + ".trr"),
            "-s", self.to_wd("md.gro"),
            "-o", self.to_wd("ss.xpm"),
            "-sc", self.to_wd("scount.xvg")
        ]

        os.environ["DSSP"] = "/usr/bin/mkdssp"
        self.run_process(
            step_no,
            "Analyze results",
            command,
            self.path_log(step_no),
            stdin_input="1"
        )
        step_no = "qw"

        command = [
            settings.g_prefix + "xpm2ps",
            "-f", self.to_wd("ss.xpm"),
            "-o", self.to_wd("plot.eps")
        ]

        self.run_process(
            step_no,
            "Convert plot",
            command,
            self.path_log(step_no)
        )

    def do_anaeig(self, step_no, egvec, egval, vecs: Tuple[int, int], output_file):
        """ Execute eigenvector analysis using 'gmx anaeig'. """
        svecs = [str(v) for v in vecs]

        command = [
            self.gromacs.anaeig,
            "-s", self.to_wd("md.gro"),
            "-f", self.to_wd(self.downsample_prefix + ".trr"),
            "-v", egvec,
            "-eig", egval,
            "-extr", output_file,
            "-first", svecs[0],
            "-last", svecs[1],
            "-nframes", "30"
        ]
        step_no += "-".join(svecs)
        self.run_process(
            step_no,
            "Covariance matrix",
            command,
            self.path_log(step_no),
            stdin_input="4\n4"
        )

    @pipeline_step
    def pca(self, step_no, step_name, arguments):
        """Run Principal Component Analysis unsing GROMACS' utilities."""
        step_no = "covar"
        SUBDIR = ["PCA"]

        EGVAL = self.to_wd("eigenval.xvg", SUBDIR)
        EGVEC = self.to_wd("eigenvec.trr", SUBDIR)

        command = [
            self.gromacs.covar,
            "-s", self.to_wd("md.gro"),
            "-f", self.to_wd(self.downsample_prefix + ".trr"),
            "-av", self.to_wd("average.pdb", SUBDIR),
            "-o", EGVAL,
            "-v", EGVEC,
            "-l", self.to_wd("covar.log", SUBDIR),
        ]

        self.run_process(
            step_no,
            "Covariance matrix",
            command,
            self.path_log(step_no),
            stdin_input="4\n4"
        )

        self.do_anaeig(
            "anaeig",
            EGVEC, EGVAL,
            (1, 2),
            self.to_wd("extreme.pdb", SUBDIR)
        )

    def nma(self):
        step_no = "RERUN_HESSIAN"

        hessian = self.to_wd("hessian.mtx", subdir=["COV"])
        command = [
            self.gromacs.mdrun,
            "-rerun", self.to_wd(self.downsample_prefix + ".trr"),
            "-s", self.to_wd("md.tpr"),
            "-mtx", hessian
        ]

        self.run_process(
            step_no,
            "Create Hessian matrix",
            command,
            self.path_log(step_no)
        )

        step_no = "ANALYZE_HESSIAN"
        command = [
            self.gromacs.nmeig,
            "-f", hessian,
            "-s", self.to_wd("md.tpr")
        ]

        self.run_process(
            step_no,
            "Analyze hessian eigenvectors",
            command,
            self.path_log(step_no)
        )


# FIXME: Unused.
def calculate_solute_nmol(solvent_nmol: int, molar_concentration: float) -> int:
    mol_mass = {
        "NaCl": 58.443,
        "H20": 18.01528
    }

    # solvent_vol = solvent_nmol * 2.989 * 10e-26
    # solute_nmol = molar_concentration * scipy.constants.Avogadro * solvent_vol

    solute_nmol = 0.0187 * molar_concentration * solvent_nmol

    return round(solute_nmol)


def read_solvent_nmol(topology_path: str) -> int:
    with open(topology_path, encoding="utf8") as f:
        content = f.read()
        number = re.findall(r"SOL +(\d+)", content)
        return int(number[0])


def fix_pme_ranks(command):
    """
    Ensure PME will be executed as a single rank.
    Both checks performed here are important, because
    GROMACS can be VERY pedantic about the mdrun arguments.
    """
    def check_pattern(pattern):
        com = ",".join(command)
        return re.findall(pattern, com)

    pme_check = check_pattern(r"-pme,gpu")

    ntmpi_pat = check_pattern(r"-ntmpi,\d+")
    ntmpi_check = ntmpi_pat and ntmpi_pat[0] not in ["-ntmpi,0", "-ntmpi,1"]

    if pme_check and ntmpi_check:
        command += ["-npme", "1"]

    return command


def get_gpu_arguments(use_gpu, is_hpc, custom_offload: str = "") -> List[str]:
    """
    Manage additional arguments for when GPUs are used,
    while also considering usual HPC constraints
    which were tested with NVIDIA Volta video cards.
    """
    known_flags = ["pme", "pmefft", "bonded", "update", "nb"]

    def set_flags(flags: [str], target: str) -> List[str]:
        _command = [
            [f"-{flag.strip()}", target]
            for flag in flags
        ]
        return list(itertools.chain.from_iterable(_command))

    if custom_offload:
        return set_flags(list(filter(None, custom_offload.split(","))), "gpu")

    if not use_gpu:
        return set_flags(known_flags, "cpu")

    command = [
            "-nb", "gpu",
            # "-update", "gpu",
    ]

    if not is_hpc:
        command += [
            "-pme", "gpu",
            "-pmefft", "gpu",
            "-bonded", "gpu",
        ]

    return command


class SessionAction(enum.Enum):
    """ Holds which pipeline part will be executed. """
    NOTHING = 0
    NEW = 1
    RESUME = 2
    POST_PROCESS_ONLY = 3
    DUMMY = 4
    ANALYSIS_ONLY = 5


def session_action_decision(arguments) -> SessionAction:
    """
    Interprete the CLI arguments to decide
    which components of the pipeline will be executed.
    """

    if not os.path.isdir(arguments.working_dir):
        if arguments.dummy:
            return SessionAction.DUMMY

        return SessionAction.NEW

    trr_present = False
    for content_file in os.listdir(arguments.working_dir):
        if content_file.endswith("md.gro"):
            if arguments.postprocess_only:
                return SessionAction.POST_PROCESS_ONLY
            if arguments.analysis_only:
                return SessionAction.ANALYSIS_ONLY
            return SessionAction.NOTHING

        if content_file.endswith("md.trr"):
            trr_present = True

    if trr_present:
        return SessionAction.RESUME

    return SessionAction.NOTHING


def backup_existing_directory(working_dir):
    """
    Copies the entire contents of
    the current working directory if it exists.
    """
    wdir = os.path.split(working_dir)

    if not wdir[1]:
        from_path = wdir[0]
    else:
        from_path = working_dir

    if os.path.isdir(from_path):
        to_path = os.path.join(os.path.basename(from_path), "BACKUP")
        if os.path.isdir(to_path):
            shutil.rmtree(to_path)
        shutil.move(from_path, to_path)


def run_pipeline(arguments):
    """ Execute the pipeline. """

    if arguments.RemoveDirectory:
        backup_existing_directory(arguments.working_dir)

    pipeline = GromacsSimulation(arguments)

    Action = session_action_decision(arguments)

    welcome()

    # DECLARE STEPS
    STEPS_GATHER = [
        pipeline.gather_files
    ]

    STEPS_PREPARE = [
        # mandatory steps;
        pipeline.pdb2gmx_coord,
        pipeline.solvate_complex,
        pipeline.solvate_box,
        pipeline.add_ions,
        pipeline.minimize_prepare,
        pipeline.minimize,
    ]

    STEPS_EXECUTE = [
        pipeline.nvt,
        pipeline.npt,
        pipeline.initmd,
        pipeline.main_md,
    ]

    STEPS_RESUME = [
        pipeline.continue_mdrun
    ]

    STEPS_POSTPROCESS = [
        pipeline.skip_frames,
        pipeline.solve_periodic_boundaries
    ]

    STEPS_ANALYSIS = [
        pipeline.secondary_structure_analysis,
        pipeline.pca
    ]

    require_pdb = False
    if Action == SessionAction.RESUME:
        steps = STEPS_RESUME + STEPS_POSTPROCESS + STEPS_ANALYSIS

    elif Action == SessionAction.NOTHING:
        print(f"""Nothing to do for working directory {arguments.working_dir}:
        .gro found.""")
        steps = []

    elif Action == SessionAction.POST_PROCESS_ONLY:
        steps = STEPS_POSTPROCESS

    elif Action == SessionAction.ANALYSIS_ONLY:
        steps = STEPS_ANALYSIS

    elif Action == SessionAction.DUMMY:
        require_pdb = True
        steps = STEPS_GATHER + STEPS_PREPARE

    elif Action == SessionAction.NEW:
        require_pdb = True
        steps = STEPS_GATHER + \
            STEPS_PREPARE + \
            STEPS_EXECUTE + \
            STEPS_POSTPROCESS + \
            STEPS_ANALYSIS

    if require_pdb:
        if not pipeline.protein_file_path:
            print("No input protein specified.")
            sys.exit(1)

    for step in steps:
        step(pipeline, arguments)

    # -- Create MDP settings summary.
    mdp_control.build_settings_summary(pipeline, arguments)


def main():
    """ Execute pipeline. """
    args = cli_arguments.parse_arguments()

    run_pipeline(args)


if __name__ == "__main__":
    main()
