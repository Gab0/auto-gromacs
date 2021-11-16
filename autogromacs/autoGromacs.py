#!/bin/python
import enum
import shutil
import argparse
import sys
import os
import re
import subprocess
import datetime
import pathlib

from typing import List, Union, Optional
from .core.messages import welcome_message
from .core import settings

from . import mdp_control

IONS_MDP = "ions.mdp"
EM_MDP = "em.mdp"
NVT_MDP = "nvt.mdp"


def message(func, message):
    def wrapper():
        print(message)
        func()


def handle_error(ERROR, step_no, log_file=None):
    if ERROR == 0:
        print('STEP%s Completed!' % step_no)
        print("")
        return
    else:
        if log_file is not None:
            READ_LOG = open(log_file).read()
            print(READ_LOG)

        if ERROR < 0:
            print(f"HEADS UP: Killed by signal {-ERROR} :(")
            sys.exit(ERROR)
        else:
            print("\n")
            print("HEADS UP: Command failed for Step %s: return code %i."
                  % (step_no, ERROR))
            sys.exit(ERROR)


class GromacsSimulation(object):
    bashlog = None

    def __init__(self, arguments):
        self.protein_file_path = arguments.protein
        if not self.protein_file_path and not arguments.resume:
            print("No input protein specified.")
            sys.exit(1)

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

        self.setup_gmx_executables()

    def setup_gmx_executables(self):
        gmx_commands = [
            "mdrun",
            "covar",
            "grompp",
            "genion",
            "do_dssp",
            "nmeig",
            "anaeig",
            "trjconv"
        ]

        for gmx_command in gmx_commands:
            setattr(self, gmx_command, settings.g_prefix + gmx_command)

    def to_wd(self, f, subdir: List[str] = []) -> str:
        basedir = os.path.join(self.working_dir, *subdir)
        pathlib.Path(basedir).mkdir(parents=True, exist_ok=True)

        return os.path.join(basedir, f)

    def path_state_file(self):
        return self.to_wd("md.cpt")

    def fix_includes(self, fpath):
        with open(fpath) as f:
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

            with open(fpath, 'w') as f:
                f.write(output)
        else:
            print("Query not found!")

    def path_log(self, n, extra: str = ""):
        m = ""
        if extra:
            m = "_" + extra
        return self.to_wd(f"step{n}{m}.log")

    @staticmethod
    def welcome():
        """
        Prints out a welcome message, license info and the version.
        """
        print(welcome_message)

    def check_file(self, filename):
        pass

    def exit_program(self):
        pass

    def run_process(
            self,
            step_no: str,
            step_name: str,
            command: Union[List[str], str],
            log_file=None,
            Input: Optional[str] = None):

        print("INFO: Attempting to execute " + step_name +
              " [STEP:" + step_no + "]")

        if isinstance(command, list):
            command = " ".join(command)

        # Assure logging is enabled;
        if ">" not in command and log_file is not None:
            command += f" >> {log_file} 2>&1"

        if self.bashlog is not None:
            self.bashlog.write("%s\n" % command)

        if self.dummy:
            return

        ret = subprocess.Popen(command, stdin=subprocess.PIPE, shell=True)

        if Input is not None:
            proc_input = Input.encode("utf-8")
            ret.communicate(proc_input)

        ret.wait()

        handle_error(ret.returncode, step_no, log_file)


    def gather_files(self):
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

    def pdb2gmx_coord(self, arguments):
        if self.pdb2gmx_proc(arguments, "protein"):
            return 1
        if self.ligand_file_path:
            if self.pdb2gmx_proc(arguments, "ligand"):
                return 1
            return self.prepare_system()

        return 0

    def pdb2gmx_proc(self, arguments, TARGET):
        assert (TARGET in ["protein", "ligand"])

        print("-> STEP 1: Initiating Procedure to generate topology for %s."
              % TARGET)
        pdb2gmx = settings.g_prefix + "pdb2gmx"
        step_no = "1"
        step_name = "Topology Generation"
        log_file = self.path_log(step_no, TARGET)

        POSRE_PATH = self.to_wd("posre.itp")
        TOPOL_PATH = self.to_wd("topol.top")
        command = [
            pdb2gmx,
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

    def prepare_system(self):
        exit()
        print("-> STEP 2: Initiating Precedure to merge two molecules.")
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

        f1 = open(self.to_wd('topol.top', 'r'))
        f2 = open(self.to_wd('topol_temp.top', 'w'))

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

    def solvate_complex(self, arguments):
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

        print(">STEP4 : Initiating Procedure to Solvate Complex")
        genbox = settings.g_prefix + "genbox"
        step_no = "4"
        step_name = "Solvating the Box"
        command = [
            genbox,
            "-cp", self.to_wd("newbox.gro"),
            "-p ", self.to_wd("topol.top"),
            "-cs spc216.gro -o ",
            self.to_wd("solv.gro")
        ]

        if self.run_process(step_no, step_name, command):
            return 1

        command = [
            "gmx solvate -cp " + self.to_wd("newbox.gro"),
            "-p", self.to_wd("topol.top"),
            "-cs spc216.gro",
            "-o", self.to_wd("solv.gro")
        ]

        return self.run_process(step_no, step_name, command)

    def add_ions(self, arguments):
        print(">STEP5 : Initiating Procedure to Add Ions & Neutralise the "
              "Complex")

        # TODO: Better name. Whats this?
        grompp = settings.g_prefix + "grompp"
        step_no = "5"
        step_name = "Check Ions "

        log_file = self.path_log(step_no)

        mdp_control.load_mdp(self, arguments, IONS_MDP)

        command = [
            grompp,
            "-f", self.to_wd(IONS_MDP),
            "-c", self.to_wd("solv.gro"),
            "-p", self.to_wd("topol.top"),
            "-o", self.to_wd("ions.tpr"),
            "-po", self.to_wd("mdout.mdp"),
            "-maxwarn", self.maxwarn,
            ">", log_file,
            "2>&1"
        ]
        command = " ".join(command)
        self.run_process(step_no, step_name, command, log_file)

        # calculating the charge of the system
        # TODO: What is this doing? word??? Better name!
        word = 'total'  # Your word
        charge = 0
        with open(self.to_wd('step5.log')) as f:
            for line in f:
                if word in line:
                    s_line = line.strip().split()
                    two_words = (s_line[s_line.index(word) + 1],
                                 s_line[s_line.index(word) + 2])
                    charge = two_words[1]
                    break

        # TODO: This charge varibale might break the code
        print("Charge of the system is %s " % charge)
        charge = float(charge)
        charge = round(charge)

        # TODO: REWORK THIS
        if charge >= 0:
            print("System has positive charge .")
            print(f"Adding {charge} CL ions to Neutralize the system")
            genion = settings.g_prefix + "genion"
            step_no = "6"
            step_name = "Adding Negative Ions "
            command = [
                genion,
                "-s", self.to_wd("ions.tpr"),
                "-o", self.to_wd("solv_ions.gro"),
                "-p", self.to_wd("topol.top"),
                "-nname CL -nn " + str(charge) + " >>",
                self.to_wd("step6.log"),
                "2>&1", "<< EOF\nSOL\nEOF"
            ]

            self.run_process(step_no, step_name, command)

        elif charge < 0:
            print("charge is negative")
            print(f"Adding {-charge} CL ions to Neutralize the system")
            genion = settings.g_prefix + "genion"
            step_no = "6"
            step_name = "Adding Positive Ions "
            command = [
                genion,
                "-s", self.to_wd("ions.tpr"),
                "-o", self.to_wd("solv_ions.gro"),
                "-p", self.to_wd("topol.top"),
                "-pname", "NA",
                "-np", str(-charge),
                "<< EOF\nSOL\nEOF"
            ]
            self.run_process(step_no, step_name, command)

        elif charge == 0:
            print("System has Neutral charge , No adjustments Required :)")
            try:
                shutil.copy(
                    self.to_wd('ions.tpr'),
                    self.to_wd("solv_ions.tpr")
                )
            except FileNotFoundError:
                pass

        print("DOUBLE CHEERS: SUCCESSFULLY PREPARED SYSTEM FOR SIMULATION")

    def minimize(self, arguments):
        print(">STEP7 : Preparing the files for Minimization")
        # grompp -f em_real.mdp -c solv_ions.gro -p topol.top -o em.tpr
        # mdrun -v -deffnm em
        grompp = settings.g_prefix + "grompp"
        step_no = "7"
        step_name = "Prepare files for Minimisation"
        # max warn 3 only for now
        command = [
            grompp + " -f " + self.to_wd(IONS_MDP),
            "-c", self.to_wd("solv_ions.gro"),
            "-p", self.to_wd("topol.top"),
            "-o", self.to_wd("em.tpr"),
            "-po", self.to_wd("mdout.mdp"),
            "-maxwarn 3"
        ]

        self.run_process(step_no, step_name, command)

        step_no = "8"
        step_name = " Minimisation"

        command = self.base_mdrun("em")

        self.run_process(step_no, step_name, command)

    def nvt(self, arguments):
        print(">STEP9 : Initiating the Procedure to Equilibrate the System")
        print("Beginging Equilibration with NVT Ensemble")
        grompp = settings.g_prefix + "grompp"
        step_no = "9"
        step_name = "Preparing files for NVT Equilibration"

        if not os.path.isfile(self.to_wd(NVT_MDP)):
            mdp_control.load_mdp(self, arguments, NVT_MDP)

        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = [
            grompp,
            "-f", self.to_wd("nvt.mdp"),
            "-c", self.to_wd("em.gro"),
            "-r", self.to_wd("em.gro"),
            "-p", self.to_wd("topol.top"),
            "-o", self.to_wd("nvt.tpr"),
            "-po", self.to_wd("mdout.mdp"),
            "-maxwarn 3"
            ]

        if not arguments.dummy:
            self.run_process(step_no, step_name, command)

        step_no = "10"
        step_name = "NVT Equilibration"
        command = self.base_mdrun("nvt")
        command += ["-deffnm", self.to_wd("nvt")]

        if not arguments.dummy:
            self.run_process(step_no, step_name, command)

    def npt(self, arguments):
        print(">STEP11 : Initiating the Procedure to Equilibrate the System")
        print("Beginging Equilibration with NPT Ensemble")
        grompp = settings.g_prefix + "grompp"
        step_no = "11"
        step_name = "Preparing files for NPT Equilibration"

        if not os.path.isfile(self.to_wd("npt.mdp")):
            mdp_control.load_mdp(self, arguments, "npt.mdp")

        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = [
            grompp,
            "-f", self.to_wd("npt.mdp"),
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
        command = self.base_mdrun("npt")

        if not arguments.dummy:
            self.run_process(step_no, step_name, command)

    def initmd(self, arguments):
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

    def base_mdrun(self, file_prefix="md"):

        mdrun_arguments_extensions = {
            "-s": ".tpr",
            "-c": ".gro",
            "-o": ".trr",
            "-e": ".edr",
            "-x": ".xtc",
            "-g": ".log",
            "-mtx": ".mtx"
        }

        arguments = [self.mdrun, "-v"]

        for arg, ext in mdrun_arguments_extensions.items():
            arguments += [arg, self.to_wd(file_prefix + ext)]

        arguments += [
            "-cpo", self.path_state_file()
        ]

        return arguments

    def md(self, arguments):
        print(">STEP14: Simulation stared for %s" % self.protein_file_path)
        step_no = "14"
        step_name = "Creating producion MD."

        command = self.base_mdrun()
        if arguments.gpu:
            command += [
                "-nb", "gpu",
                "-pme", "gpu",
                "-pmefft", "gpu",
                "-bonded", "gpu",
                # "-update", "gpu",
            ]

        command = " ".join(command)
        log_file = self.path_log(step_no)
        if not arguments.dummy:
            self.run_process(step_no, step_name, command, log_file)

    def continue_mdrun(self, arguments):
        critical_file = self.path_state_file()
        simulation_log = self.path_log("14")

        step_no = "14"
        log_file = self.path_log(step_no)

        if arguments.refresh_mdp:
            mdp_control.load_mdp(self, arguments, "md.mdp")
            GP = [
                self.grompp,
                "-f", self.to_wd("md.mdp"),
                "-c", self.to_wd("md.tpr"),
                "-p", self.to_wd("topol.top"),
                "-t", self.path_state_file()
            ]

            step_no = "14.1"
            self.run_process(
                step_no,
                "Update simulation",
                " ".join(GP),
                self.path_log(step_no)
            )

        if not os.path.isfile(critical_file):
            print(f"RESUME {critical_file} FILE NOT FOUND.")
            return None

        command = self.base_mdrun() + [
            "-cpi", critical_file,
            "-append"
        ]

        with open(simulation_log) as f:
            content = f.read()
            CMD = re.findall(r"gmx [ \-\d\w/\.]+", content)

        chosen = sorted(CMD, key=len, reverse=True)[0]

        CMD = chosen.split(" ") + command
        command = " ".join(command)

        print(command)
        self.run_process(
            step_no,
            "Resume simulation",
            command,
            log_file
        )

    def postprocess(self):
        step_no = "POST1"

        commandA = [
            self.trjconv,
            "-f", self.to_wd("md.trr"),
            "-o", self.to_wd("md5.trr"),
            "-skip", str(5)
        ]

        FinalFile = "mdf"

        commandB = [
            self.trjconv,
            "-pbc", "mol",
            "-ur", "compact",
            "-f", self.to_wd("md5.trr"),
            "-s", self.to_wd("md.tpr"),
            "-o", self.to_wd(FinalFile + ".trr")
        ]

        step_no = "POST_SKIP"
        self.run_process(
            step_no,
            "Skip frames",
            commandA,
            self.path_log(step_no)
        )

        step_no = "POST_PBC"
        self.run_process(
            step_no,
            "Resolve Periodic Boundary Conditions",
            commandB,
            self.path_log(step_no),
            Input="0"
        )

    def analysis(self):
        file_prefix = "mdf"
        step_no = "ANALYSIS"
        command = [
            self.do_dssp,
            "-f", self.to_wd(file_prefix + ".trr"),
            "-s", self.to_wd("md.gro")
        ]

        os.environ["DSSP"] = "/usr/bin/mkdssp"
        self.run_process(
            step_no,
            "Analyze results",
            command,
            self.path_log(step_no),
            Input="1"
        )
        step_no = "qw"

        command = [
            settings.g_prefix + "xpm2ps",
            "-f", self.to_wd("ss.xpm")
        ]

        self.run_process(
            step_no,
            "Convert plot",
            command,
            self.path_log(step_no)
        )

    def pca(self):
        step_no = "covar"
        SUBDIR = ["PCA"]

        EGVAL = self.to_wd("eigenval.xvg", SUBDIR)
        EGVEC = self.to_wd("eigenvec.trr", SUBDIR)

        commandCOV = [
            self.covar,
            "-s", "ref.pdb",
            "-f", "allpdb_bb.xtc",
            "-o", EGVAL,
            "-v", EGVEC,
            "-l", self.to_wd("covar.log", SUBDIR),
        ]

        commandEIG = [
            self.anaeig,
            "-s", self.to_wd("md.gro"),
            "-f", self.to_wd("mdf.trr"),
            "-v", EGVEC,
            "-eig", EGVAL,
            "-extr", "extreme1_xray.pdb",
            "-first", "1", "-last", "1", "-nframes", "30"
        ]

        self.run_process(
            step_no,
            "Covariance matrix",
            commandCOV,
            self.path_log(step_no),
            Input="44"
        )

        step_no = "anaeig"
        self.run_process(
            step_no,
            "Covariance matrix",
            commandEIG,
            self.path_log(step_no),
            Input="44"
        )

    def nma(self):
        step_no = "RERUN_HESSIAN"

        Hessian = self.to_wd("hessian.mtx", subdirs=["COV"])
        command = [
            self.mdrun,
            "-rerun", self.to_wd("mdf.trr"),
            "-s", self.to_wd("md.tpr"),
            "-mtx", Hessian
        ]

        self.run_process(
            step_no,
            "Create Hessian matrix",
            command,
            self.path_log(step_no)
        )

        step_no = "ANALYZE_HESSIAN"
        command = [
            self.nmeig,
            "-f", Hessian,
            "-s", self.to_wd("md.tpr")
        ]

        self.run_process(
            step_no,
            "Analyze hessian eigenvectors",
            command,
            self.path_log(step_no)
        )


class SessionAction(enum.Enum):
    Nothing = 0
    New = 1
    Resume = 2
    PostProcessOnly = 3
    Dummy = 4
    AnalysisOnly = 5


def session_action_decision(arguments) -> SessionAction:

    if not os.path.isdir(arguments.working_dir):
        if arguments.dummy:
            return SessionAction.Dummy

        return SessionAction.New

    TRR_Present = False
    for F in os.listdir(arguments.working_dir):
        if F.endswith("md.gro"):
            if arguments.postprocess_only:
                return SessionAction.PostProcessOnly
            if arguments.analysis_only:
                return SessionAction.AnalysisOnly
            return SessionAction.Nothing

        elif F.endswith("md.trr"):
            TRR_Present = True

    if TRR_Present:
        return SessionAction.Resume

    return SessionAction.Nothing


def backup_existing_directory(working_dir):
    W = os.path.split(working_dir)

    if not W[1]:
        From = W[0]
    else:
        From = working_dir

    if os.path.isdir(From):
        To = os.path.join(os.path.basename(From), "BACKUP")
        if os.path.isdir(To):
            shutil.rmtree(To)
        shutil.move(From, To)


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-l', '--ligand',
                        help='Input a ligand file [*.gro]')
    parser.add_argument('-i', '--itp', dest="ligand_topology",
                        help='Input a ligand topology file [*.itp]')
    parser.add_argument('-p', '--protein',
                        help='Input a protein file (default:protein.pdb)')
    parser.add_argument('-w', '--wdir', dest="working_dir",
                        help='Working Directory of project (default:work)',
                        default='work')
    parser.add_argument('-v', '--verbose', help='Loud and Noisy[default]',
                        action="store_true")
    parser.add_argument('-q', '--quiet', help='Be very quit',
                        action="store_true")

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

    mdp_control.add_option_override(parser, "MD", "dt")
    mdp_control.add_option_override(parser, "MD", "nsteps")

    return parser.parse_args()


def run_pipeline(arguments):

    if arguments.RemoveDirectory:
        backup_existing_directory(arguments.working_dir)

    obj = GromacsSimulation(arguments)

    Action = session_action_decision(arguments)

    obj.welcome()

    # DECLARE STEPS
    STEPS_GATHER = [
        obj.gather_files
    ]

    STEPS_PREPARE = [
        # mandatory steps;
        obj.pdb2gmx_coord,
        obj.solvate_complex,
        obj.add_ions,
        obj.minimize,
    ]

    STEPS_EXECUTE = [
        obj.nvt,
        obj.npt,
        obj.initmd,
        obj.md,
        obj.postprocess,
    ]

    STEPS_RESUME = [
        obj.continue_mdrun
    ]

    STEPS_POSTPROCESS = [
        obj.postprocess
    ]

    STEPS_ANALYSIS = [
        obj.analysis,
        obj.pca
    ]

    if Action == SessionAction.Resume:
        STEPS = STEPS_RESUME + STEPS_POSTPROCESS

    elif Action == SessionAction.Nothing:
        print(f"""Nothing to do for working directory {arguments.working_dir}:
        .gro found.""")
        STEPS = []

    elif Action == SessionAction.PostProcessOnly:
        STEPS = STEPS_POSTPROCESS

    elif Action == SessionAction.AnalysisOnly:
        STEPS = STEPS_ANALYSIS

    elif Action == SessionAction.Dummy:
        STEPS = STEPS_GATHER + STEPS_PREPARE

    elif Action == SessionAction.New:
        STEPS = STEPS_GATHER + \
            STEPS_PREPARE + \
            STEPS_EXECUTE + \
            STEPS_POSTPROCESS

    for STEP in STEPS:
        now = datetime.datetime.now().strftime("%H:%M:%S")
        print(f"\n\t[{now}]")
        W = 'arguments' in STEP.__code__.co_varnames
        if W:
            STEP(arguments)
        else:
            STEP()
    mdp_control.build_settings_summary(obj, arguments)


def main():
    arguments = parse_arguments()

    run_pipeline(arguments)


if __name__ == "__main__":
    main()
