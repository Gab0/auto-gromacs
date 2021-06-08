#!/bin/python
import enum
import shutil
import argparse
import sys
import os
import re
import subprocess
import shutil
import datetime
import pathlib

from typing import List, Union
from .core.messages import welcome_message, backup_folder_already_exists
from .core import settings


bashlog = None

EM_MDP = "ions.mdp"
EMW_MDP = "ions.mdp"
NVT_MDP = "nvt.mdp"


def handle_error(ERROR, step_no, log_file=None):
    if ERROR == 0:
        print('STEP%s Completed!' %  step_no)
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
            print("HEADS UP: Command failed for Step %s: return code %i." % (step_no, ERROR))
            sys.exit(ERROR)


def read_settings_from_mdp(fpath):
    pattern = r"(\w+)[\t ]*=[\t ]*([\w\d\.]+)[\t ]*;*"
    with open(fpath) as mdp:
        for line in mdp.read().splitlines():
            cat = re.findall(pattern, line)
            if cat:
                name, parameter = cat[0]
                yield name, parameter


class GromacsSimulation(object):
    def __init__(self, arguments):
        self.protein_file_path = arguments.protein
        if not self.protein_file_path and not arguments.resume:
            print("No input protein specified.")
            sys.exit(1)

        self.ligand_file_path = arguments.ligand
        if not self.ligand_file_path:
            pass
            #self.ligand_file_path = self.protein_file_path.split('.')[0] + '.gro'

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

    def to_wd(self, f):
        return os.path.join(self.working_dir, f)

    def path_state_file(self):
        return self.to_wd("md.cpt")

    def fix_includes(self, fpath):
        with open(fpath) as f:
            contents = f.read()

        BASE = '#include "'
        QUERY = BASE + self.working_dir + "/*"

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

    def path_log(self, n):
        return self.to_wd(f"step{n}.log")

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

    @staticmethod
    def run_process(step_no, step_name: str,
                    command: Union[List[str], str],
                    log_file=None, Input: str = None):
        print("INFO: Attempting to execute " + step_name + \
              " [STEP:" + step_no + "]")


        if isinstance(command, list):
            command = " ".join(command)

        # Assure logging is enabled;
        if ">" not in command and log_file is not None:
            command += f" >> {log_file} 2>&1"

        if bashlog is not None:
            bashlog.write("%s\n" % command)

        if self.dummy:
            return

        ret = subprocess.Popen(command, stdin=subprocess.PIPE, shell=True)

        if Input is not None:
            Input = Input.encode("utf-8")
            ret.communicate(Input)

        ret.wait()

        handle_error(ret.returncode, step_no, log_file)


    def build_settings_summary(self, arguments):
        settings_file = "settings.csv"

        Summary = {
            "force field": arguments.FF,
            "water model": arguments.solvent
        }

        mdp_prefixes = [
            "ions",
            "nvt",
            "npt",
            "md"
        ]

        header = ["Stage", "Parameter", "Value"]
        mdp_parameters = []

        for k in sorted(Summary.keys()):
            mdp_parameters.append(["*", k, Summary[k]])

        for mdp_prefix in mdp_prefixes:
            data = read_settings_from_mdp(
                self.to_wd(mdp_prefix + ".mdp"))
            for param, value in data:
                k = [mdp_prefix.upper(), param, value]
                mdp_parameters.append(k)

        mdp_parameters = self.compact_parameter_list(mdp_parameters)
        with open(self.to_wd(settings_file), 'w') as f:
            for parameter in [header] + mdp_parameters:
                f.write(",".join(parameter) + "\n")

    def compact_parameter_list(self, parameters):
        compat = {}

        def make_key(parameter, value):
            return f"{parameter}:{value}"

        for stage, parameter, value in parameters:
            K = make_key(parameter, value)
            if K not in compat.keys():
                compat[K] = []
            compat[K].append(stage)

        output_parameters = []
        for stage, parameter, value in parameters:
            K = make_key(parameter, value)
            if K in compat.keys():
                stages = compat[K]
                message = "+".join(stages)
                output_parameters.append([message, parameter, value])
                del compat[K]

        return output_parameters

    def gather_files(self):
        if self.ligand_file_path and not os.path.isfile(
                self.ligand_file_path):
            print('Ligand file not found at ', self.ligand_file_path)
            needStepZero = 1

        elif not self.ligand_topology_file_path or \
             not os.path.isfile(self.ligand_topology_file_path):
            print('Ligand Topology file not found at ', \
                self.ligand_topology_file_path)
            needStepZero = 1

        elif not os.path.isfile(self.protein_file_path):
            print('Protein file not found at ', self.protein_file_path)
            sys.exit()

        else:
            print('All data files found.')

        if os.path.isdir(self.working_dir):
            print("Folder '" + self.working_dir + "' Aready exist")
            if os.path.isdir("BACKUP"):
                print("WARNING: Backup folder already exists. Removing Backup folder!")
                print(backup_folder_already_exists)
                shutil.rmtree("BACKUP")


            if os.rename(self.working_dir, "BACKUP"):
                print("Old " + self.working_dir + " was moved to BACKUP/")


        pathlib.Path(self.working_dir).mkdir(parents=True, exist_ok=True)

        global bashlog
        bashlog = open(os.path.join(self.working_dir, 'bashlog'), 'w')

        print("CHEERS: Working Directory " + self.working_dir + \
              " created Successfully")
        print("Moving the files to Working Directory" + self.working_dir)

        shutil.copy2(self.protein_file_path, self.working_dir + 'protein.pdb')

        if self.ligand_file_path:
            shutil.copy2(self.ligand_file_path, self.working_dir + 'ligand.pdb')
        #shutil.copy2(self.ligand_topology_file_path,
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

        print("-> STEP 1: Initiating Procedure to generate topology for %s." % TARGET)
        pdb2gmx = settings.g_prefix + "pdb2gmx"
        step_no = "1"
        step_name = "Topology Generation"

        FFs = [
            "gromos54a7",
            "amber03"
        ]

        log_file = "%sstep1_%s.log" % (self.working_dir,TARGET)

        POSRE_PATH = self.to_wd("posre.itp")
        TOPOL_PATH = self.to_wd("topol.top")
        command = [
            pdb2gmx,
            "-f", "%s%s.pdb" % (self.working_dir, TARGET),
            "-o", "%s%s.gro" % (self.working_dir, TARGET),
            "-ignh",
            "-i", POSRE_PATH,
            "-p", TOPOL_PATH,
            "-ff", arguments.FF,
            "-water", arguments.solvent
        ]

        #assert(os.path.isfile(POSRE_PATH))
        command = " ".join(command)
        self.run_process(step_no, step_name, command, log_file)

        # FIX TOPOLOGY INCLUDE PATHS
        # Gromacs won't consider we're creating files in another folder
        # When processing 'includes' for topology files.
        self.fix_includes(TOPOL_PATH)

    # MAYBE THIS IS NOT NEEDED.
    def conjugateProteinLigand(self):
        print("-> STEP 1.5: Conjugating protein and ligand.")
        step_no = "0"
        step_name = "Complex Conjugation"

    def prepare_system(self):
        exit()
        print("-> STEP 2: Initiating Precedure to merge two molecules.")
        start_from_line = 3  # or whatever line I need to jump to

        # TODO: WHAT IS THIS?
        protein = self.working_dir + "protein.gro"
        system = self.working_dir + "system.gro"
        ligand = self.working_dir + "ligand.gro"

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

        f1 = open(self.working_dir + 'topol.top', 'r')
        f2 = open(self.working_dir + 'topol_temp.top', 'w')

        f1.close()
        f2.close()

        # swaping the files to get the original file
        f1 = open(self.working_dir + 'topol.top', 'w')
        f2 = open(self.working_dir + 'topol_temp.top', 'r')
        for line in f2:
            f1.write(line)

        #f1.write("UNK        1\n")
        f1.close()
        f2.close()
        os.unlink(self.working_dir + 'topol_temp.top')
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
            "-f", self.working_dir + "protein.gro",
            "-o", self.working_dir + "newbox.gro",
            "-bt", "cubic",
            "-d", str(arguments.box_size),
            "-c"
        ]

        self.run_process(step_no, step_name, " ".join(command), log_file)

        print(">STEP4 : Initiating Procedure to Solvate Complex")
        genbox = settings.g_prefix + "genbox"
        step_no = "4"
        step_name = "Solvating the Box"
        command = genbox + " -cp " + self.working_dir + "newbox.gro -p " + \
            self.working_dir + "topol.top -cs spc216.gro -o " + \
            self.working_dir + "solv.gro >> " + self.working_dir + \
            "step4.log 2>&1"
        if self.run_process(step_no, step_name, command):
            return 1

        command = "gmx solvate -cp " + self.working_dir + "newbox.gro -p " + \
            self.working_dir + "topol.top -cs spc216.gro -o " + \
            self.working_dir + "solv.gro >> " + self.working_dir + \
            "step4.log 2>&1"
        return self.run_process(step_no, step_name, command)

    def write_em_mdp(self):
        self.load_mdp(EM_MDP)

    def load_mdp(self, mdpname):

        if os.path.isfile(mdpname):
            print(">Using user-defined %s" % mdpname)
            Source = mdpname
        else:
            print(">Writing built-in %s" % mdpname)
            Source = os.path.join(self.module_dir, "mdp", mdpname)

        Target = os.path.join(self.working_dir, mdpname)

        shutil.copy2(Source, Target)

    def add_ions(self):
        print(">STEP5 : Initiating Procedure to Add Ions & Neutralise the " \
              "Complex")

        # TODO: Better name. Whats this?
        grompp = settings.g_prefix + "grompp"
        step_no = "5"
        step_name = "Check Ions "

        log_file = self.path_log(step_no)

        command = [
            grompp,
            "-f", self.working_dir + EM_MDP,
            "-c", self.working_dir + "solv.gro",
            "-p", self.working_dir + "topol.top",
            "-o", self.working_dir + "ions.tpr",
            "-po", self.working_dir + "mdout.mdp",
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
        with open(self.working_dir + 'step5.log') as f:
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

        if charge >= 0:
            print("System has positive charge .")
            print("Adding " + str(charge) + " CL ions to Neutralize the system")
            genion = settings.g_prefix + "genion"
            step_no = "6"
            step_name = "Adding Negative Ions "
            command = genion + " -s " + self.working_dir + "ions.tpr -o " + \
                self.working_dir + "solv_ions.gro -p " + self.working_dir + \
                "topol.top -nname CL -nn " + str(charge) + " >> " + \
                self.working_dir + "step6.log 2>&1"+\
                " << EOF\nSOL\nEOF"

            self.run_process(step_no, step_name, command)

        elif charge < 0:
            print("charge is negative")
            print("Adding " + str(-charge) + " CL ions to Neutralize the system")
            genion = settings.g_prefix + "genion"
            step_no = "6"
            step_name = "Adding Positive Ions "
            command = [
                genion,
                "-s", self.working_dir + "ions.tpr",
                "-o", self.working_dir + "solv_ions.gro",
                "-p", self.working_dir + "topol.top -pname NA -np",
                str(-charge),
                "<< EOF\nSOL\nEOF"
            ]
            self.run_process(step_no, step_name, command)
        elif charge == 0:
            print("System has Neutral charge , No adjustments Required :)")
            try:
                shutil.copy(self.working_dir + 'ions.tpr', self.working_dir + "solv_ions.tpr")
            except FileNotFoundError:
                pass

        print("DOUBLE CHEERS: SUCCESSFULLY PREPARED SYSTEM FOR SIMULATION")

    def create_em_mdp(self):
        self.load_mdp(EMW_MDP)

    def minimize(self, arguments):
        print(">STEP7 : Preparing the files for Minimisation")
        # grompp -f em_real.mdp -c solv_ions.gro -p topol.top -o em.tpr
        # mdrun -v -deffnm em
        grompp = settings.g_prefix + "grompp"
        mdrun = settings.g_prefix + "mdrun"
        step_no = "7"
        step_name = "Prepare files for Minimisation"
        # max warn 3 only for now
        command = grompp + " -f " + self.working_dir + EMW_MDP + " -c " +\
            self.working_dir + "solv_ions.gro -p " + self.working_dir\
            + "topol.top -o " + self.working_dir + "em.tpr -po " +\
            self.working_dir + "mdout.mdp -maxwarn 3 > " + self.working_dir\
            + "step7.log 2>&1"
        self.run_process(step_no, step_name, command)

        step_no = "8"
        step_name = " Minimisation"

        command = self.base_mdrun("em")

        if not arguments.dummy:
            self.run_process(step_no, step_name, command)

    def nvt(self, arguments):
        print(">STEP9 : Initiating the Procedure to Equilibrate the System")
        print("Beginging Equilibration with NVT Ensemble")
        grompp = settings.g_prefix + "grompp"
        mdrun = settings.g_prefix + "mdrun"
        step_no = "9"
        step_name = "Preparing files for NVT Equilibration"

        if not (os.path.isfile(self.working_dir + NVT_MDP)):
            self.load_mdp(NVT_MDP)

        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = [
            grompp,
            "-f", self.working_dir + "nvt.mdp",
            "-c", self.working_dir + "em.gro",
            "-r", self.working_dir + "em.gro",
            "-p", self.working_dir + "topol.top",
            "-o", self.working_dir + "nvt.tpr",
            "-po", self.working_dir + "mdout.mdp",
            "-maxwarn 3",
            ">", self.working_dir + "step9.log 2>&1"
        ]

        if not arguments.dummy:
            self.run_process(step_no, step_name, command)

        step_no = "10"
        step_name = "NVT Equilibration"
        command = self.base_mdrun("nvt")
        command +=  ["-deffnm", self.working_dir + "nvt"]

        if not arguments.dummy:
            self.run_process(step_no, step_name, command)

    def npt(self, arguments):
        print(">STEP11 : Initiating the Procedure to Equilibrate the System")
        print("Beginging Equilibration with NPT Ensemble")
        grompp = settings.g_prefix + "grompp"
        mdrun = settings.g_prefix + "mdrun"
        step_no = "11"
        step_name = "Preparing files for NPT Equilibration"

        if not (os.path.isfile(self.working_dir + "npt.mdp")):
            self.load_mdp("npt.mdp")

        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = grompp +\
                  " -f " + self.working_dir + "npt.mdp"+\
                  " -c " + self.working_dir + "nvt.gro" +\
                  " -r " + self.working_dir + "nvt.gro" +\
                  " -p " + self.working_dir + "topol.top" +\
                  " -o " + self.working_dir + "npt.tpr" +\
                  " -po "+ self.working_dir + "mdout.mdp" +\
                  " -maxwarn 3" +\
                  " > " + self.working_dir + "step11.log 2>&1"

        if not arguments.dummy:
            self.run_process(step_no, step_name, command)

        step_no = "12"
        step_name = "NPT Equilibration"
        command = self.base_mdrun("npt")

        if not arguments.dummy:
            self.run_process(step_no, step_name, command)

    def initmd(self):
        print(">STEP13 : Initiating the Production Run")
        grompp = settings.g_prefix + "grompp"
        step_no = "13"
        step_name = "Preparing files for NPT Equilibration"

        if not (os.path.isfile(self.working_dir + "md.mdp")):
            self.load_mdp("md.mdp")

        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = [
            grompp,
            "-f", self.working_dir + "md.mdp",
            "-c", self.working_dir + "npt.gro",
            "-p", self.working_dir + "topol.top",
            "-o", self.working_dir + "md.tpr",
            "-po", self.working_dir + "mdout.mdp",
            "-maxwarn", "3"
        ]

        command = " ".join(command)

        self.run_process(step_no, step_name, command, self.path_log(step_no))

    def base_mdrun(self, file_prefix="md"):
        mdrun = settings.g_prefix + "mdrun"

        mdrun_arguments_extensions = {
            "-s": ".tpr",
            "-c": ".gro",
            "-o": ".trr",
            "-e": ".edr",
            "-x": ".xtc",
            "-g": ".log"
        }

        arguments = [mdrun, "-v"]

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
                #"-update", "gpu"
            ]

        command = " ".join(command)
        log_file = self.path_log(step_no)
        if not arguments.dummy:
            self.run_process(step_no, step_name, command, log_file)

    def continue_mdrun(self, arguments):
        critical_file = self.path_state_file()
        simulation_log = self.path_log("14")

        grompp = settings.g_prefix + "grompp"
        step_no = "14"
        log_file = self.path_log(step_no)

        if arguments.refresh_mdp:
            self.load_mdp("md.mdp")
            GP = [
                grompp,
                "-f", self.to_wd("md.mdp"),
                "-c", self.to_wd("md.tpr"),
                "-p", self.to_wd("topol.top"),
                "-t", self.path_state_file()
            ]

            step_no = "14.1"
            self.run_process(step_no, "Update simulation", " ".join(GP),
                             self.path_log(step_no))

        if not os.path.isfile(critical_file):
            print(f"RESUME {critical_file} FILE NOT FOUND.")
            return None

        command = self.base_mdrun() + [
            "-cpi", critical_file,
            "-append"
        ]

        with open(simulation_log) as f:
            content = f.read()
            CMD = re.findall("gmx [ \-\d\w/\.]+", content)

        chosen = sorted(CMD, key=len, reverse=True)[0]

        CMD = chosen.split(" ") + command
        command = " ".join(command)

        print(command)
        self.run_process(
            "15",
            "Resume simulation",
            command,
            log_file
        )

    def postprocess(self):
        trjconv = settings.g_prefix +  "trjconv"
        step_no = "POST1"
        log_file = self.path_log(step_no)

        commandA = [
            trjconv,
            "-f", self.to_wd("md.trr"),
            "-o", self.to_wd("md5.trr"),
            "-skip", str(5)
        ]

        commandA = " ".join(commandA)

        FinalFile = "mdf"

        commandB = [
            trjconv,
            "-pbc", "mol",
            "-ur", "compact",
            "-f", self.to_wd("md5.trr"),
            "-s", self.to_wd("md.tpr"),
            "-o", self.to_wd("mdf.trr")
        ]

        commandB = " ".join(commandB)

        if not self.dummy:
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

        commandCOV = [
            "gmx", "covar",
            "-s", "ref.pdb",
            "-f", "allpdb_bb.xtc"
        ]

        commandEIG = [
            "gmx anaeig",
            "-s", "ref.pdb",
            "-f", "allpdb_bb.xtc",
            "-extr", "extreme1_xray.pdb",
            "-first", "1", "-last", "1", "-nframes", "30"
        ]


class SessionAction(enum.Enum):
    Nothing = 0
    New = 1
    Resume = 2
    PostProcessOnly = 3
    Dummy = 4


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
            return SessionAction.Nothing

        elif F.endswith("md.trr"):
            TRR_Present = True

    if TRR_Present:
        return SessionAction.Resume

    return SessionAction.Nothing


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

    return parser.parse_args()


def run_pipeline(arguments):

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
        obj.write_em_mdp,
        obj.add_ions,
        obj.create_em_mdp,
        obj.minimize,
    ]

    STEPS_EXECUTE = [
        obj.nvt,
        obj.npt,
        obj.initmd,
        obj.build_settings_summary,
        obj.md,
        obj.postprocess,
    ]

    STEPS_RESUME = [
        obj.continue_mdrun
    ]

    STEPS_POSTPROCESS = [
        obj.postprocess
    ]

    if Action == SessionAction.Resume:
        STEPS = STEPS_RESUME + STEPS_POSTPROCESS

    elif Action == SessionAction.Nothing:
        print(f"""Nothing to do for working directory {arguments.working_dir}:
        .gro found.""")
        STEPS = []

    elif Action == SessionAction.PostProcessOnly:
        STEPS = STEPS_POSTPROCESS

    elif Action == SessionAction.Dummy:
        STEPS = STEPS_GATHER + STEPS_PREPARE

    elif Action == SessionAction.New:
        STEPS = STEPS_GATHER + STEPS_PREPARE + STEPS_EXECUTE + STEPS_POSTPROCESS

    for STEP in STEPS:
        now = datetime.datetime.now().strftime("%H:%M:%S")
        print(f"\n\t[{now}]")
        W = 'arguments' in STEP.__code__.co_varnames
        if W:
            STEP(arguments)
        else:
            STEP()


def main():
    arguments = parse_arguments()

    # FIXME: well...
    if arguments.working_dir[-1] != "/":
        arguments.working_dir += "/"

    run_pipeline(arguments)


if __name__ == "__main__":
    main()
