#!/bin/python
import shutil
import argparse
import sys
import os
import re
import subprocess
import shutil
import datetime

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
            print("HEADS UP: Killed by signal :(", -ERROR)
            sys.exit(ERROR)
        else:
            print("\n")
            print("HEADS UP: Command failed for Step %s: return code %i." % (step_no, ERROR))
            sys.exit(ERROR)


class ProteinLigMin(object):
    def __init__(self, *args, **kwargs):
        self.protein_file_path = kwargs.pop('protein_file')
        if not self.protein_file_path:
            exit()

        self.ligand_file_path = kwargs.pop('ligand_file')
        if not self.ligand_file_path:
            pass
            #self.ligand_file_path = self.protein_file_path.split('.')[0] + '.gro'

        self.ligand_topology_file_path = kwargs.pop('ligand_topology_file')
        if not self.ligand_topology_file_path:
            self.ligand_topology_file_path = 'posre.itp'

        self.working_dir = kwargs.pop('working_dir')

        self.module_dir = os.path.dirname(__file__)

        self.verbose = kwargs.pop('verbose')
        self.quiet = kwargs.pop('quiet')

        # A user cant use both the verbose and the quiet flag together
        if self.verbose is True and self.quiet is True:
            print('Can\'t use both the verbose and quiet flags together')
            sys.exit()

    def toWD(self, f):
        return os.path.join(self.working_dir, f)

    def fix_includes(self, fpath):
        with open(fpath) as f:
            contents = f.read()

        BASE = '#include "'
        output = re.sub(
            BASE + self.working_dir + "/*",
            BASE,
            contents
        )
        with open(fpath, 'w') as f:
            f.write(output)

    def makeLogPath(self, n):
        return self.toWD(f"step{n}.log")

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
    def run_process(step_no, step_name, command, log_file=None):
        print("INFO: Attempting to execute " + step_name + \
              " [STEP:" + step_no + "]")

        # Assure logging is enabled;
        if ">" not in command and log_file is not None:
            command += f">> {log_file} 2>&1"

        bashlog.write("%s\n" % command)
        ret = subprocess.call(command, shell=True)

        handle_error(ret, step_no, log_file)

    def gather_files(self):
        if self.ligand_file_path and not os.path.isfile(self.ligand_file_path):
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

        os.mkdir(self.working_dir)

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

        POSRE_PATH = self.toWD("posre.itp")
        TOPOL_PATH = self.toWD("topol.top")
        command = [
            pdb2gmx,
            "-f", "%s%s.pdb" % (self.working_dir, TARGET),
            "-o", "%s%s.gro" % (self.working_dir, TARGET),
            "-ignh",
            "-i", POSRE_PATH,
            "-p", TOPOL_PATH,
            "-ff", arguments.FF,
            "-water spce"
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
        """
        for line in f1:
            f2.write(line.replace('; Include water topology',
                                  '; Include Ligand topology\n #include '
                                  '" ' + self.working_dir + 'ligand.itp"\n\n\n; Include water topology ')
                     )"""
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

    def solvate_complex(self):
        # editconf -f system.gro -o newbox.gro -bt cubic -d 1 -c
        # genbox -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

        print(">STEP3 : Initiating Procedure to Solvate Complex")
        editconf = settings.g_prefix + "editconf"
        step_no = "3"
        step_name = "Defining the Box"

        log_file = self.makeLogPath(step_no)

        # -f system.gro
        command = [
            editconf,
            "-f", self.working_dir + "protein.gro",
            "-o", self.working_dir + "newbox.gro",
            "-bt", "cubic",
            "-d", "1",
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

        #MODF
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

        log_file = self.makeLogPath(step_no)

        command = [
            grompp,
            "-f", self.working_dir + EM_MDP,
            "-c", self.working_dir + "solv.gro",
            "-p", self.working_dir + "topol.top",
            "-o", self.working_dir + "ions.tpr",
            "-po", self.working_dir + "mdout.mdp",
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
                "-s" + self.working_dir + "ions.tpr",
                "-o" + self.working_dir + "solv_ions.gro",
                "-p", self.working_dir + "topol.top -pname NA -np",
                str(-charge),
                "<< EOF\nSOL\nEOF"
                ]
            command = " ".join(command)
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

    def minimize(self):
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

        command = mdrun + " -v  -s " + self.working_dir + "em.tpr -c " + \
            self.working_dir + "em.gro -o " + self.working_dir + \
            "em.trr -e " + self.working_dir + "em.edr -x " + \
            self.working_dir + "em.xtc -g " + self.working_dir + \
            "em.log > " + self.working_dir + "step8.log 2>&1"
        self.run_process(step_no, step_name, command)

    def nvt(self):
        print(">STEP9 : Initiating the Procedure to Equilibrate the System")
        print("Beginging Equilibration with NVT Ensemble")
        grompp = settings.g_prefix + "grompp"
        mdrun = settings.g_prefix + "mdrun"
        step_no = "9"
        step_name = "Preparing files for NVT Equilibration"

        """
        shutil.copy2(*[
            os.path.join(self.working_dir, x)
            for x in ["#posre.itp.1#", "posre.itp"]
        ])
        """

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
        command = " ".join(command)

        self.run_process(step_no, step_name, command)

        step_no = "10"
        step_name = "NVT Equilibration"
        command = mdrun + " -v"+\
                  " -s " + self.working_dir + "nvt.tpr"+\
                  " -c " + self.working_dir + "nvt.gro"+\
                  " -o " + self.working_dir + "nvt.trr"+\
                  " -e " + self.working_dir + "nvt.edr"+\
                  " -x " + self.working_dir + "nvt.xtc"+\
                  " -g " + self.working_dir + "nvt.log"+\
                  " -deffnm " + self.working_dir + "nvt"+\
                  " > " + self.working_dir + "step10.log 2>&1"
        # command = "gmx mdrun -deffnm nvt > step10.log 2>&1"
        self.run_process(step_no, step_name, command)

    def npt(self):
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
        self.run_process(step_no, step_name, command)

        step_no = "12"
        step_name = "NPT Equilibration"
        command = mdrun + " -v  -s " + self.working_dir + "npt.tpr -c " + \
            self.working_dir + "npt.gro -o " + self.working_dir + \
            "npt.trr -e " + self.working_dir + "npt.edr -x " + \
            self.working_dir + "npt.xtc -g " + self.working_dir + "npt.log > "\
            + self.working_dir + "step12.log 2>&1"
        self.run_process(step_no, step_name, command)

    def initmd(self):
        print("CHEERS :) WE ARE CLOSE TO SUCCESS :):)")
        print(">STEP13 : Initiating the Production Run")
        grompp = settings.g_prefix + "grompp"
        step_no = "13"
        step_name = "Preparing files for NPT Equilibration"

        if not (os.path.isfile(self.working_dir + "md.mdp")):
            self.load_mdp("md.mdp")

        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = grompp +\
                  " -f " + self.working_dir + "md.mdp" +\
                  " -c " + self.working_dir + "npt.gro" +\
                  " -p " + self.working_dir + "topol.top" +\
                  " -o " + self.working_dir + "md.tpr" +\
                  " -po "+ self.working_dir + "mdout.mdp" +\
                  " -maxwarn 3 "
        self.run_process(step_no, step_name, command, self.makeLogPath(step_no))

    def md(self):
        mdrun = settings.g_prefix + "mdrun"
        step_no = "14"
        step_name = "Creating producion MD."
        command = mdrun +\
                  " -v  -s " + \
                  self.working_dir + "md.tpr"+ \
                  " -c " + \
                  self.working_dir + \
                  "md.gro -o " + \
                  self.working_dir + \
                  "md.trr -e " + \
            self.working_dir + "md.edr -x " + self.working_dir + "md.xtc -g " +\
            self.working_dir + "md.log"

        log_file = self.makeLogPath(step_no)
        self.run_process(step_no, step_name, command, log_file)


    def continue_mdrun(self):
        critical_file = self.makeLogPath(14)

        if os.path.isfile(critical_file):
            with open(critical_file) as f:
                content = f.read()
                CMD = re.findall("gmx [ \-\d\w/\.]+", content)
                chosen = sorted(CMD, key=len, reverse=True)[0]
                print(CMD)
                print(chosen)

        else:
            print(f"RESUME {critical_file} FILE NOT FOUND.")

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-l', '--ligand',
                        help='Input a ligand file [*.gro]')
    parser.add_argument('-i', '--itp',
                        help='Input a ligand topology file [*.itp]')
    parser.add_argument('-p', '--protein',
                        help='Input a protein file (default:protein.pdb)')
    parser.add_argument('-w', '--wdir',
                        help='Working Directory of project (default:work)',
                        default='work')
    parser.add_argument('-v', '--verbose', help='Loud and Noisy[default]',
                        action="store_true")
    parser.add_argument('-q', '--quiet', help='Be very quit',
                        action="store_true")

    parser.add_argument('--force-field', dest="FF",
                        default="amber03", help="Gromacs force field to use.")

    parser.add_argument('--norunmd', dest="runmd",
                        action="store_false", default=True)

    parser.add_argument('--resume',
                        action="store_true")

    return parser.parse_args()


def run_pipeline(arguments):

    # TODO: Think of a better name
    obj = ProteinLigMin(
        ligand_file=arguments.ligand,
        ligand_topology_file=arguments.itp,
        protein_file=arguments.protein,
        working_dir=arguments.wdir,
        verbose=arguments.verbose,
        quiet=arguments.quiet
    )
    if arguments.resume:
        print("RESUMING")
        obj.continue_mdrun()
        sys.exit(1)

    obj.welcome()
    obj.gather_files()


    STEPS = [
        # mandatory steps;
        obj.pdb2gmx_coord,
        obj.solvate_complex,
        obj.write_em_mdp,
        obj.add_ions,
        obj.create_em_mdp,
        obj.minimize,
    ]

    if arguments.runmd:
        STEPS += [
            # evaluation steps;
            obj.nvt,
            obj.npt,
            obj.initmd,
            obj.md
        ]


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
    if arguments.wdir[-1] != "/":
        arguments.wdir += "/"

    run_pipeline(arguments)


if __name__ == "__main__":
    main()
