import shutil
import argparse
import sys
import os
import subprocess
from core.messages import welcome_message, backup_folder_already_exists
from core import settings


class ProteinLigMin(object):
    def __init__(self, *args, **kwargs):
        self.ligand_file_path = kwargs.pop('ligand_file')
        self.ligand_topology_file_path = kwargs.pop('ligand_topology_file')
        self.protein_file_path = kwargs.pop('protein_file')
        self.working_dir = kwargs.pop('working_dir')

        self.verbose = kwargs.pop('verbose')
        self.quiet = kwargs.pop('quiet')

        # A user cant use both the verbose and the quiet flag together
        if self.verbose is True and self.quiet is True:
            print 'Can\'t use both the verbose and quiet flags together'
            sys.exit()

    @staticmethod
    def welcome():
        """
        Prints out a welcome message, license info and the version.
        """
        print welcome_message

    def check_file(self, filename):
        pass

    def exit_program(self):
        pass

    @staticmethod
    def run_process(step_no, step_name, command):
        print "INFO: Attempting to execute " + step_name + \
              " [STEP:" + step_no + "]"
        ret = subprocess.call(command, shell=True)
        if ret != 0:
            if ret < 0:
                print "HEADS UP: Killed by signal :(", -ret
                sys.exit()
            else:
                print "HEADS UP: Command failed with return code", ret
                sys.exit()
        else:
            print 'Completed!'

    def gather_files(self):
        if not os.path.isfile(self.ligand_file_path):
            print 'Ligand file not found at ', self.ligand_file_path
            sys.exit()

        elif not os.path.isfile(self.ligand_topology_file_path):
            print 'Ligand Topology file not found at ', \
                self.ligand_topology_file_path
            sys.exit()

        elif not os.path.isfile(self.protein_file_path):
            print 'Protein file not found at ', self.protein_file_path
            sys.exit()

        else:
            print 'All data files found'

        if os.path.isdir(self.working_dir):
            print "Folder '" + self.working_dir + "' Aready exist"
            if os.path.isdir("BACKUP"):
                print "ERROR: Backup folder already exists :( "
                print backup_folder_already_exists
                sys.exit()
            else:
                if os.rename(self.working_dir, "BACKUP"):
                    print "Old " + self.working_dir + " was moved to BACKUP/"

        os.mkdir(self.working_dir)
        print "CHEERS: Working Directory " + self.working_dir + \
              " created Successfully"
        print "Moving the files to Working Directory" + self.working_dir
        shutil.copy2(self.protein_file_path, self.working_dir + 'protein.pdb')
        shutil.copy2(self.ligand_file_path, self.working_dir + 'ligand.gro')
        shutil.copy2(self.ligand_topology_file_path,
                     self.working_dir + 'ligand.itp')

    def pdb2gmx_proc(self):
        print ">STEP1 : Initiating Procedure to generate topology for protein"
        pdb2gmx = settings.g_prefix + "pdb2gmx"
        step_no = "1"
        step_name = "Topology Generation"
        command = pdb2gmx + " -f " + self.working_dir + "protein.pdb -o " + \
                  self.working_dir + "protein.gro -ignh -p " + \
                  self.working_dir + "topol.top -i " + self.working_dir + \
                  "posre.itp -ff gromos53a6 -water spc >> " + \
                  self.working_dir + "step1.log 2>&1"
        self.run_process(step_no, step_name, command)

    def prepare_system(self):

        print ">STEP2 : Initiating Precedure to make system[Protein+Ligand]"
        start_from_line = 3  # or whatever line I need to jump to

        # TODO: WHAT IS THIS?
        protein = self.working_dir + "protein.gro"
        system = self.working_dir + "system.gro"
        ligand = self.working_dir + "ligand.gro"

        protein_file = open(protein, "r", 0)
        ligand_file = open(ligand, "r", 0)
        system_file = open(system, 'wa', 0)

        # get the last line of protein
        # get the count of Protein and Ligand files
        protien_lines_count = len(protein_file.readlines())
        ligand_lines_count = len(ligand_file.readlines())

        # print protien_lines_count
        # print ligand_lines_count
        # count of the system
        # TODO: Better name
        SystemCount = protien_lines_count + ligand_lines_count - 6
        # print SystemCount
        protein_file.close()
        ligand_file.close()

        # open files for reading
        protein_file = open(protein, "r", 0)
        ligand_file = open(ligand, "r", 0)

        system_file.write(
            "System.gro Designed for Simulation by [bngromacs.py]\n")
        system_file.write(str(SystemCount) + "\n")

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
        protein_file = open(protein, "r", 0)
        last_line = protein_file.readlines()[-1]
        # print last_line
        system_file.write(last_line)
        print "CHEERS: system.gro WAS GENERATED SUCCESSFULLY"

        f1 = open(self.working_dir + 'topol.top', 'r')
        f2 = open(self.working_dir + 'topol_temp.top', 'w')
        for line in f1:
            f2.write(line.replace('; Include water topology',
                                  '; Include Ligand topology\n #include "ligand.itp"\n\n\n; Include water topology '))
        f1.close()
        f2.close()
        # swaping the files to get the original file
        f1 = open(self.working_dir + 'topol.top', 'w')
        f2 = open(self.working_dir + 'topol_temp.top', 'r')
        for line in f2:
            f1.write(line)
        f1.write("UNK        1\n")
        f1.close()
        f2.close()
        os.unlink(self.working_dir + 'topol_temp.top')
        print "INFO: Topology File Updated with Ligand topology info "
        print "CHEERS: STEP[2] SUCCESSFULLY COMPLETED :)\n\n\n"

    def solvate_complex(self):
        # editconf -f system.gro -o newbox.gro -bt cubic -d 1 -c
        # genbox -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro


        print ">STEP3 : Initiating Procedure to Solvate Complex"
        editconf = settings.g_prefix + "editconf"
        step_no = "3"
        step_name = "Defining the Box"
        command = editconf + " -f " + self.working_dir + "system.gro -o " + self.working_dir + "newbox.gro -bt cubic -d 1 -c >> " + self.working_dir + "step3.log 2>&1"
        self.run_process(step_no, step_name, command)

        print ">STEP4 : Initiating Procedure to Solvate Complex"
        genbox = settings.g_prefix + "genbox"
        step_no = "4"
        step_name = "Solvating the Box"
        command = genbox + " -cp " + self.working_dir + "newbox.gro -p " + self.working_dir + "topol.top -cs spc216.gro -o " + self.working_dir + "solv.gro >> " + self.working_dir + "step4.log 2>&1"
        self.run_process(step_no, step_name, command)

    def write_em_mdp(self):
        pass

    def file_copy(self):
        pass

    def add_ions(self):
        pass

    def create_em_mdp(self):
        pass

    def minimize(self):
        pass

    def nvt(self):
        pass

    def npt(self):
        pass

    def md(self):
        pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-l', '--ligand',
                        help='Input a ligand file [*.gro]')
    parser.add_argument('-i', '--itp',
                        help='Input a ligand topology file [*.itp]')
    parser.add_argument('-p', '--protein',
                        help='Input a protein file (default:protein.pdb)')
    parser.add_argument('-w', '--wdir',
                        help='Working Directory of project (default:work)')
    parser.add_argument('-v', '--verbose', help='Loud and Noisy[default]',
                        action="store_true")
    parser.add_argument('-q', '--quiet', help='Be very quit',
                        action="store_true")

    arguments = parser.parse_args()

    # TODO: Think of a better name
    obj = ProteinLigMin(
        ligand_file=arguments.ligand,
        ligand_topology_file=arguments.itp,
        protein_file=arguments.protein,
        working_dir=arguments.wdir,
        verbose=arguments.verbose,
        quiet=arguments.quiet
    )

    obj.welcome()
    obj.gather_files()
