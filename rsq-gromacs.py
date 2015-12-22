import shutil
import argparse
import sys
import os
import subprocess
from core.messages import welcome_message, backup_folder_already_exists
from core import settings


class ProteinLigMin(object):
    def __init__(self, *args, **kwargs):
        self.ligand_file = kwargs.pop('ligand_file')
        self.ligand_topology_file = kwargs.pop('ligand_topology_file')
        self.protein_file = kwargs.pop('protein_file')
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
        if not os.path.isfile(self.ligand_file):
            print 'Ligand file not found at ', self.ligand_file
            sys.exit()

        elif not os.path.isfile(self.ligand_topology_file):
            print 'Ligand Topology file not found at ', \
                self.ligand_topology_file
            sys.exit()

        elif not os.path.isfile(self.protein_file):
            print 'Protein file not found at ', self.protein_file
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
        shutil.copy2(self.protein_file, self.working_dir + 'protein.pdb')
        shutil.copy2(self.ligand_file, self.working_dir + 'ligand.gro')
        shutil.copy2(self.ligand_topology_file, self.working_dir + 'ligand.itp')

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
        pass

    def solvate_complex(self):
        pass

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
