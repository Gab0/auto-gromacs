from core.messages import welcome_message
import argparse


class ProteinLigMin(object):
    @staticmethod
    def welcome():
        """
        Prints out a welcome message, license info and the version.
        """
        print welcome_message

    def check_file(self, fielname):
        pass

    def exit_program(self):
        pass

    def run_process(self):
        pass

    def gather_files(self):
        pass

    def pdb2gmx_proc(self):
        pass

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

    args = parser.parse_args()

    print args
