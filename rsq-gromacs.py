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
        system_count = protien_lines_count + ligand_lines_count - 6
        protein_file.close()
        ligand_file.close()

        # open files for reading
        protein_file = open(protein, "r", 0)
        ligand_file = open(ligand, "r", 0)

        system_file.write(
            "System.gro Designed for Simulation by [bngromacs.py]\n")
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
        print ">NOTE: Writing em.mdp"
        # TODO: Better name
        some_file = open(self.working_dir + "em.mdp", "w")
        data = """
        ; LINES STARTING WITH ';' ARE COMMENTS
        title         = Minimization    ; Title of run

        ; Parameters describing what to do, when to stop and what to save
        integrator    = steep        ; Algorithm (steep = steepest descent minimization)
        emtol         = 100.0      ; Stop minimization when the maximum force < 1.0 kJ/mol
        emstep        = 0.01      ; Energy step size
        nsteps        = 50000          ; Maximum number of (minimization) steps to perform
        energygrps    = system    ; Which energy group(s) to write to disk

        ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
        nstlist       = 1            ; Frequency to update the neighbor list and long range forces
        ns_type       = grid        ; Method to determine neighbor list (simple, grid)
        rlist         = 1.0        ; Cut-off for making neighbor list (short range forces)
        coulombtype   = PME        ; Treatment of long range electrostatic interactions
        rcoulomb      = 1.0        ; long range electrostatic cut-off
        rvdw          = 1.0        ; long range Van der Waals cut-off
        pbc           = xyz         ; Periodic Boundary Conditions (yes/no)
        """
        some_file.write(str(data))
        some_file.close()

    @staticmethod
    def file_copy(source, destination):
        # TODO: There must be something better in the os module?
        in_file = open(source, 'r')
        out_file = open(destination, 'w')
        temp = in_file.read()
        out_file.write(temp)
        in_file.close()
        out_file.close()

    def add_ions(self):
        print ">STEP5 : Initiating Procedure to Add Ions & Neutralise the Complex"
        # grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr
        # genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -nn X -np X
        # TODO: Better name. Whats this?
        grompp = settings.g_prefix + "grompp"
        step_no = "5"
        step_name = "Check Ions "
        command = grompp + " -f " + self.working_dir + "em.mdp -c " + self.working_dir + "solv.gro -p " + self.working_dir + "topol.top -o " + self.working_dir + "ions.tpr -po " + self.working_dir + "mdout.mdp > " + self.working_dir + "step5.log 2>&1"
        self.run_process(step_no, step_name, command)

        # calculating the charge of the system
        # TODO: What iis this doing? word??? Better name!
        word = 'total'  # Your word

        with open(self.working_dir + 'step5.log') as f:
            for line in f:
                if word in line:
                    s_line = line.strip().split()
                    two_words = (s_line[s_line.index(word) + 1],
                                 s_line[s_line.index(word) + 2])
                    charge = two_words[1]
                    break

        # TODO: This charge varibale might break the code
        print "Charge of the system is " + charge
        charge = float(charge)
        charge = round(charge)

        if charge > 0:
            print "System has positive charge ."
            print "Adding " + str(charge) + " CL ions to Neutralize the system"
            genion = settings.g_prefix + "genion"
            step_no = "6"
            step_name = "Adding Negative Ions "
            command = genion + " -s " + self.working_dir + "ions.tpr -o " + self.working_dir + "solv_ions.gro -p " + self.working_dir + "topol.top -nname CL -nn " + str(
                charge) + " -g " + self.working_dir + "step6.log 2>&1 << EOF\nSOL\nEOF"
            self.run_process(step_no, step_name, command)

        elif charge < 0:
            print "charge is negative"
            print "Adding " + str(-charge) + " CL ions to Neutralize the system"
            genion = settings.g_prefix + "genion"
            step_no = "6"
            step_name = "Adding Positive Ions "
            command = genion + " -s " + self.working_dir + "ions.tpr -o " + self.working_dir + "solv_ions.gro -p " + self.working_dir + "topol.top -pname NA -np " + str(
                -charge) + " << EOF\nSOL\nEOF"
            self.run_process(step_no, step_name, command)

        elif charge == 0:
            print "System has Neutral charge , No adjustments Required :)"
            self.file_copy('work/ions.tpr', "work/solv_ions.tpr")

        print "DOUBLE CHEERS: SUCCESFULY PREPARED SYSTEM FOR SIMULATION"

    def create_em_mdp(self):
        # TODO: better name
        some_file = open(self.working_dir + "em_real.mdp", "w")
        em_mdp = """
        ; LINES STARTING WITH ';' ARE COMMENTS
        title        = Minimization    ; Title of run

        ; Parameters describing what to do, when to stop and what to save
        integrator     = steep        ; Algorithm (steep = steepest descent minimization)
        emtol          = 100.0      ; Stop minimization when the maximum force < 1.0 kJ/mol
        emstep         = 0.01      ; Energy step size
        nsteps         = 50000          ; Maximum number of (minimization) steps to perform
        energygrps     = Protein UNK    ; Which energy group(s) to write to disk

        ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
        nstlist        = 1            ; Frequency to update the neighbor list and long range forces
        ns_type        = grid        ; Method to determine neighbor list (simple, grid)
        rlist          = 1.0        ; Cut-off for making neighbor list (short range forces)
        coulombtype    = PME        ; Treatment of long range electrostatic interactions
        rcoulomb       = 1.0        ; long range electrostatic cut-off
        rvdw           = 1.0        ; long range Van der Waals cut-off
        pbc            = xyz         ; Periodic Boundary Conditions (yes/no)
        """
        some_file.write(em_mdp)
        print "CHEERS: em_real.mdp SUCCESSFULLY GENERATED :)"

    def minimize(self):
        print "MESSAGE: "
        t = raw_input(
            "Did you check your complex !! do you wish to continue: (y/n)")

        if (t == 'y'):

            print ">STEP7 : Preparing the files for Minimisation"
            # grompp -f em_real.mdp -c solv_ions.gro -p topol.top -o em.tpr
            # mdrun -v -deffnm em
            grompp = settings.g_prefix + "grompp"
            mdrun = settings.g_prefix + "mdrun"
            step_no = "7"
            step_name = "Prepare files for Minimisation"
            # max warn 3 only for now
            command = grompp + " -f " + self.working_dir + "em_real.mdp -c " + self.working_dir + "solv_ions.gro -p " + self.working_dir + "topol.top -o " + self.working_dir + "em.tpr -po " + self.working_dir + "mdout.mdp -maxwarn 3 > " + self.working_dir + "step7.log 2>&1"
            self.run_process(step_no, step_name, command)

            step_no = "8"
            step_name = " Minimisation"
            # $command="g_mdrun -v -nt ".$nproc." -s ".$path.$tpr." -c ".$path.$out." -o ".$path.$trr." -x ".$path.$xtc ."-g ".$path.$logfile." > ".$path.$clog." 2>&1";
            command = mdrun + " -v  -s " + self.working_dir + "em.tpr -c " + self.working_dir + "em.gro -o " + self.working_dir + "em.trr -e " + self.working_dir + "em.edr -x " + self.working_dir + "em.xtc -g " + self.working_dir + "em.log > " + self.working_dir + "step8.log 2>&1"
            self.run_process(step_no, step_name, command)
        else:
            print "Exiting on user request "
            sys.exit()

    def nvt(self):
        pass
        # Never used

    def npt(self):
        print ">STEP11 : Initiating the Procedure to Equiliberate the System"
        print "Beginging Equiliberation with NPT Ensemble"
        grompp = settings.g_prefix + "grompp"
        mdrun = settings.g_prefix + "mdrun"
        step_no = "11"
        step_name = "Preparing files for NPT Equiliberation"
        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = grompp + "-f " + self.working_dir + "npt.mdp -c " + self.working_dir + "nvt.gro -p " + self.working_dir + "topol.top -o " + self.working_dir + "npt.tpr -po " + self.working_dir + "mdout.mdp -maxwarn 3 > " + self.working_dir + "step11.log 2>&1"
        self.run_process(step_no, step_name, command)

        step_no = "12"
        step_name = "NPT Equiliberation"
        command = mdrun + " -v  -s " + self.working_dir + "npt.tpr -c " + self.working_dir + "npt.gro -o " + self.working_dir + "npt.trr -e " + self.working_dir + "npt.edr -x " + self.working_dir + "npt.xtc -g " + self.working_dir + "npt.log > " + self.working_dir + "step12.log 2>&1"
        self.run_process(step_no, step_name, command)

    def md(self):
        print "CHEERS :) WE ARE CLOSE TO SUCCESS :):)"
        print ">STEP13 : Initiating the Production Run"
        grompp = settings.g_prefix + "grompp"
        mdrun = settings.g_prefix + "mdrun"
        step_no = "13"
        step_name = "Preparing files for NPT Equiliberation"
        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = grompp + "-f " + self.working_dir + "md.mdp -c " + self.working_dir + "npt.gro -p " + self.working_dir + "topol.top -o " + self.working_dir + "md.tpr -po " + self.working_dir + "mdout.mdp -maxwarn 3 > " + self.working_dir + "step13.log 2>&1"
        self.run_process(step_no, step_name, command)

        step_no = "14"
        step_name = "NPT Equiliberation"
        command = mdrun + " -v  -s " + self.working_dir + "md.tpr -c " + self.working_dir + "md.gro -o " + self.working_dir + "md.trr -e " + self.working_dir + "md.edr -x " + self.working_dir + "md.xtc -g " + self.working_dir + "md.log > " + self.working_dir + "step14.log 2>&1"
        self.run_process(step_no, step_name, command)


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
