# Developer Notes




## 1. Introduction
- This is a library/module which acts as a wrapper around the commandline tools of [Gromacs](http://www.gromacs.org).
- The main purpose of this library is to automate the tasks with least or no programming experience.
- The library makes the usage of the tool gromacs tool more interactive for the beginers. 
- This is one of the modules of the drug discovery tools automation pipeline [rsquarelabs-core](https://github.com/rsquarelabs/rsquarelabs-core) by [rsquarelabs.org](http://rsquarelabs.org).  
- The library contains the presets of the protocols based on the gromacs tutorials by written [Justin A. Lemkul, Ph.D](http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/).
- The library is programmable - well documented enough to write custom protocols, with no or least experience.  


## 2. Features
1. Automate the tasks
2. Presets of most used protocols
3. Protocol Customisation ()
4. Notification by email once the task is done.

## 2.Protocols (PT)
1. Protein Minimisation(P-M)
2. Protein Ligand Minimisation(PL-M)
3. Protein Molecular Dynamics(P-MD) - Continuation of P-M
4. Protein Ligand Molecular Dynamics(PL-MD) - Continuation of PL-M


## 3. steps and step_names

All the output files of steps 5, 6, 7 should have a namespace that gets tagged to the files creates ie., if user provide namespace "neutralise_system_", we will tag "em.log" with "neutralise_system_" and create "neutralise_system_em.log". 
Every method should have a namespace/name_of_step which will be tagged to the files.


config_file - Basic Project file (similar to npm init ie., when the script is started it will ask the user for basic questions(prompt) like "ProjectName, email, ProteinName" and saves the info to project.config  ) - project.config
state_file - Each step when started and after completed should write the status to "status_file.log"
log_file - every print statement of the current should be converted into logging ie., we should show logging on console and file. "project.log"


### 3.1 main steps
1. create_protein_topology - 
``` pdb2gmx -f  protein.pdb -o protein.gro -ignh -p topol.top -i posre.itp -ff gromos53a6 -water spc ```

2. combine_protein_ligand_topology - Step needed only for Protein Ligand Protocols ie., PT-M, PL-MD 
```
    a. system.gro = ligand.gro+protein.gro
    b. update topol.top     - add LIGAND topology (#include ligand.itp)
    c. update topol.top     - add water topology 
    d. update topol.top     - add ligand identifier at the end 
    e. update system.gro    - update system_coordinates_count
```
3. create_box - creates a box for simulation to happen
``` editconf   -f  system.gro -o  newbox.gro -bt cubic -d 1 -c ```

4. solvate_box - fills the box 
``` genbox -cp  newbox.gro -p  topol.top -cs spc216.gro -o  solv.gro ```

5. generate_mdp_file - writes a template md.mdp or em.mdp or em_real.mdp based on the protocol (P-M needs em.mdp in one step and em_real.mdp in next step, P-MD need ). We should be saving these em.mdp, em_real.mdp, md_mdp in a folder templates/

6. neutralise_system(add_ions) -
    a. Check the current charge of the system  - ```grompp -f  em.mdp -c solv.gro -p topol.top -o ions.tpr -po mdout.mdp  ```
    b. add the counter charges and neutralise the system
    b.1 adding positive charges - ``` genion -s ions.tpr -o solv_ions.gro -p topol.top -nname CL -nn str(charge_count)  -g adding_charges.log  ```
    b.2 adding negative charges - ``` genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -np str(-charge_count) -g adding_charges.log ``` 

7. run_dynamics
   a. prepare_for_md - ```grompp -f em_real.mdp -c solv_ions.gro -p topol.top -o em.tpr -po mdout.mdp -maxwarn 3 ```
   b. actual_md_run - ``` mdrun  -v  -s  em.tpr -c  em.gro -o em.trr -e em.edr -x em.xtc -g em.log ```


### 3.2 utils
1. check_gro_file - check if the file has *.gro extension (the next version will be checking the format of the internal content)
2. check_pdb_file - check if the file has *.pdb extension (the next version will be checking the format of the internal content)
3. copy_file_to_project - copies the files from the external source to ~/rsquarelabs/[project_name]
3. write_report - logger function which writes the report to report.log
4. write_status - logger function which writes the status of the each step starting and ending 
4. notify_with_email - notifies the user by email after the task is done. you get email from project.config




## 4 Code File Structure 
gromacs/__init__.py
gromacs/utils.py
