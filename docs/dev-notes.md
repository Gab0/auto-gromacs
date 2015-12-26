# Developer Notes




## 1.Protocols (PT)
1. Protein Minimisation(P-M)
2. Protein Ligand Minimisation(PL-M)
3. Protein Molecular Dynamics(P-MD) - Continuation of P-M
4. Protein Ligand Molecular Dynamics(PL-MD) - Continuation of PL-M





## 2. steps and step_names

All the output files of steps 5, 6, 7 should have a namespace that gets tagged to the files creates ie., if user provide namspace "neutralise_system_", we will tag "em.log" with "neutralise_system_" and create "neutralise_system_em.log". 
Every method should have a namespace/name_of_step which will be tagged to the files.


config_file - Basic Project file (similar to npm init ie., when the script is started it will ask the user for basic questions(prompt) like "ProjectName, email, ProteinName" and saves the info to project.config  )
state_file - Each step when started and after completed should write the status to "status_file.log"
log_file - every print statement of the current should be converted into logging ie., we should show logging on console and file. "project.log"


### 2.1 main steps
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


### 2.2 utils
1. read_gro_file
2. read_pdb_file
3. write_report
4. write_status
4. notify_with_email
