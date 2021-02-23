from .settings import version, weblink

welcome_message = """
        ********************************************
        [  RSQUARE-LABS GROMACS AUTOMATION SCRIPT  ]
        ********************************************
        Version """ + version + """
        Check out """ + weblink + """ for more information.

          ****CHEERS & HAPPY RESEARCH****

        This program is free software; you can redistribute it and/or modify it
        under the terms of the GNU General Public License version 2 (GNU GPL v2)
        Do not remove the original authors name & credits in redistributions or
        modified version.
        """
backup_folder_already_exists = """
    POSSIBLE SOLUTION: A backup folder already present.
    You may loose your work if we replace , so change your
    Working Directory or delete the folder BACKUP
"""


general_nvt_data = """
; VARIOUS PREPROCESSING OPTIONS
;title                    = NVT simulation (constant number, volume and temperature)
cpp                      = /lib/cpp

; RUN CONTROL PARAMETERS
integrator               = md
tinit                    = 0
dt                       = 0.002
nsteps                   = 1000 ; 2 ps
nstcomm                  = 0

; OUTPUT CONTROL OPTIONS
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 1
nstenergy                = 1
nstxtcout                = 0
xtc_precision            = 1000
xtc-grps                 = System
energygrps               = Protein Non-Protein

; NEIGHBORSEARCHING PARAMETERS
nstlist                  = 5
ns-type                  = Grid
pbc                      = xyz
rlist                    = 0.9

; OPTIONS FOR ELECTROSTATICS AND VDW
coulombtype              = Reaction-Field
rcoulomb                 = 1.4
epsilon_rf               = 78
epsilon_r                = 1
vdw-type                 = Cut-off
rvdw                     = 1.4

; Temperature coupling  
Tcoupl                   = Berendsen
tc-grps                  = Protein  Non-Protein
tau_t                    = 0.1      0.1
ref_t                    = 300      300
; Pressure coupling     
Pcoupl                   = No

; GENERATE VELOCITIES FOR STARTUP RUN
gen_vel                  = yes
gen_temp                 = 300.0
gen_seed                 = 173529

; OPTIONS FOR BONDS    
constraints              = all-bonds
constraint-algorithm     = Lincs
unconstrained-start      = no
lincs-order              = 4
lincs-iter               = 1
lincs-warnangle          = 30"""

 # SPECIFIC NPT DATA!! change ASAP.
general_npt_data = """
;title		= OPLS Lysozyme NPT equilibration 
include = -I..
define		= -DPOSRES	; position restrain the protein
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= 50000		; 2 * 50000 = 100 ps
dt		    = 0.002		; 2 fs
; Output control
nstxout		= 500		; save coordinates every 1.0 ps
nstvout		= 500		; save velocities every 1.0 ps
nstenergy	= 500		; save energies every 1.0 ps
nstlog		= 500		; update log file every 1.0 ps
; Bond parameters
continuation	        = yes		; Restarting after NVT 
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		    ; cubic interpolation
fourierspacing	= 0.16		; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1	        ; time constant, in ps
ref_t		= 300 	  300	        ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT
pcoupltype	        = isotropic	            ; uniform scaling of box vectors
tau_p		        = 2.0		            ; time constant, in ps
ref_p		        = 1.0		            ; reference pressure, in bar
compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com
; Periodic boundary conditions
pbc		= xyz		; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= no		; Velocity generation is off """



# SPECIFIC MD PARAMETERS! CHANGE ASAP.
general_md_data = """
;title		= OPLS Lysozyme MD simulation 
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= 500000	; 2 * 500000 = 1000 ps (1 ns)
dt		    = 0.002		; 2 fs
; Output control
nstxout		        = 5000		; save coordinates every 10.0 ps
nstvout		        = 5000		; save velocities every 10.0 ps
nstenergy	        = 5000		; save energies every 10.0 ps
nstlog		        = 5000		; update log file every 10.0 ps
nstxout-compressed  = 5000      ; save compressed coordinates every 10.0 ps
                                ; nstxout-compressed replaces nstxtcout
compressed-x-grps   = System    ; replaces xtc-grps
; Bond parameters
continuation	        = yes		; Restarting after NPT 
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		    ; cubic interpolation
fourierspacing	= 0.16		; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1	        ; time constant, in ps
ref_t		= 300 	  300	        ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT
pcoupltype	        = isotropic	            ; uniform scaling of box vectors
tau_p		        = 2.0		            ; time constant, in ps
ref_p		        = 1.0		            ; reference pressure, in bar
compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc		= xyz		; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= no		; Velocity generation is off 
"""
