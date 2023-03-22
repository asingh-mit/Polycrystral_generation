[Mesh]
[]

[Variables]
[./eta_pmi]
initial_condition = 0.0
order = FIRST
family = LAGRANGE
[../]
[./eta_gb]
initial_condition = 0.0
order = FIRST
family = LAGRANGE
[../]
[]

[Kernels]
[./PhaseField_eta_pmi]
type = Phasefield_eta_pmi
variable = eta_pmi
lc_eta_pmi =5e-08
[../]


[./PhaseField_eta_gb]
type = Phasefield_eta_gb
variable = eta_gb
lc_eta_gb =5e-08
[../]
[]

[BCs]
[./eta_pmi]
type = ADDirichletBC
variable = eta_pmi
boundary = 'PMI' 
value = 1.0  
[]

[./eta_gb]
type = ADDirichletBC
variable = eta_gb
boundary = 'GB1 GB2 GB3 GB4 GB5 GB6 GB7 GB8 GB9 GB10 GB11 GB12 GB13 GB14 GB15 GB16 GB17 GB18 GB19 GB20 GB21 GB22 GB23 GB24 GB25 GB26 GB27 GB28 GB29 GB30 GB31 GB32 GB33 GB34 GB35 GB36 GB37 GB38 GB39 GB40 GB41 GB42 GB43 GB44 GB45 GB46 GB47 GB48 GB49 GB50 GB51 GB52 GB53 GB54 GB55 GB56 GB57 GB58 GB59 GB60 GB61 GB62 GB63 GB64 GB65 GB66 GB67 GB68 GB69 GB70 GB71 GB72 GB73 GB74 GB75 GB76 GB77 GB78 GB79 GB80 GB81 GB82 GB83 '
value = 1.0  
[]
[]

[Executioner]
type = Steady
automatic_scaling = true
compute_scaling_once=false
solve_type = 'PJFNK'
petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
petsc_options_value = 'hypre boomeramg 101'
l_max_its = 50
nl_max_its = 10
nl_abs_tol = 1e-9
l_tol = 1e-04
[]

[Outputs]
exodus = true
[]

