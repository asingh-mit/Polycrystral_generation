function [Filename] = Sub_app_Section(numat,lc_frac)

Filename='Diffused_Interface_Sub_APP.i';
fileID = fopen(Filename, 'w');

fprintf(fileID,'[Mesh]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[Variables]\n');
fprintf(fileID,'[./eta_pmi]\n');
fprintf(fileID,'initial_condition = 0.0\n');
fprintf(fileID,'order = FIRST\n');
fprintf(fileID,'family = LAGRANGE\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'[./eta_gb]\n');
fprintf(fileID,'initial_condition = 0.0\n');
fprintf(fileID,'order = FIRST\n');
fprintf(fileID,'family = LAGRANGE\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[Kernels]\n');
fprintf(fileID,'[./PhaseField_eta_pmi]\n');
fprintf(fileID,'type = Phasefield_eta_pmi\n');
fprintf(fileID,'variable = eta_pmi\n');
fprintf(fileID,[strcat('lc_eta_pmi =',num2str(lc_frac/2)) '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');
fprintf(fileID,'\n');
fprintf(fileID,'[./PhaseField_eta_gb]\n');
fprintf(fileID,'type = Phasefield_eta_gb\n');
fprintf(fileID,'variable = eta_gb\n');
fprintf(fileID,[strcat('lc_eta_gb =',num2str(lc_frac/2)) '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[BCs]\n');
fprintf(fileID,'[./eta_pmi]\n');
fprintf(fileID,'type = ADDirichletBC\n');
fprintf(fileID,'variable = eta_pmi\n');
fprintf(fileID,'boundary = ''PMI'' \n');
fprintf(fileID,'value = 1.0  \n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./eta_gb]\n');
fprintf(fileID,'type = ADDirichletBC\n');
fprintf(fileID,'variable = eta_gb\n');
fprintf(fileID,'boundary = ''');
for i = 1:1:numat-1
    fprintf(fileID,strcat('GB',[num2str(i)]));
    fprintf(fileID,' ');
end
fprintf(fileID,'''\n');
fprintf(fileID,'value = 1.0  \n');
fprintf(fileID,'[]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[Executioner]\n');
fprintf(fileID,'type = Steady\n');
fprintf(fileID,'automatic_scaling = true\n');
fprintf(fileID,'compute_scaling_once=false\n');
fprintf(fileID,'solve_type = ''PJFNK''\n');
fprintf(fileID,'petsc_options_iname = ''-pc_type -pc_hypre_type -ksp_gmres_restart''\n');
fprintf(fileID,'petsc_options_value = ''hypre boomeramg 101''\n');
fprintf(fileID,'l_max_its = 50\n');
fprintf(fileID,'nl_max_its = 10\n');
fprintf(fileID,'nl_abs_tol = 1e-9\n');
fprintf(fileID,'l_tol = 1e-04\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[Outputs]\n');
fprintf(fileID,'exodus = true\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');