function [Filename,l] = Moose_Physics_Section(mattype,Scale)


cmax = 4.8976e4;
c_ref = cmax;
ki = 1.0;

gc_prop_b = 5;
gc_prop_gb = 2;
gc_prop_pmi = 2;
gc_prop_mat = 10;
l = 0.1e-6;
visco = 2.0e-6;

Density = 4750000;
Particle_diameter = 10e-6;

In_Out_Switch=-1;
C_Rate=1;
Faraday=96485.3415;
Th_Cap=300;
Density=4750000;
R = 8.314;
T = 300;

Diffusivity_particle = [1.0e-13 0 0 0 1.0e-14 0 0 0 0];
Diffusivity_matrix = [1.0e-12 0 0 0 1.0e-12 0 0 0 0];

C1_ijkl_particle = [235.0e9 127.0e9 66.0e9 235.0e9 66.0e9 193.0e9 39.5e9 39.5e9 39.5e9];
C0_ijkl_particle = [284.0e9 70.0e9 5.8e9 284.0e9 5.8e9 12.5e9 1.6e9 1.6e9 1.6e9];

C_ijkl_matrix = [1.0962e10 7.3077e+09];

Periodic_translation_left_right = [Scale 0 0];
Periodic_translation_bottom_top = [0 Scale 0];

Eigen_vector_particle = [4.44112967954821e-07 0 0 0 -3.85150691418097e-07 0 0 0 0];
Eigen_vector_matrix = [0 0 0 0 0 0 0 0 0];

Filename='PHYSICS.i';
fileID = fopen(Filename, 'w');

%==========================================================================
%================================Variables=================================
%==========================================================================
fprintf(fileID,'[GlobalParams]\n');
fprintf(fileID,'displacements = ''disp_x disp_y''\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[Variables]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./diffused]\n');
fprintf(fileID,'order = FIRST\n');
fprintf(fileID,'family = LAGRANGE\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./damage]\n');
fprintf(fileID,'order = FIRST\n');
fprintf(fileID,'family = LAGRANGE\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./disp_x]\n');
fprintf(fileID,'order = FIRST\n');
fprintf(fileID,'family = LAGRANGE\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./disp_y]\n');
fprintf(fileID,'order = FIRST\n');
fprintf(fileID,'family = LAGRANGE\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[]\n');
fprintf(fileID,'\n');
%==========================================================================
%==========================================================================


%==========================================================================
%==========================Initial condition for diffusion=================
%==========================================================================
fprintf(fileID,'[ICs]\n');
fprintf(fileID,'[./ic_particle]\n');
fprintf(fileID,'type = ConstantIC\n');
fprintf(fileID,'variable = diffused\n');
fprintf(fileID,[strcat('value=',num2str(cmax)) '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');
%==========================================================================
%==========================================================================

%==========================================================================
%=======================Aux variables======================================
%==========================================================================
fprintf(fileID,'[AuxVariables]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./eta_gb]\n');
fprintf(fileID,'order = FIRST\n');
fprintf(fileID,'family = LAGRANGE\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./eta_pmi]\n');
fprintf(fileID,'order = FIRST\n');
fprintf(fileID,'family = LAGRANGE\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./bounds_dummy_d]\n');
fprintf(fileID,'order = FIRST\n');
fprintf(fileID,'family = LAGRANGE\n');
fprintf(fileID,'[../]\n');


fprintf(fileID,'[./sigma_h]\n');
fprintf(fileID,'family = MONOMIAL\n');
fprintf(fileID,'order = FIRST\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./von_mises]\n');
fprintf(fileID,'family = MONOMIAL\n');
fprintf(fileID,'order = FIRST\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');
%==========================================================================
%==========================================================================



%==========================================================================
%=============================Kernels======================================
%==========================================================================
fprintf(fileID,'[Kernels]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./diffusion_particle]\n');
fprintf(fileID,'type = Concentration\n'); 
fprintf(fileID,'variable = diffused\n');
fprintf(fileID,[strcat('cmax=',num2str(cmax)) '\n']);
fprintf(fileID,[strcat('ki=',num2str(ki)) '\n']);
fprintf(fileID,[strcat('block=','''',num2str(0:1:mattype-2)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./diffusion_matrix]\n');
fprintf(fileID,'type = Concentration \n');
fprintf(fileID,'variable = diffused\n');
fprintf(fileID,[strcat('cmax=',num2str(cmax)) '\n']);
fprintf(fileID,[strcat('ki=',num2str(0)) '\n']);
fprintf(fileID,[strcat('block=','''',num2str(mattype-1)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./euler]\n');
fprintf(fileID,'type = TimeDerivative\n');
fprintf(fileID,'variable = diffused\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./ACbulk]\n');
fprintf(fileID,'type = AllenCahn\n');
fprintf(fileID,'variable = damage\n');
fprintf(fileID,'f_name = F\n');
fprintf(fileID,'args = ''diffused''\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./ACInterfaceCleavageFracture]\n');
fprintf(fileID,'type = ACInterface # Anisotropic_PF_Fracture\n');
fprintf(fileID,'variable = damage\n');
fprintf(fileID,'kappa_name = kappa_op\n');
%fprintf(fileID,'beta_penalty = 0\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./dcdt]\n');
fprintf(fileID,'type = TimeDerivative\n');
fprintf(fileID,'variable = damage\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./solid_x]\n');
fprintf(fileID,'type = PhaseFieldFractureMechanicsOffDiag\n');
fprintf(fileID,'variable = disp_x\n');
fprintf(fileID,'component = 0\n');
fprintf(fileID,'c = damage\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./solid_y]\n');
fprintf(fileID,'type = PhaseFieldFractureMechanicsOffDiag\n');
fprintf(fileID,'variable = disp_y\n');
fprintf(fileID,'component = 1\n');
fprintf(fileID,'c = damage\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./off_disp]\n');
fprintf(fileID,'type = AllenCahnElasticEnergyOffDiag\n');
fprintf(fileID,'variable = damage\n');
fprintf(fileID,'displacements = ''disp_x disp_y''\n');
fprintf(fileID,'mob_name = L\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');
%==========================================================================
%==========================================================================

%==========================================================================
%========================Tensor Mechanics==================================
%==========================================================================
fprintf(fileID,'[Modules/TensorMechanics/Master]\n');
fprintf(fileID,'[all]\n');
fprintf(fileID,'add_variables = true\n');
fprintf(fileID,'strain = SMALL\n');
fprintf(fileID,'automatic_eigenstrain_names = true\n');
fprintf(fileID,'generate_output = ''vonmises_stress strain_xx strain_yy stress_xx stress_yy''\n');
fprintf(fileID,'eigenstrain_names = eigenstrain\n');
fprintf(fileID,'planar_formulation = PLANE_STRAIN\n');
%fprintf(fileID,'incremental = true\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');
%==========================================================================
%==========================================================================


%==========================================================================
%=============================User Objectss================================
%==========================================================================
fprintf(fileID,'[UserObjects]\n');
fprintf(fileID,'[./euler_angle_read]\n');
fprintf(fileID,'type = PropertyReadFile\n');
fprintf(fileID,'prop_file_name = ''Crystal_Orientation.txt''\n');
fprintf(fileID,'nprop = 3\n');
fprintf(fileID,'read_type = block\n');
fprintf(fileID,[strcat('nblock=',num2str(mattype)) '\n']);
fprintf(fileID,'[../]\n');
%fprintf(fileID,'[./cleavage_angle_read]\n');
%fprintf(fileID,'type = ElementPropertyReadFile\n');
%fprintf(fileID,'prop_file_name = ''cleavage_planes.txt''\n');
%fprintf(fileID,'nprop = 3\n');
%fprintf(fileID,'read_type = block\n');
%fprintf(fileID,'nblock= 56\n');
%fprintf(fileID,'[../]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');
%==========================================================================
%==========================================================================


fprintf(fileID,'[Materials]\n');
fprintf(fileID,'\n');

%fprintf(fileID,'[./Cleavage_planes]\n');
%fprintf(fileID,'type = Anisotropic_PF\n');
%fprintf(fileID,'read_prop_user_object = cleavage_angle_read\n');
%fprintf(fileID,'[../]\n');
%fprintf(fileID,'\n');

fprintf(fileID,'[./pfbulkmat_particle]\n');
fprintf(fileID,'type = GenericConstantMaterial\n');
fprintf(fileID,'prop_names = ''gc_prop_b gc_prop_gb gc_prop_pmi l visco''\n');
fprintf(fileID,[strcat('prop_values=','''',num2str([gc_prop_b gc_prop_gb gc_prop_pmi l visco])),'''' '\n']);
fprintf(fileID,[strcat('block=','''',num2str(0:1:mattype-2)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./pfbulkmat_matrix]\n');
fprintf(fileID,'type = GenericConstantMaterial\n');
fprintf(fileID,'prop_names = ''gc_prop l visco''\n');
fprintf(fileID,[strcat('prop_values=','''',num2str([gc_prop_mat l visco])),'''' '\n']);
fprintf(fileID,[strcat('block=','''',num2str(mattype-1)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./define_gc_prop_particle]\n');
fprintf(fileID,'type = ParsedMaterial\n');
fprintf(fileID,'material_property_names = ''gc_prop_b gc_prop_gb gc_prop_pmi''\n');
fprintf(fileID,'f_name = gc_prop\n');
fprintf(fileID,'args = ''eta_gb eta_pmi''\n');
fprintf(fileID,'function = ''(1-eta_gb)*(1-eta_gb)*(1-eta_pmi)*(1-eta_pmi)*gc_prop_b + (1-(1-eta_gb)*(1-eta_gb))*gc_prop_gb + (1-(1-eta_pmi)*(1-eta_pmi))*gc_prop_pmi''\n');
fprintf(fileID,[strcat('block=','''',num2str(0:1:mattype-2)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./diffusivity_tensor_with_Euler_particle]\n');
fprintf(fileID,'type = AnisoDiffusivity\n');
fprintf(fileID,'read_prop_user_object = euler_angle_read\n');
fprintf(fileID,[strcat('diffusivity=','''',num2str(Diffusivity_particle)),'''' '\n']);
fprintf(fileID,[strcat('block=','''',num2str(0:1:mattype-2)),'''' '\n']);
fprintf(fileID,'d = damage\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./diffusivity_tensor_with_Euler_matrix]\n');
fprintf(fileID,'type = AnisoDiffusivity\n');
fprintf(fileID,'read_prop_user_object = euler_angle_read\n');
fprintf(fileID,[strcat('diffusivity=','''',num2str(Diffusivity_matrix)),'''' '\n']);
fprintf(fileID,[strcat('block=','''',num2str(mattype-1)),'''' '\n']);
fprintf(fileID,'d = damage\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./elasticity_tensor_with_Euler_particle]\n');
fprintf(fileID,'type = Aniso_Elasticity_Tensor\n');
fprintf(fileID,'c = diffused\n');
fprintf(fileID,'cmax =  4.8976e4\n');
fprintf(fileID,[strcat('C1_ijkl=','''',num2str(C1_ijkl_particle)),'''' '\n']);
fprintf(fileID,[strcat('C0_ijkl=','''',num2str(C0_ijkl_particle)),'''' '\n']);
fprintf(fileID,'fill_method1 = symmetric9\n');
fprintf(fileID,'fill_method0 = symmetric9\n');
fprintf(fileID,'read_prop_user_object = euler_angle_read\n');
fprintf(fileID,[strcat('block=','''',num2str(0:1:mattype-2)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');
  
fprintf(fileID,'[./elasticity_tensor_with_Euler_matrix]\n');
fprintf(fileID,'type = ComputeElasticityTensor\n');
fprintf(fileID,[strcat('C_ijkl=','''',num2str(C_ijkl_matrix)),'''' '\n']);
fprintf(fileID,'fill_method = symmetric_isotropic\n');
fprintf(fileID,[strcat('block=','''',num2str(mattype-1)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./eigen_strain_prefactor_particle]\n');
fprintf(fileID,'type = DerivativeParsedMaterial \n');
fprintf(fileID,'args = diffused\n');
fprintf(fileID,'f_name = eigen_strain_prefactor\n');
fprintf(fileID,'constant_names = ''c_ref''\n');
fprintf(fileID,[strcat('constant_expressions = ','''',num2str(c_ref)),'''' '\n']);
fprintf(fileID,'function = (diffused-c_ref)\n');
fprintf(fileID,[strcat('block=','''',num2str(0:1:mattype-2)),'''' '\n']);
fprintf(fileID,'  [../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./eigen_strain_prefactor_matrix]\n');
fprintf(fileID,'type = DerivativeParsedMaterial\n');
fprintf(fileID,'args = diffused\n');
fprintf(fileID,'f_name = eigen_strain_prefactor\n');
fprintf(fileID,'constant_names = ''c_ref''\n');
fprintf(fileID,[strcat('constant_expressions = ','''',num2str(c_ref)),'''' '\n']);
fprintf(fileID,'function = (diffused-c_ref)\n');
fprintf(fileID,[strcat('block=','''',num2str(mattype-1)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./eigenstrain_particle]\n');
fprintf(fileID,'type = Aniso_Eigen_Strain\n');
fprintf(fileID,[strcat('eigen_base = ','''',num2str(Eigen_vector_particle)),'''' '\n']);
fprintf(fileID,'prefactor = eigen_strain_prefactor\n');
fprintf(fileID,'eigenstrain_name = eigenstrain\n');
fprintf(fileID,'read_prop_user_object = euler_angle_read\n');
fprintf(fileID,[strcat('block=','''',num2str(0:1:mattype-2)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./eigenstrain_matrix]\n');
fprintf(fileID,'type = ComputeEigenstrain\n');
fprintf(fileID,[strcat('eigen_base = ','''',num2str(Eigen_vector_matrix)),'''' '\n']);
fprintf(fileID,'prefactor = eigen_strain_prefactor\n');
fprintf(fileID,'eigenstrain_name = eigenstrain\n');
fprintf(fileID,[strcat('block=','''',num2str(mattype-1)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./define_mobility]\n');
fprintf(fileID,'type = ParsedMaterial\n');
fprintf(fileID,'material_property_names = ''gc_prop visco''\n');
fprintf(fileID,'f_name = L\n');
fprintf(fileID,'function = ''1.0/(gc_prop * visco)''\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./define_kappa]\n');
fprintf(fileID,'type = ParsedMaterial\n');
fprintf(fileID,'material_property_names = ''gc_prop l''\n');
fprintf(fileID,'f_name = kappa_op\n');
fprintf(fileID,'function = ''gc_prop * l''\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./damage_stress]\n');
fprintf(fileID,'type = ComputeLinearElasticPFFractureStress\n');
fprintf(fileID,'c = damage\n');
fprintf(fileID,'E_name = ''elastic_energy''\n');
fprintf(fileID,'D_name = ''degradation''\n');
fprintf(fileID,'F_name = ''local_fracture_energy''\n');
fprintf(fileID,'decomposition_type = stress_spectral\n');  
fprintf(fileID,'use_snes_vi_solver = true\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./degradation]\n');
fprintf(fileID,'type = DerivativeParsedMaterial\n');
fprintf(fileID,'f_name = degradation\n');
fprintf(fileID,'args = ''damage''\n');
fprintf(fileID,'function = ''(1.0-damage)^2*(1.0 - eta) + eta''\n');
fprintf(fileID,'constant_names       = ''eta''\n');
fprintf(fileID,'constant_expressions = ''1.0e-5''\n');
fprintf(fileID,'derivative_order = 2\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./pfbulkmat]\n');
fprintf(fileID,'type = GenericConstantMaterial\n');
fprintf(fileID,'prop_names = ''R T cmax ki''\n');
fprintf(fileID,[strcat('prop_values = ','''',num2str([R T cmax ki])),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./chemical_energy]\n');
fprintf(fileID,'type = DerivativeParsedMaterial\n');
fprintf(fileID,'f_name = chemical_energy\n');
fprintf(fileID,'args = ''diffused''\n');
fprintf(fileID,'material_property_names = ''R T cmax ki''\n');
fprintf(fileID,'function = ''R*T*cmax*((diffused/cmax)*log(1-diffused/cmax) + (1-diffused/cmax)*log(1-diffused/cmax) + ki*(diffused/cmax)*(1.0-diffused/cmax))''\n');
fprintf(fileID,'derivative_order = 2\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./local_fracture_energy]\n');
fprintf(fileID,'type = DerivativeParsedMaterial\n');
fprintf(fileID,'f_name = local_fracture_energy\n');
fprintf(fileID,'args = ''damage''\n');
fprintf(fileID,'material_property_names = ''gc_prop l''\n');
fprintf(fileID,'function = ''damage^2 * gc_prop / 2 / l''\n');
fprintf(fileID,'derivative_order = 2\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./fracture_driving_energy]\n');
fprintf(fileID,'type = DerivativeSumMaterial\n');
fprintf(fileID,'args = ''damage''\n');
fprintf(fileID,'sum_materials = ''elastic_energy local_fracture_energy''\n');
fprintf(fileID,'derivative_order = 2\n');
fprintf(fileID,'f_name = F\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[AuxKernels]\n');
fprintf(fileID,'[./gethystress]\n');
fprintf(fileID,'type = RankTwoScalarAux\n');
fprintf(fileID,'rank_two_tensor = stress\n');
fprintf(fileID,'variable = sigma_h\n');
fprintf(fileID,'scalar_type = Hydrostatic\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./vonstress]\n');
fprintf(fileID,'type = RankTwoScalarAux\n');
fprintf(fileID,'rank_two_tensor = stress\n');
fprintf(fileID,'variable = von_mises\n');
fprintf(fileID,'scalar_type = VonMisesStress\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[Functions]\n');
fprintf(fileID,'[./Particle_mass_surf_area]\n');
fprintf(fileID,'type = ParsedFunction\n');
fprintf(fileID,'value = (0.25*pi*Particle_diameter*Particle_diameter*1)*Density/(pi*Particle_diameter)  # mass (grams)/fluxarea (m^2)\n');
fprintf(fileID,'vars = ''Density Particle_diameter''\n');
fprintf(fileID,[strcat('vals=','''',num2str([Density Particle_diameter])),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./FluxVal_temp]\n');
fprintf(fileID,'type = ParsedFunction\n');
fprintf(fileID,'value = In_Out_Switch*(Th_Cap/1000)*C_Rate/Faraday\n');
fprintf(fileID,'vars = ''In_Out_Switch C_Rate Faraday Th_Cap Density''   # theoretical_capacity (Ah/g) *crate (1/h)/F(As/mol)\n');
fprintf(fileID,[strcat('vals=','''',num2str([In_Out_Switch C_Rate Faraday Th_Cap Density])),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./FluxVal]\n');
fprintf(fileID,'type = CompositeFunction\n');
fprintf(fileID,'functions = ''Particle_mass_surf_area FluxVal_temp''\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./dts]\n');
fprintf(fileID,'type = PiecewiseLinear\n');
fprintf(fileID,'x = ''0 3000''\n');
fprintf(fileID,'y = ''2 2''\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[BCs]\n');

fprintf(fileID,'\n');
fprintf(fileID,'[./flux_bc]\n');
fprintf(fileID,'type = FunctionNeumannBC\n');
fprintf(fileID,'variable = diffused\n');
fprintf(fileID,'boundary = ''left right top bottom''\n');
fprintf(fileID,'function = FluxVal\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[corner_x]\n');
fprintf(fileID,'type = DirichletBC\n');
fprintf(fileID,'variable = disp_x\n');
fprintf(fileID,'boundary = ''left_bottom left_top right_bottom right_top''\n');
fprintf(fileID,'value = 0\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[corner_y]\n');
fprintf(fileID,'type = DirichletBC\n');
fprintf(fileID,'variable = disp_y\n');
fprintf(fileID,'boundary = ''left_bottom left_top right_bottom right_top''\n');
fprintf(fileID,'value = 0\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./Periodic]\n');
fprintf(fileID,'[./bottom_top_disp_xy_periodic]\n');
fprintf(fileID,'variable = ''disp_x disp_y''\n');
fprintf(fileID,'primary =''bottom''\n');
fprintf(fileID,'secondary = ''top''\n');
fprintf(fileID,[strcat('translation=','''',num2str(Periodic_translation_bottom_top)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'[./left_right_disp_xy_periodic]\n');
fprintf(fileID,'variable = ''disp_x disp_y''\n');
fprintf(fileID,'primary =''left''\n');
fprintf(fileID,'secondary = ''right''\n');
fprintf(fileID,[strcat('translation=','''',num2str(Periodic_translation_left_right)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[Bounds]\n');
fprintf(fileID,'[./d_upper_bound]\n');
fprintf(fileID,'type = ConstantBoundsAux\n');
fprintf(fileID,'variable = bounds_dummy_d\n');
fprintf(fileID,'bounded_variable = damage\n');
fprintf(fileID,'bound_type = upper\n');
fprintf(fileID,'bound_value = 1.0\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'[./d_lower_bound]\n');
fprintf(fileID,'type = VariableOldValueBoundsAux\n');
fprintf(fileID,'variable = bounds_dummy_d\n');
fprintf(fileID,'bounded_variable = damage\n');
fprintf(fileID,'bound_type = lower\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[Preconditioning]\n');
fprintf(fileID,'[./smp]\n');
fprintf(fileID,'type = SMP\n');
fprintf(fileID,'full = true\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[Postprocessors]\n');
fprintf(fileID,'[./avg_surface_concentration]\n');
fprintf(fileID,'type = SideAverageValue\n');
fprintf(fileID,'variable = diffused\n');
fprintf(fileID,'boundary = ''PMI''\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./ave_damage_particle]\n');
fprintf(fileID,'type = ElementAverageValue  # Volumetric_average\n');
fprintf(fileID,'variable = damage\n');
fprintf(fileID,[strcat('block=','''',num2str(0:1:mattype-2)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./ave_damage_matrix]\n');
fprintf(fileID,'type = ElementAverageValue  # Volumetric_average\n');
fprintf(fileID,'variable = damage\n');
fprintf(fileID,[strcat('block=','''',num2str(mattype-1)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./ave_damage_all]\n');
fprintf(fileID,'type = ElementAverageValue  # Volumetric_average\n');
fprintf(fileID,'variable = damage\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./ave_conc_particle]\n');
fprintf(fileID,'type = ElementAverageValue  # Volumetric_average\n');
fprintf(fileID,'variable = diffused\n');
fprintf(fileID,[strcat('block=','''',num2str(0:1:mattype-2)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./ave_conc_matrix]\n');
fprintf(fileID,'type = ElementAverageValue  # Volumetric_average\n');
fprintf(fileID,'variable = diffused\n');
fprintf(fileID,[strcat('block=','''',num2str(mattype-1)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./ave_conc_all]\n');
fprintf(fileID,'type = ElementAverageValue  # Volumetric_average\n');
fprintf(fileID,'variable = diffused\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./max_damage_value_particle]\n');
fprintf(fileID,'type = NodalExtremeValue\n');
fprintf(fileID,'variable = damage\n');
fprintf(fileID,[strcat('block=','''',num2str(0:1:mattype-2)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./max_damage_value_matrix]\n');
fprintf(fileID,'type = NodalExtremeValue\n');
fprintf(fileID,'variable = damage\n');
fprintf(fileID,[strcat('block=','''',num2str(mattype-1)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./max_damage_value_all]\n');
fprintf(fileID,'type = NodalExtremeValue\n');
fprintf(fileID,'variable = damage\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./avg_pmi_damage]\n');
fprintf(fileID,'type = SideAverageValue\n');
fprintf(fileID,'variable = damage\n');
fprintf(fileID,'boundary = ''PMI''\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');


fprintf(fileID,'[./ave_vonmises_particle]\n');
fprintf(fileID,'type = ElementAverageValue  # Volumetric_average\n');
fprintf(fileID,'variable = von_mises\n');
fprintf(fileID,[strcat('block=','''',num2str(0:1:mattype-2)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./ave_von_mises_matrix]\n');
fprintf(fileID,'type = ElementAverageValue  # Volumetric_average\n');
fprintf(fileID,'variable = von_mises\n');
fprintf(fileID,[strcat('block=','''',num2str(mattype-1)),'''' '\n']);
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./ave_von_mises_all]\n');
fprintf(fileID,'type = ElementAverageValue  # Volumetric_average\n');
fprintf(fileID,'variable = von_mises\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[Executioner]\n');
fprintf(fileID,'type = Transient\n');
fprintf(fileID,'scheme = bdf2\n');
fprintf(fileID,'\n');

fprintf(fileID,'solve_type = PJFNK\n');
fprintf(fileID,'petsc_options_iname = ''-pc_type -pc_factor_mat_solver_package -snes_type''\n');
fprintf(fileID,'petsc_options_value = ''lu       superlu_dist                  vinewtonrsls''\n');
fprintf(fileID,'automatic_scaling = true\n');
fprintf(fileID,'\n');

fprintf(fileID,'l_max_its = 100\n');
fprintf(fileID,'nl_max_its = 30\n');
fprintf(fileID,'nl_abs_tol = 5e-3\n');
fprintf(fileID,'l_tol = 1e-04\n');
fprintf(fileID,'num_steps = 1500\n');
fprintf(fileID,'\n');
%fprintf(fileID,'dt = 2\n');

fprintf(fileID,'[./TimeStepper]\n');
fprintf(fileID,'type = FunctionDT\n');
fprintf(fileID,'function = dts\n');
fprintf(fileID,'min_dt = 0.5\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[Outputs]\n');
fprintf(fileID,'execute_on = ''initial timestep_end''\n');
fprintf(fileID,'interval = 5\n');
fprintf(fileID,'exodus = true\n');
fprintf(fileID,'csv = true\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./my_checkpoint]\n');
fprintf(fileID,'type = Checkpoint\n');
fprintf(fileID,'num_files = 5\n');
fprintf(fileID,'interval = 25\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[MultiApps]\n');
fprintf(fileID,'[full_solve]\n');
fprintf(fileID,'type = FullSolveMultiApp\n');
fprintf(fileID,'execute_on = initial\n');
fprintf(fileID,'positions = ''0 0 0''\n');
fprintf(fileID,'input_files = Diffused_Interface_Sub_APP.i\n');
fprintf(fileID,'clone_master_mesh = true\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'[]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[Transfers]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./from_eta_gb]\n');
fprintf(fileID,'type = MultiAppCopyTransfer\n');
fprintf(fileID,'source_variable = eta_gb\n');
fprintf(fileID,'variable = eta_gb\n');
fprintf(fileID,'from_multi_app = full_solve\n');
fprintf(fileID,'[../]\n');
fprintf(fileID,'\n');

fprintf(fileID,'[./from_eta_pmi]\n');
fprintf(fileID,'type = MultiAppCopyTransfer\n');
fprintf(fileID,'source_variable = eta_pmi\n');
fprintf(fileID,'variable = eta_pmi\n');
fprintf(fileID,'from_multi_app = full_solve\n');
fprintf(fileID,'[../]\n');

fprintf(fileID,'\n');
fprintf(fileID,'[]\n');

fprintf(fileID,'\n');

