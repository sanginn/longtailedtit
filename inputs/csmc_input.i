[Mesh]
  type = FileMesh # Can generate simple lines, rectangles and rectangular prisms
  file = mesh/mesh-q9.msh
  boundary_id = '1 2 3 4 5 6 7 8 9' # Assign names to boundaries to make things clearer
  boundary_name = 'b1 b2 b3 b4 b5 b6 b7 b8 b9'
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
    save_in_disp_x = force_x
    save_in_disp_y = force_y
    use_displaced_mesh = false
  [../]
  [./gravity_y]
    type = Gravity
    variable = disp_y
    value = -9.81
  [../]
[]

[BCs]
  [./bottom_x]
    type = PresetBC
    variable = disp_x
    boundary = 'b1 b2'
    value = 0
  [../]
  [./bottom_y]
    type = PresetBC
    variable = disp_y
    boundary = 'b1 b2'
    value = 0
  [../]
  [./right_x]
    type = PresetBC
    variable = disp_x
    boundary = 'b3 b4'
    value = 0
  [../]
  [./left_x]
    type = PresetBC
    variable = disp_x
    boundary = 'b7 b8'
    value = 0
  [../]
  [./push_x]
    type = PresetBC
    variable = disp_x
    boundary = b6
    value = 0
  [../]
  [./push_y]
    type = FunctionPresetBC
    variable = disp_y
    boundary = b6
    function = disp_footing
  [../]
[]

[AuxVariables]
  [./force_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./force_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./stress_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_mises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pstrain_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pstrain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pstrain_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pstrain_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./mc_int]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./yield_fcn]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./void_ratio]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dilatancy_angle]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./stress_x]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_x
    index_i = 0
    index_j = 0
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  [../]
  [./stress_y]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_y
    index_i = 1
    index_j = 1
  [../]
  [./stress_z]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_z
    index_i = 2
    index_j = 2
  [../]
  [./stress_mises_kernel]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = stress_mises
    scalar_type = VonMisesStress
  [../]
  [./pstrain_x]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = pstrain_x
    index_i = 0
    index_j = 0
  [../]
  [./pstrain_xy]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = pstrain_xy
    index_i = 0
    index_j = 1
  [../]
  [./pstrain_y]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = pstrain_y
    index_i = 1
    index_j = 1
  [../]
  [./pstrain_z]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = pstrain_z
    index_i = 2
    index_j = 2
  [../]
  [./mc_int_auxk]
    type = MaterialStdVectorAux
    index = 0
    property = plastic_internal_parameter
    variable = mc_int
  [../]
  [./yield_fcn_auxk]
    type = MaterialStdVectorAux
    index = 0
    property = plastic_yield_function
    variable = yield_fcn
  [../]
  [./void_ratio_auxk]
    type = MaterialRealAux
    index = 0
    property = mc_void_ratio
    variable = void_ratio
  [../]
  [./dilatancy_angle_auxk]
    type = MaterialRealAux
    index = 0
    property = mc_dilatancy_angle
    variable = dilatancy_angle
  [../]
[]

[Functions]
  [./sig_ini_v]
    type = ParsedFunction
    value = 15*y # initial stress that should result from the weight force
  [../]
  [./sig_ini_h]
    type = ParsedFunction
    value = 7.5*y # some arbitrary xx and yy stress that should not affect the result
  [../]
  [./disp_footing]
    type = ParsedFunction
    value = -0.01*t
  [../]
[]

[Materials]
  [./mc]
    type = FiniteStrainCriticalStateMohrCoulomb
    block = 0
    initial_void_ratio = 0.7331
    gamma = 0.934
    lambda = 0.019
    xi = 0.7
    p_ref = 100
    psi_0 = 2.7
    fixed_elastic_moduli = false
    c_g = 125
    poisson_ratio = 0.25
    disp_x = disp_x
    disp_y = disp_y
    fill_method = symmetric_isotropic
    C_ijkl = '3.2E4 3.2E4'
    mc_cohesion = 0
    mc_friction_angle = 31.5
    mc_dilation_angle = 3.15
    mc_tip_smoother = 0
    mc_edge_smoother = 25
    yield_function_tolerance = 1E-8
    ep_plastic_tolerance = 1
    internal_constraint_tolerance = 1
    initial_stress = 'sig_ini_h 0 0  0 sig_ini_v 0  0 0 sig_ini_h'
  [../]
  [./density]
    type = GenericConstantMaterial
    block = 0
    prop_names = density
    prop_values = 1.52905199 # in ton/m^3
  [../]
[]

[Postprocessors]
  [./react_y_bottom]
    type = NodalSum
    variable = force_y
    boundary = b6
  [../]
[]

[Preconditioning]
  [./SMP]
    # Creates the entire Jacobian, for the Newton solve
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  end_time = 8.0
  solve_type = LINEAR
  nl_abs_tol = 1e-3
  dt = 0.0002
  dtmax = 0.0002
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 0.0002
  [../]
  [./TimeIntegrator]
    type = ExplicitMidpoint
  [../]
  [./Quadrature]
    type = GAUSS
    order = FIRST
  [../]
[]

[Outputs]
  interval = 5
  console = true
  exodus = true
  tecplot = true
  csv = true
[]
