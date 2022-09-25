Include "cable_data_group16.geo";

DefineConstant[
  Flag_AnalysisType = {0,
    Choices{
      0="Electric",
      1="Magnetic",
      2="Magneto-thermal (linear)",
      3="Magneto-thermal (nonlinear)"
    },
    Name "{00Parameters/00Type of analysis", Highlight "Blue",
    ServerAction Str["Reset","GetDP/1ResolutionChoices"]
  }

  Flag_sigma_funcT = (Flag_AnalysisType==3)?1:0, // Only useful in the magneto-thermal coupled case

  nb_iter = 20, // Maximum number of nonlinear iterations (You may adapt)
  relaxation_factor = 1, // value in [0,1]; if 1, there is no relaxation; if <1, you used the solution of previous iteration for helping convergence
  stop_criterion = 1e-6, // prescribed tolerance, iterative process stops when the difference between two consecutive iterations is smaller than this value

  // You can predefine the default Resolution, ComputeCommand and Operation (-solve, -pos)
  // In these files:
  // Resolution depends on Flag_AnalysisType => the Default depends on the formulation file that is included
  // PostOperations are called in the Resolution => No need to indicate anything else
  r_ = {"", Name "GetDP/1ResolutionChoices", Visible 1}
  c_ = {"-solve -v2", Name "GetDP/9ComputeCommand", Visible 1},
  p_ = {"", Name "GetDP/2PostOperationChoices", Visible 1}
];
Group {

  semiconductorin= Region[{SEMI_IN}];
  xlpe= Region[{ XLPE }]; // The defect becomes Air
  defect = Region[{ DEFECT }];
  semiconductorout= Region[{SEMI_OUT}];
  tpsoutpolypropilene= Region[{TPS_OUT_PE}];
  tpsoutaluminium= Region[{TPS_OUT_AL}];

  polypropilene= Region[{POLYPROPILENE }];
  aluminiumarmourin= Region[{ALUMINIUM_ARMOUR_IN}];
  polyethylenesheatout= Region[{POLYETHYLENE_SHEAT_OUT}];

  water_EM= Region[{WATER_EM}];
  water_TH= Region[{WATER_TH}];
  soil_EM= Region[{SOIL_EM}];
  soil_TH= Region[{SOIL_TH}];
  For k In {1:3}
    Ind~{k} = Region[{(WIRE+k-1) }];
    Inds += Region[{(WIRE+k-1) }];
  EndFor

  Ind_1 = Region[{1000}];
  Ind_2 = Region[{1001}];
  Ind_3 = Region[{1002}];

  Cable = Region[{Ind_1, Ind_2, Ind_3,semiconductorin,xlpe,semiconductorout,tpsoutpolypropilene,tpsoutaluminium,polypropilene,aluminiumarmourin,polyethylenesheatout,defect}];
  //Cable += Region [{}]; in case we use a defect


  SurfaceGe0 = Region[{OUTBND_WATER_EM,OUTBNB_SOIL_EM,INTERFACE_WATER_SOIL}];

  Sur_Dirichlet_Ele = Region[{SurfaceGe0}]; // outer boundary for the electrodynamic computation
  Sur_Dirichlet_Mag = Region[{SurfaceGe0}]; // outer boundary for the magnetic computation

  Domain_Ele = Region[{Cable, water_EM, soil_EM}]; // total domain for electrodynamic computation

  DomainCC_Mag = Region[{Ind_1, Ind_2, Ind_3}]; // The rest of the cable.
  DomainCC_Mag += Region[{xlpe,tpsoutpolypropilene,polypropilene,defect,polyethylenesheatout}];
  DomainC_Mag = Region[{tpsoutaluminium,aluminiumarmourin,semiconductorin,semiconductorout}]; // Metallic parts (Aluminum)

  DomainS0_Mag = Region[{}]; // Not used because we are using Current_2D constraints
  DomainS_Mag = Region[{Ind_1, Ind_2, Ind_3}]; // The actual conductors using Current_2D constraints

  DomainCWithI_Mag = Region[{}];
  Domain_Mag = Region[ {DomainCC_Mag, DomainC_Mag} ];


  // Thermal domain
  Vol_Thermal          = Region[{Cable,water_EM, soil_EM,water_TH, soil_TH}]; // The cable
  Vol_QSource_Thermal  = Region[{DomainCC_Mag}]; //The actual conductors using Current_2D constraints
  Vol_QSource0_Thermal = Region[{DomainS0_Mag}]; // Not used, as S0_Mag is not necessary for this model
  Vol_QSourceB_Thermal = Region[{}]; // The actual conductors using Current_2D constraints
  Sur_Convection_Thermal =  Region[{water_TH, soil_TH,water_EM, soil_EM}]; // The surroundings
  Sur_Dirichlet_Thermal  = Region[{water_TH, soil_TH}]; // The boundary conditionm in this case, the region around the cable.
  //Sur_Dirichlet_Thermal  = Region[{OUTBND_WATER_TH, OUTBNB_SOIL_TH}]; // The boundary conditionm in this case, the region around the cable.
  Domain_Thermal = Region[{Vol_Thermal}];



DomainDummy = Region[{12345}];

}
Function {
  mu0 = 4.e-7 * Pi;
    eps0 = 8.854187818e-12;

    // TO DEFINE FOR ALL MATERIALS
    // nu[], sigma[], epsilon[]... are piecewise defined functions
    // - if no Region is indicated, the program assumes the same value all over: e.g. nu[] = 1/mu0;
    // - if we have only one region for a particular value, you can put write: nu[MyRegion]=  1/(mu0*mur_steel);
    // - if you have more than one region with the same characteristic: nu[Region[{MyRegion1, MyRegion2}]] = 1/(mu0*mur_steel);
    // ATTENTION: You can't define the function twice for a Physical Region, you will get a runtime error
    epsilon[Region[{1000,1001,1002}]]  = eps0; // Copper cables
    epsilon[Region[{SEMI_IN}]]  = eps0*espr_semiconductor; // Semiconductor cable in
    epsilon[Region[{XLPE}]]  = eps0*epsr_xlpe; // Isolating material
    epsilon[Region[{SEMI_OUT}]]  = eps0*espr_semiconductor; // semiconductor Cable out
    epsilon[Region[{TPS_OUT_PE}]]  = eps0*epsr_polyethylene; //Polyethilene sheath Bundle
    epsilon[Region[{TPS_OUT_AL}]]  = eps0; // Alumiepsilonm shield
    epsilon[Region[{POLYPROPILENE}]]  = eps0*epsr_polypropilene;
    epsilon[Region[{ALUMINIUM_ARMOUR_IN}]]  = eps0;
    epsilon[Region[{POLYETHYLENE_SHEAT_OUT}]]  = eps0*epsr_polyethylene;
    epsilon[Region[{WATER_EM}]]  = eps0*epsr_water;
    epsilon[Region[{WATER_TH}]]  = eps0*epsr_water;
    epsilon[Region[{SOIL_EM}]]  = eps0*epsr_soil;
    epsilon[Region[{SOIL_TH}]]  = eps0*epsr_soil;

  //Defect Characterization

    If(Flag_defect_in_XLPE)
      epsilon[Region[{DEFECT}]]  = eps0; // The defect becomes Air
    Else
      epsilon[Region[{DEFECT}]]  = eps0*epsr_xlpe; // The defect becomes XLPE: No defect
    EndIf

    //conductivity sigma
    //sigma[Region[{1000,1001,1002}]]  = sigma_cu; // Copper cables
    sigma[Region[{SEMI_IN}]]  = sigma_semiconductor; // Semiconductor cable in
    sigma[Region[{XLPE}]]  = sigma_xlpe; // Isolating material
    sigma[Region[{SEMI_OUT}]]  = sigma_semiconductor; // semiconductor Cable out
    sigma[Region[{TPS_OUT_PE}]]  = sigma_polyethylene; //Polyethilene sheath Bundle
    //sigma[Region[{TPS_OUT_AL}]]  = sigma_al; // Aluminum shield
    sigma[Region[{POLYPROPILENE}]]  = sigma_polypropilene;
    //sigma[Region[{ALUMINIUM_ARMOUR_IN}]]  = sigma_al;
    sigma[Region[{POLYETHYLENE_SHEAT_OUT}]]  = sigma_polyethylene;
    sigma[Region[{WATER_EM}]]  = sigma_water;
    sigma[Region[{WATER_TH}]]  = sigma_water;
    sigma[Region[{SOIL_EM}]]  = sigma_soil_sea;
    sigma[Region[{SOIL_TH}]]  = sigma_soil_sea;

    //Defect Characterization

    If(Flag_defect_in_XLPE)
      sigma[Region[{DEFECT}]]  = sigma_air; // The defect becomes Air
    Else
      sigma[Region[{DEFECT}]]  = sigma_xlpe; // The defect becomes XLPE: No defect
    EndIf
    // Examples of nonlinear functions for the sigma dependence with the temperature
    // Attention alpha_cu, alpha_al, Tref have to be defined somewhere before this
    // $1 is the argument, you have to provide it when calling it; in this case the temperature you compute {T} (see magneto-thermal formulation)

    alpha_al = 0.00390; // K^-1 for Aluminum
    alpha_cu = 0.00386; // K^-1 for copper
    Tref = 276.15; // 3degC Average Temperature Underwater in the north sea

    fT_cu[] = (1+alpha_cu*($1-Tref)); // $1 is current temperature in [K], alpha in [1/K]
    fT_al[] = (1+alpha_al*($1-Tref));


    If (!Flag_sigma_funcT)
      sigma[Region[{1000,1001,1002}]] = sigma_cu;
      sigma[Region[{TPS_OUT_AL,ALUMINIUM_ARMOUR_IN}]] = sigma_al;
    Else
      sigma[Region[{1000,1001,1002}]] = sigma_cu/fT_cu[$1];
      sigma[Region[{TPS_OUT_AL,ALUMINIUM_ARMOUR_IN}]] = sigma_al/fT_al[$1];

    EndIf

    Freq = 50; // Adapt if needed
    Omega = 2*Pi*Freq;

    // Example for a three phase system
    Pa = 0.; Pb = -120./180.*Pi; Pc = -240./180.*Pi;
    I = 406; // maximum value current in data sheet
    DefineFunction[js0]; // Only needed if source in DomainS0_Mag
    // Example for a three-phase cable ***Not used in our case***
    //js0[Ind_1] = Vector[0,0,1] * I / SurfaceArea[] * F_Cos_wt_p[]{Omega, Pa};
    //js0[Ind_2] = Vector[0,0,1] * I / SurfaceArea[] * F_Cos_wt_p[]{Omega, Pb};
    //js0[Ind_3] = Vector[0,0,1] * I / SurfaceArea[] * F_Cos_wt_p[]{Omega, Pc};

    Ns[]= 1; // Number of strands in your source domain
    Sc[]= SurfaceArea[]; // Section area of your source domain

    // Using second order hierarchical basis functions if set to 2
    Flag_Degree_a = 1;
    Flag_Degree_v = 1;

    // thermal parameters
    Tamb = Tref; // Average temperature underwater.
    Tambient[] = Tamb; // [K]

    // thermal conductivities [W/(m K)]
    // also piecewise defined
    k[Region[{1000,1001,1002}]]  = kappa_cu; // Copper cables
    k[Region[{SEMI_IN}]]  = kappa_semiconductor; // Semiconductor cable in
    k[Region[{XLPE}]]  = kappa_xlpe; // Isolating material
    k[Region[{SEMI_OUT}]]  = kappa_semiconductor; // semiconductor Cable out
    k[Region[{TPS_OUT_PE}]]  = kappa_polyethilene; //Polyethilene sheath Bundle
    k[Region[{TPS_OUT_AL}]]  = kappa_al; // Aluminum shield
    k[Region[{POLYPROPILENE}]]  = kappa_polypropilene;
    k[Region[{ALUMINIUM_ARMOUR_IN}]]  = kappa_al;
    k[Region[{POLYETHYLENE_SHEAT_OUT}]]  = kappa_polyethilene;
    k[Region[{WATER_EM}]]  = kappa_water;
    k[Region[{WATER_TH}]]  = kappa_water;
    k[Region[{SOIL_EM}]]  = kappa_soil_sea;
    k[Region[{SOIL_TH}]]  = kappa_soil_sea;

    //Defect Characterization

    If(Flag_defect_in_XLPE)
      k[Region[{DEFECT}]]  = kappa_air; // The defect becomes Air
    Else
      k[Region[{DEFECT}]]  = kappa_xlpe; // The defect becomes XLPE: No defect
    EndIf
    // * heat conduction mechanism is the main heat transfer mechanism for an underground cable system
    // * all materials have constant thermal properties, including the thermal resistivity of the soil
    // * radiation and convection are not considered

    h_water = 100; // Based on values on the internet. https://www.sciencedirect.com/topics/engineering/convection-heat-transfer-coefficient
    h[] = h_water;//7.371 + 6.43*v_wind^0.75; // 1, 10 ... Convective coefficient [W/(m^2 K)]

//Electromagnetic part, definition of MU

    nu[Region[{1000,1001,1002}]] = 1/mu0*mur_other_metals;
    nu[Region[{TPS_OUT_AL}]] = 1/mu0*mur_other_metals;
    nu[Region[{ALUMINIUM_ARMOUR_IN}]] = 1/mu0*mur_other_metals;
    nu[Region[{semiconductorin,xlpe,semiconductorout,tpsoutpolypropilene,polypropilene,polyethylenesheatout,defect}]] = 1/mu0;

}


Constraint {
  // All the constraint hereafter must be adapted to your problem. Commented definitions are kept as example.

  // Electrical constraints
  { Name ElectricScalarPotential;
    Case {
      { Region Ind_1; Value V0; TimeFunction F_Cos_wt_p[]{2*Pi*Freq, Pa}; }
      { Region Ind_2; Value V0; TimeFunction F_Cos_wt_p[]{2*Pi*Freq, Pb}; }
      { Region Ind_3; Value V0; TimeFunction F_Cos_wt_p[]{2*Pi*Freq, Pc}; }

      { Region SurfaceGe0; Value 0; }
    }
  }
  { Name MagneticVectorPotential_2D;
    Case {
      { Region Sur_Dirichlet_Mag; Value 0.; }
    }
  }
  { Name Voltage_2D;
    Case {
    }
  }
  { Name Current_2D;
    Case {
      // constraint used if Inds in DomainS_Mag
      // example for a three-phase cable
       { Region Ind_1; Value I; TimeFunction F_Cos_wt_p[]{2*Pi*Freq, Pa}; }
       { Region Ind_2; Value I; TimeFunction F_Cos_wt_p[]{2*Pi*Freq, Pb}; }
       { Region Ind_3; Value I; TimeFunction F_Cos_wt_p[]{2*Pi*Freq, Pc}; }
    }
  }

  { Name DirichletTemp ;
    Case {
      { Type Assign; Region Sur_Dirichlet_Thermal ; Value Tambient[]; }
    }
  }

}

Include "Jacobian_Integration.pro"; // Normally no modification is needed

// The following files contain: basis functions, formulations, resolution, post-processing, post-operation
// Some adaptations may be needed
If (Flag_AnalysisType ==0)
  Include "electrodynamic_formulation.pro";
EndIf
If (Flag_AnalysisType ==1)
  Include "darwin_formulation.pro";
EndIf
If (Flag_AnalysisType > 2)
  Include "magneto-thermal_formulation.pro";
EndIf
