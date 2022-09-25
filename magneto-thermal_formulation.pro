FunctionSpace {

   { Name Hcurl_a_Mag_2D; Type Form1P;
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_PerpendicularEdge;
        Support Domain_Mag; Entity NodesOf[ All ]; }
      If (Flag_Degree_a == 2)
        { Name se2; NameOfCoef ae2; Function BF_PerpendicularEdge_2E;
          Support Domain_Mag; Entity EdgesOf[ All ]; }
      EndIf
    }
    Constraint {
      { NameOfCoef ae;
        EntityType NodesOf; NameOfConstraint MagneticVectorPotential_2D; }
      If (Flag_Degree_a == 2)
	{ NameOfCoef ae2; // Only OK if homogeneous BC, otherwise specify zero-BC
          EntityType EdgesOf; NameOfConstraint MagneticVectorPotential_2D; }
      EndIf
    }
  }

  { Name Hregion_i_2D ; Type Vector ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_RegionZ ;
        Support DomainS_Mag ; Entity DomainS_Mag ; }
    }
    GlobalQuantity {
      { Name Is ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Us ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Us ; EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef Is ; EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }

  // Temperature is discretised with this BF
  { Name Hgrad_Thermal; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef t; Function BF_Node;
        Support Domain_Thermal; Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef t; EntityType NodesOf ; NameOfConstraint DirichletTemp; }
    }
  }

}



Formulation {

  { Name Darwin_a_2D; Type FemEquation; // Magnetodynamics + displacement current, no coupling
    Quantity {
      { Name a;  Type Local; NameOfSpace Hcurl_a_Mag_2D; }

      { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
      { Name Us ; Type Global ; NameOfSpace Hregion_i_2D[Us] ; }
      { Name Is ; Type Global ; NameOfSpace Hregion_i_2D[Is] ; }

      { Name T; Type Local ; NameOfSpace Hgrad_Thermal; }
    }
    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      Galerkin { DtDof [ sigma[{T}] * Dof{a} , {a} ];
        In DomainC_Mag; Jacobian Vol; Integration I1; }
      Galerkin { DtDtDof [ epsilon[] * Dof{a} , {a} ]; // Added term => Darwin approximation
        In DomainC_Mag; Jacobian Vol; Integration I1; }

      // Either you impose directly the function js0[]
      Galerkin { [ -js0[] , {a} ];
        In DomainS0_Mag; Jacobian Vol; Integration I1; }

      // or you use the constraints => allows accounting for sigma[{T}]
      Galerkin { [ -Ns[]/Sc[] * Dof{ir}, {a} ] ;
        In DomainS_Mag ; Jacobian Vol ; Integration I1 ; }
      Galerkin { DtDof [ Ns[]/Sc[] * Dof{a}, {ir} ] ;
        In DomainS_Mag ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ Ns[]/Sc[] / sigma[{T}] * Ns[]/Sc[]* Dof{ir} , {ir} ] ; // resistance term
        In DomainS_Mag ; Jacobian Vol ; Integration I1 ; }
      //GlobalTerm { [ Rdc * Dof{Is} , {Is} ] ; In DomainS ; } // OR this resitance term
      GlobalTerm { [ Dof{Us}, {Is} ] ; In DomainS_Mag ; }

    }
  }

  { Name ThermalSta ; Type FemEquation;
    Quantity {
      { Name T; Type Local ; NameOfSpace Hgrad_Thermal; }

      // quantities from previous formulations, not unknowns
      { Name a; Type Local ; NameOfSpace Hcurl_a_Mag_2D ; }
      { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
    }
    Equation {
      Galerkin { [ k[] * Dof{d T} , {d T} ];
	In Vol_Thermal; Integration I1; Jacobian Vol;  }

      // Thermal source = Joule losses from EM computation
      // Use <a> when mixing frequency domain and time domain computations
      // {a} is complex, losses are computed taking that into account
      Galerkin { [ -0.5*sigma[{T}] * <a>[ SquNorm[Dt[{a}]] ], {T} ];
        In Vol_QSource_Thermal; Integration I1; Jacobian Vol;  }
      Galerkin { [ -0.5/sigma[{T}] * <a>[ SquNorm[js0[]] ], {T} ];
        In Vol_QSource0_Thermal; Integration I1; Jacobian Vol;  }
      Galerkin { [ -0.5/sigma[{T}] * <ir>[ SquNorm[Ns[]/Sc[]*{ir}] ], {T} ];
        In Vol_QSourceB_Thermal; Integration I1; Jacobian Vol;  }

      // Convection boundary condition
      Galerkin { [ h[] * Dof{T} , {T} ] ;
	In Sur_Convection_Thermal; Jacobian Sur ; Integration I1 ; }
      Galerkin { [-h[] * Tambient[] , {T} ] ;
        In Sur_Convection_Thermal ; Jacobian Sur ; Integration I1 ; }
    }
  }

}

Resolution {

  { Name Magneto_thermal;
    System {
      { Name Sys_Mag; NameOfFormulation Darwin_a_2D;
        Type Complex; Frequency Freq; }
      { Name Sys_The; NameOfFormulation ThermalSta; }
    }
    Operation {
      CreateDir["res"];
      InitSolution[Sys_Mag];
      InitSolution[Sys_The];

      If(!Flag_sigma_funcT)
        Generate[Sys_Mag]; Solve[Sys_Mag];
        Generate[Sys_The]; Solve[Sys_The];
      Else
        IterativeLoop[ nb_iter, stop_criterion, relaxation_factor]{
          GenerateJac[Sys_Mag]; SolveJac[Sys_Mag];
          GenerateJac[Sys_The]; SolveJac[Sys_The];
        }
      EndIf
      SaveSolution[Sys_Mag];
      SaveSolution[Sys_The];

      PostOperation[Mag_Maps];
      PostOperation[Mag_Global];
      PostOperation[The_Maps];
    }
  }

}

PostProcessing {

  { Name Darwin_a_2D; NameOfFormulation Darwin_a_2D;
    PostQuantity {
      { Name a; Value { Term { [ {a} ]; In Domain_Mag; Jacobian Vol; } } }
      { Name az; Value { Term { [ CompZ[{a}] ]; In Domain_Mag; Jacobian Vol; } } }
      { Name b; Value { Term { [ {d a} ]; In Domain_Mag; Jacobian Vol; } } }
      { Name bm; Value { Term { [ Norm[{d a}] ]; In Domain_Mag; Jacobian Vol; } } }
      { Name j; Value {
          Term { [ -sigma[{T}]*Dt[{a}] ]; In DomainC_Mag;  Jacobian Vol; }
          Term { [ js0[] ]              ; In DomainS0_Mag; Jacobian Vol; }
          Term { [ Ns[]/Sc[]*{ir} ]     ; In DomainS_Mag;  Jacobian Vol; }
        } }
      { Name jz; Value {
          Term { [ CompZ[-sigma[{T}]*Dt[{a}]] ]; In DomainC_Mag; Jacobian Vol; }
          Term { [ CompZ[ js0[] ] ];             In DomainS0_Mag; Jacobian Vol; }
          Term { [ CompZ[ Ns[]/Sc[]*{ir} ]];     In DomainS_Mag; Jacobian Vol; }
        } }
      { Name jm;
        Value { Term { [ Norm[-sigma[{T}]*Dt[{a}]] ]; In DomainC_Mag; Jacobian Vol; }
         Term { [ Norm[ js0[] ] ];                    In DomainS0_Mag; Jacobian Vol; }
          Term { [ Norm[ Ns[]/Sc[]*{ir} ]];           In DomainS_Mag; Jacobian Vol; }
        } }

      { Name d; Value { Term { [ epsilon[] * Dt[Dt[{a}]] ]; In DomainC_Mag; Jacobian Vol; } } }
      { Name dz; Value { Term { [ CompZ[epsilon[] *  Dt[Dt[{a}]] ] ]; In DomainC_Mag; Jacobian Vol; } } }
      { Name dm; Value { Term { [ Norm[epsilon[]  *  Dt[Dt[{a}]] ] ]; In DomainC_Mag; Jacobian Vol; } } }


      { Name rhoj2; Value { // local losses
          Term { [ 0.5*sigma[{T}]*SquNorm[Dt[{a}]] ]; In DomainC_Mag; Jacobian Vol; }
          Term { [ 0.5/sigma[{T}]*SquNorm[js0[]] ]; In DomainS0_Mag ; Jacobian Vol; }
          Term { [ 0.5/sigma[{T}]*SquNorm[Ns[]/Sc[]*{ir}] ]; In DomainS_Mag; Jacobian Vol; }
        }
      }
      { Name JouleLosses; Value { // global losses
          Integral { [ 0.5*sigma[{T}]*SquNorm[Dt[{a}]] ]        ; In DomainC_Mag  ; Jacobian Vol ; Integration I1 ; }
          Integral { [ 0.5/sigma[{T}]*SquNorm[js0[]] ]          ; In DomainS0_Mag ; Jacobian Vol ; Integration I1 ; }
          Integral { [ 0.5/sigma[{T}]*SquNorm[Ns[]/Sc[]*{ir}] ] ; In DomainS_Mag  ; Jacobian Vol ; Integration I1 ; }
        }
      }

      { Name U ; Value {
          Term { [ {Us} ] ; In DomainS_Mag ; }
        }
      }

      { Name I ; Value {
          Term { [ {Is} ] ; In DomainS_Mag ; }
        }
      }

      { Name S ; Value {
          Term { [ {Us}*Conj[{Is}] ] ; In DomainS_Mag ; }
        }
      }

      { Name R ; Value {
          Term { [ -Re[{Us}/{Is}] ] ; In DomainS_Mag ; }
        }
      }

      { Name L ; Value {
          Term { [ -Im[{Us}/{Is}]/(2*Pi*Freq) ] ; In DomainS_Mag ; }
        }
      }

    }
  }



  { Name ThermalSta; NameOfFormulation ThermalSta; //NameOfSystem T;
    Quantity {
      { Name T; Value{ Local{ [ {T} ] ; In Vol_Thermal; Jacobian Vol; } } }
      { Name TinC; Value{ Local{ [ {T}-273.15 ] ; In Vol_Thermal; Jacobian Vol; } } }
      { Name q; Value{ Local{ [ -k[]*{d T} ] ; In Vol_Thermal; Jacobian Vol; } } } // heat flux [W/m^2]
      
    }
  }

}


PostOperation{
  // Magnetic
  //-------------------------------
 { Name Mag_Maps; NameOfPostProcessing Darwin_a_2D;
    Operation {
      // You may want to see the result only in part of the domain => adapt domain that follows OnElementsOf
      // Name is not compulsory, it can be modified
      Print[ bm , OnElementsOf Domain_Mag,
        Name "|B| [T]", File "res/bm.pos" ];
      Print[ jm , OnElementsOf DomainC_Mag,
        Name "|j| [A/m²] Al sheath", File "res/jm.pos" ];
      Print[ dm , OnElementsOf DomainC_Mag,
        Name "|D| [C/m²]", File "res/dm.pos" ];
    }
  }

  // This two definitions refer to the GUI => organisation of the results
  po = "{01Losses/";
  po2 = "{02PU-parameters/";
  { Name Mag_Global; NameOfPostProcessing Darwin_a_2D;
    Operation {
      // You may restrict DomainC_Mag to part of it
      // You may also computed the losses in DomainS_Mag
      Print[ JouleLosses[DomainC_Mag], OnGlobal, Format Table,
        SendToServer StrCat[po,"All conducting domain"],
        Units "W/m", File "res/losses_total.dat" ];

      // Region to adapt according to your cable
      Print[ R, OnRegion Ind_1, Format Table,
        SendToServer StrCat[po2,"0R"],
        Units "Ω", File "res/Rinds.dat" ];
      Print[ L, OnRegion Ind_1, Format Table,
        SendToServer StrCat[po2,"1L"],
        Units "H", File "res/Linds.dat" ];

    }
  }

  // Thermal
  // -------------------------------
  { Name The_Maps; NameOfPostProcessing ThermalSta;
    Operation {
      Print[ TinC , OnElementsOf Region[{Vol_Thermal,-Cable}], // all the thermal domain but the cable
        Name "T [°C] around cable", File "res/T.pos" ];
      Print[ TinC , OnElementsOf Cable,
        Name "T [°C] cable", File "res/Tcable.pos" ];
      Print[ q , OnElementsOf Region[{Vol_Thermal,-Cable}], // all the thermal domain but the cable
        Name "heat flux [W/m²] around cable", File "res/q.pos" ];
      Print[ q , OnElementsOf Cable,
        Name "heat flux [W/m²] cable", File "res/qcable.pos" ];

    }
  }

}
