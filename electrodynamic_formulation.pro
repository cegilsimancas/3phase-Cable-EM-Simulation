FunctionSpace {

  { Name Hgrad_v_Ele; Type Form0;
    BasisFunction {
      // v = \sum_n v_n  s_n,  for all nodes
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Domain_Ele; Entity NodesOf[ All ]; }
      If (Flag_Degree_v == 2)
        { Name sn2; NameOfCoef vn2; Function BF_Node_2E;
          Support Domain_Ele; Entity EdgesOf[ All ]; }
      EndIf
    }

    Constraint {
      { NameOfCoef vn; EntityType NodesOf;
        NameOfConstraint ElectricScalarPotential; }
      If (Flag_Degree_v == 2)
        { NameOfCoef vn2;
          EntityType EdgesOf; NameOfConstraint ZeroElectricScalarPotential; }
      EndIf
    }
  }

}

Formulation {

  { Name Electrodynamics_v; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
    }
    Equation {
      Galerkin { [ sigma[] * Dof{d v} , {d v} ] ;
        In Domain_Ele; Jacobian Vol ; Integration I1 ; }
      Galerkin { DtDof[ epsilon[] * Dof{d v} , {d v} ];
        In Domain_Ele; Jacobian Vol; Integration I1; }
    }
  }

}

Resolution {

  { Name Electrodynamics;
    System {
      { Name Sys_Ele; NameOfFormulation Electrodynamics_v;
        Type Complex; Frequency Freq; }
    }
    Operation {
      CreateDir["res"];

      Generate[Sys_Ele]; Solve[Sys_Ele]; SaveSolution[Sys_Ele];
      PostOperation[Ele_Maps];
      // PostOperation[Ele_Cuts]; // To adapt for your cable
    }
  }

}


PostProcessing {

  { Name EleDyn_v; NameOfFormulation Electrodynamics_v;
    Quantity {
      { Name v; Value { Term { [ {v} ]; In Domain_Ele; Jacobian Vol; } } }
      { Name e; Value { Term { [ -{d v} ]; In Domain_Ele; Jacobian Vol; } } }
      { Name em; Value { Term { [ Norm[-{d v}] ]; In Domain_Ele; Jacobian Vol; } } }

      { Name d; Value { Term { [ -epsilon[] * {d v} ]; In Domain_Ele; Jacobian Vol; } } }
      { Name dm; Value { Term { [ Norm[-epsilon[] * {d v}] ]; In Domain_Ele; Jacobian Vol; } } }

      { Name j ; Value { Term { [ -sigma[] * {d v} ] ; In Domain_Ele ; Jacobian Vol; } } }
      { Name jm ; Value { Term { [ Norm[-sigma[] * {d v}] ] ; In Domain_Ele ; Jacobian Vol; } } }

      { Name jtot ; Value {
          Term { [ -sigma[] * {d v} ] ;       In Domain_Ele ; Jacobian Vol; }
          Term { [ -epsilon[] * Dt[{d v}] ] ; In Domain_Ele ; Jacobian Vol; }
        } }

      { Name ElectricEnergy; Value {
          Integral {
            [ 0.5 * epsilon[] * SquNorm[{d v}] ];
            In Domain_Ele; Jacobian Vol; Integration I1;
          }
	}
      }

      { Name V0 ; Value {
          // For recovering the imposed voltage in post-pro
          // Most likely you will need to adapt for your cable
          // The default hereafter is for a three-phase cable
          Term { Type Global ; [ V0 * F_Cos_wt_p[]{2*Pi*Freq, Pa}] ; In Ind_1 ; }
          Term { Type Global ; [ V0 * F_Cos_wt_p[]{2*Pi*Freq, Pb}] ; In Ind_2 ; }
          Term { Type Global ; [ V0 * F_Cos_wt_p[]{2*Pi*Freq, Pc}] ; In Ind_3 ; }
        } }

      { Name C_from_Energy ; Value { Term { Type Global; [ 2*$We/SquNorm[$voltage] ] ; In DomainDummy ; } } }
    }
  }

}

PostOperation{

  // Electric
  //-------------------------------

  po0 = "{01Capacitance/"; // Useful only for the GUI

  { Name Ele_Maps; NameOfPostProcessing EleDyn_v;
    Operation {
      Print[ v,  OnElementsOf Domain_Ele, File "res/v.pos" ];
      Print[ em, OnElementsOf Domain_Ele, Name "|E| [V/m]",  File "res/em.pos" ]; // Name is not compulsory, it may be adapted
      Print[ dm, OnElementsOf Domain_Ele, Name "|D| [C/m²]", File "res/dm.pos" ];
      Print[ ElectricEnergy[Domain_Ele], OnGlobal, Format Table, StoreInVariable $We,
        SendToServer StrCat[po0,"0Electric energy"], File "res/energy.dat" ];

      Print[ V0, OnRegion Ind_1, Format Table, StoreInVariable $voltage,
        SendToServer StrCat[po0,"0U1"], Units "V", File "res/U.dat" ];
      Print[ C_from_Energy, OnRegion DomainDummy, Format Table, StoreInVariable $C1,
        SendToServer StrCat[po0,"1Cpha"], Units "F/m", File "res/C.dat" ];
    }
  }

  /*
  // To adapt for your cable
  dist_cab = dc + 2*(ti+txlpe+to+tapl)+tps;
  h = dist_cab * Sin[Pi/3]; // height of equilateral triangle
  x0 = 0; y0 = 2*h/3;
  x1 = -dist_cab/2; y1 = -h/3;
  x2 =  dist_cab/2; y2 = -h/3;


  // Default visualisation (3D) is not handy for cuts
  // Change the plot type to 2D: Options->View[]->General->Plot type
  { Name Ele_Cuts; NameOfPostProcessing EleDyn_v;
    Operation {
      Print[ em , OnLine { {x2,y2,0} {x2+dc/2+ti+txlpe+to+tapl,y2,0} } {100},
        Name "|E| [V/m] cut in phase 2", File "res/em_cut.pos"];
    }
  }
  */
}
