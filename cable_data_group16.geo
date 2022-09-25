// Geometrical data for cable of group 16

mm=1e-3; // Millimeters

//D_i = 33.7*mm;// Diameter of conductors
D_i = 25*mm;// Diameter of conductors


t_i = 1.2*mm;// Thickeness of inner semiconductor of first layer
t_xlpe = 14*mm; // Thickness of insulating material
t_o = 1*mm; //Thickness of insulation screen
//t_apl = 2.0*mm; // Thickness of APL sheath *******************
t_swa = 2*mm; // Thickness of the aluminum wire armour
t_ps = 1*mm; // Thickness of polyethilene sheath
t_pc = 1*mm; // Thickness of the polyethilene recovering
t_steel_outer = 4*mm;
D_tot = 70.5*mm; //Outer Diameter
D_bundle = 2*D_tot + 4*t_steel_outer + 4*t_pc + 2*2*mm; //Overall diameter of the cable
D_inf = 7*D_tot; // Radius of the EM domain
depth_cable=2.5; // [m] Laying depth of the cable
dinf_th=3; //thermal analysis
dinf_th_ground=3;
/// Material properties

epsr_polyethylene = 2.25;
espr_semiconductor = 2.25;
epsr_xlpe = 2.5;
epsr_metals = 1;

//Variables Added for Function 03122021
epsr_water = 70; // Normally 80, reduced because salt concentration. https://chemistry.stackexchange.com/questions/16434/salt-concentration-and-electrical-permittivity-of-water
epsr_soil = 30;
epsr_polypropilene = 2.36;

// Relative permeability
mur_steel = 4;
mur_other_metals = 1;

// Electrical conductivity S/m

sigma_cu = 5.99e7 ;
sigma_steel = 4.7e6;
sigma_al = 3.77e7;
sigma_polyethylene = 1.0e-18;
sigma_semiconductor = 2;
sigma_xlpe = 1.0e-18;
sigma_soil = 50e0;

//Variables Added for Function 03122021
sigma_water = 1.630; // https://aip.scitation.org/doi/pdf/10.1063/5.0018105
sigma_soil_sea = 1;
sigma_polypropilene = 1.0e-18;
sigma_air = 1.0e-15;

// Thermal conductivity

kappa_cu = 400;
kappa_steel = 50.2;
kappa_al = 237;
kappa_polyethilene = 0.46;
kappa_semiconductor = 10;
kappa_xlpe = 0.46;
kappa_soil = 0.4;
kappa_air = 0.025499;

//Variables Added for Function 03122021
kappa_water = 0.593; // https://aip.scitation.org/doi/pdf/10.1063/5.0018105
kappa_soil_sea = 2;
kappa_polypropilene = 0.11;//2.3e-4; /// Can be 0.11 according to https://www.intechopen.com/chapters/65584


V0 = 150000*1.4142  ;
//mesh properties
DefineConstant [  s = {1., Name " Parameters/Global mesh size factor"}];

// Defect check box
DefineConstant[
  Flag_defect_in_XLPE = {0,  Choices{0,1},
  Name "{00Parameters/10Defect in XLPE (phase 2)", Highlight "Green"}
  //dd = {2, Choices{2,5},Name "Parameters/11activation of defect area",
  //Visible Flag_defect_in_XLPE}
];

DefineConstant[
  Flag_sigma_funcT = {0,  Choices{0,1},
  Name "{01Parameters/10Temperature Dependency of Sigma", Highlight "Blue"}
  //dd = {2, Choices{2,5},Name "Parameters/11activation of defect area",
  //Visible Flag_defect_in_XLPE}
];

//physical number

AIR_IN = 900;
AIR_OUT = 901;
WIRE = 1000;

SEMI_IN  = 2000;
XLPE = 3000; DEFECT = 3001 ;
SEMI_OUT = 4000;
TPS_OUT_PE = 5000;
TPS_OUT_AL = 5001;

POLYPROPILENE = 5002;
POLYETHYLENE_SHEAT_OUT =6000;
ALUMINIUM_ARMOUR_IN = 7000;

WATER_EM = 8000;
WATER_TH = 9000;
SOIL_EM = 10000;
SOIL_TH = 11000;


OUTBND_WATER_EM = 1111;
OUTBND_WATER_TH = 2222;
OUTBNB_SOIL_EM = 3333;
OUTBNB_SOIL_TH = 4444;
INTERFACE_WATER_SOIL = 5555;
