Include "cable_data_group16.geo";

//Mesh.Algorithm = 6;

//Mesh.ElementOrder = 1;
Mesh.MinimumCircleNodes = 7*5;
Mesh.MinimumCurveNodes = 7;

SetFactory("OpenCASCADE");

dist_cab = D_i + 2*(t_i+t_xlpe+t_o+t_ps);
h = dist_cab * Sin(Pi/3);

// Points of centers of cables

x0 = 0.0;         y0 = 2*h/3;
x1 = -dist_cab/2; y1 = -h/3;
x2 =  dist_cab/2; y2 = -h/3;

//Conductor surface and isolator
sur_wire() = {};
  sur_wire(0) = news; Disk(news) = {x0, y0, 0., D_i/2};
  sur_wire(1) = news; Disk(news) = {x1, y1, 0., D_i/2};
  sur_wire(2) = news; Disk(news) = {x2, y2, 0., D_i/2};

sur_semi_in() = {};
  sur_semi_in(0) = news; Disk(news) = {x0, y0, 0., D_i/2+t_i};
  sur_semi_in(1) = news; Disk(news) = {x1, y1, 0., D_i/2+t_i};
  sur_semi_in(2) = news; Disk(news) = {x2, y2, 0., D_i/2+t_i};

sur_xlpe() = {};
  sur_xlpe(0) = news; Disk(news) = {x0, y0, 0., D_i/2+t_i+t_xlpe};
  sur_xlpe(1) = news; Disk(news) = {x1, y1, 0., D_i/2+t_i+t_xlpe};
  sur_xlpe(2) = news; Disk(news) = {x2, y2, 0., D_i/2+t_i+t_xlpe};

sur_semi_out() = {};
  sur_semi_out(0) = news; Disk(news) = {x0, y0, 0., D_i/2+t_i+t_xlpe+t_o};
  sur_semi_out(1) = news; Disk(news) = {x1, y1, 0., D_i/2+t_i+t_xlpe+t_o};
  sur_semi_out(2) = news; Disk(news) = {x2, y2, 0., D_i/2+t_i+t_xlpe+t_o};

sur_tps_out() = {};
  sur_tps_out(0) = news; Disk(news) = {x0, y0, 0., D_i/2+t_i+t_xlpe+t_o+t_ps};
  sur_tps_out(1) = news; Disk(news) = {x1, y1, 0., D_i/2+t_i+t_xlpe+t_o+t_ps};
  sur_tps_out(2) = news; Disk(news) = {x2, y2, 0., D_i/2+t_i+t_xlpe+t_o+t_ps};


// Outer triangles rounded

R0 = D_i/2 + t_i + t_xlpe + t_o;
R1 = R0 + t_ps;
R2 = R1 + t_swa;
arc_angle = 2*Pi/3;

arc_aux() += newl; Circle(newl) = {x0, y0, 0, R0, Pi/6, Pi/6+arc_angle}; // Arc from Pi/6 to Pi/6 + arc_angle
arc_aux() += newl; Circle(newl) = {x1, y1, 0, R0, 5*Pi/6, 5*Pi/6+arc_angle}; // Arc from Pi/6 to Pi/6 + arc_angle
arc_aux() += newl; Circle(newl) = {x2, y2, 0, R0, 9*Pi/6, 9*Pi/6+arc_angle}; // Arc from Pi/6 to Pi/6 + arc_angle

lin_aux() += newl; Line(newl) = {21,18};
lin_aux() += newl; Line(newl) = {17,20};
lin_aux() += newl; Line(newl) = {19,16};
Curve Loop(newll) = {21, 16, 20, 18, 19, 17};
sur_tri() += news; Plane Surface(news) = {newll-1};

arc_aux() += newl; Circle(newl) = {x0, y0, 0, R1, Pi/6, Pi/6+arc_angle}; // Arc from Pi/6 to Pi/6 + arc_angle
arc_aux() += newl; Circle(newl) = {x1, y1, 0, R1, 5*Pi/6, 5*Pi/6+arc_angle}; // Arc from Pi/6 to Pi/6 + arc_angle
arc_aux() += newl; Circle(newl) = {x2, y2, 0, R1, 9*Pi/6, 9*Pi/6+arc_angle}; // Arc from Pi/6 to Pi/6 + arc_angle

lin_aux() += newl; Line(newl) = {25,22};
lin_aux() += newl; Line(newl) = {23,26};
lin_aux() += newl; Line(newl) = {27,24};
Curve Loop(newll) = {27, 24, 28, 26, 29, 25};
sur_tri() += news; Plane Surface(news) = {newll-1};


arc_aux() += newl; Circle(newl) = {x0, y0, 0, R2, Pi/6, Pi/6+arc_angle}; // Arc from Pi/6 to Pi/6 + arc_angle
arc_aux() += newl; Circle(newl) = {x1, y1, 0, R2, 5*Pi/6, 5*Pi/6+arc_angle}; // Arc from Pi/6 to Pi/6 + arc_angle
arc_aux() += newl; Circle(newl) = {x2, y2, 0, R2, 9*Pi/6, 9*Pi/6+arc_angle}; // Arc from Pi/6 to Pi/6 + arc_angle

// Lines for the triangle

lin_aux() += newl; Line(newl) = {29,32};
lin_aux() += newl; Line(newl) = {33,30};
lin_aux() += newl; Line(newl) = {31,28};
Curve Loop(newll) = {37, 32, 35, 34, 36, 33};
sur_tri() += news; Plane Surface(news) = {newll-1};

// Outer Sheell

sur_covering() = {};
  sur_covering(0) = news; Disk(news) = {0., 0., 0., D_bundle/2};
  sur_covering(1) = news; Disk(news) = {0., 0., 0., D_bundle/2-t_steel_outer};
  sur_covering(2) = news; Disk(news) = {0., 0., 0., D_bundle/2-t_steel_outer - t_pc};


// Split areas
/*

*/
//Per Cable specific surface Arrays

sur_wire_cover0() = {sur_wire(0), sur_semi_in(0), sur_xlpe(0), sur_semi_out(0), sur_tps_out(0)};
sur_wire_cover1() = {sur_wire(1), sur_semi_in(1), sur_xlpe(1), sur_semi_out(1), sur_tps_out(1)};
sur_wire_cover2() = {sur_wire(2), sur_semi_in(2), sur_xlpe(2), sur_semi_out(2), sur_tps_out(2)};
/* Maybe we do not need this
sur_ps2(0) = news; BooleanDifference(news) = {Surface{sur_tps_out(0)}; Delete;}{Surface{sur_wire_cover0()}; };
sur_ps2(1) = news; BooleanDifference(news) = {Surface{sur_tps_out(1)}; Delete;}{Surface{sur_wire_cover1()}; };
sur_ps2(2) = news; BooleanDifference(news) = {Surface{sur_tps_out(2)}; Delete;}{Surface{sur_wire_cover2()}; };
*/

BooleanFragments{
  Surface{:}; Delete;
}{}


sur_wire() = {1,2,3};
sur_semi_in() = {4,5,6};
sur_xlpe() = {7,8};
sur_semi_out() = {10,11,12};
sur_tps_out() = {13,14,15,16,17,18,23,24,25,26};
sur_covering() = {27,28};



// EM Domain Definition: Disk of radius D_inf and make a BooleanDifference to remove the surfaces of the cables.
all_sur_cable() = Surface{:};

sur_EMdom() = news; Disk(news) = {0,0,0., D_inf/2};
//BooleanDifference(news) = {Surface{sur_EMdom};Delete;}{Surface{all_sur_cable()}; };
//sur_EMdom() = news-1;
//Printf("",sur_EMdom());
//Create the Rectangles of water and Soil (seabed)

sur_sea_waterTH_EM=news; Rectangle(news)={-dinf_th,(-D_bundle-0.01)/2,0,2*dinf_th,dinf_th_ground};
//BooleanDifference(news) = {Surface{sur_sea_waterTH_EM};Delete;}{Surface{sur_EMdom(),all_sur_cable()}; };
//sur_sea_waterTH_EM() = news-1;


sur_soil=news; Rectangle(news)={-dinf_th,(-D_bundle-0.01)/2-dinf_th_ground,0,2*dinf_th,dinf_th_ground,0};

//Surface of defect in xlpe
sur_defect = news; Disk(news) = {x2-D_i/2-4*t_i, y2, 0., 2*t_i};

BooleanFragments{
  Surface{:}; Delete;
}{}

sur_EMdom()= {71,72};
sur_defect = 69;
sur_xlpe() = {7,8,70};

//Thermal domain

//ajusting the charasteristic mesh size
cl = D_bundle/s;

//Characteristic Length {Point(:)} = cl;

Characteristic Length { PointsOf { Surface {sur_covering()}; } } = cl/16 ;
Characteristic Length { PointsOf { Surface {sur_wire(), sur_semi_in(),sur_xlpe(),sur_semi_out(), sur_tps_out(), sur_defect}; } } = cl/32;
Characteristic Length { PointsOf { Surface {sur_EMdom()}; } } = cl/16 ;
Characteristic Length { PointsOf { Surface {73,74}; } } = cl*3 ;
Characteristic Length { PointsOf { Line {62,63,64}; } } = cl/4 ;
//Characteristic Length { PointsOf { Line {63}; } } = cl/16 ;


Physical Surface ("wire 1 ", WIRE + 0 ) = 1;
Physical Surface ("wire 2 ", WIRE + 1 ) = 2;
Physical Surface ("wire 3 ", WIRE + 2 ) = 3;

Physical Surface ("inner semiconductor ", SEMI_IN) = sur_semi_in();
Physical Surface ("XLPE ", XLPE ) = sur_xlpe();
Physical Surface ("Defect ", DEFECT ) = sur_defect;
Physical Surface ("outer semiconductor ", SEMI_OUT) = sur_semi_out();
Physical Surface ("tps out polyethylene ", TPS_OUT_PE ) = {13,14,15,16,17,18,23,24,25};
Physical Surface ("tps out aluminium ", TPS_OUT_AL ) = 26;

Physical Surface ("Polypropilene filling",POLYPROPILENE ) = {19,20,21,22,29};
Physical Surface ("inner aluminium armour ", ALUMINIUM_ARMOUR_IN ) = 28;
Physical Surface ("outer polyethylene sheat ",POLYETHYLENE_SHEAT_OUT ) = 27;

Physical Surface ("water (EM) ", WATER_EM ) = 71;
Physical Surface ("water (TH) ", WATER_TH ) = 73;
Physical Surface ("soil under the sea (EM) ", SOIL_EM ) = 72;
Physical Surface ("soil under the sea(TH) ", SOIL_TH ) = 74;


Physical Line ("outer boundary water (EM) ", OUTBND_WATER_EM ) = {64,62};
Physical Line ("outer boundary water (TH) ", OUTBND_WATER_TH ) = {67,68,69};
Physical Line ("Outer boundary soil (EM) ", OUTBNB_SOIL_EM ) = {66};
Physical Line ("Outer boundary soil(TH) ", OUTBNB_SOIL_TH ) = {72,73,74};
Physical Line ("interface water/soil", INTERFACE_WATER_SOIL ) = {71,63,70};


/*
*/
