

# include "TrackerAlignment/inc/AlignmentDerivatives.hh"
# include <math.h>
# include <vector>

double CosmicTrack_DCA(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = 1.0*R5*R6;
    double R8 = sin(plane_b);
    double R9 = sin(panel_g);
    double R10 = cos(panel_a);
    double R11 = R10*R9;
    double R12 = sin(panel_b);
    double R13 = cos(panel_g);
    double R14 = R13*R5;
    double R15 = -R11 + R12*R14;
    double R16 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R17 = 1.0/R16;
    double R18 = R17*panel_straw0x;
    double R19 = R10*R13;
    double R20 = R5*R9;
    double R21 = R12*R20 + R19;
    double R22 = R17*panel_straw0y;
    double R23 = 1.0*R15*R18 - 1.0*R21*R22;
    double R24 = 1.0*R15*R22 + 1.0*R18*R21;
    double R25 = sin(plane_a);
    double R26 = R25*R3;
    double R27 = -R23*R8 + R24*R26 + R4*R7;
    double R28 = sin(plane_g);
    double R29 = R25*R28;
    double R30 = cos(plane_g);
    double R31 = R2*R30;
    double R32 = R29 + R31*R8;
    double R33 = R3*R30;
    double R34 = R2*R28;
    double R35 = R25*R30;
    double R36 = -R34 + R35*R8;
    double R37 = R23*R33 + R24*R36 + R32*R7;
    double R38 = R34*R8 - R35;
    double R39 = R28*R3;
    double R40 = R29*R8 + R31;
    double R41 = R23*R39 + R24*R40 + R38*R7;
    double R42 = pow(pow(R27, 2) + pow(R37, 2) + pow(R41, 2), -1.0/2.0);
    double R43 = R27*R42;
    double R44 = R41*R42;
    double R45 = R0*R44;
    double R46 = R37*R42;
    double R47 = R0*a1;
    double R48 = R1*R43 - R45 + R46*R47;
    double R49 = 1.0/(1.0 - pow(R48, 2));
    double R50 = -plane_z;
    double R51 = R50 + panel_dz + panel_straw0z;
    double R52 = R18*panel_dx - R22*panel_dy + panel_straw0x;
    double R53 = R18*panel_dy + R22*panel_dx + panel_straw0y;
    double R54 = -R16 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R55 = R18*R6;
    double R56 = R22*R6;
    double R57 = R13*R55 - R56*R9;
    double R58 = R13*R56 + R55*R9;
    double R59 = -panel_straw0z + wire_z;
    double R60 = R10*R6;
    double R61 = R12*R19 + R20;
    double R62 = R11*R12 - R14;
    double R63 = R18*R61 - R22*R62;
    double R64 = R18*R62 + R22*R61;
    double R65 = -R26*R53 - R4*R51 + R50 + R52*R8 - R54*(-R12*R4 + R26*R58 - R57*R8) - R59*(R26*R64 + R4*R60 - R63*R8) + b0 - plane_dz;
    double R66 = R1*R65;
    double R67 = -R38*R51 - R39*R52 - R40*R53 - R54*(-R12*R38 + R39*R57 + R40*R58) - R59*(R38*R60 + R39*R63 + R40*R64) - plane_dy;
    double R68 = R0*R67;
    double R69 = -R32*R51 - R33*R52 - R36*R53 - R54*(-R12*R32 + R33*R57 + R36*R58) - R59*(R32*R60 + R33*R63 + R36*R64) + a0 - plane_dx;
    double R70 = R47*R69;
    double R71 = R43*R65 + R44*R67 + R46*R69;
    double R72 = R49*(R48*R71 - R66 + R68 - R70);
    double R73 = R49*(-R48*(R66 - R68 + R70) + R71);
    double R74 = sqrt(pow(-R0*R72 - R44*R73 + R67, 2) + pow(R1*R72 - R43*R73 + R65, 2) + pow(-R46*R73 + R47*R72 + R69, 2));
    double R75 = R0*R46 + R45*a1;
    double R76 = R1*R46 - R43*R47;
    double R77 = -R0*R43 - R45*b1;
    double R78 = pow(pow(R75, 2) + pow(R76, 2) + pow(R77, 2), -1.0/2.0);
    double result = ((R65*R75*R78 + R67*R76*R78 + R69*R77*R78 > 0) ? (
   -R74
)
: (
   R74
));
    return result;
}


double CosmicTrack_DCAalignpos_x(double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z)
{
    double R0 = sin(plane_a);
    double R1 = sin(plane_g);
    double R2 = cos(plane_a);
    double R3 = cos(plane_g);
    double R4 = R3*sin(plane_b);
    double R5 = R0*R1 + R2*R4;
    double R6 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R7 = 1.0/R6;
    double R8 = R7*panel_straw0x;
    double R9 = R7*panel_straw0y;
    double R10 = R3*cos(plane_b);
    double R11 = R0*R4 - R1*R2;
    double R12 = sin(panel_b);
    double R13 = cos(panel_g);
    double R14 = cos(panel_b);
    double R15 = R14*R8;
    double R16 = sin(panel_g);
    double R17 = R14*R9;
    double R18 = cos(panel_a);
    double R19 = sin(panel_a);
    double R20 = R12*R18;
    double R21 = R13*R20 + R16*R19;
    double R22 = -R13*R19 + R16*R20;
    double result = R10*(R8*panel_dx - R9*panel_dy + panel_straw0x) + R11*(R8*panel_dy + R9*panel_dx + panel_straw0y) + R5*(panel_dz + panel_straw0z - plane_z) + plane_dx + (-R6 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)))*(R10*(R13*R15 - R16*R17) + R11*(R13*R17 + R15*R16) - R12*R5) + (-panel_straw0z + wire_z)*(R10*(R21*R8 - R22*R9) + R11*(R21*R9 + R22*R8) + R14*R18*R5);
    return result;
}


double CosmicTrack_DCAalignpos_y(double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z)
{
    double R0 = sin(plane_a);
    double R1 = cos(plane_g);
    double R2 = cos(plane_a);
    double R3 = sin(plane_g);
    double R4 = R3*sin(plane_b);
    double R5 = -R0*R1 + R2*R4;
    double R6 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R7 = 1.0/R6;
    double R8 = R7*panel_straw0x;
    double R9 = R7*panel_straw0y;
    double R10 = R3*cos(plane_b);
    double R11 = R0*R4 + R1*R2;
    double R12 = sin(panel_b);
    double R13 = cos(panel_g);
    double R14 = cos(panel_b);
    double R15 = R14*R8;
    double R16 = sin(panel_g);
    double R17 = R14*R9;
    double R18 = cos(panel_a);
    double R19 = sin(panel_a);
    double R20 = R12*R18;
    double R21 = R13*R20 + R16*R19;
    double R22 = -R13*R19 + R16*R20;
    double result = R10*(R8*panel_dx - R9*panel_dy + panel_straw0x) + R11*(R8*panel_dy + R9*panel_dx + panel_straw0y) + R5*(panel_dz + panel_straw0z - plane_z) + plane_dy + (-R6 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)))*(R10*(R13*R15 - R16*R17) + R11*(R13*R17 + R15*R16) - R12*R5) + (-panel_straw0z + wire_z)*(R10*(R21*R8 - R22*R9) + R11*(R21*R9 + R22*R8) + R14*R18*R5);
    return result;
}


double CosmicTrack_DCAalignpos_z(double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z)
{
    double R0 = cos(plane_b);
    double R1 = R0*cos(plane_a);
    double R2 = sin(plane_b);
    double R3 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R4 = 1.0/R3;
    double R5 = R4*panel_straw0x;
    double R6 = R4*panel_straw0y;
    double R7 = R0*sin(plane_a);
    double R8 = sin(panel_b);
    double R9 = cos(panel_g);
    double R10 = cos(panel_b);
    double R11 = R10*R5;
    double R12 = sin(panel_g);
    double R13 = R10*R6;
    double R14 = cos(panel_a);
    double R15 = sin(panel_a);
    double R16 = R14*R8;
    double R17 = R12*R15 + R16*R9;
    double R18 = R12*R16 - R15*R9;
    double result = R1*(panel_dz + panel_straw0z - plane_z) - R2*(R5*panel_dx - R6*panel_dy + panel_straw0x) + R7*(R5*panel_dy + R6*panel_dx + panel_straw0y) + plane_dz + plane_z + (-R3 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)))*(-R1*R8 - R2*(R11*R9 - R12*R13) + R7*(R11*R12 + R13*R9)) + (-panel_straw0z + wire_z)*(R1*R10*R14 - R2*(R17*R5 - R18*R6) + R7*(R17*R6 + R18*R5));
    return result;
}


double CosmicTrack_DCAaligndir_x(double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z)
{
    double R0 = sin(panel_a);
    double R1 = sin(plane_a);
    double R2 = sin(plane_g);
    double R3 = cos(plane_a);
    double R4 = cos(plane_g);
    double R5 = R4*sin(plane_b);
    double R6 = sin(panel_g);
    double R7 = cos(panel_a);
    double R8 = cos(panel_g);
    double R9 = R0*sin(panel_b);
    double R10 = -R6*R7 + R8*R9;
    double R11 = pow(pow(panel_straw0x, 2) + pow(panel_straw0y, 2), -1.0/2.0);
    double R12 = R11*panel_straw0x;
    double R13 = R6*R9 + R7*R8;
    double R14 = R11*panel_straw0y;
    double result = 1.0*R0*(R1*R2 + R3*R5)*cos(panel_b) + 1.0*R4*(R10*R12 - R13*R14)*cos(plane_b) + 1.0*(R1*R5 - R2*R3)*(R10*R14 + R12*R13);
    return result;
}


double CosmicTrack_DCAaligndir_y(double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z)
{
    double R0 = sin(panel_a);
    double R1 = sin(plane_a);
    double R2 = cos(plane_g);
    double R3 = cos(plane_a);
    double R4 = sin(plane_g);
    double R5 = R4*sin(plane_b);
    double R6 = sin(panel_g);
    double R7 = cos(panel_a);
    double R8 = cos(panel_g);
    double R9 = R0*sin(panel_b);
    double R10 = -R6*R7 + R8*R9;
    double R11 = pow(pow(panel_straw0x, 2) + pow(panel_straw0y, 2), -1.0/2.0);
    double R12 = R11*panel_straw0x;
    double R13 = R6*R9 + R7*R8;
    double R14 = R11*panel_straw0y;
    double result = 1.0*R0*(-R1*R2 + R3*R5)*cos(panel_b) + 1.0*R4*(R10*R12 - R13*R14)*cos(plane_b) + 1.0*(R1*R5 + R2*R3)*(R10*R14 + R12*R13);
    return result;
}


double CosmicTrack_DCAaligndir_z(double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z)
{
    double R0 = sin(panel_a);
    double R1 = 1.0*cos(plane_b);
    double R2 = sin(panel_g);
    double R3 = cos(panel_a);
    double R4 = cos(panel_g);
    double R5 = R0*sin(panel_b);
    double R6 = -R2*R3 + R4*R5;
    double R7 = pow(pow(panel_straw0x, 2) + pow(panel_straw0y, 2), -1.0/2.0);
    double R8 = R7*panel_straw0x;
    double R9 = R2*R5 + R3*R4;
    double R10 = R7*panel_straw0y;
    double result = R0*R1*cos(panel_b)*cos(plane_a) + R1*(R10*R6 + R8*R9)*sin(plane_a) - 1.0*(-R10*R9 + R6*R8)*sin(plane_b);
    return result;
}


double CosmicTrack_DCA_Deriv_a0(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = 1.0*R5*R6;
    double R8 = sin(plane_b);
    double R9 = sin(panel_g);
    double R10 = cos(panel_a);
    double R11 = R10*R9;
    double R12 = sin(panel_b);
    double R13 = cos(panel_g);
    double R14 = R13*R5;
    double R15 = -R11 + R12*R14;
    double R16 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R17 = 1.0/R16;
    double R18 = R17*panel_straw0x;
    double R19 = R10*R13;
    double R20 = R5*R9;
    double R21 = R12*R20 + R19;
    double R22 = R17*panel_straw0y;
    double R23 = 1.0*R15*R18 - 1.0*R21*R22;
    double R24 = 1.0*R15*R22 + 1.0*R18*R21;
    double R25 = sin(plane_a);
    double R26 = R25*R3;
    double R27 = -R23*R8 + R24*R26 + R4*R7;
    double R28 = sin(plane_g);
    double R29 = R25*R28;
    double R30 = cos(plane_g);
    double R31 = R2*R30;
    double R32 = R29 + R31*R8;
    double R33 = R3*R30;
    double R34 = R2*R28;
    double R35 = R25*R30;
    double R36 = -R34 + R35*R8;
    double R37 = R23*R33 + R24*R36 + R32*R7;
    double R38 = R34*R8 - R35;
    double R39 = R28*R3;
    double R40 = R29*R8 + R31;
    double R41 = R23*R39 + R24*R40 + R38*R7;
    double R42 = pow(pow(R27, 2) + pow(R37, 2) + pow(R41, 2), -1.0/2.0);
    double R43 = R27*R42;
    double R44 = R41*R42;
    double R45 = R0*R44;
    double R46 = R0*a1;
    double R47 = R37*R42;
    double R48 = R1*R43 - R45 + R46*R47;
    double R49 = 1.0/(1.0 - pow(R48, 2));
    double R50 = -plane_z;
    double R51 = R50 + panel_dz + panel_straw0z;
    double R52 = R18*panel_dx - R22*panel_dy + panel_straw0x;
    double R53 = R18*panel_dy + R22*panel_dx + panel_straw0y;
    double R54 = -R16 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R55 = R18*R6;
    double R56 = R22*R6;
    double R57 = R13*R55 - R56*R9;
    double R58 = R13*R56 + R55*R9;
    double R59 = -panel_straw0z + wire_z;
    double R60 = R10*R6;
    double R61 = R12*R19 + R20;
    double R62 = R11*R12 - R14;
    double R63 = R18*R61 - R22*R62;
    double R64 = R18*R62 + R22*R61;
    double R65 = -R26*R53 - R4*R51 + R50 + R52*R8 - R54*(-R12*R4 + R26*R58 - R57*R8) - R59*(R26*R64 + R4*R60 - R63*R8) + b0 - plane_dz;
    double R66 = R1*R65;
    double R67 = -R38*R51 - R39*R52 - R40*R53 - R54*(-R12*R38 + R39*R57 + R40*R58) - R59*(R38*R60 + R39*R63 + R40*R64) - plane_dy;
    double R68 = R0*R67;
    double R69 = -R32*R51 - R33*R52 - R36*R53 - R54*(-R12*R32 + R33*R57 + R36*R58) - R59*(R32*R60 + R33*R63 + R36*R64) + a0 - plane_dx;
    double R70 = R46*R69;
    double R71 = R43*R65 + R44*R67 + R47*R69;
    double R72 = R49*(R48*R71 - R66 + R68 - R70);
    double R73 = R49*(-R48*(R66 - R68 + R70) + R71);
    double R74 = R1*R72 - R43*R73 + R65;
    double R75 = -R0*R72 - R44*R73 + R67;
    double R76 = R46*R72 - R47*R73 + R69;
    double R77 = 2*R49;
    double R78 = R77*(-R46 + R47*R48);
    double R79 = R77*(-R46*R48 + R47);
    double R80 = ((1.0/2.0)*R74*(R1*R78 - R43*R79) + (1.0/2.0)*R75*(-R0*R78 - R44*R79) + (1.0/2.0)*R76*(R46*R78 - R47*R79 + 2))/sqrt(pow(R74, 2) + pow(R75, 2) + pow(R76, 2));
    double R81 = R0*R47 + R44*R46;
    double R82 = R1*R47 - R43*R46;
    double R83 = -R0*R43 - R45*b1;
    double R84 = pow(pow(R81, 2) + pow(R82, 2) + pow(R83, 2), -1.0/2.0);
    double result = ((R65*R81*R84 + R67*R82*R84 + R69*R83*R84 > 0) ? (
   -R80
)
: (
   R80
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_b0(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = 1.0*R5*R6;
    double R8 = sin(plane_b);
    double R9 = sin(panel_g);
    double R10 = cos(panel_a);
    double R11 = R10*R9;
    double R12 = sin(panel_b);
    double R13 = cos(panel_g);
    double R14 = R13*R5;
    double R15 = -R11 + R12*R14;
    double R16 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R17 = 1.0/R16;
    double R18 = R17*panel_straw0x;
    double R19 = R10*R13;
    double R20 = R5*R9;
    double R21 = R12*R20 + R19;
    double R22 = R17*panel_straw0y;
    double R23 = 1.0*R15*R18 - 1.0*R21*R22;
    double R24 = 1.0*R15*R22 + 1.0*R18*R21;
    double R25 = sin(plane_a);
    double R26 = R25*R3;
    double R27 = -R23*R8 + R24*R26 + R4*R7;
    double R28 = sin(plane_g);
    double R29 = R25*R28;
    double R30 = cos(plane_g);
    double R31 = R2*R30;
    double R32 = R29 + R31*R8;
    double R33 = R3*R30;
    double R34 = R2*R28;
    double R35 = R25*R30;
    double R36 = -R34 + R35*R8;
    double R37 = R23*R33 + R24*R36 + R32*R7;
    double R38 = R34*R8 - R35;
    double R39 = R28*R3;
    double R40 = R29*R8 + R31;
    double R41 = R23*R39 + R24*R40 + R38*R7;
    double R42 = pow(pow(R27, 2) + pow(R37, 2) + pow(R41, 2), -1.0/2.0);
    double R43 = R27*R42;
    double R44 = R41*R42;
    double R45 = R0*R44;
    double R46 = R37*R42;
    double R47 = R0*a1;
    double R48 = R1*R43 - R45 + R46*R47;
    double R49 = 1.0/(1.0 - pow(R48, 2));
    double R50 = -plane_z;
    double R51 = R50 + panel_dz + panel_straw0z;
    double R52 = R18*panel_dx - R22*panel_dy + panel_straw0x;
    double R53 = R18*panel_dy + R22*panel_dx + panel_straw0y;
    double R54 = -R16 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R55 = R18*R6;
    double R56 = R22*R6;
    double R57 = R13*R55 - R56*R9;
    double R58 = R13*R56 + R55*R9;
    double R59 = -panel_straw0z + wire_z;
    double R60 = R10*R6;
    double R61 = R12*R19 + R20;
    double R62 = R11*R12 - R14;
    double R63 = R18*R61 - R22*R62;
    double R64 = R18*R62 + R22*R61;
    double R65 = -R26*R53 - R4*R51 + R50 + R52*R8 - R54*(-R12*R4 + R26*R58 - R57*R8) - R59*(R26*R64 + R4*R60 - R63*R8) + b0 - plane_dz;
    double R66 = R1*R65;
    double R67 = -R38*R51 - R39*R52 - R40*R53 - R54*(-R12*R38 + R39*R57 + R40*R58) - R59*(R38*R60 + R39*R63 + R40*R64) - plane_dy;
    double R68 = R0*R67;
    double R69 = -R32*R51 - R33*R52 - R36*R53 - R54*(-R12*R32 + R33*R57 + R36*R58) - R59*(R32*R60 + R33*R63 + R36*R64) + a0 - plane_dx;
    double R70 = R47*R69;
    double R71 = R43*R65 + R44*R67 + R46*R69;
    double R72 = R49*(R48*R71 - R66 + R68 - R70);
    double R73 = R49*(-R48*(R66 - R68 + R70) + R71);
    double R74 = R1*R72 - R43*R73 + R65;
    double R75 = -R0*R72 - R44*R73 + R67;
    double R76 = -R46*R73 + R47*R72 + R69;
    double R77 = 2*R49;
    double R78 = R77*(-R1 + R43*R48);
    double R79 = R77*(-R1*R48 + R43);
    double R80 = ((1.0/2.0)*R74*(R1*R78 - R43*R79 + 2) + (1.0/2.0)*R75*(-R0*R78 - R44*R79) + (1.0/2.0)*R76*(-R46*R79 + R47*R78))/sqrt(pow(R74, 2) + pow(R75, 2) + pow(R76, 2));
    double R81 = R0*R46 + R45*a1;
    double R82 = R1*R46 - R43*R47;
    double R83 = -R0*R43 - R1*R44;
    double R84 = pow(pow(R81, 2) + pow(R82, 2) + pow(R83, 2), -1.0/2.0);
    double result = ((R65*R81*R84 + R67*R82*R84 + R69*R83*R84 > 0) ? (
   -R80
)
: (
   R80
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_a1(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(a1, 2);
    double R1 = R0 + pow(b1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = cos(plane_a);
    double R5 = cos(plane_b);
    double R6 = R4*R5;
    double R7 = sin(panel_a);
    double R8 = cos(panel_b);
    double R9 = 1.0*R7*R8;
    double R10 = sin(plane_b);
    double R11 = sin(panel_g);
    double R12 = cos(panel_a);
    double R13 = R11*R12;
    double R14 = sin(panel_b);
    double R15 = cos(panel_g);
    double R16 = R15*R7;
    double R17 = -R13 + R14*R16;
    double R18 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R19 = 1.0/R18;
    double R20 = R19*panel_straw0x;
    double R21 = R12*R15;
    double R22 = R11*R7;
    double R23 = R14*R22 + R21;
    double R24 = R19*panel_straw0y;
    double R25 = 1.0*R17*R20 - 1.0*R23*R24;
    double R26 = 1.0*R17*R24 + 1.0*R20*R23;
    double R27 = sin(plane_a);
    double R28 = R27*R5;
    double R29 = -R10*R25 + R26*R28 + R6*R9;
    double R30 = sin(plane_g);
    double R31 = R27*R30;
    double R32 = cos(plane_g);
    double R33 = R32*R4;
    double R34 = R10*R33 + R31;
    double R35 = R32*R5;
    double R36 = R30*R4;
    double R37 = R27*R32;
    double R38 = R10*R37 - R36;
    double R39 = R25*R35 + R26*R38 + R34*R9;
    double R40 = R10*R36 - R37;
    double R41 = R30*R5;
    double R42 = R10*R31 + R33;
    double R43 = R25*R41 + R26*R42 + R40*R9;
    double R44 = pow(pow(R29, 2) + pow(R39, 2) + pow(R43, 2), -1.0/2.0);
    double R45 = R29*R44;
    double R46 = R43*R44;
    double R47 = R2*R46;
    double R48 = R39*R44;
    double R49 = R2*R48;
    double R50 = R3*R45 - R47 + R49*a1;
    double R51 = 1.0 - pow(R50, 2);
    double R52 = 1.0/R51;
    double R53 = -plane_z;
    double R54 = R53 + panel_dz + panel_straw0z;
    double R55 = R20*panel_dx - R24*panel_dy + panel_straw0x;
    double R56 = R20*panel_dy + R24*panel_dx + panel_straw0y;
    double R57 = -R18 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R58 = R20*R8;
    double R59 = R24*R8;
    double R60 = -R11*R59 + R15*R58;
    double R61 = R11*R58 + R15*R59;
    double R62 = -panel_straw0z + wire_z;
    double R63 = R12*R8;
    double R64 = R14*R21 + R22;
    double R65 = R13*R14 - R16;
    double R66 = R20*R64 - R24*R65;
    double R67 = R20*R65 + R24*R64;
    double R68 = R10*R55 - R28*R56 + R53 - R54*R6 - R57*(-R10*R60 - R14*R6 + R28*R61) - R62*(-R10*R66 + R28*R67 + R6*R63) + b0 - plane_dz;
    double R69 = -R40*R54 - R41*R55 - R42*R56 - R57*(-R14*R40 + R41*R60 + R42*R61) - R62*(R40*R63 + R41*R66 + R42*R67) - plane_dy;
    double R70 = -R34*R54 - R35*R55 - R38*R56 - R57*(-R14*R34 + R35*R60 + R38*R61) - R62*(R34*R63 + R35*R66 + R38*R67) + a0 - plane_dx;
    double R71 = R45*R68 + R46*R69 + R48*R70;
    double R72 = R3*R68;
    double R73 = R2*R69;
    double R74 = R2*R70;
    double R75 = R74*a1;
    double R76 = -R72 + R73 - R75;
    double R77 = R50*R71 + R76;
    double R78 = R52*R77;
    double R79 = -R50*(R72 - R73 + R75) + R71;
    double R80 = R52*R79;
    double R81 = R3*R78 - R45*R80 + R68;
    double R82 = R2*R78;
    double R83 = -R46*R80 + R69 - R82;
    double R84 = -R48*R80 + R70 + R82*a1;
    double R85 = pow(R1, -3.0/2.0);
    double R86 = R85*a1;
    double R87 = R86*b1;
    double R88 = 2*R78;
    double R89 = R45*R87;
    double R90 = R46*R86;
    double R91 = R0*R85;
    double R92 = R48*R91;
    double R93 = R49 - R89 + R90 - R92;
    double R94 = R68*R87 - R69*R86 + R70*R91 - R74;
    double R95 = 2*R52;
    double R96 = R95*(R71*R93 + R94);
    double R97 = R95*(R50*R94 + R76*R93);
    double R98 = 2*R50*(2*R49 - 2*R89 + 2*R90 - 2*R92)/pow(R51, 2);
    double R99 = R77*R98;
    double R100 = R79*R98;
    double R101 = R2*R96;
    double R102 = R2*R99;
    double R103 = ((1.0/2.0)*R81*(-R100*R45 + R3*R96 + R3*R99 - R45*R97 - R87*R88) + (1.0/2.0)*R83*(-R100*R46 - R101 - R102 - R46*R97 + R86*R88) + (1.0/2.0)*R84*(-R100*R48 + R101*a1 + R102*a1 - R48*R97 + 2*R82 - R88*R91))/sqrt(pow(R81, 2) + pow(R83, 2) + pow(R84, 2));
    double R104 = R47*a1 + R49;
    double R105 = R2*R45;
    double R106 = -R105*a1 + R3*R48;
    double R107 = -R105 - R47*b1;
    double R108 = pow(pow(R104, 2) + pow(R106, 2) + pow(R107, 2), -1.0/2.0);
    double result = ((R104*R108*R68 + R106*R108*R69 + R107*R108*R70 > 0) ? (
   -R103
)
: (
   R103
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_b1(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = cos(plane_a);
    double R4 = cos(plane_b);
    double R5 = R3*R4;
    double R6 = sin(panel_a);
    double R7 = cos(panel_b);
    double R8 = 1.0*R6*R7;
    double R9 = sin(plane_b);
    double R10 = sin(panel_g);
    double R11 = cos(panel_a);
    double R12 = R10*R11;
    double R13 = sin(panel_b);
    double R14 = cos(panel_g);
    double R15 = R14*R6;
    double R16 = -R12 + R13*R15;
    double R17 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R18 = 1.0/R17;
    double R19 = R18*panel_straw0x;
    double R20 = R11*R14;
    double R21 = R10*R6;
    double R22 = R13*R21 + R20;
    double R23 = R18*panel_straw0y;
    double R24 = 1.0*R16*R19 - 1.0*R22*R23;
    double R25 = 1.0*R16*R23 + 1.0*R19*R22;
    double R26 = sin(plane_a);
    double R27 = R26*R4;
    double R28 = -R24*R9 + R25*R27 + R5*R8;
    double R29 = sin(plane_g);
    double R30 = R26*R29;
    double R31 = cos(plane_g);
    double R32 = R3*R31;
    double R33 = R30 + R32*R9;
    double R34 = R31*R4;
    double R35 = R29*R3;
    double R36 = R26*R31;
    double R37 = -R35 + R36*R9;
    double R38 = R24*R34 + R25*R37 + R33*R8;
    double R39 = R35*R9 - R36;
    double R40 = R29*R4;
    double R41 = R30*R9 + R32;
    double R42 = R24*R40 + R25*R41 + R39*R8;
    double R43 = pow(pow(R28, 2) + pow(R38, 2) + pow(R42, 2), -1.0/2.0);
    double R44 = R28*R43;
    double R45 = R2*R44;
    double R46 = R42*R43;
    double R47 = R2*R46;
    double R48 = R38*R43;
    double R49 = R2*a1;
    double R50 = R45*b1 - R47 + R48*R49;
    double R51 = 1.0 - pow(R50, 2);
    double R52 = 1.0/R51;
    double R53 = -plane_z;
    double R54 = R53 + panel_dz + panel_straw0z;
    double R55 = R19*panel_dx - R23*panel_dy + panel_straw0x;
    double R56 = R19*panel_dy + R23*panel_dx + panel_straw0y;
    double R57 = -R17 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R58 = R19*R7;
    double R59 = R23*R7;
    double R60 = -R10*R59 + R14*R58;
    double R61 = R10*R58 + R14*R59;
    double R62 = -panel_straw0z + wire_z;
    double R63 = R11*R7;
    double R64 = R13*R20 + R21;
    double R65 = R12*R13 - R15;
    double R66 = R19*R64 - R23*R65;
    double R67 = R19*R65 + R23*R64;
    double R68 = -R27*R56 - R5*R54 + R53 + R55*R9 - R57*(-R13*R5 + R27*R61 - R60*R9) - R62*(R27*R67 + R5*R63 - R66*R9) + b0 - plane_dz;
    double R69 = -R39*R54 - R40*R55 - R41*R56 - R57*(-R13*R39 + R40*R60 + R41*R61) - R62*(R39*R63 + R40*R66 + R41*R67) - plane_dy;
    double R70 = -R33*R54 - R34*R55 - R37*R56 - R57*(-R13*R33 + R34*R60 + R37*R61) - R62*(R33*R63 + R34*R66 + R37*R67) + a0 - plane_dx;
    double R71 = R44*R68 + R46*R69 + R48*R70;
    double R72 = R2*R68;
    double R73 = R72*b1;
    double R74 = R2*R69;
    double R75 = R49*R70;
    double R76 = -R73 + R74 - R75;
    double R77 = R50*R71 + R76;
    double R78 = R52*R77;
    double R79 = R2*R78;
    double R80 = -R50*(R73 - R74 + R75) + R71;
    double R81 = R52*R80;
    double R82 = -R44*R81 + R68 + R79*b1;
    double R83 = -R46*R81 + R69 - R79;
    double R84 = -R48*R81 + R49*R78 + R70;
    double R85 = pow(R1, -3.0/2.0);
    double R86 = R85*b1;
    double R87 = 2*R78;
    double R88 = R0*R85;
    double R89 = R44*R88;
    double R90 = R46*R86;
    double R91 = R86*a1;
    double R92 = R48*R91;
    double R93 = R45 - R89 + R90 - R92;
    double R94 = R68*R88 - R69*R86 + R70*R91 - R72;
    double R95 = 2*R52;
    double R96 = R95*(R71*R93 + R94);
    double R97 = R2*R96;
    double R98 = R95*(R50*R94 + R76*R93);
    double R99 = 2*R50*(2*R45 - 2*R89 + 2*R90 - 2*R92)/pow(R51, 2);
    double R100 = R77*R99;
    double R101 = R100*R2;
    double R102 = R80*R99;
    double R103 = ((1.0/2.0)*R82*(R101*b1 - R102*R44 - R44*R98 + 2*R79 - R87*R88 + R97*b1) + (1.0/2.0)*R83*(-R101 - R102*R46 - R46*R98 + R86*R87 - R97) + (1.0/2.0)*R84*(R100*R49 - R102*R48 - R48*R98 + R49*R96 - R87*R91))/sqrt(pow(R82, 2) + pow(R83, 2) + pow(R84, 2));
    double R104 = R2*R48;
    double R105 = R104 + R47*a1;
    double R106 = R104*b1 - R44*R49;
    double R107 = -R45 - R47*b1;
    double R108 = pow(pow(R105, 2) + pow(R106, 2) + pow(R107, 2), -1.0/2.0);
    double result = ((R105*R108*R68 + R106*R108*R69 + R107*R108*R70 > 0) ? (
   -R103
)
: (
   R103
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_t0(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double result = -1;
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dx(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = 1.0*R5*R6;
    double R8 = sin(plane_b);
    double R9 = sin(panel_g);
    double R10 = cos(panel_a);
    double R11 = R10*R9;
    double R12 = sin(panel_b);
    double R13 = cos(panel_g);
    double R14 = R13*R5;
    double R15 = -R11 + R12*R14;
    double R16 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R17 = 1.0/R16;
    double R18 = R17*panel_straw0x;
    double R19 = R10*R13;
    double R20 = R5*R9;
    double R21 = R12*R20 + R19;
    double R22 = R17*panel_straw0y;
    double R23 = 1.0*R15*R18 - 1.0*R21*R22;
    double R24 = 1.0*R15*R22 + 1.0*R18*R21;
    double R25 = sin(plane_a);
    double R26 = R25*R3;
    double R27 = -R23*R8 + R24*R26 + R4*R7;
    double R28 = sin(plane_g);
    double R29 = R25*R28;
    double R30 = cos(plane_g);
    double R31 = R2*R30;
    double R32 = R29 + R31*R8;
    double R33 = R3*R30;
    double R34 = R2*R28;
    double R35 = R25*R30;
    double R36 = -R34 + R35*R8;
    double R37 = R23*R33 + R24*R36 + R32*R7;
    double R38 = R34*R8 - R35;
    double R39 = R28*R3;
    double R40 = R29*R8 + R31;
    double R41 = R23*R39 + R24*R40 + R38*R7;
    double R42 = pow(pow(R27, 2) + pow(R37, 2) + pow(R41, 2), -1.0/2.0);
    double R43 = R27*R42;
    double R44 = R41*R42;
    double R45 = R0*R44;
    double R46 = R0*a1;
    double R47 = R37*R42;
    double R48 = R1*R43 - R45 + R46*R47;
    double R49 = 1.0/(1.0 - pow(R48, 2));
    double R50 = -plane_z;
    double R51 = R50 + panel_dz + panel_straw0z;
    double R52 = R18*panel_dx - R22*panel_dy + panel_straw0x;
    double R53 = R18*panel_dy + R22*panel_dx + panel_straw0y;
    double R54 = -R16 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R55 = R18*R6;
    double R56 = R22*R6;
    double R57 = R13*R55 - R56*R9;
    double R58 = R13*R56 + R55*R9;
    double R59 = -panel_straw0z + wire_z;
    double R60 = R10*R6;
    double R61 = R12*R19 + R20;
    double R62 = R11*R12 - R14;
    double R63 = R18*R61 - R22*R62;
    double R64 = R18*R62 + R22*R61;
    double R65 = -R26*R53 - R4*R51 + R50 + R52*R8 - R54*(-R12*R4 + R26*R58 - R57*R8) - R59*(R26*R64 + R4*R60 - R63*R8) + b0 - plane_dz;
    double R66 = R1*R65;
    double R67 = -R38*R51 - R39*R52 - R40*R53 - R54*(-R12*R38 + R39*R57 + R40*R58) - R59*(R38*R60 + R39*R63 + R40*R64) - plane_dy;
    double R68 = R0*R67;
    double R69 = -R32*R51 - R33*R52 - R36*R53 - R54*(-R12*R32 + R33*R57 + R36*R58) - R59*(R32*R60 + R33*R63 + R36*R64) + a0 - plane_dx;
    double R70 = R46*R69;
    double R71 = R43*R65 + R44*R67 + R47*R69;
    double R72 = R49*(R48*R71 - R66 + R68 - R70);
    double R73 = R49*(-R48*(R66 - R68 + R70) + R71);
    double R74 = R1*R72 - R43*R73 + R65;
    double R75 = -R0*R72 - R44*R73 + R67;
    double R76 = R46*R72 - R47*R73 + R69;
    double R77 = 2*R49;
    double R78 = R77*(R46 - R47*R48);
    double R79 = R77*(R46*R48 - R47);
    double R80 = ((1.0/2.0)*R74*(R1*R78 - R43*R79) + (1.0/2.0)*R75*(-R0*R78 - R44*R79) + (1.0/2.0)*R76*(R46*R78 - R47*R79 - 2))/sqrt(pow(R74, 2) + pow(R75, 2) + pow(R76, 2));
    double R81 = R0*R47 + R44*R46;
    double R82 = R1*R47 - R43*R46;
    double R83 = -R0*R43 - R45*b1;
    double R84 = pow(pow(R81, 2) + pow(R82, 2) + pow(R83, 2), -1.0/2.0);
    double result = ((R65*R81*R84 + R67*R82*R84 + R69*R83*R84 > 0) ? (
   -R80
)
: (
   R80
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dy(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = 1.0*R5*R6;
    double R8 = sin(plane_b);
    double R9 = sin(panel_g);
    double R10 = cos(panel_a);
    double R11 = R10*R9;
    double R12 = sin(panel_b);
    double R13 = cos(panel_g);
    double R14 = R13*R5;
    double R15 = -R11 + R12*R14;
    double R16 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R17 = 1.0/R16;
    double R18 = R17*panel_straw0x;
    double R19 = R10*R13;
    double R20 = R5*R9;
    double R21 = R12*R20 + R19;
    double R22 = R17*panel_straw0y;
    double R23 = 1.0*R15*R18 - 1.0*R21*R22;
    double R24 = 1.0*R15*R22 + 1.0*R18*R21;
    double R25 = sin(plane_a);
    double R26 = R25*R3;
    double R27 = -R23*R8 + R24*R26 + R4*R7;
    double R28 = sin(plane_g);
    double R29 = R25*R28;
    double R30 = cos(plane_g);
    double R31 = R2*R30;
    double R32 = R29 + R31*R8;
    double R33 = R3*R30;
    double R34 = R2*R28;
    double R35 = R25*R30;
    double R36 = -R34 + R35*R8;
    double R37 = R23*R33 + R24*R36 + R32*R7;
    double R38 = R34*R8 - R35;
    double R39 = R28*R3;
    double R40 = R29*R8 + R31;
    double R41 = R23*R39 + R24*R40 + R38*R7;
    double R42 = pow(pow(R27, 2) + pow(R37, 2) + pow(R41, 2), -1.0/2.0);
    double R43 = R27*R42;
    double R44 = R41*R42;
    double R45 = R0*R44;
    double R46 = R37*R42;
    double R47 = R0*a1;
    double R48 = R1*R43 - R45 + R46*R47;
    double R49 = 1.0/(1.0 - pow(R48, 2));
    double R50 = -plane_z;
    double R51 = R50 + panel_dz + panel_straw0z;
    double R52 = R18*panel_dx - R22*panel_dy + panel_straw0x;
    double R53 = R18*panel_dy + R22*panel_dx + panel_straw0y;
    double R54 = -R16 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R55 = R18*R6;
    double R56 = R22*R6;
    double R57 = R13*R55 - R56*R9;
    double R58 = R13*R56 + R55*R9;
    double R59 = -panel_straw0z + wire_z;
    double R60 = R10*R6;
    double R61 = R12*R19 + R20;
    double R62 = R11*R12 - R14;
    double R63 = R18*R61 - R22*R62;
    double R64 = R18*R62 + R22*R61;
    double R65 = -R26*R53 - R4*R51 + R50 + R52*R8 - R54*(-R12*R4 + R26*R58 - R57*R8) - R59*(R26*R64 + R4*R60 - R63*R8) + b0 - plane_dz;
    double R66 = R1*R65;
    double R67 = -R38*R51 - R39*R52 - R40*R53 - R54*(-R12*R38 + R39*R57 + R40*R58) - R59*(R38*R60 + R39*R63 + R40*R64) - plane_dy;
    double R68 = R0*R67;
    double R69 = -R32*R51 - R33*R52 - R36*R53 - R54*(-R12*R32 + R33*R57 + R36*R58) - R59*(R32*R60 + R33*R63 + R36*R64) + a0 - plane_dx;
    double R70 = R47*R69;
    double R71 = R43*R65 + R44*R67 + R46*R69;
    double R72 = R49*(R48*R71 - R66 + R68 - R70);
    double R73 = R49*(-R48*(R66 - R68 + R70) + R71);
    double R74 = R1*R72 - R43*R73 + R65;
    double R75 = -R0*R72 - R44*R73 + R67;
    double R76 = -R46*R73 + R47*R72 + R69;
    double R77 = 2*R49;
    double R78 = R77*(-R0 - R44*R48);
    double R79 = R77*(-R0*R48 - R44);
    double R80 = ((1.0/2.0)*R74*(R1*R78 - R43*R79) + (1.0/2.0)*R75*(-R0*R78 - R44*R79 - 2) + (1.0/2.0)*R76*(-R46*R79 + R47*R78))/sqrt(pow(R74, 2) + pow(R75, 2) + pow(R76, 2));
    double R81 = R0*R46 + R45*a1;
    double R82 = R1*R46 - R43*R47;
    double R83 = -R0*R43 - R45*b1;
    double R84 = pow(pow(R81, 2) + pow(R82, 2) + pow(R83, 2), -1.0/2.0);
    double result = ((R65*R81*R84 + R67*R82*R84 + R69*R83*R84 > 0) ? (
   -R80
)
: (
   R80
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dz(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = 1.0*R5*R6;
    double R8 = sin(plane_b);
    double R9 = sin(panel_g);
    double R10 = cos(panel_a);
    double R11 = R10*R9;
    double R12 = sin(panel_b);
    double R13 = cos(panel_g);
    double R14 = R13*R5;
    double R15 = -R11 + R12*R14;
    double R16 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R17 = 1.0/R16;
    double R18 = R17*panel_straw0x;
    double R19 = R10*R13;
    double R20 = R5*R9;
    double R21 = R12*R20 + R19;
    double R22 = R17*panel_straw0y;
    double R23 = 1.0*R15*R18 - 1.0*R21*R22;
    double R24 = 1.0*R15*R22 + 1.0*R18*R21;
    double R25 = sin(plane_a);
    double R26 = R25*R3;
    double R27 = -R23*R8 + R24*R26 + R4*R7;
    double R28 = sin(plane_g);
    double R29 = R25*R28;
    double R30 = cos(plane_g);
    double R31 = R2*R30;
    double R32 = R29 + R31*R8;
    double R33 = R3*R30;
    double R34 = R2*R28;
    double R35 = R25*R30;
    double R36 = -R34 + R35*R8;
    double R37 = R23*R33 + R24*R36 + R32*R7;
    double R38 = R34*R8 - R35;
    double R39 = R28*R3;
    double R40 = R29*R8 + R31;
    double R41 = R23*R39 + R24*R40 + R38*R7;
    double R42 = pow(pow(R27, 2) + pow(R37, 2) + pow(R41, 2), -1.0/2.0);
    double R43 = R27*R42;
    double R44 = R41*R42;
    double R45 = R0*R44;
    double R46 = R37*R42;
    double R47 = R0*a1;
    double R48 = R1*R43 - R45 + R46*R47;
    double R49 = 1.0/(1.0 - pow(R48, 2));
    double R50 = -plane_z;
    double R51 = R50 + panel_dz + panel_straw0z;
    double R52 = R18*panel_dx - R22*panel_dy + panel_straw0x;
    double R53 = R18*panel_dy + R22*panel_dx + panel_straw0y;
    double R54 = -R16 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R55 = R18*R6;
    double R56 = R22*R6;
    double R57 = R13*R55 - R56*R9;
    double R58 = R13*R56 + R55*R9;
    double R59 = -panel_straw0z + wire_z;
    double R60 = R10*R6;
    double R61 = R12*R19 + R20;
    double R62 = R11*R12 - R14;
    double R63 = R18*R61 - R22*R62;
    double R64 = R18*R62 + R22*R61;
    double R65 = -R26*R53 - R4*R51 + R50 + R52*R8 - R54*(-R12*R4 + R26*R58 - R57*R8) - R59*(R26*R64 + R4*R60 - R63*R8) + b0 - plane_dz;
    double R66 = R1*R65;
    double R67 = -R38*R51 - R39*R52 - R40*R53 - R54*(-R12*R38 + R39*R57 + R40*R58) - R59*(R38*R60 + R39*R63 + R40*R64) - plane_dy;
    double R68 = R0*R67;
    double R69 = -R32*R51 - R33*R52 - R36*R53 - R54*(-R12*R32 + R33*R57 + R36*R58) - R59*(R32*R60 + R33*R63 + R36*R64) + a0 - plane_dx;
    double R70 = R47*R69;
    double R71 = R43*R65 + R44*R67 + R46*R69;
    double R72 = R49*(R48*R71 - R66 + R68 - R70);
    double R73 = R49*(-R48*(R66 - R68 + R70) + R71);
    double R74 = R1*R72 - R43*R73 + R65;
    double R75 = -R0*R72 - R44*R73 + R67;
    double R76 = -R46*R73 + R47*R72 + R69;
    double R77 = 2*R49;
    double R78 = R77*(R1 - R43*R48);
    double R79 = R77*(R1*R48 - R43);
    double R80 = ((1.0/2.0)*R74*(R1*R78 - R43*R79 - 2) + (1.0/2.0)*R75*(-R0*R78 - R44*R79) + (1.0/2.0)*R76*(-R46*R79 + R47*R78))/sqrt(pow(R74, 2) + pow(R75, 2) + pow(R76, 2));
    double R81 = R0*R46 + R45*a1;
    double R82 = R1*R46 - R43*R47;
    double R83 = -R0*R43 - R1*R44;
    double R84 = pow(pow(R81, 2) + pow(R82, 2) + pow(R83, 2), -1.0/2.0);
    double result = ((R65*R81*R84 + R67*R82*R84 + R69*R83*R84 > 0) ? (
   -R80
)
: (
   R80
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_plane_a(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = R5*R6;
    double R8 = 1.0*R7;
    double R9 = sin(plane_b);
    double R10 = sin(panel_g);
    double R11 = cos(panel_a);
    double R12 = R10*R11;
    double R13 = sin(panel_b);
    double R14 = cos(panel_g);
    double R15 = R14*R5;
    double R16 = -R12 + R13*R15;
    double R17 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R18 = 1.0/R17;
    double R19 = R18*panel_straw0x;
    double R20 = R11*R14;
    double R21 = R10*R5;
    double R22 = R13*R21 + R20;
    double R23 = R18*panel_straw0y;
    double R24 = 1.0*R16*R19 - 1.0*R22*R23;
    double R25 = R19*R22;
    double R26 = R16*R23;
    double R27 = R25 + R26;
    double R28 = 1.0*R27;
    double R29 = sin(plane_a);
    double R30 = R29*R3;
    double R31 = -R24*R9 + R28*R30 + R4*R8;
    double R32 = sin(plane_g);
    double R33 = R29*R32;
    double R34 = cos(plane_g);
    double R35 = R2*R34;
    double R36 = R35*R9;
    double R37 = R33 + R36;
    double R38 = R3*R34;
    double R39 = R2*R32;
    double R40 = R29*R34;
    double R41 = R40*R9;
    double R42 = -R39 + R41;
    double R43 = R24*R38 + R28*R42 + R37*R8;
    double R44 = R39*R9;
    double R45 = -R40 + R44;
    double R46 = R3*R32;
    double R47 = R33*R9;
    double R48 = R35 + R47;
    double R49 = R24*R46 + R28*R48 + R45*R8;
    double R50 = pow(R31, 2) + pow(R43, 2) + pow(R49, 2);
    double R51 = pow(R50, -1.0/2.0);
    double R52 = R31*R51;
    double R53 = R49*R51;
    double R54 = R0*R53;
    double R55 = R0*a1;
    double R56 = R43*R51;
    double R57 = R1*R52 - R54 + R55*R56;
    double R58 = 1.0 - pow(R57, 2);
    double R59 = 1.0/R58;
    double R60 = -plane_z;
    double R61 = R60 + panel_dz + panel_straw0z;
    double R62 = R19*panel_dx - R23*panel_dy + panel_straw0x;
    double R63 = R23*panel_dx;
    double R64 = R19*panel_dy;
    double R65 = R63 + R64 + panel_straw0y;
    double R66 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R67 = -R17 + R66;
    double R68 = R19*R6;
    double R69 = R23*R6;
    double R70 = -R10*R69 + R14*R68;
    double R71 = R10*R68 + R14*R69;
    double R72 = -panel_straw0z;
    double R73 = R72 + wire_z;
    double R74 = R11*R6;
    double R75 = R13*R20 + R21;
    double R76 = R12*R13 - R15;
    double R77 = R19*R75 - R23*R76;
    double R78 = R19*R76 + R23*R75;
    double R79 = -R30*R65 - R4*R61 + R60 + R62*R9 - R67*(-R13*R4 + R30*R71 - R70*R9) - R73*(R30*R78 + R4*R74 - R77*R9) + b0 - plane_dz;
    double R80 = R51*R79;
    double R81 = -R45*R61 - R46*R62 - R48*R65 - R67*(-R13*R45 + R46*R70 + R48*R71) - R73*(R45*R74 + R46*R77 + R48*R78) - plane_dy;
    double R82 = -R37*R61 - R38*R62 - R42*R65 - R67*(-R13*R37 + R38*R70 + R42*R71) - R73*(R37*R74 + R38*R77 + R42*R78) + a0 - plane_dx;
    double R83 = R51*R82;
    double R84 = R31*R80 + R43*R83 + R53*R81;
    double R85 = R1*R79;
    double R86 = R0*R81;
    double R87 = R55*R82;
    double R88 = -R85 + R86 - R87;
    double R89 = R57*R84 + R88;
    double R90 = R59*R89;
    double R91 = -R57*(R85 - R86 + R87) + R84;
    double R92 = R59*R91;
    double R93 = R1*R90 - R52*R92 + R79;
    double R94 = -R0*R90 - R53*R92 + R81;
    double R95 = R55*R90 - R56*R92 + R82;
    double R96 = R30*(R72 - panel_dz + plane_z);
    double R97 = R4*(-R63 - R64 - panel_straw0y);
    double R98 = R17 - R66;
    double R99 = R98*(R13*R30 + R4*R71);
    double R100 = panel_straw0z - wire_z;
    double R101 = R100*(-R30*R74 + R4*R78);
    double R102 = R4*(1.0*R25 + 1.0*R26);
    double R103 = R102 - R30*R8;
    double R104 = R103*R51;
    double R105 = 2*R92;
    double R106 = R7*(1.0*R39 - 1.0*R41);
    double R107 = R27*(1.0*R33 + 1.0*R36);
    double R108 = R7*(-1.0*R35 - 1.0*R47);
    double R109 = R27*(-1.0*R40 + 1.0*R44);
    double R110 = (-1.0/2.0*R31*(2*R102 - 2.0*R30*R7) - 1.0/2.0*R43*(2*R106 + 2*R107) - 1.0/2.0*R49*(2*R108 + 2*R109))/pow(R50, 3.0/2.0);
    double R111 = R110*R31;
    double R112 = R1*R104;
    double R113 = R51*(R108 + R109);
    double R114 = R0*R113;
    double R115 = R106 + R107;
    double R116 = R115*R51;
    double R117 = R116*R55;
    double R118 = R1*R111;
    double R119 = R110*R49;
    double R120 = R0*R119;
    double R121 = R110*R43;
    double R122 = R121*R55;
    double R123 = 2*R57*(2*R112 - 2*R114 + 2*R117 + 2*R118 - 2*R120 + 2*R122)/pow(R58, 2);
    double R124 = R123*R89;
    double R125 = R123*R91;
    double R126 = R101 - R96 + R97 + R99;
    double R127 = R48*R61;
    double R128 = R65*(R40 - R44);
    double R129 = R98*(R13*R48 + R45*R71);
    double R130 = R100*(R45*R78 + R74*(-R35 - R47));
    double R131 = R127 + R128 + R129 + R130;
    double R132 = R42*R61;
    double R133 = R65*(-R33 - R36);
    double R134 = R98*(R13*R42 + R37*R71);
    double R135 = R100*(R37*R78 + R74*(R39 - R41));
    double R136 = R132 + R133 + R134 + R135;
    double R137 = R0*R131 - R1*R126 - R136*R55;
    double R138 = R112 - R114 + R117 + R118 - R120 + R122;
    double R139 = R103*R80 + R111*R79 + R113*R81 + R115*R83 + R119*R81 + R121*R82 + R126*R52 + R131*R53 + R136*R56;
    double R140 = 2*R59;
    double R141 = R140*(R137*R57 + R138*R88 + R139);
    double R142 = R140*(R137 + R138*R84 + R139*R57);
    double R143 = ((1.0/2.0)*R93*(R1*R124 + R1*R142 + 2*R101 - R104*R105 - R105*R111 - R125*R52 - R141*R52 - 2*R96 + 2*R97 + 2*R99) + (1.0/2.0)*R94*(-R0*R124 - R0*R142 - R105*R113 - R105*R119 - R125*R53 + 2*R127 + 2*R128 + 2*R129 + 2*R130 - R141*R53) + (1.0/2.0)*R95*(-R105*R116 - R105*R121 + R124*R55 - R125*R56 + 2*R132 + 2*R133 + 2*R134 + 2*R135 - R141*R56 + R142*R55))/sqrt(pow(R93, 2) + pow(R94, 2) + pow(R95, 2));
    double R144 = R0*R56 + R54*a1;
    double R145 = R1*R56 - R52*R55;
    double R146 = -R0*R52 - R54*b1;
    double R147 = pow(pow(R144, 2) + pow(R145, 2) + pow(R146, 2), -1.0/2.0);
    double result = ((R144*R147*R79 + R145*R147*R81 + R146*R147*R82 > 0) ? (
   -R143
)
: (
   R143
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_plane_b(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = R5*R6;
    double R8 = 1.0*R7;
    double R9 = sin(plane_b);
    double R10 = sin(panel_g);
    double R11 = cos(panel_a);
    double R12 = R10*R11;
    double R13 = sin(panel_b);
    double R14 = cos(panel_g);
    double R15 = R14*R5;
    double R16 = -R12 + R13*R15;
    double R17 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R18 = 1.0/R17;
    double R19 = R18*panel_straw0x;
    double R20 = R16*R19;
    double R21 = R11*R14;
    double R22 = R10*R5;
    double R23 = R13*R22 + R21;
    double R24 = R18*panel_straw0y;
    double R25 = R23*R24;
    double R26 = 1.0*R20 - 1.0*R25;
    double R27 = R19*R23;
    double R28 = R16*R24;
    double R29 = R27 + R28;
    double R30 = 1.0*R29;
    double R31 = sin(plane_a);
    double R32 = R3*R31;
    double R33 = -R26*R9 + R30*R32 + R4*R8;
    double R34 = sin(plane_g);
    double R35 = R31*R34;
    double R36 = cos(plane_g);
    double R37 = R2*R36;
    double R38 = R35 + R37*R9;
    double R39 = R3*R36;
    double R40 = R2*R34;
    double R41 = R31*R36;
    double R42 = -R40 + R41*R9;
    double R43 = R26*R39 + R30*R42 + R38*R8;
    double R44 = R40*R9 - R41;
    double R45 = R3*R34;
    double R46 = R35*R9 + R37;
    double R47 = R26*R45 + R30*R46 + R44*R8;
    double R48 = pow(R33, 2) + pow(R43, 2) + pow(R47, 2);
    double R49 = pow(R48, -1.0/2.0);
    double R50 = R33*R49;
    double R51 = R47*R49;
    double R52 = R0*R51;
    double R53 = R43*R49;
    double R54 = R0*a1;
    double R55 = R1*R50 - R52 + R53*R54;
    double R56 = 1.0 - pow(R55, 2);
    double R57 = 1.0/R56;
    double R58 = -plane_z;
    double R59 = R58 + panel_dz + panel_straw0z;
    double R60 = R19*panel_dx;
    double R61 = R24*panel_dy;
    double R62 = R60 - R61 + panel_straw0x;
    double R63 = R24*panel_dx;
    double R64 = R19*panel_dy;
    double R65 = R63 + R64 + panel_straw0y;
    double R66 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R67 = -R17 + R66;
    double R68 = R19*R6;
    double R69 = R14*R68;
    double R70 = R24*R6;
    double R71 = R10*R70;
    double R72 = R69 - R71;
    double R73 = R72*R9;
    double R74 = R10*R68 + R14*R70;
    double R75 = -panel_straw0z;
    double R76 = R75 + wire_z;
    double R77 = R11*R6;
    double R78 = R13*R21 + R22;
    double R79 = R19*R78;
    double R80 = R12*R13 - R15;
    double R81 = R24*R80;
    double R82 = R79 - R81;
    double R83 = R82*R9;
    double R84 = R19*R80 + R24*R78;
    double R85 = -R32*R65 - R4*R59 + R58 + R62*R9 - R67*(-R13*R4 + R32*R74 - R73) - R76*(R32*R84 + R4*R77 - R83) + b0 - plane_dz;
    double R86 = R49*R85;
    double R87 = R3*R62;
    double R88 = -R34*R87 - R44*R59 - R46*R65 - R67*(-R13*R44 + R45*R72 + R46*R74) - R76*(R44*R77 + R45*R82 + R46*R84) - plane_dy;
    double R89 = -R36*R87 - R38*R59 - R42*R65 - R67*(-R13*R38 + R39*R72 + R42*R74) - R76*(R38*R77 + R39*R82 + R42*R84) + a0 - plane_dx;
    double R90 = R33*R86 + R51*R88 + R53*R89;
    double R91 = R1*R85;
    double R92 = R0*R88;
    double R93 = R54*R89;
    double R94 = -R91 + R92 - R93;
    double R95 = R55*R90 + R94;
    double R96 = R57*R95;
    double R97 = -R55*(R91 - R92 + R93) + R90;
    double R98 = R57*R97;
    double R99 = R1*R96 - R50*R98 + R85;
    double R100 = -R0*R96 - R51*R98 + R88;
    double R101 = -R53*R98 + R54*R96 + R89;
    double R102 = R2*R9;
    double R103 = R102*(R75 - panel_dz + plane_z);
    double R104 = R31*R9;
    double R105 = R104*(-R63 - R64 - panel_straw0y);
    double R106 = R17 - R66;
    double R107 = R106*(R102*R13 - R104*R74 + R3*(-R69 + R71));
    double R108 = panel_straw0z - wire_z;
    double R109 = R108*(-R102*R77 - R104*R84 + R3*(-R79 + R81));
    double R110 = 1.0*R25;
    double R111 = 1.0*R20;
    double R112 = R3*(R110 - R111);
    double R113 = R104*(1.0*R27 + 1.0*R28);
    double R114 = -R102*R8 + R112 - R113;
    double R115 = R114*R49;
    double R116 = 2*R98;
    double R117 = 2.0*R7;
    double R118 = R3*R37;
    double R119 = -R110 + R111;
    double R120 = R36*R9;
    double R121 = R119*R120;
    double R122 = R3*R41;
    double R123 = 2.0*R29;
    double R124 = R3*R40;
    double R125 = R34*R9;
    double R126 = R119*R125;
    double R127 = R3*R35;
    double R128 = (-1.0/2.0*R33*(-R102*R117 + 2*R112 - 2*R113) - 1.0/2.0*R43*(R117*R118 - 2*R121 + R122*R123) - 1.0/2.0*R47*(R117*R124 + R123*R127 - 2*R126))/pow(R48, 3.0/2.0);
    double R129 = R128*R33;
    double R130 = R1*R115;
    double R131 = R49*(R124*R8 - R126 + R127*R30);
    double R132 = R0*R131;
    double R133 = R49*(R118*R8 - R121 + R122*R30);
    double R134 = R133*R54;
    double R135 = R1*R129;
    double R136 = R128*R47;
    double R137 = R0*R136;
    double R138 = R128*R43;
    double R139 = R138*R54;
    double R140 = 2*R55*(2*R130 - 2*R132 + 2*R134 + 2*R135 - 2*R137 + 2*R139)/pow(R56, 2);
    double R141 = R140*R95;
    double R142 = R140*R97;
    double R143 = -R103 - R105 + R107 + R109 + R87;
    double R144 = R124*R59;
    double R145 = R127*R65;
    double R146 = -R60 + R61 - panel_straw0x;
    double R147 = R125*R146;
    double R148 = R106*(-R124*R13 + R127*R74 - R34*R73);
    double R149 = R108*(R124*R77 + R127*R84 - R34*R83);
    double R150 = -R144 - R145 - R147 + R148 + R149;
    double R151 = R118*R59;
    double R152 = R122*R65;
    double R153 = R120*R146;
    double R154 = R106*(-R118*R13 + R122*R74 - R36*R73);
    double R155 = R108*(R118*R77 + R122*R84 - R36*R83);
    double R156 = -R151 - R152 - R153 + R154 + R155;
    double R157 = R0*R150 - R1*R143 - R156*R54;
    double R158 = R130 - R132 + R134 + R135 - R137 + R139;
    double R159 = R114*R86 + R129*R85 + R131*R88 + R133*R89 + R136*R88 + R138*R89 + R143*R50 + R150*R51 + R156*R53;
    double R160 = 2*R57;
    double R161 = R160*(R157*R55 + R158*R94 + R159);
    double R162 = R160*(R157 + R158*R90 + R159*R55);
    double R163 = ((1.0/2.0)*R100*(-R0*R141 - R0*R162 - R116*R131 - R116*R136 - R142*R51 - 2*R144 - 2*R145 - 2*R147 + 2*R148 + 2*R149 - R161*R51) + (1.0/2.0)*R101*(-R116*R133 - R116*R138 + R141*R54 - R142*R53 - 2*R151 - 2*R152 - 2*R153 + 2*R154 + 2*R155 - R161*R53 + R162*R54) + (1.0/2.0)*R99*(R1*R141 + R1*R162 - 2*R103 - 2*R105 + 2*R107 + 2*R109 - R115*R116 - R116*R129 - R142*R50 - R161*R50 + 2*R87))/sqrt(pow(R100, 2) + pow(R101, 2) + pow(R99, 2));
    double R164 = R0*R53 + R52*a1;
    double R165 = R1*R53 - R50*R54;
    double R166 = -R0*R50 - R52*b1;
    double R167 = pow(pow(R164, 2) + pow(R165, 2) + pow(R166, 2), -1.0/2.0);
    double result = ((R164*R167*R85 + R165*R167*R88 + R166*R167*R89 > 0) ? (
   -R163
)
: (
   R163
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_plane_g(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = R5*R6;
    double R8 = 1.0*R7;
    double R9 = sin(plane_b);
    double R10 = sin(panel_g);
    double R11 = cos(panel_a);
    double R12 = R10*R11;
    double R13 = sin(panel_b);
    double R14 = cos(panel_g);
    double R15 = R14*R5;
    double R16 = -R12 + R13*R15;
    double R17 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R18 = 1.0/R17;
    double R19 = R18*panel_straw0x;
    double R20 = R16*R19;
    double R21 = R11*R14;
    double R22 = R10*R5;
    double R23 = R13*R22 + R21;
    double R24 = R18*panel_straw0y;
    double R25 = R23*R24;
    double R26 = 1.0*R20 - 1.0*R25;
    double R27 = R16*R24 + R19*R23;
    double R28 = 1.0*R27;
    double R29 = sin(plane_a);
    double R30 = R29*R3;
    double R31 = -R26*R9 + R28*R30 + R4*R8;
    double R32 = sin(plane_g);
    double R33 = R29*R32;
    double R34 = cos(plane_g);
    double R35 = R2*R34;
    double R36 = R35*R9;
    double R37 = R33 + R36;
    double R38 = R3*R34;
    double R39 = R2*R32;
    double R40 = R29*R34;
    double R41 = R40*R9;
    double R42 = -R39 + R41;
    double R43 = R26*R38 + R28*R42 + R37*R8;
    double R44 = R39*R9;
    double R45 = -R40 + R44;
    double R46 = R3*R32;
    double R47 = R33*R9;
    double R48 = R35 + R47;
    double R49 = R26*R46 + R28*R48 + R45*R8;
    double R50 = pow(R31, 2) + pow(R43, 2) + pow(R49, 2);
    double R51 = pow(R50, -1.0/2.0);
    double R52 = R31*R51;
    double R53 = R49*R51;
    double R54 = R0*R53;
    double R55 = R0*a1;
    double R56 = R43*R51;
    double R57 = R1*R52 - R54 + R55*R56;
    double R58 = 1.0 - pow(R57, 2);
    double R59 = 1.0/R58;
    double R60 = -plane_z;
    double R61 = R60 + panel_dz + panel_straw0z;
    double R62 = R19*panel_dx;
    double R63 = R24*panel_dy;
    double R64 = R62 - R63 + panel_straw0x;
    double R65 = R19*panel_dy + R24*panel_dx + panel_straw0y;
    double R66 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R67 = -R17 + R66;
    double R68 = R19*R6;
    double R69 = R24*R6;
    double R70 = -R10*R69 + R14*R68;
    double R71 = R10*R68 + R14*R69;
    double R72 = -panel_straw0z + wire_z;
    double R73 = R11*R6;
    double R74 = R13*R21 + R22;
    double R75 = R12*R13 - R15;
    double R76 = R19*R74 - R24*R75;
    double R77 = R19*R75 + R24*R74;
    double R78 = -R30*R65 - R4*R61 + R60 + R64*R9 - R67*(-R13*R4 + R30*R71 - R70*R9) - R72*(R30*R77 + R4*R73 - R76*R9) + b0 - plane_dz;
    double R79 = R45*R61;
    double R80 = R48*R65;
    double R81 = R13*R45;
    double R82 = R46*R70;
    double R83 = R46*R76;
    double R84 = -R46*R64 - R67*(R48*R71 - R81 + R82) - R72*(R45*R73 + R48*R77 + R83) - R79 - R80 - plane_dy;
    double R85 = R38*R70 + R42*R71;
    double R86 = R37*R73 + R38*R76 + R42*R77;
    double R87 = -R37*R61 - R38*R64 - R42*R65 - R67*(-R13*R37 + R85) - R72*R86 + a0 - plane_dx;
    double R88 = R51*R87;
    double R89 = R43*R88 + R52*R78 + R53*R84;
    double R90 = R1*R78;
    double R91 = R0*R84;
    double R92 = R55*R87;
    double R93 = -R90 + R91 - R92;
    double R94 = R57*R89 + R93;
    double R95 = R59*R94;
    double R96 = -R57*(R90 - R91 + R92) + R89;
    double R97 = R59*R96;
    double R98 = R1*R95 - R52*R97 + R78;
    double R99 = -R0*R95 - R53*R97 + R84;
    double R100 = R55*R95 - R56*R97 + R87;
    double R101 = R7*(1.0*R40 - 1.0*R44);
    double R102 = 1.0*R20 - 1.0*R25;
    double R103 = R102*R46;
    double R104 = R27*(-1.0*R35 - 1.0*R47);
    double R105 = R7*(1.0*R33 + 1.0*R36);
    double R106 = R102*R38;
    double R107 = R27*(-1.0*R39 + 1.0*R41);
    double R108 = (-1.0/2.0*R43*(2*R101 - 2*R103 + 2*R104) - 1.0/2.0*R49*(2*R105 + 2*R106 + 2*R107))/pow(R50, 3.0/2.0);
    double R109 = R108*R31;
    double R110 = 2*R97;
    double R111 = R51*(R105 + R106 + R107);
    double R112 = R0*R111;
    double R113 = R101 - R103 + R104;
    double R114 = R113*R51;
    double R115 = R114*R55;
    double R116 = R1*R109;
    double R117 = R108*R49;
    double R118 = R0*R117;
    double R119 = R108*R43;
    double R120 = R119*R55;
    double R121 = 2*R57*(-2*R112 + 2*R115 + 2*R116 - 2*R118 + 2*R120)/pow(R58, 2);
    double R122 = R121*R94;
    double R123 = R121*R96;
    double R124 = -R33 - R36;
    double R125 = R124*R61;
    double R126 = -R62 + R63 - panel_straw0x;
    double R127 = R126*R38;
    double R128 = R65*(R39 - R41);
    double R129 = R17 - R66;
    double R130 = R129*(R124*R13 + R85);
    double R131 = panel_straw0z - wire_z;
    double R132 = R131*R86;
    double R133 = R125 + R127 + R128 + R130 + R132;
    double R134 = R126*R46;
    double R135 = -R35 - R47;
    double R136 = R129*(R135*R71 + R81 - R82);
    double R137 = R131*(R135*R77 + R73*(R40 - R44) - R83);
    double R138 = -R134 + R136 + R137 + R79 + R80;
    double R139 = R0*R133 - R138*R55;
    double R140 = -R112 + R115 + R116 - R118 + R120;
    double R141 = R109*R78 + R111*R84 + R113*R88 + R117*R84 + R119*R87 + R133*R53 + R138*R56;
    double R142 = 2*R59;
    double R143 = R142*(R139*R57 + R140*R93 + R141);
    double R144 = R142*(R139 + R140*R89 + R141*R57);
    double R145 = ((1.0/2.0)*R100*(-R110*R114 - R110*R119 + R122*R55 - R123*R56 - 2*R134 + 2*R136 + 2*R137 - R143*R56 + R144*R55 + 2*R79 + 2*R80) + (1.0/2.0)*R98*(R1*R122 + R1*R144 - R109*R110 - R123*R52 - R143*R52) + (1.0/2.0)*R99*(-R0*R122 - R0*R144 - R110*R111 - R110*R117 - R123*R53 + 2*R125 + 2*R127 + 2*R128 + 2*R130 + 2*R132 - R143*R53))/sqrt(pow(R100, 2) + pow(R98, 2) + pow(R99, 2));
    double R146 = R0*R56 + R54*a1;
    double R147 = R1*R56 - R52*R55;
    double R148 = -R0*R52 - R54*b1;
    double R149 = pow(pow(R146, 2) + pow(R147, 2) + pow(R148, 2), -1.0/2.0);
    double result = ((R146*R149*R78 + R147*R149*R84 + R148*R149*R87 > 0) ? (
   -R145
)
: (
   R145
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dx(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = 1.0*R5*R6;
    double R8 = sin(plane_b);
    double R9 = sin(panel_g);
    double R10 = cos(panel_a);
    double R11 = R10*R9;
    double R12 = sin(panel_b);
    double R13 = cos(panel_g);
    double R14 = R13*R5;
    double R15 = -R11 + R12*R14;
    double R16 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R17 = 1.0/R16;
    double R18 = R17*panel_straw0x;
    double R19 = R10*R13;
    double R20 = R5*R9;
    double R21 = R12*R20 + R19;
    double R22 = R17*panel_straw0y;
    double R23 = 1.0*R15*R18 - 1.0*R21*R22;
    double R24 = 1.0*R15*R22 + 1.0*R18*R21;
    double R25 = sin(plane_a);
    double R26 = R25*R3;
    double R27 = -R23*R8 + R24*R26 + R4*R7;
    double R28 = sin(plane_g);
    double R29 = R25*R28;
    double R30 = cos(plane_g);
    double R31 = R2*R30;
    double R32 = R29 + R31*R8;
    double R33 = R3*R30;
    double R34 = R2*R28;
    double R35 = R25*R30;
    double R36 = R35*R8;
    double R37 = -R34 + R36;
    double R38 = R23*R33 + R24*R37 + R32*R7;
    double R39 = R34*R8 - R35;
    double R40 = R28*R3;
    double R41 = R29*R8;
    double R42 = R31 + R41;
    double R43 = R23*R40 + R24*R42 + R39*R7;
    double R44 = pow(pow(R27, 2) + pow(R38, 2) + pow(R43, 2), -1.0/2.0);
    double R45 = R27*R44;
    double R46 = R43*R44;
    double R47 = R0*R46;
    double R48 = R38*R44;
    double R49 = R0*a1;
    double R50 = R1*R45 - R47 + R48*R49;
    double R51 = 1.0/(1.0 - pow(R50, 2));
    double R52 = -plane_z;
    double R53 = R52 + panel_dz + panel_straw0z;
    double R54 = R18*panel_dx - R22*panel_dy + panel_straw0x;
    double R55 = R18*panel_dy + R22*panel_dx + panel_straw0y;
    double R56 = -R16 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R57 = R18*R6;
    double R58 = R22*R6;
    double R59 = R13*R57 - R58*R9;
    double R60 = R13*R58 + R57*R9;
    double R61 = -panel_straw0z + wire_z;
    double R62 = R10*R6;
    double R63 = R12*R19 + R20;
    double R64 = R11*R12 - R14;
    double R65 = R18*R63 - R22*R64;
    double R66 = R18*R64 + R22*R63;
    double R67 = -R26*R55 - R4*R53 + R52 + R54*R8 - R56*(-R12*R4 + R26*R60 - R59*R8) - R61*(R26*R66 + R4*R62 - R65*R8) + b0 - plane_dz;
    double R68 = R1*R67;
    double R69 = -R39*R53 - R40*R54 - R42*R55 - R56*(-R12*R39 + R40*R59 + R42*R60) - R61*(R39*R62 + R40*R65 + R42*R66) - plane_dy;
    double R70 = R0*R69;
    double R71 = -R32*R53 - R33*R54 - R37*R55 - R56*(-R12*R32 + R33*R59 + R37*R60) - R61*(R32*R62 + R33*R65 + R37*R66) + a0 - plane_dx;
    double R72 = R49*R71;
    double R73 = R45*R67 + R46*R69 + R48*R71;
    double R74 = R51*(R50*R73 - R68 + R70 - R72);
    double R75 = R51*(-R50*(R68 - R70 + R72) + R73);
    double R76 = R1*R74 - R45*R75 + R67;
    double R77 = -R0*R74 - R46*R75 + R69;
    double R78 = -R48*R75 + R49*R74 + R71;
    double R79 = R18*R8;
    double R80 = R22*R26;
    double R81 = R79 - R80;
    double R82 = R18*R33;
    double R83 = R22*(R34 - R36);
    double R84 = -R82 + R83;
    double R85 = R18*R40;
    double R86 = R22*(-R31 - R41);
    double R87 = -R85 + R86;
    double R88 = R45*R81 + R46*R87 + R48*R84;
    double R89 = R0*R87 - R1*R81 - R49*R84;
    double R90 = 2*R51;
    double R91 = R90*(R50*R88 + R89);
    double R92 = R90*(R50*R89 + R88);
    double R93 = ((1.0/2.0)*R76*(R1*R91 - R45*R92 + 2*R79 - 2*R80) + (1.0/2.0)*R77*(-R0*R91 - R46*R92 - 2*R85 + 2*R86) + (1.0/2.0)*R78*(-R48*R92 + R49*R91 - 2*R82 + 2*R83))/sqrt(pow(R76, 2) + pow(R77, 2) + pow(R78, 2));
    double R94 = R0*R48 + R47*a1;
    double R95 = R1*R48 - R45*R49;
    double R96 = -R0*R45 - R47*b1;
    double R97 = pow(pow(R94, 2) + pow(R95, 2) + pow(R96, 2), -1.0/2.0);
    double result = ((R67*R94*R97 + R69*R95*R97 + R71*R96*R97 > 0) ? (
   -R93
)
: (
   R93
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dy(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = 1.0*R5*R6;
    double R8 = sin(plane_b);
    double R9 = sin(panel_g);
    double R10 = cos(panel_a);
    double R11 = R10*R9;
    double R12 = sin(panel_b);
    double R13 = cos(panel_g);
    double R14 = R13*R5;
    double R15 = -R11 + R12*R14;
    double R16 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R17 = 1.0/R16;
    double R18 = R17*panel_straw0x;
    double R19 = R10*R13;
    double R20 = R5*R9;
    double R21 = R12*R20 + R19;
    double R22 = R17*panel_straw0y;
    double R23 = 1.0*R15*R18 - 1.0*R21*R22;
    double R24 = 1.0*R15*R22 + 1.0*R18*R21;
    double R25 = sin(plane_a);
    double R26 = R25*R3;
    double R27 = -R23*R8 + R24*R26 + R4*R7;
    double R28 = sin(plane_g);
    double R29 = R25*R28;
    double R30 = cos(plane_g);
    double R31 = R2*R30;
    double R32 = R29 + R31*R8;
    double R33 = R3*R30;
    double R34 = R2*R28;
    double R35 = R25*R30;
    double R36 = R35*R8;
    double R37 = -R34 + R36;
    double R38 = R23*R33 + R24*R37 + R32*R7;
    double R39 = R34*R8 - R35;
    double R40 = R28*R3;
    double R41 = R29*R8;
    double R42 = R31 + R41;
    double R43 = R23*R40 + R24*R42 + R39*R7;
    double R44 = pow(pow(R27, 2) + pow(R38, 2) + pow(R43, 2), -1.0/2.0);
    double R45 = R27*R44;
    double R46 = R43*R44;
    double R47 = R0*R46;
    double R48 = R38*R44;
    double R49 = R0*a1;
    double R50 = R1*R45 - R47 + R48*R49;
    double R51 = 1.0/(1.0 - pow(R50, 2));
    double R52 = -plane_z;
    double R53 = R52 + panel_dz + panel_straw0z;
    double R54 = R18*panel_dx - R22*panel_dy + panel_straw0x;
    double R55 = R18*panel_dy + R22*panel_dx + panel_straw0y;
    double R56 = -R16 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R57 = R18*R6;
    double R58 = R22*R6;
    double R59 = R13*R57 - R58*R9;
    double R60 = R13*R58 + R57*R9;
    double R61 = -panel_straw0z + wire_z;
    double R62 = R10*R6;
    double R63 = R12*R19 + R20;
    double R64 = R11*R12 - R14;
    double R65 = R18*R63 - R22*R64;
    double R66 = R18*R64 + R22*R63;
    double R67 = -R26*R55 - R4*R53 + R52 + R54*R8 - R56*(-R12*R4 + R26*R60 - R59*R8) - R61*(R26*R66 + R4*R62 - R65*R8) + b0 - plane_dz;
    double R68 = R1*R67;
    double R69 = -R39*R53 - R40*R54 - R42*R55 - R56*(-R12*R39 + R40*R59 + R42*R60) - R61*(R39*R62 + R40*R65 + R42*R66) - plane_dy;
    double R70 = R0*R69;
    double R71 = -R32*R53 - R33*R54 - R37*R55 - R56*(-R12*R32 + R33*R59 + R37*R60) - R61*(R32*R62 + R33*R65 + R37*R66) + a0 - plane_dx;
    double R72 = R49*R71;
    double R73 = R45*R67 + R46*R69 + R48*R71;
    double R74 = R51*(R50*R73 - R68 + R70 - R72);
    double R75 = R51*(-R50*(R68 - R70 + R72) + R73);
    double R76 = R1*R74 - R45*R75 + R67;
    double R77 = -R0*R74 - R46*R75 + R69;
    double R78 = -R48*R75 + R49*R74 + R71;
    double R79 = R22*R8;
    double R80 = R18*R26;
    double R81 = -R79 - R80;
    double R82 = R22*R33;
    double R83 = R18*(R34 - R36);
    double R84 = R82 + R83;
    double R85 = R22*R40;
    double R86 = R18*(-R31 - R41);
    double R87 = R85 + R86;
    double R88 = R45*R81 + R46*R87 + R48*R84;
    double R89 = R0*R87 - R1*R81 - R49*R84;
    double R90 = 2*R51;
    double R91 = R90*(R50*R88 + R89);
    double R92 = R90*(R50*R89 + R88);
    double R93 = ((1.0/2.0)*R76*(R1*R91 - R45*R92 - 2*R79 - 2*R80) + (1.0/2.0)*R77*(-R0*R91 - R46*R92 + 2*R85 + 2*R86) + (1.0/2.0)*R78*(-R48*R92 + R49*R91 + 2*R82 + 2*R83))/sqrt(pow(R76, 2) + pow(R77, 2) + pow(R78, 2));
    double R94 = R0*R48 + R47*a1;
    double R95 = R1*R48 - R45*R49;
    double R96 = -R0*R45 - R47*b1;
    double R97 = pow(pow(R94, 2) + pow(R95, 2) + pow(R96, 2), -1.0/2.0);
    double result = ((R67*R94*R97 + R69*R95*R97 + R71*R96*R97 > 0) ? (
   -R93
)
: (
   R93
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dz(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = 1.0*R5*R6;
    double R8 = sin(plane_b);
    double R9 = sin(panel_g);
    double R10 = cos(panel_a);
    double R11 = R10*R9;
    double R12 = sin(panel_b);
    double R13 = cos(panel_g);
    double R14 = R13*R5;
    double R15 = -R11 + R12*R14;
    double R16 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R17 = 1.0/R16;
    double R18 = R17*panel_straw0x;
    double R19 = R10*R13;
    double R20 = R5*R9;
    double R21 = R12*R20 + R19;
    double R22 = R17*panel_straw0y;
    double R23 = 1.0*R15*R18 - 1.0*R21*R22;
    double R24 = 1.0*R15*R22 + 1.0*R18*R21;
    double R25 = sin(plane_a);
    double R26 = R25*R3;
    double R27 = -R23*R8 + R24*R26 + R4*R7;
    double R28 = sin(plane_g);
    double R29 = R25*R28;
    double R30 = cos(plane_g);
    double R31 = R2*R30;
    double R32 = R31*R8;
    double R33 = R29 + R32;
    double R34 = R3*R30;
    double R35 = R2*R28;
    double R36 = R25*R30;
    double R37 = -R35 + R36*R8;
    double R38 = R23*R34 + R24*R37 + R33*R7;
    double R39 = R35*R8;
    double R40 = -R36 + R39;
    double R41 = R28*R3;
    double R42 = R29*R8 + R31;
    double R43 = R23*R41 + R24*R42 + R40*R7;
    double R44 = pow(pow(R27, 2) + pow(R38, 2) + pow(R43, 2), -1.0/2.0);
    double R45 = R27*R44;
    double R46 = R43*R44;
    double R47 = R0*R46;
    double R48 = R38*R44;
    double R49 = R0*a1;
    double R50 = R1*R45 - R47 + R48*R49;
    double R51 = 1.0/(1.0 - pow(R50, 2));
    double R52 = -plane_z;
    double R53 = R52 + panel_dz + panel_straw0z;
    double R54 = R18*panel_dx - R22*panel_dy + panel_straw0x;
    double R55 = R18*panel_dy + R22*panel_dx + panel_straw0y;
    double R56 = -R16 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R57 = R18*R6;
    double R58 = R22*R6;
    double R59 = R13*R57 - R58*R9;
    double R60 = R13*R58 + R57*R9;
    double R61 = -panel_straw0z + wire_z;
    double R62 = R10*R6;
    double R63 = R12*R19 + R20;
    double R64 = R11*R12 - R14;
    double R65 = R18*R63 - R22*R64;
    double R66 = R18*R64 + R22*R63;
    double R67 = -R26*R55 - R4*R53 + R52 + R54*R8 - R56*(-R12*R4 + R26*R60 - R59*R8) - R61*(R26*R66 + R4*R62 - R65*R8) + b0 - plane_dz;
    double R68 = R1*R67;
    double R69 = -R40*R53 - R41*R54 - R42*R55 - R56*(-R12*R40 + R41*R59 + R42*R60) - R61*(R40*R62 + R41*R65 + R42*R66) - plane_dy;
    double R70 = R0*R69;
    double R71 = -R33*R53 - R34*R54 - R37*R55 - R56*(-R12*R33 + R34*R59 + R37*R60) - R61*(R33*R62 + R34*R65 + R37*R66) + a0 - plane_dx;
    double R72 = R49*R71;
    double R73 = R45*R67 + R46*R69 + R48*R71;
    double R74 = R51*(R50*R73 - R68 + R70 - R72);
    double R75 = R51*(-R50*(R68 - R70 + R72) + R73);
    double R76 = R1*R74 - R45*R75 + R67;
    double R77 = -R0*R74 - R46*R75 + R69;
    double R78 = -R48*R75 + R49*R74 + R71;
    double R79 = R36 - R39;
    double R80 = -R29 - R32;
    double R81 = -R4*R45 + R46*R79 + R48*R80;
    double R82 = R0*R79 + R1*R4 - R49*R80;
    double R83 = 2*R51;
    double R84 = R83*(R50*R81 + R82);
    double R85 = R83*(R50*R82 + R81);
    double R86 = ((1.0/2.0)*R76*(R1*R84 - 2*R4 - R45*R85) + (1.0/2.0)*R77*(-R0*R84 + 2*R36 - 2*R39 - R46*R85) + (1.0/2.0)*R78*(-2*R29 - 2*R32 - R48*R85 + R49*R84))/sqrt(pow(R76, 2) + pow(R77, 2) + pow(R78, 2));
    double R87 = R0*R48 + R47*a1;
    double R88 = R1*R48 - R45*R49;
    double R89 = -R0*R45 - R47*b1;
    double R90 = pow(pow(R87, 2) + pow(R88, 2) + pow(R89, 2), -1.0/2.0);
    double result = ((R67*R87*R90 + R69*R88*R90 + R71*R89*R90 > 0) ? (
   -R86
)
: (
   R86
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_panel_a(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = R5*R6;
    double R8 = R4*R7;
    double R9 = sin(plane_b);
    double R10 = sin(panel_g);
    double R11 = cos(panel_a);
    double R12 = R10*R11;
    double R13 = sin(panel_b);
    double R14 = cos(panel_g);
    double R15 = R14*R5;
    double R16 = R13*R15;
    double R17 = -R12 + R16;
    double R18 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R19 = 1.0/R18;
    double R20 = R19*panel_straw0x;
    double R21 = R11*R14;
    double R22 = R10*R5;
    double R23 = R13*R22;
    double R24 = R21 + R23;
    double R25 = R19*panel_straw0y;
    double R26 = 1.0*R17*R20 - 1.0*R24*R25;
    double R27 = 1.0*R17*R25 + 1.0*R20*R24;
    double R28 = sin(plane_a);
    double R29 = R28*R3;
    double R30 = -R26*R9 + R27*R29 + 1.0*R8;
    double R31 = sin(plane_g);
    double R32 = R28*R31;
    double R33 = cos(plane_g);
    double R34 = R2*R33;
    double R35 = R34*R9;
    double R36 = R32 + R35;
    double R37 = R36*R7;
    double R38 = R3*R33;
    double R39 = R2*R31;
    double R40 = R28*R33;
    double R41 = R40*R9;
    double R42 = -R39 + R41;
    double R43 = R26*R38 + R27*R42 + 1.0*R37;
    double R44 = R39*R9;
    double R45 = -R40 + R44;
    double R46 = R45*R7;
    double R47 = R3*R31;
    double R48 = R32*R9;
    double R49 = R34 + R48;
    double R50 = R26*R47 + R27*R49 + 1.0*R46;
    double R51 = pow(R30, 2) + pow(R43, 2) + pow(R50, 2);
    double R52 = pow(R51, -1.0/2.0);
    double R53 = R30*R52;
    double R54 = R50*R52;
    double R55 = R0*R54;
    double R56 = R0*a1;
    double R57 = R43*R52;
    double R58 = R1*R53 - R55 + R56*R57;
    double R59 = 1.0 - pow(R58, 2);
    double R60 = 1.0/R59;
    double R61 = -plane_z;
    double R62 = R61 + panel_dz + panel_straw0z;
    double R63 = R20*panel_dx - R25*panel_dy + panel_straw0x;
    double R64 = R20*panel_dy + R25*panel_dx + panel_straw0y;
    double R65 = -R18 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R66 = R20*R6;
    double R67 = R25*R6;
    double R68 = -R10*R67 + R14*R66;
    double R69 = R10*R66 + R14*R67;
    double R70 = -panel_straw0z + wire_z;
    double R71 = R11*R6;
    double R72 = R4*R71;
    double R73 = R13*R21 + R22;
    double R74 = R20*R73;
    double R75 = R12*R13 - R15;
    double R76 = R25*R75;
    double R77 = R74 - R76;
    double R78 = R25*R73;
    double R79 = R20*R75;
    double R80 = R78 + R79;
    double R81 = -R29*R64 - R4*R62 + R61 + R63*R9 - R65*(-R13*R4 + R29*R69 - R68*R9) - R70*(R29*R80 + R72 - R77*R9) + b0 - plane_dz;
    double R82 = R52*R81;
    double R83 = -R45*R62 - R47*R63 - R49*R64 - R65*(-R13*R45 + R47*R68 + R49*R69) - R70*(R45*R71 + R47*R77 + R49*R80) - plane_dy;
    double R84 = -R36*R62 - R38*R63 - R42*R64 - R65*(-R13*R36 + R38*R68 + R42*R69) - R70*(R36*R71 + R38*R77 + R42*R80) + a0 - plane_dx;
    double R85 = R52*R84;
    double R86 = R30*R82 + R43*R85 + R54*R83;
    double R87 = R1*R81;
    double R88 = R0*R83;
    double R89 = R56*R84;
    double R90 = -R87 + R88 - R89;
    double R91 = R58*R86 + R90;
    double R92 = R60*R91;
    double R93 = -R58*(R87 - R88 + R89) + R86;
    double R94 = R60*R93;
    double R95 = R1*R92 - R53*R94 + R81;
    double R96 = -R0*R92 - R54*R94 + R83;
    double R97 = R56*R92 - R57*R94 + R84;
    double R98 = R12 - R16;
    double R99 = R20*R98;
    double R100 = -R21 - R23;
    double R101 = R100*R25;
    double R102 = R100*R20 + R25*R98;
    double R103 = R102*R29 - R8 + R9*(R101 - R99);
    double R104 = panel_straw0z - wire_z;
    double R105 = 2*R104;
    double R106 = 1.0*R74;
    double R107 = 1.0*R76;
    double R108 = R9*(-R106 + R107);
    double R109 = R29*(1.0*R78 + 1.0*R79);
    double R110 = R108 + R109 + 1.0*R72;
    double R111 = R110*R52;
    double R112 = 2*R94;
    double R113 = R71*(1.0*R32 + 1.0*R35);
    double R114 = R106 - R107;
    double R115 = R114*R38;
    double R116 = R80*(-1.0*R39 + 1.0*R41);
    double R117 = R71*(-1.0*R40 + 1.0*R44);
    double R118 = R114*R47;
    double R119 = R80*(1.0*R34 + 1.0*R48);
    double R120 = (-1.0/2.0*R30*(2*R108 + 2*R109 + 2.0*R72) - 1.0/2.0*R43*(2*R113 + 2*R115 + 2*R116) - 1.0/2.0*R50*(2*R117 + 2*R118 + 2*R119))/pow(R51, 3.0/2.0);
    double R121 = R120*R30;
    double R122 = R1*R111;
    double R123 = R52*(R117 + R118 + R119);
    double R124 = R0*R123;
    double R125 = R113 + R115 + R116;
    double R126 = R125*R52;
    double R127 = R126*R56;
    double R128 = R1*R121;
    double R129 = R120*R50;
    double R130 = R0*R129;
    double R131 = R120*R43;
    double R132 = R131*R56;
    double R133 = 2*R58*(2*R122 - 2*R124 + 2*R127 + 2*R128 - 2*R130 + 2*R132)/pow(R59, 2);
    double R134 = R133*R91;
    double R135 = R133*R93;
    double R136 = R122 - R124 + R127 + R128 - R130 + R132;
    double R137 = R103*R104;
    double R138 = -R101 + R99;
    double R139 = R104*(R102*R49 + R138*R47 - R46);
    double R140 = R102*R42 + R138*R38 - R37;
    double R141 = R104*R140;
    double R142 = R0*R139 - R1*R137 - R141*R56;
    double R143 = R110*R82 + R121*R81 + R123*R83 + R125*R85 + R129*R83 + R131*R84 + R137*R53 + R139*R54 + R141*R57;
    double R144 = 2*R60;
    double R145 = R144*(R136*R90 + R142*R58 + R143);
    double R146 = R144*(R136*R86 + R142 + R143*R58);
    double R147 = ((1.0/2.0)*R95*(R1*R134 + R1*R146 + R103*R105 - R111*R112 - R112*R121 - R135*R53 - R145*R53) + (1.0/2.0)*R96*(-R0*R134 - R0*R146 - R112*R123 - R112*R129 - R135*R54 + 2*R139 - R145*R54) + (1.0/2.0)*R97*(R105*R140 - R112*R126 - R112*R131 + R134*R56 - R135*R57 - R145*R57 + R146*R56))/sqrt(pow(R95, 2) + pow(R96, 2) + pow(R97, 2));
    double R148 = R0*R52;
    double R149 = R148*R43 + R55*a1;
    double R150 = R1*R57 - R53*R56;
    double R151 = -R148*R30 - R55*b1;
    double R152 = pow(pow(R149, 2) + pow(R150, 2) + pow(R151, 2), -1.0/2.0);
    double result = ((R149*R152*R81 + R150*R152*R83 + R151*R152*R84 > 0) ? (
   -R147
)
: (
   R147
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_panel_b(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(panel_b);
    double R3 = cos(plane_a);
    double R4 = cos(plane_b);
    double R5 = R3*R4;
    double R6 = R2*R5;
    double R7 = sin(panel_a);
    double R8 = 1.0*R7;
    double R9 = sin(plane_b);
    double R10 = sin(panel_g);
    double R11 = cos(panel_a);
    double R12 = R10*R11;
    double R13 = sin(panel_b);
    double R14 = cos(panel_g);
    double R15 = R14*R7;
    double R16 = -R12 + R13*R15;
    double R17 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R18 = 1.0/R17;
    double R19 = R18*panel_straw0x;
    double R20 = R11*R14;
    double R21 = R10*R7;
    double R22 = R13*R21 + R20;
    double R23 = R18*panel_straw0y;
    double R24 = 1.0*R16*R19 - 1.0*R22*R23;
    double R25 = 1.0*R16*R23 + 1.0*R19*R22;
    double R26 = sin(plane_a);
    double R27 = R26*R4;
    double R28 = -R24*R9 + R25*R27 + R6*R8;
    double R29 = sin(plane_g);
    double R30 = R26*R29;
    double R31 = cos(plane_g);
    double R32 = R3*R31;
    double R33 = R32*R9;
    double R34 = R30 + R33;
    double R35 = R2*R8;
    double R36 = R31*R4;
    double R37 = R29*R3;
    double R38 = R26*R31;
    double R39 = R38*R9;
    double R40 = -R37 + R39;
    double R41 = R24*R36 + R25*R40 + R34*R35;
    double R42 = R37*R9;
    double R43 = -R38 + R42;
    double R44 = R29*R4;
    double R45 = R30*R9;
    double R46 = R32 + R45;
    double R47 = R24*R44 + R25*R46 + R35*R43;
    double R48 = pow(R28, 2) + pow(R41, 2) + pow(R47, 2);
    double R49 = pow(R48, -1.0/2.0);
    double R50 = R28*R49;
    double R51 = R47*R49;
    double R52 = R0*R51;
    double R53 = R41*R49;
    double R54 = R0*a1;
    double R55 = R1*R50 - R52 + R53*R54;
    double R56 = 1.0 - pow(R55, 2);
    double R57 = 1.0/R56;
    double R58 = -plane_z;
    double R59 = R58 + panel_dz + panel_straw0z;
    double R60 = R19*panel_dx - R23*panel_dy + panel_straw0x;
    double R61 = R19*panel_dy + R23*panel_dx + panel_straw0y;
    double R62 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R63 = -R17 + R62;
    double R64 = R13*R5;
    double R65 = R19*R2;
    double R66 = R2*R23;
    double R67 = -R10*R66 + R14*R65;
    double R68 = R10*R65 + R14*R66;
    double R69 = -panel_straw0z + wire_z;
    double R70 = R11*R2;
    double R71 = R13*R20 + R21;
    double R72 = R12*R13 - R15;
    double R73 = R19*R71 - R23*R72;
    double R74 = R19*R72 + R23*R71;
    double R75 = -R27*R61 - R5*R59 + R58 + R60*R9 - R63*(R27*R68 - R64 - R67*R9) - R69*(R27*R74 + R5*R70 - R73*R9) + b0 - plane_dz;
    double R76 = R13*R43;
    double R77 = -R43*R59 - R44*R60 - R46*R61 - R63*(R44*R67 + R46*R68 - R76) - R69*(R43*R70 + R44*R73 + R46*R74) - plane_dy;
    double R78 = R13*R34;
    double R79 = -R34*R59 - R36*R60 - R40*R61 - R63*(R36*R67 + R40*R68 - R78) - R69*(R34*R70 + R36*R73 + R40*R74) + a0 - plane_dx;
    double R80 = R50*R75 + R51*R77 + R53*R79;
    double R81 = R1*R75;
    double R82 = R0*R77;
    double R83 = R54*R79;
    double R84 = -R81 + R82 - R83;
    double R85 = R55*R80 + R84;
    double R86 = R57*R85;
    double R87 = -R55*(R81 - R82 + R83) + R80;
    double R88 = R57*R87;
    double R89 = R1*R86 - R50*R88 + R75;
    double R90 = -R0*R86 - R51*R88 + R77;
    double R91 = -R53*R88 + R54*R86 + R79;
    double R92 = panel_straw0z - wire_z;
    double R93 = R12*R66;
    double R94 = R20*R65;
    double R95 = R12*R65 + R20*R66;
    double R96 = R92*(-R11*R64 + R27*R95 + R9*(R93 - R94));
    double R97 = R17 - R62;
    double R98 = R13*R14;
    double R99 = R19*R98;
    double R100 = R10*R13;
    double R101 = R100*R23;
    double R102 = -R100*R19 - R23*R98;
    double R103 = R97*(R102*R27 - R6 + R9*(-R101 + R99));
    double R104 = 1.0*R15*R65;
    double R105 = 1.0*R21*R66;
    double R106 = R9*(-R104 + R105);
    double R107 = R21*R65;
    double R108 = R15*R66;
    double R109 = R27*(1.0*R107 + 1.0*R108);
    double R110 = R49*(R106 + R109 - R64*R8);
    double R111 = 2*R88;
    double R112 = R13*R7;
    double R113 = R112*(1.0*R30 + 1.0*R33);
    double R114 = R104 - R105;
    double R115 = R114*R36;
    double R116 = R107 + R108;
    double R117 = R116*(-1.0*R37 + 1.0*R39);
    double R118 = R112*(-1.0*R38 + 1.0*R42);
    double R119 = R114*R44;
    double R120 = R116*(1.0*R32 + 1.0*R45);
    double R121 = (-1.0/2.0*R28*(2*R106 + 2*R109 - 2.0*R64*R7) - 1.0/2.0*R41*(-2*R113 + 2*R115 + 2*R117) - 1.0/2.0*R47*(-2*R118 + 2*R119 + 2*R120))/pow(R48, 3.0/2.0);
    double R122 = R121*R28;
    double R123 = R1*R110;
    double R124 = R49*(-R118 + R119 + R120);
    double R125 = R0*R124;
    double R126 = R49*(-R113 + R115 + R117);
    double R127 = R126*R54;
    double R128 = R1*R122;
    double R129 = R121*R47;
    double R130 = R0*R129;
    double R131 = R121*R41;
    double R132 = R131*R54;
    double R133 = 2*R55*(2*R123 - 2*R125 + 2*R127 + 2*R128 - 2*R130 + 2*R132)/pow(R56, 2);
    double R134 = R133*R85;
    double R135 = R133*R87;
    double R136 = R103 + R96;
    double R137 = -R93 + R94;
    double R138 = R92*(-R11*R76 + R137*R44 + R46*R95);
    double R139 = R101 - R99;
    double R140 = R97*(R102*R46 + R139*R44 + R2*(R38 - R42));
    double R141 = R138 + R140;
    double R142 = R92*(-R11*R78 + R137*R36 + R40*R95);
    double R143 = R97*(R102*R40 + R139*R36 + R2*(-R30 - R33));
    double R144 = R142 + R143;
    double R145 = R0*R141 - R1*R136 - R144*R54;
    double R146 = R123 - R125 + R127 + R128 - R130 + R132;
    double R147 = R110*R75 + R122*R75 + R124*R77 + R126*R79 + R129*R77 + R131*R79 + R136*R50 + R141*R51 + R144*R53;
    double R148 = 2*R57;
    double R149 = R148*(R145*R55 + R146*R84 + R147);
    double R150 = R148*(R145 + R146*R80 + R147*R55);
    double R151 = ((1.0/2.0)*R89*(R1*R134 + R1*R150 + 2*R103 - R110*R111 - R111*R122 - R135*R50 - R149*R50 + 2*R96) + (1.0/2.0)*R90*(-R0*R134 - R0*R150 - R111*R124 - R111*R129 - R135*R51 + 2*R138 + 2*R140 - R149*R51) + (1.0/2.0)*R91*(-R111*R126 - R111*R131 + R134*R54 - R135*R53 + 2*R142 + 2*R143 - R149*R53 + R150*R54))/sqrt(pow(R89, 2) + pow(R90, 2) + pow(R91, 2));
    double R152 = R0*R53 + R52*a1;
    double R153 = R1*R53 - R50*R54;
    double R154 = -R0*R50 - R52*b1;
    double R155 = pow(pow(R152, 2) + pow(R153, 2) + pow(R154, 2), -1.0/2.0);
    double result = ((R152*R155*R75 + R153*R155*R77 + R154*R155*R79 > 0) ? (
   -R151
)
: (
   R151
))/driftvel;
    return result;
}


double CosmicTrack_DCA_Deriv_panel_g(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = cos(plane_a);
    double R3 = cos(plane_b);
    double R4 = R2*R3;
    double R5 = sin(panel_a);
    double R6 = cos(panel_b);
    double R7 = 1.0*R5*R6;
    double R8 = sin(plane_b);
    double R9 = sin(panel_g);
    double R10 = cos(panel_a);
    double R11 = R10*R9;
    double R12 = sin(panel_b);
    double R13 = cos(panel_g);
    double R14 = R13*R5;
    double R15 = -R11 + R12*R14;
    double R16 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R17 = 1.0/R16;
    double R18 = R17*panel_straw0x;
    double R19 = R15*R18;
    double R20 = R10*R13;
    double R21 = R5*R9;
    double R22 = R12*R21;
    double R23 = R20 + R22;
    double R24 = R17*panel_straw0y;
    double R25 = 1.0*R19 - 1.0*R23*R24;
    double R26 = R15*R24;
    double R27 = 1.0*R18*R23 + 1.0*R26;
    double R28 = sin(plane_a);
    double R29 = R28*R3;
    double R30 = -R25*R8 + R27*R29 + R4*R7;
    double R31 = sin(plane_g);
    double R32 = R28*R31;
    double R33 = cos(plane_g);
    double R34 = R2*R33;
    double R35 = R32 + R34*R8;
    double R36 = R3*R33;
    double R37 = R2*R31;
    double R38 = R28*R33;
    double R39 = R38*R8;
    double R40 = -R37 + R39;
    double R41 = R25*R36 + R27*R40 + R35*R7;
    double R42 = R37*R8 - R38;
    double R43 = R3*R31;
    double R44 = R32*R8;
    double R45 = R34 + R44;
    double R46 = R25*R43 + R27*R45 + R42*R7;
    double R47 = pow(R30, 2) + pow(R41, 2) + pow(R46, 2);
    double R48 = pow(R47, -1.0/2.0);
    double R49 = R30*R48;
    double R50 = R46*R48;
    double R51 = R0*R50;
    double R52 = R41*R48;
    double R53 = R0*a1;
    double R54 = R1*R49 - R51 + R52*R53;
    double R55 = 1.0 - pow(R54, 2);
    double R56 = 1.0/R55;
    double R57 = -plane_z;
    double R58 = R57 + panel_dz + panel_straw0z;
    double R59 = R18*panel_dx - R24*panel_dy + panel_straw0x;
    double R60 = R18*panel_dy + R24*panel_dx + panel_straw0y;
    double R61 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R62 = -R16 + R61;
    double R63 = R18*R6;
    double R64 = R24*R6;
    double R65 = R13*R63 - R64*R9;
    double R66 = R63*R9;
    double R67 = R13*R64;
    double R68 = R66 + R67;
    double R69 = -panel_straw0z + wire_z;
    double R70 = R10*R6;
    double R71 = R12*R20 + R21;
    double R72 = R18*R71;
    double R73 = R11*R12;
    double R74 = -R14 + R73;
    double R75 = -R24*R74 + R72;
    double R76 = R24*R71;
    double R77 = R18*R74 + R76;
    double R78 = -R29*R60 - R4*R58 + R57 + R59*R8 - R62*(-R12*R4 + R29*R68 - R65*R8) - R69*(R29*R77 + R4*R70 - R75*R8) + b0 - plane_dz;
    double R79 = -R42*R58 - R43*R59 - R45*R60 - R62*(-R12*R42 + R43*R65 + R45*R68) - R69*(R42*R70 + R43*R75 + R45*R77) - plane_dy;
    double R80 = -R35*R58 - R36*R59 - R40*R60 - R62*(-R12*R35 + R36*R65 + R40*R68) - R69*(R35*R70 + R36*R75 + R40*R77) + a0 - plane_dx;
    double R81 = R49*R78 + R50*R79 + R52*R80;
    double R82 = R1*R78;
    double R83 = R0*R79;
    double R84 = R53*R80;
    double R85 = -R82 + R83 - R84;
    double R86 = R54*R81 + R85;
    double R87 = R56*R86;
    double R88 = -R54*(R82 - R83 + R84) + R81;
    double R89 = R56*R88;
    double R90 = R1*R87 - R49*R89 + R78;
    double R91 = -R0*R87 - R50*R89 + R79;
    double R92 = -R52*R89 + R53*R87 + R80;
    double R93 = R16 - R61;
    double R94 = R93*(R29*R65 + R68*R8);
    double R95 = panel_straw0z - wire_z;
    double R96 = R14 - R73;
    double R97 = R18*R96;
    double R98 = R24*R96 + R72;
    double R99 = R95*(R29*R98 + R8*(R76 - R97));
    double R100 = 1.0*R26;
    double R101 = -R20 - R22;
    double R102 = 1.0*R101*R18;
    double R103 = R8*(R100 - R102);
    double R104 = R101*R24;
    double R105 = R29*(1.0*R104 + 1.0*R19);
    double R106 = R48*(R103 + R105);
    double R107 = 2*R89;
    double R108 = -R100 + R102;
    double R109 = R108*R43;
    double R110 = R104 + R19;
    double R111 = R110*(1.0*R34 + 1.0*R44);
    double R112 = R108*R36;
    double R113 = R110*(-1.0*R37 + 1.0*R39);
    double R114 = (-1.0/2.0*R30*(2*R103 + 2*R105) - 1.0/2.0*R41*(2*R112 + 2*R113) - 1.0/2.0*R46*(2*R109 + 2*R111))/pow(R47, 3.0/2.0);
    double R115 = R114*R30;
    double R116 = R1*R106;
    double R117 = R48*(R109 + R111);
    double R118 = R0*R117;
    double R119 = R48*(R112 + R113);
    double R120 = R119*R53;
    double R121 = R1*R115;
    double R122 = R114*R46;
    double R123 = R0*R122;
    double R124 = R114*R41;
    double R125 = R124*R53;
    double R126 = 2*R54*(2*R116 - 2*R118 + 2*R120 + 2*R121 - 2*R123 + 2*R125)/pow(R55, 2);
    double R127 = R126*R86;
    double R128 = R126*R88;
    double R129 = R94 + R99;
    double R130 = -R66 - R67;
    double R131 = R93*(R130*R43 + R45*R65);
    double R132 = -R76 + R97;
    double R133 = R95*(R132*R43 + R45*R98);
    double R134 = R131 + R133;
    double R135 = R93*(R130*R36 + R40*R65);
    double R136 = R95*(R132*R36 + R40*R98);
    double R137 = R135 + R136;
    double R138 = R0*R134 - R1*R129 - R137*R53;
    double R139 = R116 - R118 + R120 + R121 - R123 + R125;
    double R140 = R106*R78 + R115*R78 + R117*R79 + R119*R80 + R122*R79 + R124*R80 + R129*R49 + R134*R50 + R137*R52;
    double R141 = 2*R56;
    double R142 = R141*(R138*R54 + R139*R85 + R140);
    double R143 = R141*(R138 + R139*R81 + R140*R54);
    double R144 = ((1.0/2.0)*R90*(R1*R127 + R1*R143 - R106*R107 - R107*R115 - R128*R49 - R142*R49 + 2*R94 + 2*R99) + (1.0/2.0)*R91*(-R0*R127 - R0*R143 - R107*R117 - R107*R122 - R128*R50 + 2*R131 + 2*R133 - R142*R50) + (1.0/2.0)*R92*(-R107*R119 - R107*R124 + R127*R53 - R128*R52 + 2*R135 + 2*R136 - R142*R52 + R143*R53))/sqrt(pow(R90, 2) + pow(R91, 2) + pow(R92, 2));
    double R145 = R0*R52 + R51*a1;
    double R146 = R1*R52 - R49*R53;
    double R147 = -R0*R49 - R51*b1;
    double R148 = pow(pow(R145, 2) + pow(R146, 2) + pow(R147, 2), -1.0/2.0);
    double result = ((R145*R148*R78 + R146*R148*R79 + R147*R148*R80 > 0) ? (
   -R144
)
: (
   R144
))/driftvel;
    return result;
}


std::vector<double> CosmicTrack_DCA_LocalDeriv(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
        std::vector<double> result = {CosmicTrack_DCA_Deriv_a0(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_b0(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_a1(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_b1(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_t0(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel)};
return result;
}

std::vector<double> CosmicTrack_DCA_GlobalDeriv(double const& a0, double const& b0, double const& a1, double const& b1, double const& t0, double const& plane_dx, double const& plane_dy, double const& plane_dz, double const& plane_a, double const& plane_b, double const& plane_g, double const& panel_dx, double const& panel_dy, double const& panel_dz, double const& panel_a, double const& panel_b, double const& panel_g, double const& wire_x, double const& wire_y, double const& wire_z, double const& wdir_x, double const& wdir_y, double const& wdir_z, double const& plane_x, double const& plane_y, double const& plane_z, double const& panel_straw0x, double const& panel_straw0y, double const& panel_straw0z, double const& driftvel)
{
        std::vector<double> result = {CosmicTrack_DCA_Deriv_plane_dx(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_plane_dy(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_plane_dz(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_plane_a(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_plane_b(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_plane_g(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_panel_dx(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_panel_dy(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_panel_dz(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_panel_a(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_panel_b(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel),
CosmicTrack_DCA_Deriv_panel_g(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z,driftvel)};
return result;
}

