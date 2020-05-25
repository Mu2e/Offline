

# include "TrackerAlignment/inc/AlignmentDerivatives.hh"
# include <math.h>
# include <vector>

double CosmicTrack_DCA(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = 1.0*R3*R4;
    double R6 = sin(plane_b);
    double R7 = sin(panel_g);
    double R8 = cos(panel_a);
    double R9 = R7*R8;
    double R10 = sin(panel_b);
    double R11 = cos(panel_g);
    double R12 = R11*R3;
    double R13 = R10*R12 - R9;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = 1.0/R14;
    double R16 = R15*panel_straw0x;
    double R17 = R11*R8;
    double R18 = R3*R7;
    double R19 = R10*R18 + R17;
    double R20 = R15*panel_straw0y;
    double R21 = 1.0*R13*R16 - 1.0*R19*R20;
    double R22 = 1.0*R13*R20 + 1.0*R16*R19;
    double R23 = sin(plane_a);
    double R24 = R1*R23;
    double R25 = R2*R5 - R21*R6 + R22*R24;
    double R26 = sin(plane_g);
    double R27 = R23*R26;
    double R28 = cos(plane_g);
    double R29 = R0*R28;
    double R30 = R27 + R29*R6;
    double R31 = R1*R28;
    double R32 = R0*R26;
    double R33 = R23*R28;
    double R34 = -R32 + R33*R6;
    double R35 = R21*R31 + R22*R34 + R30*R5;
    double R36 = R32*R6 - R33;
    double R37 = R1*R26;
    double R38 = R27*R6 + R29;
    double R39 = R21*R37 + R22*R38 + R36*R5;
    double R40 = pow(pow(R25, 2) + pow(R35, 2) + pow(R39, 2), -1.0/2.0);
    double R41 = R25*R40;
    double R42 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R43 = R42*b1;
    double R44 = R39*R40;
    double R45 = R35*R40;
    double R46 = R42*a1;
    double R47 = R41*R43 - R42*R44 + R45*R46;
    double R48 = 1.0/(1.0 - pow(R47, 2));
    double R49 = panel_dz + panel_straw0z - plane_z;
    double R50 = R16*panel_dx - R20*panel_dy + panel_straw0x;
    double R51 = R16*panel_dy + R20*panel_dx + panel_straw0y;
    double R52 = -R14 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R53 = R16*R4;
    double R54 = R20*R4;
    double R55 = R11*R53 - R54*R7;
    double R56 = R11*R54 + R53*R7;
    double R57 = -panel_straw0z + wire_z;
    double R58 = R4*R8;
    double R59 = R10*R17 + R18;
    double R60 = R10*R9 - R12;
    double R61 = R16*R59 - R20*R60;
    double R62 = R16*R60 + R20*R59;
    double R63 = R2*R49 + R24*R51 - R50*R6 + R52*(-R10*R2 + R24*R56 - R55*R6) + R57*(R2*R58 + R24*R62 - R6*R61) - b0 + plane_dz + plane_z;
    double R64 = R41*R63;
    double R65 = R36*R49 + R37*R50 + R38*R51 + R52*(-R10*R36 + R37*R55 + R38*R56) + R57*(R36*R58 + R37*R61 + R38*R62) + plane_dy;
    double R66 = R44*R65;
    double R67 = R30*R49 + R31*R50 + R34*R51 + R52*(-R10*R30 + R31*R55 + R34*R56) + R57*(R30*R58 + R31*R61 + R34*R62) - a0 + plane_dx;
    double R68 = R45*R67;
    double R69 = -R42*R65 + R43*R63 + R46*R67;
    double R70 = R48*(R47*R69 - R64 - R66 - R68);
    double R71 = R48*(-R47*(R64 + R66 + R68) + R69);
    double R72 = sqrt(pow(R41*R70 - R43*R71 + R63, 2) + pow(R42*R71 + R44*R70 + R65, 2) + pow(R45*R70 - R46*R71 + R67, 2));
    double result = ((R71 > 0) ? (
   R72
)
: (
   -1.0*R72
));
    return result;
}


double CosmicTrack_DCAalignpos_x(double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
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


double CosmicTrack_DCAalignpos_y(double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
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


double CosmicTrack_DCAalignpos_z(double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
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


double CosmicTrack_DCAaligndir_x(double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
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


double CosmicTrack_DCAaligndir_y(double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
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


double CosmicTrack_DCAaligndir_z(double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
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


double CosmicTrack_DCA_Deriv_a0(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = 1.0*R3*R4;
    double R6 = sin(plane_b);
    double R7 = sin(panel_g);
    double R8 = cos(panel_a);
    double R9 = R7*R8;
    double R10 = sin(panel_b);
    double R11 = cos(panel_g);
    double R12 = R11*R3;
    double R13 = R10*R12 - R9;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = 1.0/R14;
    double R16 = R15*panel_straw0x;
    double R17 = R11*R8;
    double R18 = R3*R7;
    double R19 = R10*R18 + R17;
    double R20 = R15*panel_straw0y;
    double R21 = 1.0*R13*R16 - 1.0*R19*R20;
    double R22 = 1.0*R13*R20 + 1.0*R16*R19;
    double R23 = sin(plane_a);
    double R24 = R1*R23;
    double R25 = R2*R5 - R21*R6 + R22*R24;
    double R26 = sin(plane_g);
    double R27 = R23*R26;
    double R28 = cos(plane_g);
    double R29 = R0*R28;
    double R30 = R27 + R29*R6;
    double R31 = R1*R28;
    double R32 = R0*R26;
    double R33 = R23*R28;
    double R34 = -R32 + R33*R6;
    double R35 = R21*R31 + R22*R34 + R30*R5;
    double R36 = R32*R6 - R33;
    double R37 = R1*R26;
    double R38 = R27*R6 + R29;
    double R39 = R21*R37 + R22*R38 + R36*R5;
    double R40 = pow(pow(R25, 2) + pow(R35, 2) + pow(R39, 2), -1.0/2.0);
    double R41 = R25*R40;
    double R42 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R43 = R42*b1;
    double R44 = R39*R40;
    double R45 = R42*a1;
    double R46 = R35*R40;
    double R47 = R41*R43 - R42*R44 + R45*R46;
    double R48 = 1.0/(1.0 - pow(R47, 2));
    double R49 = panel_dz + panel_straw0z - plane_z;
    double R50 = R16*panel_dx - R20*panel_dy + panel_straw0x;
    double R51 = R16*panel_dy + R20*panel_dx + panel_straw0y;
    double R52 = -R14 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R53 = R16*R4;
    double R54 = R20*R4;
    double R55 = R11*R53 - R54*R7;
    double R56 = R11*R54 + R53*R7;
    double R57 = -panel_straw0z + wire_z;
    double R58 = R4*R8;
    double R59 = R10*R17 + R18;
    double R60 = R10*R9 - R12;
    double R61 = R16*R59 - R20*R60;
    double R62 = R16*R60 + R20*R59;
    double R63 = R2*R49 + R24*R51 - R50*R6 + R52*(-R10*R2 + R24*R56 - R55*R6) + R57*(R2*R58 + R24*R62 - R6*R61) - b0 + plane_dz + plane_z;
    double R64 = R41*R63;
    double R65 = R36*R49 + R37*R50 + R38*R51 + R52*(-R10*R36 + R37*R55 + R38*R56) + R57*(R36*R58 + R37*R61 + R38*R62) + plane_dy;
    double R66 = R44*R65;
    double R67 = R30*R49 + R31*R50 + R34*R51 + R52*(-R10*R30 + R31*R55 + R34*R56) + R57*(R30*R58 + R31*R61 + R34*R62) - a0 + plane_dx;
    double R68 = R46*R67;
    double R69 = -R42*R65 + R43*R63 + R45*R67;
    double R70 = R48*(R47*R69 - R64 - R66 - R68);
    double R71 = R48*(-R47*(R64 + R66 + R68) + R69);
    double R72 = R41*R70 - R43*R71 + R63;
    double R73 = R42*R71 + R44*R70 + R65;
    double R74 = -R45*R71 + R46*R70 + R67;
    double R75 = 2*R48;
    double R76 = R75*(-R45 + R46*R47);
    double R77 = R75*(-R45*R47 + R46);
    double R78 = ((1.0/2.0)*R72*(R41*R77 - R43*R76) + (1.0/2.0)*R73*(R42*R76 + R44*R77) + (1.0/2.0)*R74*(-R45*R76 + R46*R77 - 2))/sqrt(pow(R72, 2) + pow(R73, 2) + pow(R74, 2));
    double result = 16.0*((R71 > 0) ? (
   R78
)
: (
   -1.0*R78
));
    return result;
}


double CosmicTrack_DCA_Deriv_b0(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = 1.0*R3*R4;
    double R6 = sin(plane_b);
    double R7 = sin(panel_g);
    double R8 = cos(panel_a);
    double R9 = R7*R8;
    double R10 = sin(panel_b);
    double R11 = cos(panel_g);
    double R12 = R11*R3;
    double R13 = R10*R12 - R9;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = 1.0/R14;
    double R16 = R15*panel_straw0x;
    double R17 = R11*R8;
    double R18 = R3*R7;
    double R19 = R10*R18 + R17;
    double R20 = R15*panel_straw0y;
    double R21 = 1.0*R13*R16 - 1.0*R19*R20;
    double R22 = 1.0*R13*R20 + 1.0*R16*R19;
    double R23 = sin(plane_a);
    double R24 = R1*R23;
    double R25 = R2*R5 - R21*R6 + R22*R24;
    double R26 = sin(plane_g);
    double R27 = R23*R26;
    double R28 = cos(plane_g);
    double R29 = R0*R28;
    double R30 = R27 + R29*R6;
    double R31 = R1*R28;
    double R32 = R0*R26;
    double R33 = R23*R28;
    double R34 = -R32 + R33*R6;
    double R35 = R21*R31 + R22*R34 + R30*R5;
    double R36 = R32*R6 - R33;
    double R37 = R1*R26;
    double R38 = R27*R6 + R29;
    double R39 = R21*R37 + R22*R38 + R36*R5;
    double R40 = pow(pow(R25, 2) + pow(R35, 2) + pow(R39, 2), -1.0/2.0);
    double R41 = R25*R40;
    double R42 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R43 = R42*b1;
    double R44 = R39*R40;
    double R45 = R35*R40;
    double R46 = R42*a1;
    double R47 = R41*R43 - R42*R44 + R45*R46;
    double R48 = 1.0/(1.0 - pow(R47, 2));
    double R49 = panel_dz + panel_straw0z - plane_z;
    double R50 = R16*panel_dx - R20*panel_dy + panel_straw0x;
    double R51 = R16*panel_dy + R20*panel_dx + panel_straw0y;
    double R52 = -R14 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R53 = R16*R4;
    double R54 = R20*R4;
    double R55 = R11*R53 - R54*R7;
    double R56 = R11*R54 + R53*R7;
    double R57 = -panel_straw0z + wire_z;
    double R58 = R4*R8;
    double R59 = R10*R17 + R18;
    double R60 = R10*R9 - R12;
    double R61 = R16*R59 - R20*R60;
    double R62 = R16*R60 + R20*R59;
    double R63 = R2*R49 + R24*R51 - R50*R6 + R52*(-R10*R2 + R24*R56 - R55*R6) + R57*(R2*R58 + R24*R62 - R6*R61) - b0 + plane_dz + plane_z;
    double R64 = R41*R63;
    double R65 = R36*R49 + R37*R50 + R38*R51 + R52*(-R10*R36 + R37*R55 + R38*R56) + R57*(R36*R58 + R37*R61 + R38*R62) + plane_dy;
    double R66 = R44*R65;
    double R67 = R30*R49 + R31*R50 + R34*R51 + R52*(-R10*R30 + R31*R55 + R34*R56) + R57*(R30*R58 + R31*R61 + R34*R62) - a0 + plane_dx;
    double R68 = R45*R67;
    double R69 = -R42*R65 + R43*R63 + R46*R67;
    double R70 = R48*(R47*R69 - R64 - R66 - R68);
    double R71 = R48*(-R47*(R64 + R66 + R68) + R69);
    double R72 = R41*R70 - R43*R71 + R63;
    double R73 = R42*R71 + R44*R70 + R65;
    double R74 = R45*R70 - R46*R71 + R67;
    double R75 = 2*R48;
    double R76 = R75*(R41*R47 - R43);
    double R77 = R75*(R41 - R43*R47);
    double R78 = ((1.0/2.0)*R72*(R41*R77 - R43*R76 - 2) + (1.0/2.0)*R73*(R42*R76 + R44*R77) + (1.0/2.0)*R74*(R45*R77 - R46*R76))/sqrt(pow(R72, 2) + pow(R73, 2) + pow(R74, 2));
    double result = 16.0*((R71 > 0) ? (
   R78
)
: (
   -1.0*R78
));
    return result;
}


double CosmicTrack_DCA_Deriv_a1(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = 1.0*R3*R4;
    double R6 = sin(plane_b);
    double R7 = sin(panel_g);
    double R8 = cos(panel_a);
    double R9 = R7*R8;
    double R10 = sin(panel_b);
    double R11 = cos(panel_g);
    double R12 = R11*R3;
    double R13 = R10*R12 - R9;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = 1.0/R14;
    double R16 = R15*panel_straw0x;
    double R17 = R11*R8;
    double R18 = R3*R7;
    double R19 = R10*R18 + R17;
    double R20 = R15*panel_straw0y;
    double R21 = 1.0*R13*R16 - 1.0*R19*R20;
    double R22 = 1.0*R13*R20 + 1.0*R16*R19;
    double R23 = sin(plane_a);
    double R24 = R1*R23;
    double R25 = R2*R5 - R21*R6 + R22*R24;
    double R26 = sin(plane_g);
    double R27 = R23*R26;
    double R28 = cos(plane_g);
    double R29 = R0*R28;
    double R30 = R27 + R29*R6;
    double R31 = R1*R28;
    double R32 = R0*R26;
    double R33 = R23*R28;
    double R34 = -R32 + R33*R6;
    double R35 = R21*R31 + R22*R34 + R30*R5;
    double R36 = R32*R6 - R33;
    double R37 = R1*R26;
    double R38 = R27*R6 + R29;
    double R39 = R21*R37 + R22*R38 + R36*R5;
    double R40 = pow(pow(R25, 2) + pow(R35, 2) + pow(R39, 2), -1.0/2.0);
    double R41 = R25*R40;
    double R42 = pow(a1, 2);
    double R43 = R42 + pow(b1, 2) + 1;
    double R44 = pow(R43, -1.0/2.0);
    double R45 = R44*b1;
    double R46 = R39*R40;
    double R47 = R35*R40;
    double R48 = R44*R47;
    double R49 = R41*R45 - R44*R46 + R48*a1;
    double R50 = 1.0 - pow(R49, 2);
    double R51 = 1.0/R50;
    double R52 = panel_dz + panel_straw0z - plane_z;
    double R53 = R16*panel_dx - R20*panel_dy + panel_straw0x;
    double R54 = R16*panel_dy + R20*panel_dx + panel_straw0y;
    double R55 = -R14 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R56 = R16*R4;
    double R57 = R20*R4;
    double R58 = R11*R56 - R57*R7;
    double R59 = R11*R57 + R56*R7;
    double R60 = -panel_straw0z + wire_z;
    double R61 = R4*R8;
    double R62 = R10*R17 + R18;
    double R63 = R10*R9 - R12;
    double R64 = R16*R62 - R20*R63;
    double R65 = R16*R63 + R20*R62;
    double R66 = R2*R52 + R24*R54 - R53*R6 + R55*(-R10*R2 + R24*R59 - R58*R6) + R60*(R2*R61 + R24*R65 - R6*R64) - b0 + plane_dz + plane_z;
    double R67 = R36*R52 + R37*R53 + R38*R54 + R55*(-R10*R36 + R37*R58 + R38*R59) + R60*(R36*R61 + R37*R64 + R38*R65) + plane_dy;
    double R68 = R30*R52 + R31*R53 + R34*R54 + R55*(-R10*R30 + R31*R58 + R34*R59) + R60*(R30*R61 + R31*R64 + R34*R65) - a0 + plane_dx;
    double R69 = R44*R68;
    double R70 = -R44*R67 + R45*R66 + R69*a1;
    double R71 = R41*R66;
    double R72 = R46*R67;
    double R73 = R47*R68;
    double R74 = -R71 - R72 - R73;
    double R75 = R49*R70 + R74;
    double R76 = R51*R75;
    double R77 = -R49*(R71 + R72 + R73) + R70;
    double R78 = R51*R77;
    double R79 = R41*R76 - R45*R78 + R66;
    double R80 = R44*R78;
    double R81 = R46*R76 + R67 + R80;
    double R82 = R47*R76 + R68 - R80*a1;
    double R83 = pow(R43, -3.0/2.0);
    double R84 = R83*a1;
    double R85 = R84*b1;
    double R86 = 2*R78;
    double R87 = R41*R85;
    double R88 = R46*R84;
    double R89 = R42*R83;
    double R90 = R47*R89;
    double R91 = R48 - R87 + R88 - R90;
    double R92 = -R66*R85 + R67*R84 - R68*R89 + R69;
    double R93 = 2*R51;
    double R94 = R93*(R74*R91 + R92);
    double R95 = R93*(R49*R92 + R70*R91);
    double R96 = 2*R49*(2*R48 - 2*R87 + 2*R88 - 2*R90)/pow(R50, 2);
    double R97 = R77*R96;
    double R98 = R75*R96;
    double R99 = R44*R94;
    double R100 = R44*R97;
    double R101 = ((1.0/2.0)*R79*(R41*R95 + R41*R98 - R45*R94 - R45*R97 + R85*R86) + (1.0/2.0)*R81*(R100 + R46*R95 + R46*R98 - R84*R86 + R99) + (1.0/2.0)*R82*(-R100*a1 + R47*R95 + R47*R98 - 2*R80 + R86*R89 - R99*a1))/sqrt(pow(R79, 2) + pow(R81, 2) + pow(R82, 2));
    double result = 16.0*((R78 > 0) ? (
   R101
)
: (
   -1.0*R101
));
    return result;
}


double CosmicTrack_DCA_Deriv_b1(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = 1.0*R3*R4;
    double R6 = sin(plane_b);
    double R7 = sin(panel_g);
    double R8 = cos(panel_a);
    double R9 = R7*R8;
    double R10 = sin(panel_b);
    double R11 = cos(panel_g);
    double R12 = R11*R3;
    double R13 = R10*R12 - R9;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = 1.0/R14;
    double R16 = R15*panel_straw0x;
    double R17 = R11*R8;
    double R18 = R3*R7;
    double R19 = R10*R18 + R17;
    double R20 = R15*panel_straw0y;
    double R21 = 1.0*R13*R16 - 1.0*R19*R20;
    double R22 = 1.0*R13*R20 + 1.0*R16*R19;
    double R23 = sin(plane_a);
    double R24 = R1*R23;
    double R25 = R2*R5 - R21*R6 + R22*R24;
    double R26 = sin(plane_g);
    double R27 = R23*R26;
    double R28 = cos(plane_g);
    double R29 = R0*R28;
    double R30 = R27 + R29*R6;
    double R31 = R1*R28;
    double R32 = R0*R26;
    double R33 = R23*R28;
    double R34 = -R32 + R33*R6;
    double R35 = R21*R31 + R22*R34 + R30*R5;
    double R36 = R32*R6 - R33;
    double R37 = R1*R26;
    double R38 = R27*R6 + R29;
    double R39 = R21*R37 + R22*R38 + R36*R5;
    double R40 = pow(pow(R25, 2) + pow(R35, 2) + pow(R39, 2), -1.0/2.0);
    double R41 = R25*R40;
    double R42 = pow(b1, 2);
    double R43 = R42 + pow(a1, 2) + 1;
    double R44 = pow(R43, -1.0/2.0);
    double R45 = R41*R44;
    double R46 = R39*R40;
    double R47 = R35*R40;
    double R48 = R44*a1;
    double R49 = -R44*R46 + R45*b1 + R47*R48;
    double R50 = 1.0 - pow(R49, 2);
    double R51 = 1.0/R50;
    double R52 = panel_dz + panel_straw0z - plane_z;
    double R53 = R16*panel_dx - R20*panel_dy + panel_straw0x;
    double R54 = R16*panel_dy + R20*panel_dx + panel_straw0y;
    double R55 = -R14 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R56 = R16*R4;
    double R57 = R20*R4;
    double R58 = R11*R56 - R57*R7;
    double R59 = R11*R57 + R56*R7;
    double R60 = -panel_straw0z + wire_z;
    double R61 = R4*R8;
    double R62 = R10*R17 + R18;
    double R63 = R10*R9 - R12;
    double R64 = R16*R62 - R20*R63;
    double R65 = R16*R63 + R20*R62;
    double R66 = R2*R52 + R24*R54 - R53*R6 + R55*(-R10*R2 + R24*R59 - R58*R6) + R60*(R2*R61 + R24*R65 - R6*R64) - b0 + plane_dz + plane_z;
    double R67 = R44*R66;
    double R68 = R36*R52 + R37*R53 + R38*R54 + R55*(-R10*R36 + R37*R58 + R38*R59) + R60*(R36*R61 + R37*R64 + R38*R65) + plane_dy;
    double R69 = R30*R52 + R31*R53 + R34*R54 + R55*(-R10*R30 + R31*R58 + R34*R59) + R60*(R30*R61 + R31*R64 + R34*R65) - a0 + plane_dx;
    double R70 = -R44*R68 + R48*R69 + R67*b1;
    double R71 = R41*R66;
    double R72 = R46*R68;
    double R73 = R47*R69;
    double R74 = -R71 - R72 - R73;
    double R75 = R49*R70 + R74;
    double R76 = R51*R75;
    double R77 = -R49*(R71 + R72 + R73) + R70;
    double R78 = R51*R77;
    double R79 = R44*R78;
    double R80 = R41*R76 + R66 - R79*b1;
    double R81 = R46*R76 + R68 + R79;
    double R82 = R47*R76 - R48*R78 + R69;
    double R83 = pow(R43, -3.0/2.0);
    double R84 = R83*b1;
    double R85 = 2*R78;
    double R86 = R42*R83;
    double R87 = R41*R86;
    double R88 = R46*R84;
    double R89 = R84*a1;
    double R90 = R47*R89;
    double R91 = R45 - R87 + R88 - R90;
    double R92 = -R66*R86 + R67 + R68*R84 - R69*R89;
    double R93 = 2*R51;
    double R94 = R93*(R74*R91 + R92);
    double R95 = R44*R94;
    double R96 = R93*(R49*R92 + R70*R91);
    double R97 = 2*R49*(2*R45 - 2*R87 + 2*R88 - 2*R90)/pow(R50, 2);
    double R98 = R77*R97;
    double R99 = R44*R98;
    double R100 = R75*R97;
    double R101 = ((1.0/2.0)*R80*(R100*R41 + R41*R96 - 2*R79 + R85*R86 - R95*b1 - R99*b1) + (1.0/2.0)*R81*(R100*R46 + R46*R96 - R84*R85 + R95 + R99) + (1.0/2.0)*R82*(R100*R47 + R47*R96 - R48*R94 - R48*R98 + R85*R89))/sqrt(pow(R80, 2) + pow(R81, 2) + pow(R82, 2));
    double result = 16.0*((R78 > 0) ? (
   R101
)
: (
   -1.0*R101
));
    return result;
}


double CosmicTrack_DCA_Deriv_t0(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double result = -1;
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dx(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = 1.0*R3*R4;
    double R6 = sin(plane_b);
    double R7 = sin(panel_g);
    double R8 = cos(panel_a);
    double R9 = R7*R8;
    double R10 = sin(panel_b);
    double R11 = cos(panel_g);
    double R12 = R11*R3;
    double R13 = R10*R12 - R9;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = 1.0/R14;
    double R16 = R15*panel_straw0x;
    double R17 = R11*R8;
    double R18 = R3*R7;
    double R19 = R10*R18 + R17;
    double R20 = R15*panel_straw0y;
    double R21 = 1.0*R13*R16 - 1.0*R19*R20;
    double R22 = 1.0*R13*R20 + 1.0*R16*R19;
    double R23 = sin(plane_a);
    double R24 = R1*R23;
    double R25 = R2*R5 - R21*R6 + R22*R24;
    double R26 = sin(plane_g);
    double R27 = R23*R26;
    double R28 = cos(plane_g);
    double R29 = R0*R28;
    double R30 = R27 + R29*R6;
    double R31 = R1*R28;
    double R32 = R0*R26;
    double R33 = R23*R28;
    double R34 = -R32 + R33*R6;
    double R35 = R21*R31 + R22*R34 + R30*R5;
    double R36 = R32*R6 - R33;
    double R37 = R1*R26;
    double R38 = R27*R6 + R29;
    double R39 = R21*R37 + R22*R38 + R36*R5;
    double R40 = pow(pow(R25, 2) + pow(R35, 2) + pow(R39, 2), -1.0/2.0);
    double R41 = R25*R40;
    double R42 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R43 = R42*b1;
    double R44 = R39*R40;
    double R45 = R42*a1;
    double R46 = R35*R40;
    double R47 = R41*R43 - R42*R44 + R45*R46;
    double R48 = 1.0/(1.0 - pow(R47, 2));
    double R49 = panel_dz + panel_straw0z - plane_z;
    double R50 = R16*panel_dx - R20*panel_dy + panel_straw0x;
    double R51 = R16*panel_dy + R20*panel_dx + panel_straw0y;
    double R52 = -R14 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R53 = R16*R4;
    double R54 = R20*R4;
    double R55 = R11*R53 - R54*R7;
    double R56 = R11*R54 + R53*R7;
    double R57 = -panel_straw0z + wire_z;
    double R58 = R4*R8;
    double R59 = R10*R17 + R18;
    double R60 = R10*R9 - R12;
    double R61 = R16*R59 - R20*R60;
    double R62 = R16*R60 + R20*R59;
    double R63 = R2*R49 + R24*R51 - R50*R6 + R52*(-R10*R2 + R24*R56 - R55*R6) + R57*(R2*R58 + R24*R62 - R6*R61) - b0 + plane_dz + plane_z;
    double R64 = R41*R63;
    double R65 = R36*R49 + R37*R50 + R38*R51 + R52*(-R10*R36 + R37*R55 + R38*R56) + R57*(R36*R58 + R37*R61 + R38*R62) + plane_dy;
    double R66 = R44*R65;
    double R67 = R30*R49 + R31*R50 + R34*R51 + R52*(-R10*R30 + R31*R55 + R34*R56) + R57*(R30*R58 + R31*R61 + R34*R62) - a0 + plane_dx;
    double R68 = R46*R67;
    double R69 = -R42*R65 + R43*R63 + R45*R67;
    double R70 = R48*(R47*R69 - R64 - R66 - R68);
    double R71 = R48*(-R47*(R64 + R66 + R68) + R69);
    double R72 = R41*R70 - R43*R71 + R63;
    double R73 = R42*R71 + R44*R70 + R65;
    double R74 = -R45*R71 + R46*R70 + R67;
    double R75 = 2*R48;
    double R76 = R75*(R45 - R46*R47);
    double R77 = R75*(R45*R47 - R46);
    double R78 = ((1.0/2.0)*R72*(R41*R77 - R43*R76) + (1.0/2.0)*R73*(R42*R76 + R44*R77) + (1.0/2.0)*R74*(-R45*R76 + R46*R77 + 2))/sqrt(pow(R72, 2) + pow(R73, 2) + pow(R74, 2));
    double result = 16.0*((R71 > 0) ? (
   R78
)
: (
   -1.0*R78
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dy(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = 1.0*R3*R4;
    double R6 = sin(plane_b);
    double R7 = sin(panel_g);
    double R8 = cos(panel_a);
    double R9 = R7*R8;
    double R10 = sin(panel_b);
    double R11 = cos(panel_g);
    double R12 = R11*R3;
    double R13 = R10*R12 - R9;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = 1.0/R14;
    double R16 = R15*panel_straw0x;
    double R17 = R11*R8;
    double R18 = R3*R7;
    double R19 = R10*R18 + R17;
    double R20 = R15*panel_straw0y;
    double R21 = 1.0*R13*R16 - 1.0*R19*R20;
    double R22 = 1.0*R13*R20 + 1.0*R16*R19;
    double R23 = sin(plane_a);
    double R24 = R1*R23;
    double R25 = R2*R5 - R21*R6 + R22*R24;
    double R26 = sin(plane_g);
    double R27 = R23*R26;
    double R28 = cos(plane_g);
    double R29 = R0*R28;
    double R30 = R27 + R29*R6;
    double R31 = R1*R28;
    double R32 = R0*R26;
    double R33 = R23*R28;
    double R34 = -R32 + R33*R6;
    double R35 = R21*R31 + R22*R34 + R30*R5;
    double R36 = R32*R6 - R33;
    double R37 = R1*R26;
    double R38 = R27*R6 + R29;
    double R39 = R21*R37 + R22*R38 + R36*R5;
    double R40 = pow(pow(R25, 2) + pow(R35, 2) + pow(R39, 2), -1.0/2.0);
    double R41 = R25*R40;
    double R42 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R43 = R42*b1;
    double R44 = R39*R40;
    double R45 = R35*R40;
    double R46 = R42*a1;
    double R47 = R41*R43 - R42*R44 + R45*R46;
    double R48 = 1.0/(1.0 - pow(R47, 2));
    double R49 = panel_dz + panel_straw0z - plane_z;
    double R50 = R16*panel_dx - R20*panel_dy + panel_straw0x;
    double R51 = R16*panel_dy + R20*panel_dx + panel_straw0y;
    double R52 = -R14 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R53 = R16*R4;
    double R54 = R20*R4;
    double R55 = R11*R53 - R54*R7;
    double R56 = R11*R54 + R53*R7;
    double R57 = -panel_straw0z + wire_z;
    double R58 = R4*R8;
    double R59 = R10*R17 + R18;
    double R60 = R10*R9 - R12;
    double R61 = R16*R59 - R20*R60;
    double R62 = R16*R60 + R20*R59;
    double R63 = R2*R49 + R24*R51 - R50*R6 + R52*(-R10*R2 + R24*R56 - R55*R6) + R57*(R2*R58 + R24*R62 - R6*R61) - b0 + plane_dz + plane_z;
    double R64 = R41*R63;
    double R65 = R36*R49 + R37*R50 + R38*R51 + R52*(-R10*R36 + R37*R55 + R38*R56) + R57*(R36*R58 + R37*R61 + R38*R62) + plane_dy;
    double R66 = R44*R65;
    double R67 = R30*R49 + R31*R50 + R34*R51 + R52*(-R10*R30 + R31*R55 + R34*R56) + R57*(R30*R58 + R31*R61 + R34*R62) - a0 + plane_dx;
    double R68 = R45*R67;
    double R69 = -R42*R65 + R43*R63 + R46*R67;
    double R70 = R48*(R47*R69 - R64 - R66 - R68);
    double R71 = R48*(-R47*(R64 + R66 + R68) + R69);
    double R72 = R41*R70 - R43*R71 + R63;
    double R73 = R42*R71 + R44*R70 + R65;
    double R74 = R45*R70 - R46*R71 + R67;
    double R75 = 2*R48;
    double R76 = R75*(-R42 - R44*R47);
    double R77 = R75*(-R42*R47 - R44);
    double R78 = ((1.0/2.0)*R72*(R41*R77 - R43*R76) + (1.0/2.0)*R73*(R42*R76 + R44*R77 + 2) + (1.0/2.0)*R74*(R45*R77 - R46*R76))/sqrt(pow(R72, 2) + pow(R73, 2) + pow(R74, 2));
    double result = 16.0*((R71 > 0) ? (
   R78
)
: (
   -1.0*R78
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dz(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = 1.0*R3*R4;
    double R6 = sin(plane_b);
    double R7 = sin(panel_g);
    double R8 = cos(panel_a);
    double R9 = R7*R8;
    double R10 = sin(panel_b);
    double R11 = cos(panel_g);
    double R12 = R11*R3;
    double R13 = R10*R12 - R9;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = 1.0/R14;
    double R16 = R15*panel_straw0x;
    double R17 = R11*R8;
    double R18 = R3*R7;
    double R19 = R10*R18 + R17;
    double R20 = R15*panel_straw0y;
    double R21 = 1.0*R13*R16 - 1.0*R19*R20;
    double R22 = 1.0*R13*R20 + 1.0*R16*R19;
    double R23 = sin(plane_a);
    double R24 = R1*R23;
    double R25 = R2*R5 - R21*R6 + R22*R24;
    double R26 = sin(plane_g);
    double R27 = R23*R26;
    double R28 = cos(plane_g);
    double R29 = R0*R28;
    double R30 = R27 + R29*R6;
    double R31 = R1*R28;
    double R32 = R0*R26;
    double R33 = R23*R28;
    double R34 = -R32 + R33*R6;
    double R35 = R21*R31 + R22*R34 + R30*R5;
    double R36 = R32*R6 - R33;
    double R37 = R1*R26;
    double R38 = R27*R6 + R29;
    double R39 = R21*R37 + R22*R38 + R36*R5;
    double R40 = pow(pow(R25, 2) + pow(R35, 2) + pow(R39, 2), -1.0/2.0);
    double R41 = R25*R40;
    double R42 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R43 = R42*b1;
    double R44 = R39*R40;
    double R45 = R35*R40;
    double R46 = R42*a1;
    double R47 = R41*R43 - R42*R44 + R45*R46;
    double R48 = 1.0/(1.0 - pow(R47, 2));
    double R49 = panel_dz + panel_straw0z - plane_z;
    double R50 = R16*panel_dx - R20*panel_dy + panel_straw0x;
    double R51 = R16*panel_dy + R20*panel_dx + panel_straw0y;
    double R52 = -R14 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R53 = R16*R4;
    double R54 = R20*R4;
    double R55 = R11*R53 - R54*R7;
    double R56 = R11*R54 + R53*R7;
    double R57 = -panel_straw0z + wire_z;
    double R58 = R4*R8;
    double R59 = R10*R17 + R18;
    double R60 = R10*R9 - R12;
    double R61 = R16*R59 - R20*R60;
    double R62 = R16*R60 + R20*R59;
    double R63 = R2*R49 + R24*R51 - R50*R6 + R52*(-R10*R2 + R24*R56 - R55*R6) + R57*(R2*R58 + R24*R62 - R6*R61) - b0 + plane_dz + plane_z;
    double R64 = R41*R63;
    double R65 = R36*R49 + R37*R50 + R38*R51 + R52*(-R10*R36 + R37*R55 + R38*R56) + R57*(R36*R58 + R37*R61 + R38*R62) + plane_dy;
    double R66 = R44*R65;
    double R67 = R30*R49 + R31*R50 + R34*R51 + R52*(-R10*R30 + R31*R55 + R34*R56) + R57*(R30*R58 + R31*R61 + R34*R62) - a0 + plane_dx;
    double R68 = R45*R67;
    double R69 = -R42*R65 + R43*R63 + R46*R67;
    double R70 = R48*(R47*R69 - R64 - R66 - R68);
    double R71 = R48*(-R47*(R64 + R66 + R68) + R69);
    double R72 = R41*R70 - R43*R71 + R63;
    double R73 = R42*R71 + R44*R70 + R65;
    double R74 = R45*R70 - R46*R71 + R67;
    double R75 = 2*R48;
    double R76 = R75*(-R41*R47 + R43);
    double R77 = R75*(-R41 + R43*R47);
    double R78 = ((1.0/2.0)*R72*(R41*R77 - R43*R76 + 2) + (1.0/2.0)*R73*(R42*R76 + R44*R77) + (1.0/2.0)*R74*(R45*R77 - R46*R76))/sqrt(pow(R72, 2) + pow(R73, 2) + pow(R74, 2));
    double result = 16.0*((R71 > 0) ? (
   R78
)
: (
   -1.0*R78
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_a(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = R3*R4;
    double R6 = 1.0*R5;
    double R7 = sin(plane_b);
    double R8 = sin(panel_g);
    double R9 = cos(panel_a);
    double R10 = R8*R9;
    double R11 = sin(panel_b);
    double R12 = cos(panel_g);
    double R13 = R12*R3;
    double R14 = -R10 + R11*R13;
    double R15 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R16 = 1.0/R15;
    double R17 = R16*panel_straw0x;
    double R18 = R12*R9;
    double R19 = R3*R8;
    double R20 = R11*R19 + R18;
    double R21 = R16*panel_straw0y;
    double R22 = 1.0*R14*R17 - 1.0*R20*R21;
    double R23 = R17*R20;
    double R24 = R14*R21;
    double R25 = R23 + R24;
    double R26 = 1.0*R25;
    double R27 = sin(plane_a);
    double R28 = R1*R27;
    double R29 = R2*R6 - R22*R7 + R26*R28;
    double R30 = sin(plane_g);
    double R31 = R27*R30;
    double R32 = cos(plane_g);
    double R33 = R0*R32;
    double R34 = R33*R7;
    double R35 = R31 + R34;
    double R36 = R1*R32;
    double R37 = R0*R30;
    double R38 = R27*R32;
    double R39 = R38*R7;
    double R40 = -R37 + R39;
    double R41 = R22*R36 + R26*R40 + R35*R6;
    double R42 = R37*R7;
    double R43 = -R38 + R42;
    double R44 = R1*R30;
    double R45 = R31*R7;
    double R46 = R33 + R45;
    double R47 = R22*R44 + R26*R46 + R43*R6;
    double R48 = pow(R29, 2) + pow(R41, 2) + pow(R47, 2);
    double R49 = pow(R48, -1.0/2.0);
    double R50 = R29*R49;
    double R51 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R52 = R51*b1;
    double R53 = R47*R49;
    double R54 = R51*a1;
    double R55 = R41*R49;
    double R56 = R50*R52 - R51*R53 + R54*R55;
    double R57 = 1.0 - pow(R56, 2);
    double R58 = 1.0/R57;
    double R59 = panel_dz + panel_straw0z - plane_z;
    double R60 = R17*panel_dx - R21*panel_dy + panel_straw0x;
    double R61 = R17*panel_dy + R21*panel_dx + panel_straw0y;
    double R62 = -R15 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R63 = R17*R4;
    double R64 = R21*R4;
    double R65 = R12*R63 - R64*R8;
    double R66 = R12*R64 + R63*R8;
    double R67 = -panel_straw0z + wire_z;
    double R68 = R4*R9;
    double R69 = R11*R18 + R19;
    double R70 = R10*R11 - R13;
    double R71 = R17*R69 - R21*R70;
    double R72 = R17*R70 + R21*R69;
    double R73 = R2*R59 + R28*R61 - R60*R7 + R62*(-R11*R2 + R28*R66 - R65*R7) + R67*(R2*R68 + R28*R72 - R7*R71) - b0 + plane_dz + plane_z;
    double R74 = R43*R59 + R44*R60 + R46*R61 + R62*(-R11*R43 + R44*R65 + R46*R66) + R67*(R43*R68 + R44*R71 + R46*R72) + plane_dy;
    double R75 = R35*R59 + R36*R60 + R40*R61 + R62*(-R11*R35 + R36*R65 + R40*R66) + R67*(R35*R68 + R36*R71 + R40*R72) - a0 + plane_dx;
    double R76 = -R51*R74 + R52*R73 + R54*R75;
    double R77 = R49*R73;
    double R78 = R29*R77;
    double R79 = R53*R74;
    double R80 = R49*R75;
    double R81 = R41*R80;
    double R82 = -R78 - R79 - R81;
    double R83 = R56*R76 + R82;
    double R84 = R58*R83;
    double R85 = -R56*(R78 + R79 + R81) + R76;
    double R86 = R58*R85;
    double R87 = R50*R84 - R52*R86 + R73;
    double R88 = R51*R86 + R53*R84 + R74;
    double R89 = -R54*R86 + R55*R84 + R75;
    double R90 = R28*R59;
    double R91 = R2*R61;
    double R92 = R62*(R11*R28 + R2*R66);
    double R93 = R67*(R2*R72 - R28*R68);
    double R94 = R2*(1.0*R23 + 1.0*R24);
    double R95 = -R28*R6 + R94;
    double R96 = R49*R95;
    double R97 = 2*R84;
    double R98 = R5*(1.0*R37 - 1.0*R39);
    double R99 = R25*(1.0*R31 + 1.0*R34);
    double R100 = R5*(-1.0*R33 - 1.0*R45);
    double R101 = R25*(-1.0*R38 + 1.0*R42);
    double R102 = (-1.0/2.0*R29*(-2.0*R28*R5 + 2*R94) - 1.0/2.0*R41*(2*R98 + 2*R99) - 1.0/2.0*R47*(2*R100 + 2*R101))/pow(R48, 3.0/2.0);
    double R103 = R102*R29;
    double R104 = R52*R96;
    double R105 = R49*(R100 + R101);
    double R106 = R105*R51;
    double R107 = R98 + R99;
    double R108 = R107*R49;
    double R109 = R108*R54;
    double R110 = R103*R52;
    double R111 = R102*R47;
    double R112 = R111*R51;
    double R113 = R102*R41;
    double R114 = R113*R54;
    double R115 = 2*R56*(2*R104 - 2*R106 + 2*R109 + 2*R110 - 2*R112 + 2*R114)/pow(R57, 2);
    double R116 = R115*R85;
    double R117 = R115*R83;
    double R118 = -R90 + R91 + R92 + R93;
    double R119 = R37 - R39;
    double R120 = R119*R59;
    double R121 = R35*R61;
    double R122 = R62*(R11*R40 + R35*R66);
    double R123 = R67*(R119*R68 + R35*R72);
    double R124 = R120 + R121 + R122 + R123;
    double R125 = -R33 - R45;
    double R126 = R125*R59;
    double R127 = R43*R61;
    double R128 = R62*(R11*R46 + R43*R66);
    double R129 = R67*(R125*R68 + R43*R72);
    double R130 = R126 + R127 + R128 + R129;
    double R131 = R118*R52 + R124*R54 - R130*R51;
    double R132 = R104 - R106 + R109 + R110 - R112 + R114;
    double R133 = -R103*R73 - R105*R74 - R107*R80 - R111*R74 - R113*R75 - R118*R50 - R124*R55 - R130*R53 - R77*R95;
    double R134 = 2*R58;
    double R135 = R134*(R131*R56 + R132*R76 + R133);
    double R136 = R134*(R131 + R132*R82 + R133*R56);
    double R137 = ((1.0/2.0)*R87*(R103*R97 - R116*R52 + R117*R50 + R135*R50 - R136*R52 - 2*R90 + 2*R91 + 2*R92 + 2*R93 + R96*R97) + (1.0/2.0)*R88*(R105*R97 + R111*R97 + R116*R51 + R117*R53 + 2*R126 + 2*R127 + 2*R128 + 2*R129 + R135*R53 + R136*R51) + (1.0/2.0)*R89*(R108*R97 + R113*R97 - R116*R54 + R117*R55 + 2*R120 + 2*R121 + 2*R122 + 2*R123 + R135*R55 - R136*R54))/sqrt(pow(R87, 2) + pow(R88, 2) + pow(R89, 2));
    double result = 16.0*((R86 > 0) ? (
   R137
)
: (
   -1.0*R137
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_b(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = R3*R4;
    double R6 = 1.0*R5;
    double R7 = sin(plane_b);
    double R8 = sin(panel_g);
    double R9 = cos(panel_a);
    double R10 = R8*R9;
    double R11 = sin(panel_b);
    double R12 = cos(panel_g);
    double R13 = R12*R3;
    double R14 = -R10 + R11*R13;
    double R15 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R16 = 1.0/R15;
    double R17 = R16*panel_straw0x;
    double R18 = R14*R17;
    double R19 = R12*R9;
    double R20 = R3*R8;
    double R21 = R11*R20 + R19;
    double R22 = R16*panel_straw0y;
    double R23 = R21*R22;
    double R24 = 1.0*R18 - 1.0*R23;
    double R25 = R17*R21;
    double R26 = R14*R22;
    double R27 = R25 + R26;
    double R28 = 1.0*R27;
    double R29 = sin(plane_a);
    double R30 = R1*R29;
    double R31 = R2*R6 - R24*R7 + R28*R30;
    double R32 = sin(plane_g);
    double R33 = R29*R32;
    double R34 = cos(plane_g);
    double R35 = R0*R34;
    double R36 = R33 + R35*R7;
    double R37 = R1*R34;
    double R38 = R0*R32;
    double R39 = R29*R34;
    double R40 = -R38 + R39*R7;
    double R41 = R24*R37 + R28*R40 + R36*R6;
    double R42 = R38*R7 - R39;
    double R43 = R1*R32;
    double R44 = R33*R7 + R35;
    double R45 = R24*R43 + R28*R44 + R42*R6;
    double R46 = pow(R31, 2) + pow(R41, 2) + pow(R45, 2);
    double R47 = pow(R46, -1.0/2.0);
    double R48 = R31*R47;
    double R49 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R50 = R49*b1;
    double R51 = R45*R47;
    double R52 = R41*R47;
    double R53 = R49*a1;
    double R54 = R48*R50 - R49*R51 + R52*R53;
    double R55 = 1.0 - pow(R54, 2);
    double R56 = 1.0/R55;
    double R57 = panel_dz + panel_straw0z - plane_z;
    double R58 = R17*panel_dx;
    double R59 = R22*panel_dy;
    double R60 = R58 - R59 + panel_straw0x;
    double R61 = R60*R7;
    double R62 = R17*panel_dy + R22*panel_dx + panel_straw0y;
    double R63 = -R15 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R64 = R17*R4;
    double R65 = R12*R64;
    double R66 = R22*R4;
    double R67 = R66*R8;
    double R68 = R65 - R67;
    double R69 = R68*R7;
    double R70 = R12*R66 + R64*R8;
    double R71 = -panel_straw0z + wire_z;
    double R72 = R4*R9;
    double R73 = R11*R19 + R20;
    double R74 = R17*R73;
    double R75 = R10*R11 - R13;
    double R76 = R22*R75;
    double R77 = R74 - R76;
    double R78 = R7*R77;
    double R79 = R17*R75 + R22*R73;
    double R80 = R2*R57 + R30*R62 - R61 + R63*(-R11*R2 + R30*R70 - R69) + R71*(R2*R72 + R30*R79 - R78) - b0 + plane_dz + plane_z;
    double R81 = R42*R57 + R43*R60 + R44*R62 + R63*(-R11*R42 + R43*R68 + R44*R70) + R71*(R42*R72 + R43*R77 + R44*R79) + plane_dy;
    double R82 = R36*R57 + R37*R60 + R40*R62 + R63*(-R11*R36 + R37*R68 + R40*R70) + R71*(R36*R72 + R37*R77 + R40*R79) - a0 + plane_dx;
    double R83 = -R49*R81 + R50*R80 + R53*R82;
    double R84 = R47*R80;
    double R85 = R31*R84;
    double R86 = R51*R81;
    double R87 = R52*R82;
    double R88 = -R85 - R86 - R87;
    double R89 = R54*R83 + R88;
    double R90 = R56*R89;
    double R91 = -R54*(R85 + R86 + R87) + R83;
    double R92 = R56*R91;
    double R93 = R48*R90 - R50*R92 + R80;
    double R94 = R49*R92 + R51*R90 + R81;
    double R95 = R52*R90 - R53*R92 + R82;
    double R96 = R0*R7;
    double R97 = R57*R96;
    double R98 = R29*R7;
    double R99 = R62*R98;
    double R100 = R1*(-R58 + R59 - panel_straw0x);
    double R101 = R63*(R1*(-R65 + R67) + R11*R96 - R70*R98);
    double R102 = R71*(R1*(-R74 + R76) - R72*R96 - R79*R98);
    double R103 = 1.0*R23;
    double R104 = 1.0*R18;
    double R105 = R1*(R103 - R104);
    double R106 = R98*(1.0*R25 + 1.0*R26);
    double R107 = R105 - R106 - R6*R96;
    double R108 = R107*R47;
    double R109 = 2*R90;
    double R110 = 2.0*R5;
    double R111 = R1*R35;
    double R112 = R7*(-R103 + R104);
    double R113 = R112*R34;
    double R114 = R1*R39;
    double R115 = 2.0*R27;
    double R116 = R1*R38;
    double R117 = R112*R32;
    double R118 = R1*R33;
    double R119 = (-1.0/2.0*R31*(2*R105 - 2*R106 - R110*R96) - 1.0/2.0*R41*(R110*R111 - 2*R113 + R114*R115) - 1.0/2.0*R45*(R110*R116 + R115*R118 - 2*R117))/pow(R46, 3.0/2.0);
    double R120 = R119*R31;
    double R121 = R108*R50;
    double R122 = R47*(R116*R6 - R117 + R118*R28);
    double R123 = R122*R49;
    double R124 = R47*(R111*R6 - R113 + R114*R28);
    double R125 = R124*R53;
    double R126 = R120*R50;
    double R127 = R119*R45;
    double R128 = R127*R49;
    double R129 = R119*R41;
    double R130 = R129*R53;
    double R131 = 2*R54*(2*R121 - 2*R123 + 2*R125 + 2*R126 - 2*R128 + 2*R130)/pow(R55, 2);
    double R132 = R131*R91;
    double R133 = R131*R89;
    double R134 = R100 + R101 + R102 - R97 - R99;
    double R135 = R116*R57;
    double R136 = R32*R61;
    double R137 = R118*R62;
    double R138 = R63*(-R11*R116 + R118*R70 - R32*R69);
    double R139 = R71*(R116*R72 + R118*R79 - R32*R78);
    double R140 = R135 - R136 + R137 + R138 + R139;
    double R141 = R111*R57;
    double R142 = R34*R61;
    double R143 = R114*R62;
    double R144 = R63*(-R11*R111 + R114*R70 - R34*R69);
    double R145 = R71*(R111*R72 + R114*R79 - R34*R78);
    double R146 = R141 - R142 + R143 + R144 + R145;
    double R147 = R134*R50 - R140*R49 + R146*R53;
    double R148 = R121 - R123 + R125 + R126 - R128 + R130;
    double R149 = -R107*R84 - R120*R80 - R122*R81 - R124*R82 - R127*R81 - R129*R82 - R134*R48 - R140*R51 - R146*R52;
    double R150 = 2*R56;
    double R151 = R150*(R147*R54 + R148*R83 + R149);
    double R152 = R150*(R147 + R148*R88 + R149*R54);
    double R153 = ((1.0/2.0)*R93*(2*R100 + 2*R101 + 2*R102 + R108*R109 + R109*R120 - R132*R50 + R133*R48 + R151*R48 - R152*R50 - 2*R97 - 2*R99) + (1.0/2.0)*R94*(R109*R122 + R109*R127 + R132*R49 + R133*R51 + 2*R135 - 2*R136 + 2*R137 + 2*R138 + 2*R139 + R151*R51 + R152*R49) + (1.0/2.0)*R95*(R109*R124 + R109*R129 - R132*R53 + R133*R52 + 2*R141 - 2*R142 + 2*R143 + 2*R144 + 2*R145 + R151*R52 - R152*R53))/sqrt(pow(R93, 2) + pow(R94, 2) + pow(R95, 2));
    double result = 16.0*((R92 > 0) ? (
   R153
)
: (
   -1.0*R153
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_g(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = R3*R4;
    double R6 = 1.0*R5;
    double R7 = sin(plane_b);
    double R8 = sin(panel_g);
    double R9 = cos(panel_a);
    double R10 = R8*R9;
    double R11 = sin(panel_b);
    double R12 = cos(panel_g);
    double R13 = R12*R3;
    double R14 = -R10 + R11*R13;
    double R15 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R16 = 1.0/R15;
    double R17 = R16*panel_straw0x;
    double R18 = R14*R17;
    double R19 = R12*R9;
    double R20 = R3*R8;
    double R21 = R11*R20 + R19;
    double R22 = R16*panel_straw0y;
    double R23 = R21*R22;
    double R24 = 1.0*R18 - 1.0*R23;
    double R25 = R14*R22 + R17*R21;
    double R26 = 1.0*R25;
    double R27 = sin(plane_a);
    double R28 = R1*R27;
    double R29 = R2*R6 - R24*R7 + R26*R28;
    double R30 = sin(plane_g);
    double R31 = R27*R30;
    double R32 = cos(plane_g);
    double R33 = R0*R32;
    double R34 = R33*R7;
    double R35 = R31 + R34;
    double R36 = R1*R32;
    double R37 = R0*R30;
    double R38 = R27*R32;
    double R39 = R38*R7;
    double R40 = -R37 + R39;
    double R41 = R24*R36 + R26*R40 + R35*R6;
    double R42 = R37*R7;
    double R43 = -R38 + R42;
    double R44 = R1*R30;
    double R45 = R31*R7;
    double R46 = R33 + R45;
    double R47 = R24*R44 + R26*R46 + R43*R6;
    double R48 = pow(R29, 2) + pow(R41, 2) + pow(R47, 2);
    double R49 = pow(R48, -1.0/2.0);
    double R50 = R29*R49;
    double R51 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R52 = R51*b1;
    double R53 = R47*R49;
    double R54 = R51*a1;
    double R55 = R41*R49;
    double R56 = R50*R52 - R51*R53 + R54*R55;
    double R57 = 1.0 - pow(R56, 2);
    double R58 = 1.0/R57;
    double R59 = panel_dz + panel_straw0z - plane_z;
    double R60 = R17*panel_dx - R22*panel_dy + panel_straw0x;
    double R61 = R17*panel_dy + R22*panel_dx + panel_straw0y;
    double R62 = -R15 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R63 = R17*R4;
    double R64 = R22*R4;
    double R65 = R12*R63 - R64*R8;
    double R66 = R12*R64 + R63*R8;
    double R67 = -panel_straw0z + wire_z;
    double R68 = R4*R9;
    double R69 = R11*R19 + R20;
    double R70 = R10*R11 - R13;
    double R71 = R17*R69 - R22*R70;
    double R72 = R17*R70 + R22*R69;
    double R73 = R2*R59 + R28*R61 - R60*R7 + R62*(-R11*R2 + R28*R66 - R65*R7) + R67*(R2*R68 + R28*R72 - R7*R71) - b0 + plane_dz + plane_z;
    double R74 = R44*R60;
    double R75 = R11*R43;
    double R76 = R44*R65;
    double R77 = R44*R71;
    double R78 = R43*R59 + R46*R61 + R62*(R46*R66 - R75 + R76) + R67*(R43*R68 + R46*R72 + R77) + R74 + plane_dy;
    double R79 = R36*R65 + R40*R66;
    double R80 = R67*(R35*R68 + R36*R71 + R40*R72);
    double R81 = R35*R59;
    double R82 = R40*R61;
    double R83 = R36*R60;
    double R84 = R80 + R81 + R82 + R83;
    double R85 = R62*(-R11*R35 + R79) + R84 - a0 + plane_dx;
    double R86 = -R51*R78 + R52*R73 + R54*R85;
    double R87 = R50*R73;
    double R88 = R53*R78;
    double R89 = R49*R85;
    double R90 = R41*R89;
    double R91 = -R87 - R88 - R90;
    double R92 = R56*R86 + R91;
    double R93 = R58*R92;
    double R94 = -R56*(R87 + R88 + R90) + R86;
    double R95 = R58*R94;
    double R96 = R50*R93 - R52*R95 + R73;
    double R97 = R51*R95 + R53*R93 + R78;
    double R98 = -R54*R95 + R55*R93 + R85;
    double R99 = R5*(1.0*R38 - 1.0*R42);
    double R100 = 1.0*R18 - 1.0*R23;
    double R101 = R100*R44;
    double R102 = R25*(-1.0*R33 - 1.0*R45);
    double R103 = R5*(1.0*R31 + 1.0*R34);
    double R104 = R100*R36;
    double R105 = R25*(-1.0*R37 + 1.0*R39);
    double R106 = (-1.0/2.0*R41*(-2*R101 + 2*R102 + 2*R99) - 1.0/2.0*R47*(2*R103 + 2*R104 + 2*R105))/pow(R48, 3.0/2.0);
    double R107 = R106*R29;
    double R108 = 2*R93;
    double R109 = R49*(R103 + R104 + R105);
    double R110 = R109*R51;
    double R111 = -R101 + R102 + R99;
    double R112 = R111*R49;
    double R113 = R112*R54;
    double R114 = R107*R52;
    double R115 = R106*R47;
    double R116 = R115*R51;
    double R117 = R106*R41;
    double R118 = R117*R54;
    double R119 = 2*R56*(-2*R110 + 2*R113 + 2*R114 - 2*R116 + 2*R118)/pow(R57, 2);
    double R120 = R119*R94;
    double R121 = R119*R92;
    double R122 = R62*(R11*(-R31 - R34) + R79);
    double R123 = R122 + R84;
    double R124 = R38 - R42;
    double R125 = R124*R59;
    double R126 = -R33 - R45;
    double R127 = R126*R61;
    double R128 = R62*(R126*R66 + R75 - R76);
    double R129 = R67*(R124*R68 + R126*R72 - R77);
    double R130 = R125 + R127 + R128 + R129 - R74;
    double R131 = -R123*R51 + R130*R54;
    double R132 = -R110 + R113 + R114 - R116 + R118;
    double R133 = -R107*R73 - R109*R78 - R111*R89 - R115*R78 - R117*R85 - R123*R53 - R130*R55;
    double R134 = 2*R58;
    double R135 = R134*(R131*R56 + R132*R86 + R133);
    double R136 = R134*(R131 + R132*R91 + R133*R56);
    double R137 = ((1.0/2.0)*R96*(R107*R108 - R120*R52 + R121*R50 + R135*R50 - R136*R52) + (1.0/2.0)*R97*(R108*R109 + R108*R115 + R120*R51 + R121*R53 + 2*R122 + R135*R53 + R136*R51 + 2*R80 + 2*R81 + 2*R82 + 2*R83) + (1.0/2.0)*R98*(R108*R112 + R108*R117 - R120*R54 + R121*R55 + 2*R125 + 2*R127 + 2*R128 + 2*R129 + R135*R55 - R136*R54 - 2*R74))/sqrt(pow(R96, 2) + pow(R97, 2) + pow(R98, 2));
    double result = 16.0*((R95 > 0) ? (
   R137
)
: (
   -1.0*R137
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dx(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = 1.0*R3*R4;
    double R6 = sin(plane_b);
    double R7 = sin(panel_g);
    double R8 = cos(panel_a);
    double R9 = R7*R8;
    double R10 = sin(panel_b);
    double R11 = cos(panel_g);
    double R12 = R11*R3;
    double R13 = R10*R12 - R9;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = 1.0/R14;
    double R16 = R15*panel_straw0x;
    double R17 = R11*R8;
    double R18 = R3*R7;
    double R19 = R10*R18 + R17;
    double R20 = R15*panel_straw0y;
    double R21 = 1.0*R13*R16 - 1.0*R19*R20;
    double R22 = 1.0*R13*R20 + 1.0*R16*R19;
    double R23 = sin(plane_a);
    double R24 = R1*R23;
    double R25 = R2*R5 - R21*R6 + R22*R24;
    double R26 = sin(plane_g);
    double R27 = R23*R26;
    double R28 = cos(plane_g);
    double R29 = R0*R28;
    double R30 = R27 + R29*R6;
    double R31 = R1*R28;
    double R32 = R0*R26;
    double R33 = R23*R28;
    double R34 = -R32 + R33*R6;
    double R35 = R21*R31 + R22*R34 + R30*R5;
    double R36 = R32*R6 - R33;
    double R37 = R1*R26;
    double R38 = R27*R6 + R29;
    double R39 = R21*R37 + R22*R38 + R36*R5;
    double R40 = pow(pow(R25, 2) + pow(R35, 2) + pow(R39, 2), -1.0/2.0);
    double R41 = R25*R40;
    double R42 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R43 = R42*b1;
    double R44 = R39*R40;
    double R45 = R35*R40;
    double R46 = R42*a1;
    double R47 = R41*R43 - R42*R44 + R45*R46;
    double R48 = 1.0/(1.0 - pow(R47, 2));
    double R49 = panel_dz + panel_straw0z - plane_z;
    double R50 = R16*panel_dx - R20*panel_dy + panel_straw0x;
    double R51 = R16*panel_dy + R20*panel_dx + panel_straw0y;
    double R52 = -R14 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R53 = R16*R4;
    double R54 = R20*R4;
    double R55 = R11*R53 - R54*R7;
    double R56 = R11*R54 + R53*R7;
    double R57 = -panel_straw0z + wire_z;
    double R58 = R4*R8;
    double R59 = R10*R17 + R18;
    double R60 = R10*R9 - R12;
    double R61 = R16*R59 - R20*R60;
    double R62 = R16*R60 + R20*R59;
    double R63 = R2*R49 + R24*R51 - R50*R6 + R52*(-R10*R2 + R24*R56 - R55*R6) + R57*(R2*R58 + R24*R62 - R6*R61) - b0 + plane_dz + plane_z;
    double R64 = R41*R63;
    double R65 = R36*R49 + R37*R50 + R38*R51 + R52*(-R10*R36 + R37*R55 + R38*R56) + R57*(R36*R58 + R37*R61 + R38*R62) + plane_dy;
    double R66 = R44*R65;
    double R67 = R30*R49 + R31*R50 + R34*R51 + R52*(-R10*R30 + R31*R55 + R34*R56) + R57*(R30*R58 + R31*R61 + R34*R62) - a0 + plane_dx;
    double R68 = R45*R67;
    double R69 = -R42*R65 + R43*R63 + R46*R67;
    double R70 = R48*(R47*R69 - R64 - R66 - R68);
    double R71 = R48*(-R47*(R64 + R66 + R68) + R69);
    double R72 = R41*R70 - R43*R71 + R63;
    double R73 = R42*R71 + R44*R70 + R65;
    double R74 = R45*R70 - R46*R71 + R67;
    double R75 = R16*R6;
    double R76 = R20*R24;
    double R77 = -R75 + R76;
    double R78 = R16*R37;
    double R79 = R20*R38;
    double R80 = R78 + R79;
    double R81 = R16*R31;
    double R82 = R20*R34;
    double R83 = R81 + R82;
    double R84 = -R41*R77 - R44*R80 - R45*R83;
    double R85 = -R42*R80 + R43*R77 + R46*R83;
    double R86 = 2*R48;
    double R87 = R86*(R47*R84 + R85);
    double R88 = R86*(R47*R85 + R84);
    double R89 = ((1.0/2.0)*R72*(R41*R88 - R43*R87 - 2*R75 + 2*R76) + (1.0/2.0)*R73*(R42*R87 + R44*R88 + 2*R78 + 2*R79) + (1.0/2.0)*R74*(R45*R88 - R46*R87 + 2*R81 + 2*R82))/sqrt(pow(R72, 2) + pow(R73, 2) + pow(R74, 2));
    double result = 16.0*((R71 > 0) ? (
   R89
)
: (
   -1.0*R89
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dy(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = 1.0*R3*R4;
    double R6 = sin(plane_b);
    double R7 = sin(panel_g);
    double R8 = cos(panel_a);
    double R9 = R7*R8;
    double R10 = sin(panel_b);
    double R11 = cos(panel_g);
    double R12 = R11*R3;
    double R13 = R10*R12 - R9;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = 1.0/R14;
    double R16 = R15*panel_straw0x;
    double R17 = R11*R8;
    double R18 = R3*R7;
    double R19 = R10*R18 + R17;
    double R20 = R15*panel_straw0y;
    double R21 = 1.0*R13*R16 - 1.0*R19*R20;
    double R22 = 1.0*R13*R20 + 1.0*R16*R19;
    double R23 = sin(plane_a);
    double R24 = R1*R23;
    double R25 = R2*R5 - R21*R6 + R22*R24;
    double R26 = sin(plane_g);
    double R27 = R23*R26;
    double R28 = cos(plane_g);
    double R29 = R0*R28;
    double R30 = R27 + R29*R6;
    double R31 = R1*R28;
    double R32 = R0*R26;
    double R33 = R23*R28;
    double R34 = -R32 + R33*R6;
    double R35 = R21*R31 + R22*R34 + R30*R5;
    double R36 = R32*R6 - R33;
    double R37 = R1*R26;
    double R38 = R27*R6 + R29;
    double R39 = R21*R37 + R22*R38 + R36*R5;
    double R40 = pow(pow(R25, 2) + pow(R35, 2) + pow(R39, 2), -1.0/2.0);
    double R41 = R25*R40;
    double R42 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R43 = R42*b1;
    double R44 = R39*R40;
    double R45 = R35*R40;
    double R46 = R42*a1;
    double R47 = R41*R43 - R42*R44 + R45*R46;
    double R48 = 1.0/(1.0 - pow(R47, 2));
    double R49 = panel_dz + panel_straw0z - plane_z;
    double R50 = R16*panel_dx - R20*panel_dy + panel_straw0x;
    double R51 = R16*panel_dy + R20*panel_dx + panel_straw0y;
    double R52 = -R14 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R53 = R16*R4;
    double R54 = R20*R4;
    double R55 = R11*R53 - R54*R7;
    double R56 = R11*R54 + R53*R7;
    double R57 = -panel_straw0z + wire_z;
    double R58 = R4*R8;
    double R59 = R10*R17 + R18;
    double R60 = R10*R9 - R12;
    double R61 = R16*R59 - R20*R60;
    double R62 = R16*R60 + R20*R59;
    double R63 = R2*R49 + R24*R51 - R50*R6 + R52*(-R10*R2 + R24*R56 - R55*R6) + R57*(R2*R58 + R24*R62 - R6*R61) - b0 + plane_dz + plane_z;
    double R64 = R41*R63;
    double R65 = R36*R49 + R37*R50 + R38*R51 + R52*(-R10*R36 + R37*R55 + R38*R56) + R57*(R36*R58 + R37*R61 + R38*R62) + plane_dy;
    double R66 = R44*R65;
    double R67 = R30*R49 + R31*R50 + R34*R51 + R52*(-R10*R30 + R31*R55 + R34*R56) + R57*(R30*R58 + R31*R61 + R34*R62) - a0 + plane_dx;
    double R68 = R45*R67;
    double R69 = -R42*R65 + R43*R63 + R46*R67;
    double R70 = R48*(R47*R69 - R64 - R66 - R68);
    double R71 = R48*(-R47*(R64 + R66 + R68) + R69);
    double R72 = R41*R70 - R43*R71 + R63;
    double R73 = R42*R71 + R44*R70 + R65;
    double R74 = R45*R70 - R46*R71 + R67;
    double R75 = R20*R6;
    double R76 = R16*R24;
    double R77 = R75 + R76;
    double R78 = R20*R37;
    double R79 = R16*R38;
    double R80 = -R78 + R79;
    double R81 = R20*R31;
    double R82 = R16*R34;
    double R83 = -R81 + R82;
    double R84 = -R41*R77 - R44*R80 - R45*R83;
    double R85 = -R42*R80 + R43*R77 + R46*R83;
    double R86 = 2*R48;
    double R87 = R86*(R47*R84 + R85);
    double R88 = R86*(R47*R85 + R84);
    double R89 = ((1.0/2.0)*R72*(R41*R88 - R43*R87 + 2*R75 + 2*R76) + (1.0/2.0)*R73*(R42*R87 + R44*R88 - 2*R78 + 2*R79) + (1.0/2.0)*R74*(R45*R88 - R46*R87 - 2*R81 + 2*R82))/sqrt(pow(R72, 2) + pow(R73, 2) + pow(R74, 2));
    double result = 16.0*((R71 > 0) ? (
   R89
)
: (
   -1.0*R89
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dz(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = 1.0*R3*R4;
    double R6 = sin(plane_b);
    double R7 = sin(panel_g);
    double R8 = cos(panel_a);
    double R9 = R7*R8;
    double R10 = sin(panel_b);
    double R11 = cos(panel_g);
    double R12 = R11*R3;
    double R13 = R10*R12 - R9;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = 1.0/R14;
    double R16 = R15*panel_straw0x;
    double R17 = R11*R8;
    double R18 = R3*R7;
    double R19 = R10*R18 + R17;
    double R20 = R15*panel_straw0y;
    double R21 = 1.0*R13*R16 - 1.0*R19*R20;
    double R22 = 1.0*R13*R20 + 1.0*R16*R19;
    double R23 = sin(plane_a);
    double R24 = R1*R23;
    double R25 = R2*R5 - R21*R6 + R22*R24;
    double R26 = sin(plane_g);
    double R27 = R23*R26;
    double R28 = cos(plane_g);
    double R29 = R0*R28;
    double R30 = R29*R6;
    double R31 = R27 + R30;
    double R32 = R1*R28;
    double R33 = R0*R26;
    double R34 = R23*R28;
    double R35 = -R33 + R34*R6;
    double R36 = R21*R32 + R22*R35 + R31*R5;
    double R37 = R33*R6;
    double R38 = -R34 + R37;
    double R39 = R1*R26;
    double R40 = R27*R6 + R29;
    double R41 = R21*R39 + R22*R40 + R38*R5;
    double R42 = pow(pow(R25, 2) + pow(R36, 2) + pow(R41, 2), -1.0/2.0);
    double R43 = R25*R42;
    double R44 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R45 = R44*b1;
    double R46 = R41*R42;
    double R47 = R36*R42;
    double R48 = R44*a1;
    double R49 = R43*R45 - R44*R46 + R47*R48;
    double R50 = 1.0/(1.0 - pow(R49, 2));
    double R51 = panel_dz + panel_straw0z - plane_z;
    double R52 = R16*panel_dx - R20*panel_dy + panel_straw0x;
    double R53 = R16*panel_dy + R20*panel_dx + panel_straw0y;
    double R54 = -R14 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R55 = R16*R4;
    double R56 = R20*R4;
    double R57 = R11*R55 - R56*R7;
    double R58 = R11*R56 + R55*R7;
    double R59 = -panel_straw0z + wire_z;
    double R60 = R4*R8;
    double R61 = R10*R17 + R18;
    double R62 = R10*R9 - R12;
    double R63 = R16*R61 - R20*R62;
    double R64 = R16*R62 + R20*R61;
    double R65 = R2*R51 + R24*R53 - R52*R6 + R54*(-R10*R2 + R24*R58 - R57*R6) + R59*(R2*R60 + R24*R64 - R6*R63) - b0 + plane_dz + plane_z;
    double R66 = R43*R65;
    double R67 = R38*R51 + R39*R52 + R40*R53 + R54*(-R10*R38 + R39*R57 + R40*R58) + R59*(R38*R60 + R39*R63 + R40*R64) + plane_dy;
    double R68 = R46*R67;
    double R69 = R31*R51 + R32*R52 + R35*R53 + R54*(-R10*R31 + R32*R57 + R35*R58) + R59*(R31*R60 + R32*R63 + R35*R64) - a0 + plane_dx;
    double R70 = R47*R69;
    double R71 = -R44*R67 + R45*R65 + R48*R69;
    double R72 = R50*(R49*R71 - R66 - R68 - R70);
    double R73 = R50*(-R49*(R66 + R68 + R70) + R71);
    double R74 = R43*R72 - R45*R73 + R65;
    double R75 = R44*R73 + R46*R72 + R67;
    double R76 = R47*R72 - R48*R73 + R69;
    double R77 = -R2*R43 - R31*R47 - R38*R46;
    double R78 = R2*R45 + R31*R48 - R38*R44;
    double R79 = 2*R50;
    double R80 = R79*(R49*R77 + R78);
    double R81 = R79*(R49*R78 + R77);
    double R82 = ((1.0/2.0)*R74*(2*R2 + R43*R81 - R45*R80) + (1.0/2.0)*R75*(-2*R34 + 2*R37 + R44*R80 + R46*R81) + (1.0/2.0)*R76*(2*R27 + 2*R30 + R47*R81 - R48*R80))/sqrt(pow(R74, 2) + pow(R75, 2) + pow(R76, 2));
    double result = 16.0*((R73 > 0) ? (
   R82
)
: (
   -1.0*R82
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_a(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = R3*R4;
    double R6 = R2*R5;
    double R7 = sin(plane_b);
    double R8 = sin(panel_g);
    double R9 = cos(panel_a);
    double R10 = R8*R9;
    double R11 = sin(panel_b);
    double R12 = cos(panel_g);
    double R13 = R12*R3;
    double R14 = R11*R13;
    double R15 = -R10 + R14;
    double R16 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R17 = 1.0/R16;
    double R18 = R17*panel_straw0x;
    double R19 = R12*R9;
    double R20 = R3*R8;
    double R21 = R11*R20;
    double R22 = R19 + R21;
    double R23 = R17*panel_straw0y;
    double R24 = 1.0*R15*R18 - 1.0*R22*R23;
    double R25 = 1.0*R15*R23 + 1.0*R18*R22;
    double R26 = sin(plane_a);
    double R27 = R1*R26;
    double R28 = -R24*R7 + R25*R27 + 1.0*R6;
    double R29 = sin(plane_g);
    double R30 = R26*R29;
    double R31 = cos(plane_g);
    double R32 = R0*R31;
    double R33 = R32*R7;
    double R34 = R30 + R33;
    double R35 = R34*R5;
    double R36 = R1*R31;
    double R37 = R0*R29;
    double R38 = R26*R31;
    double R39 = R38*R7;
    double R40 = -R37 + R39;
    double R41 = R24*R36 + R25*R40 + 1.0*R35;
    double R42 = R37*R7;
    double R43 = -R38 + R42;
    double R44 = R43*R5;
    double R45 = R1*R29;
    double R46 = R30*R7;
    double R47 = R32 + R46;
    double R48 = R24*R45 + R25*R47 + 1.0*R44;
    double R49 = pow(R28, 2) + pow(R41, 2) + pow(R48, 2);
    double R50 = pow(R49, -1.0/2.0);
    double R51 = R28*R50;
    double R52 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R53 = R52*b1;
    double R54 = R48*R50;
    double R55 = R52*a1;
    double R56 = R41*R50;
    double R57 = R51*R53 - R52*R54 + R55*R56;
    double R58 = 1.0 - pow(R57, 2);
    double R59 = 1.0/R58;
    double R60 = panel_dz + panel_straw0z - plane_z;
    double R61 = R18*panel_dx - R23*panel_dy + panel_straw0x;
    double R62 = R18*panel_dy + R23*panel_dx + panel_straw0y;
    double R63 = -R16 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R64 = R18*R4;
    double R65 = R23*R4;
    double R66 = R12*R64 - R65*R8;
    double R67 = R12*R65 + R64*R8;
    double R68 = -panel_straw0z + wire_z;
    double R69 = R4*R9;
    double R70 = R2*R69;
    double R71 = R11*R19 + R20;
    double R72 = R18*R71;
    double R73 = R10*R11 - R13;
    double R74 = R23*R73;
    double R75 = R72 - R74;
    double R76 = R23*R71;
    double R77 = R18*R73;
    double R78 = R76 + R77;
    double R79 = R2*R60 + R27*R62 - R61*R7 + R63*(-R11*R2 + R27*R67 - R66*R7) + R68*(R27*R78 - R7*R75 + R70) - b0 + plane_dz + plane_z;
    double R80 = R43*R60 + R45*R61 + R47*R62 + R63*(-R11*R43 + R45*R66 + R47*R67) + R68*(R43*R69 + R45*R75 + R47*R78) + plane_dy;
    double R81 = R34*R60 + R36*R61 + R40*R62 + R63*(-R11*R34 + R36*R66 + R40*R67) + R68*(R34*R69 + R36*R75 + R40*R78) - a0 + plane_dx;
    double R82 = -R52*R80 + R53*R79 + R55*R81;
    double R83 = R50*R79;
    double R84 = R28*R83;
    double R85 = R54*R80;
    double R86 = R50*R81;
    double R87 = R41*R86;
    double R88 = -R84 - R85 - R87;
    double R89 = R57*R82 + R88;
    double R90 = R59*R89;
    double R91 = -R57*(R84 + R85 + R87) + R82;
    double R92 = R59*R91;
    double R93 = R51*R90 - R53*R92 + R79;
    double R94 = R52*R92 + R54*R90 + R80;
    double R95 = -R55*R92 + R56*R90 + R81;
    double R96 = R10 - R14;
    double R97 = R18*R96;
    double R98 = -R19 - R21;
    double R99 = R23*R98;
    double R100 = R18*R98 + R23*R96;
    double R101 = R100*R27 - R6 + R7*(-R97 + R99);
    double R102 = 2*R68;
    double R103 = 1.0*R72;
    double R104 = 1.0*R74;
    double R105 = R7*(-R103 + R104);
    double R106 = R27*(1.0*R76 + 1.0*R77);
    double R107 = R105 + R106 + 1.0*R70;
    double R108 = R107*R50;
    double R109 = 2*R90;
    double R110 = R69*(1.0*R30 + 1.0*R33);
    double R111 = R103 - R104;
    double R112 = R111*R36;
    double R113 = R78*(-1.0*R37 + 1.0*R39);
    double R114 = R69*(-1.0*R38 + 1.0*R42);
    double R115 = R111*R45;
    double R116 = R78*(1.0*R32 + 1.0*R46);
    double R117 = (-1.0/2.0*R28*(2*R105 + 2*R106 + 2.0*R70) - 1.0/2.0*R41*(2*R110 + 2*R112 + 2*R113) - 1.0/2.0*R48*(2*R114 + 2*R115 + 2*R116))/pow(R49, 3.0/2.0);
    double R118 = R117*R28;
    double R119 = R108*R53;
    double R120 = R50*(R114 + R115 + R116);
    double R121 = R120*R52;
    double R122 = R110 + R112 + R113;
    double R123 = R122*R50;
    double R124 = R123*R55;
    double R125 = R118*R53;
    double R126 = R117*R48;
    double R127 = R126*R52;
    double R128 = R117*R41;
    double R129 = R128*R55;
    double R130 = 2*R57*(2*R119 - 2*R121 + 2*R124 + 2*R125 - 2*R127 + 2*R129)/pow(R58, 2);
    double R131 = R130*R91;
    double R132 = R130*R89;
    double R133 = R119 - R121 + R124 + R125 - R127 + R129;
    double R134 = R101*R68;
    double R135 = R97 - R99;
    double R136 = R68*(R100*R47 + R135*R45 - R44);
    double R137 = R100*R40 + R135*R36 - R35;
    double R138 = R137*R68;
    double R139 = R134*R53 - R136*R52 + R138*R55;
    double R140 = -R107*R83 - R118*R79 - R120*R80 - R122*R86 - R126*R80 - R128*R81 - R134*R51 - R136*R54 - R138*R56;
    double R141 = 2*R59;
    double R142 = R141*(R133*R82 + R139*R57 + R140);
    double R143 = R141*(R133*R88 + R139 + R140*R57);
    double R144 = ((1.0/2.0)*R93*(R101*R102 + R108*R109 + R109*R118 - R131*R53 + R132*R51 + R142*R51 - R143*R53) + (1.0/2.0)*R94*(R109*R120 + R109*R126 + R131*R52 + R132*R54 + 2*R136 + R142*R54 + R143*R52) + (1.0/2.0)*R95*(R102*R137 + R109*R123 + R109*R128 - R131*R55 + R132*R56 + R142*R56 - R143*R55))/sqrt(pow(R93, 2) + pow(R94, 2) + pow(R95, 2));
    double result = 16.0*((R92 > 0) ? (
   R144
)
: (
   -1.0*R144
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_b(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(panel_b);
    double R1 = cos(plane_a);
    double R2 = cos(plane_b);
    double R3 = R1*R2;
    double R4 = R0*R3;
    double R5 = sin(panel_a);
    double R6 = 1.0*R5;
    double R7 = sin(plane_b);
    double R8 = sin(panel_g);
    double R9 = cos(panel_a);
    double R10 = R8*R9;
    double R11 = sin(panel_b);
    double R12 = cos(panel_g);
    double R13 = R12*R5;
    double R14 = -R10 + R11*R13;
    double R15 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R16 = 1.0/R15;
    double R17 = R16*panel_straw0x;
    double R18 = R12*R9;
    double R19 = R5*R8;
    double R20 = R11*R19 + R18;
    double R21 = R16*panel_straw0y;
    double R22 = 1.0*R14*R17 - 1.0*R20*R21;
    double R23 = 1.0*R14*R21 + 1.0*R17*R20;
    double R24 = sin(plane_a);
    double R25 = R2*R24;
    double R26 = -R22*R7 + R23*R25 + R4*R6;
    double R27 = sin(plane_g);
    double R28 = R24*R27;
    double R29 = cos(plane_g);
    double R30 = R1*R29;
    double R31 = R30*R7;
    double R32 = R28 + R31;
    double R33 = R0*R6;
    double R34 = R2*R29;
    double R35 = R1*R27;
    double R36 = R24*R29;
    double R37 = R36*R7;
    double R38 = -R35 + R37;
    double R39 = R22*R34 + R23*R38 + R32*R33;
    double R40 = R35*R7;
    double R41 = -R36 + R40;
    double R42 = R2*R27;
    double R43 = R28*R7;
    double R44 = R30 + R43;
    double R45 = R22*R42 + R23*R44 + R33*R41;
    double R46 = pow(R26, 2) + pow(R39, 2) + pow(R45, 2);
    double R47 = pow(R46, -1.0/2.0);
    double R48 = R26*R47;
    double R49 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R50 = R49*b1;
    double R51 = R45*R47;
    double R52 = R39*R47;
    double R53 = R49*a1;
    double R54 = R48*R50 - R49*R51 + R52*R53;
    double R55 = 1.0 - pow(R54, 2);
    double R56 = 1.0/R55;
    double R57 = panel_dz + panel_straw0z - plane_z;
    double R58 = R17*panel_dx - R21*panel_dy + panel_straw0x;
    double R59 = R17*panel_dy + R21*panel_dx + panel_straw0y;
    double R60 = -R15 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R61 = R11*R3;
    double R62 = R0*R17;
    double R63 = R0*R21;
    double R64 = R12*R62 - R63*R8;
    double R65 = R12*R63 + R62*R8;
    double R66 = -panel_straw0z + wire_z;
    double R67 = R0*R9;
    double R68 = R11*R18 + R19;
    double R69 = R10*R11 - R13;
    double R70 = R17*R68 - R21*R69;
    double R71 = R17*R69 + R21*R68;
    double R72 = R25*R59 + R3*R57 - R58*R7 + R60*(R25*R65 - R61 - R64*R7) + R66*(R25*R71 + R3*R67 - R7*R70) - b0 + plane_dz + plane_z;
    double R73 = R11*R41;
    double R74 = R41*R57 + R42*R58 + R44*R59 + R60*(R42*R64 + R44*R65 - R73) + R66*(R41*R67 + R42*R70 + R44*R71) + plane_dy;
    double R75 = R11*R32;
    double R76 = R32*R57 + R34*R58 + R38*R59 + R60*(R34*R64 + R38*R65 - R75) + R66*(R32*R67 + R34*R70 + R38*R71) - a0 + plane_dx;
    double R77 = -R49*R74 + R50*R72 + R53*R76;
    double R78 = R48*R72;
    double R79 = R51*R74;
    double R80 = R52*R76;
    double R81 = -R78 - R79 - R80;
    double R82 = R54*R77 + R81;
    double R83 = R56*R82;
    double R84 = -R54*(R78 + R79 + R80) + R77;
    double R85 = R56*R84;
    double R86 = R48*R83 - R50*R85 + R72;
    double R87 = R49*R85 + R51*R83 + R74;
    double R88 = R52*R83 - R53*R85 + R76;
    double R89 = R10*R63;
    double R90 = R18*R62;
    double R91 = R10*R62 + R18*R63;
    double R92 = R66*(R25*R91 - R61*R9 + R7*(R89 - R90));
    double R93 = R11*R12;
    double R94 = R17*R93;
    double R95 = R11*R8;
    double R96 = R21*R95;
    double R97 = -R17*R95 - R21*R93;
    double R98 = R60*(R25*R97 - R4 + R7*(R94 - R96));
    double R99 = 1.0*R13*R62;
    double R100 = 1.0*R19*R63;
    double R101 = R7*(R100 - R99);
    double R102 = R19*R62;
    double R103 = R13*R63;
    double R104 = R25*(1.0*R102 + 1.0*R103);
    double R105 = R47*(R101 + R104 - R6*R61);
    double R106 = 2*R83;
    double R107 = R11*R5;
    double R108 = R107*(1.0*R28 + 1.0*R31);
    double R109 = -R100 + R99;
    double R110 = R109*R34;
    double R111 = R102 + R103;
    double R112 = R111*(-1.0*R35 + 1.0*R37);
    double R113 = R107*(-1.0*R36 + 1.0*R40);
    double R114 = R109*R42;
    double R115 = R111*(1.0*R30 + 1.0*R43);
    double R116 = (-1.0/2.0*R26*(2*R101 + 2*R104 - 2.0*R5*R61) - 1.0/2.0*R39*(-2*R108 + 2*R110 + 2*R112) - 1.0/2.0*R45*(-2*R113 + 2*R114 + 2*R115))/pow(R46, 3.0/2.0);
    double R117 = R116*R26;
    double R118 = R105*R50;
    double R119 = R47*(-R113 + R114 + R115);
    double R120 = R119*R49;
    double R121 = R47*(-R108 + R110 + R112);
    double R122 = R121*R53;
    double R123 = R117*R50;
    double R124 = R116*R45;
    double R125 = R124*R49;
    double R126 = R116*R39;
    double R127 = R126*R53;
    double R128 = 2*R54*(2*R118 - 2*R120 + 2*R122 + 2*R123 - 2*R125 + 2*R127)/pow(R55, 2);
    double R129 = R128*R84;
    double R130 = R128*R82;
    double R131 = R92 + R98;
    double R132 = -R89 + R90;
    double R133 = R66*(R132*R42 + R44*R91 - R73*R9);
    double R134 = -R94 + R96;
    double R135 = R60*(R0*(R36 - R40) + R134*R42 + R44*R97);
    double R136 = R133 + R135;
    double R137 = R66*(R132*R34 + R38*R91 - R75*R9);
    double R138 = R60*(R0*(-R28 - R31) + R134*R34 + R38*R97);
    double R139 = R137 + R138;
    double R140 = R131*R50 - R136*R49 + R139*R53;
    double R141 = R118 - R120 + R122 + R123 - R125 + R127;
    double R142 = -R105*R72 - R117*R72 - R119*R74 - R121*R76 - R124*R74 - R126*R76 - R131*R48 - R136*R51 - R139*R52;
    double R143 = 2*R56;
    double R144 = R143*(R140*R54 + R141*R77 + R142);
    double R145 = R143*(R140 + R141*R81 + R142*R54);
    double R146 = ((1.0/2.0)*R86*(R105*R106 + R106*R117 - R129*R50 + R130*R48 + R144*R48 - R145*R50 + 2*R92 + 2*R98) + (1.0/2.0)*R87*(R106*R119 + R106*R124 + R129*R49 + R130*R51 + 2*R133 + 2*R135 + R144*R51 + R145*R49) + (1.0/2.0)*R88*(R106*R121 + R106*R126 - R129*R53 + R130*R52 + 2*R137 + 2*R138 + R144*R52 - R145*R53))/sqrt(pow(R86, 2) + pow(R87, 2) + pow(R88, 2));
    double result = 16.0*((R85 > 0) ? (
   R146
)
: (
   -1.0*R146
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_g(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = cos(plane_a);
    double R1 = cos(plane_b);
    double R2 = R0*R1;
    double R3 = sin(panel_a);
    double R4 = cos(panel_b);
    double R5 = 1.0*R3*R4;
    double R6 = sin(plane_b);
    double R7 = sin(panel_g);
    double R8 = cos(panel_a);
    double R9 = R7*R8;
    double R10 = sin(panel_b);
    double R11 = cos(panel_g);
    double R12 = R11*R3;
    double R13 = R10*R12 - R9;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = 1.0/R14;
    double R16 = R15*panel_straw0x;
    double R17 = R13*R16;
    double R18 = R11*R8;
    double R19 = R3*R7;
    double R20 = R10*R19;
    double R21 = R18 + R20;
    double R22 = R15*panel_straw0y;
    double R23 = 1.0*R17 - 1.0*R21*R22;
    double R24 = R13*R22;
    double R25 = 1.0*R16*R21 + 1.0*R24;
    double R26 = sin(plane_a);
    double R27 = R1*R26;
    double R28 = R2*R5 - R23*R6 + R25*R27;
    double R29 = sin(plane_g);
    double R30 = R26*R29;
    double R31 = cos(plane_g);
    double R32 = R0*R31;
    double R33 = R30 + R32*R6;
    double R34 = R1*R31;
    double R35 = R0*R29;
    double R36 = R26*R31;
    double R37 = R36*R6;
    double R38 = -R35 + R37;
    double R39 = R23*R34 + R25*R38 + R33*R5;
    double R40 = R35*R6 - R36;
    double R41 = R1*R29;
    double R42 = R30*R6;
    double R43 = R32 + R42;
    double R44 = R23*R41 + R25*R43 + R40*R5;
    double R45 = pow(R28, 2) + pow(R39, 2) + pow(R44, 2);
    double R46 = pow(R45, -1.0/2.0);
    double R47 = R28*R46;
    double R48 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R49 = R48*b1;
    double R50 = R44*R46;
    double R51 = R39*R46;
    double R52 = R48*a1;
    double R53 = R47*R49 - R48*R50 + R51*R52;
    double R54 = 1.0 - pow(R53, 2);
    double R55 = 1.0/R54;
    double R56 = panel_dz + panel_straw0z - plane_z;
    double R57 = R16*panel_dx - R22*panel_dy + panel_straw0x;
    double R58 = R16*panel_dy + R22*panel_dx + panel_straw0y;
    double R59 = -R14 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R60 = R16*R4;
    double R61 = R22*R4;
    double R62 = R11*R60 - R61*R7;
    double R63 = R60*R7;
    double R64 = R11*R61;
    double R65 = R63 + R64;
    double R66 = -panel_straw0z + wire_z;
    double R67 = R4*R8;
    double R68 = R10*R18 + R19;
    double R69 = R16*R68;
    double R70 = R10*R9;
    double R71 = -R12 + R70;
    double R72 = -R22*R71 + R69;
    double R73 = R22*R68;
    double R74 = R16*R71 + R73;
    double R75 = R2*R56 + R27*R58 - R57*R6 + R59*(-R10*R2 + R27*R65 - R6*R62) + R66*(R2*R67 + R27*R74 - R6*R72) - b0 + plane_dz + plane_z;
    double R76 = R40*R56 + R41*R57 + R43*R58 + R59*(-R10*R40 + R41*R62 + R43*R65) + R66*(R40*R67 + R41*R72 + R43*R74) + plane_dy;
    double R77 = R33*R56 + R34*R57 + R38*R58 + R59*(-R10*R33 + R34*R62 + R38*R65) + R66*(R33*R67 + R34*R72 + R38*R74) - a0 + plane_dx;
    double R78 = -R48*R76 + R49*R75 + R52*R77;
    double R79 = R47*R75;
    double R80 = R50*R76;
    double R81 = R51*R77;
    double R82 = -R79 - R80 - R81;
    double R83 = R53*R78 + R82;
    double R84 = R55*R83;
    double R85 = -R53*(R79 + R80 + R81) + R78;
    double R86 = R55*R85;
    double R87 = R47*R84 - R49*R86 + R75;
    double R88 = R48*R86 + R50*R84 + R76;
    double R89 = R51*R84 - R52*R86 + R77;
    double R90 = R59*(R27*R62 + R6*R65);
    double R91 = R12 - R70;
    double R92 = R16*R91;
    double R93 = R22*R91 + R69;
    double R94 = R66*(R27*R93 + R6*(R73 - R92));
    double R95 = 1.0*R24;
    double R96 = -R18 - R20;
    double R97 = 1.0*R16*R96;
    double R98 = R6*(R95 - R97);
    double R99 = R22*R96;
    double R100 = R27*(1.0*R17 + 1.0*R99);
    double R101 = R46*(R100 + R98);
    double R102 = 2*R84;
    double R103 = -R95 + R97;
    double R104 = R103*R41;
    double R105 = R17 + R99;
    double R106 = R105*(1.0*R32 + 1.0*R42);
    double R107 = R103*R34;
    double R108 = R105*(-1.0*R35 + 1.0*R37);
    double R109 = (-1.0/2.0*R28*(2*R100 + 2*R98) - 1.0/2.0*R39*(2*R107 + 2*R108) - 1.0/2.0*R44*(2*R104 + 2*R106))/pow(R45, 3.0/2.0);
    double R110 = R109*R28;
    double R111 = R101*R49;
    double R112 = R46*(R104 + R106);
    double R113 = R112*R48;
    double R114 = R46*(R107 + R108);
    double R115 = R114*R52;
    double R116 = R110*R49;
    double R117 = R109*R44;
    double R118 = R117*R48;
    double R119 = R109*R39;
    double R120 = R119*R52;
    double R121 = 2*R53*(2*R111 - 2*R113 + 2*R115 + 2*R116 - 2*R118 + 2*R120)/pow(R54, 2);
    double R122 = R121*R85;
    double R123 = R121*R83;
    double R124 = R90 + R94;
    double R125 = -R63 - R64;
    double R126 = R59*(R125*R41 + R43*R62);
    double R127 = -R73 + R92;
    double R128 = R66*(R127*R41 + R43*R93);
    double R129 = R126 + R128;
    double R130 = R59*(R125*R34 + R38*R62);
    double R131 = R66*(R127*R34 + R38*R93);
    double R132 = R130 + R131;
    double R133 = R124*R49 - R129*R48 + R132*R52;
    double R134 = R111 - R113 + R115 + R116 - R118 + R120;
    double R135 = -R101*R75 - R110*R75 - R112*R76 - R114*R77 - R117*R76 - R119*R77 - R124*R47 - R129*R50 - R132*R51;
    double R136 = 2*R55;
    double R137 = R136*(R133*R53 + R134*R78 + R135);
    double R138 = R136*(R133 + R134*R82 + R135*R53);
    double R139 = ((1.0/2.0)*R87*(R101*R102 + R102*R110 - R122*R49 + R123*R47 + R137*R47 - R138*R49 + 2*R90 + 2*R94) + (1.0/2.0)*R88*(R102*R112 + R102*R117 + R122*R48 + R123*R50 + 2*R126 + 2*R128 + R137*R50 + R138*R48) + (1.0/2.0)*R89*(R102*R114 + R102*R119 - R122*R52 + R123*R51 + 2*R130 + 2*R131 + R137*R51 - R138*R52))/sqrt(pow(R87, 2) + pow(R88, 2) + pow(R89, 2));
    double result = 16.0*((R86 > 0) ? (
   R139
)
: (
   -1.0*R139
));
    return result;
}


std::vector<double> CosmicTrack_DCA_LocalDeriv(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
        std::vector<double> result = {CosmicTrack_DCA_Deriv_a0(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_b0(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_a1(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_b1(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_t0(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z)};
return result;
}

std::vector<double> CosmicTrack_DCA_GlobalDeriv(double a0, double b0, double a1, double b1, double t0, double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
        std::vector<double> result = {CosmicTrack_DCA_Deriv_plane_dx(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_plane_dy(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_plane_dz(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_plane_a(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_plane_b(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_plane_g(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_panel_dx(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_panel_dy(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_panel_dz(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_panel_a(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_panel_b(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_panel_g(a0,b0,a1,b1,t0,plane_dx,plane_dy,plane_dz,plane_a,plane_b,plane_g,panel_dx,panel_dy,panel_dz,panel_a,panel_b,panel_g,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z)};
return result;
}

