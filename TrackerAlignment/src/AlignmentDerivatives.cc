

# include "TrackerAlignment/inc/AlignmentDerivatives.hh"
# include <math.h>
# include <vector>

double CosmicTrack_DCA(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = pow(pow(wdir_x, 2) + pow(wdir_y, 2) + pow(wdir_z, 2), -1.0/2.0);
    double R2 = R1*wdir_y;
    double R3 = R1*wdir_x;
    double R4 = R0*a1;
    double R5 = R1*wdir_z;
    double R6 = R0*b1;
    double R7 = -R0*R2 + R3*R4 + R5*R6;
    double R8 = 1.0/(1 - pow(R7, 2));
    double R9 = R0*wire_y;
    double R10 = a0 - wire_x;
    double R11 = R10*R4;
    double R12 = b0 - wire_z;
    double R13 = R12*R6;
    double R14 = R10*R3 + R12*R5 - R2*wire_y;
    double R15 = R8*(-R11 - R13 + R14*R7 - R9);
    double R16 = R8*(R14 - R7*(R11 + R13 + R9));
    double R17 = sqrt(pow(R10 + R15*R4 - R16*R3, 2) + pow(R12 + R15*R6 - R16*R5, 2) + pow(-R0*R15 - R16*R2 - wire_y, 2));
    double result = ((R16 > 0) ? (
   R17
)
: (
   -R17
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


double CosmicTrack_DCA_Deriv_a0(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = pow(panel_straw0x, 2);
    double R3 = pow(panel_straw0y, 2);
    double R4 = R2 + R3;
    double R5 = sqrt(R4);
    double R6 = 1.0/R5;
    double R7 = 1.0/R4;
    double R8 = R6/sqrt(R2*R7 + R3*R7);
    double R9 = 1.0*R8;
    double R10 = R9*panel_straw0x;
    double R11 = R0*a1;
    double R12 = R9*panel_straw0y;
    double R13 = -R0*R10 - R11*R12;
    double R14 = 1.0/(1 - pow(R13, 2));
    double R15 = b0 - wire_z;
    double R16 = R1*R15;
    double R17 = R6*(-R5 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R18 = -R17*panel_straw0y - panel_straw0y;
    double R19 = R0*R18;
    double R20 = -R17*panel_straw0x + a0 - panel_straw0x;
    double R21 = R11*R20;
    double R22 = R10*R18 - R12*R20;
    double R23 = R14*(R13*R22 - R16 + R19 - R21);
    double R24 = R1*R23 + R15;
    double R25 = R14*(-R13*(R16 - R19 + R21) + R22);
    double R26 = -R0*R23 - R10*R25 + R18;
    double R27 = R11*R23 + R12*R25 + R20;
    double R28 = R14*(-R11 - R12*R13);
    double R29 = 2*R28;
    double R30 = 2.0*R14*R8*(-R11*R13 - R12);
    double R31 = (R1*R24*R28 + (1.0/2.0)*R26*(-R0*R29 - R30*panel_straw0x) + (1.0/2.0)*R27*(R11*R29 + R30*panel_straw0y + 2))/sqrt(pow(R24, 2) + pow(R26, 2) + pow(R27, 2));
    double result = 16.0*((R25 > 0) ? (
   R31
)
: (
   -R31
));
    return result;
}


double CosmicTrack_DCA_Deriv_b0(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = pow(panel_straw0x, 2);
    double R5 = pow(panel_straw0y, 2);
    double R6 = R4 + R5;
    double R7 = sqrt(R6);
    double R8 = 1.0/R7;
    double R9 = 1.0/R6;
    double R10 = R8/sqrt(R4*R9 + R5*R9);
    double R11 = 1.0*R10;
    double R12 = R11*panel_straw0x;
    double R13 = R2*a1;
    double R14 = R11*panel_straw0y;
    double R15 = -R12*R2 - R13*R14;
    double R16 = 1.0/(1 - pow(R15, 2));
    double R17 = b0 - wire_z;
    double R18 = R17*R3;
    double R19 = R8*(-R7 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R20 = -R19*panel_straw0y - panel_straw0y;
    double R21 = R2*R20;
    double R22 = -R19*panel_straw0x + a0 - panel_straw0x;
    double R23 = R13*R22;
    double R24 = R12*R20 - R14*R22;
    double R25 = R16*(R15*R24 - R18 + R21 - R23);
    double R26 = R17 + R25*R3;
    double R27 = R16*(-R15*(R18 - R21 + R23) + R24);
    double R28 = -R12*R27 - R2*R25 + R20;
    double R29 = R13*R25 + R14*R27 + R22;
    double R30 = 2*R16/R1;
    double R31 = R30*b1;
    double R32 = 2.0*R10*R15*R16*R3;
    double R33 = ((1.0/2.0)*R26*(-R0*R30 + 2) + (1.0/2.0)*R28*(R31 + R32*panel_straw0x) + (1.0/2.0)*R29*(-R31*a1 - R32*panel_straw0y))/sqrt(pow(R26, 2) + pow(R28, 2) + pow(R29, 2));
    double result = 16.0*((R27 > 0) ? (
   R33
)
: (
   -R33
));
    return result;
}


double CosmicTrack_DCA_Deriv_a1(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(a1, 2);
    double R1 = R0 + pow(b1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = pow(panel_straw0x, 2);
    double R5 = pow(panel_straw0y, 2);
    double R6 = R4 + R5;
    double R7 = sqrt(R6);
    double R8 = 1.0/R7;
    double R9 = 1.0/R6;
    double R10 = R8/sqrt(R4*R9 + R5*R9);
    double R11 = 1.0*R10;
    double R12 = R11*panel_straw0x;
    double R13 = R11*panel_straw0y;
    double R14 = R13*R2;
    double R15 = -R12*R2 - R14*a1;
    double R16 = 1 - pow(R15, 2);
    double R17 = 1.0/R16;
    double R18 = R8*(-R7 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R19 = -R18*panel_straw0y - panel_straw0y;
    double R20 = -R18*panel_straw0x + a0 - panel_straw0x;
    double R21 = R12*R19 - R13*R20;
    double R22 = b0 - wire_z;
    double R23 = R22*R3;
    double R24 = R19*R2;
    double R25 = R2*R20;
    double R26 = R25*a1;
    double R27 = -R23 + R24 - R26;
    double R28 = R15*R21 + R27;
    double R29 = R17*R28;
    double R30 = R22 + R29*R3;
    double R31 = R2*R29;
    double R32 = -R15*(R23 - R24 + R26) + R21;
    double R33 = R17*R32;
    double R34 = -R12*R33 + R19 - R31;
    double R35 = R13*R33 + R20 + R31*a1;
    double R36 = pow(R1, -3.0/2.0);
    double R37 = R36*a1;
    double R38 = R37*b1;
    double R39 = 2*R29;
    double R40 = R0*R36;
    double R41 = R12*R37 + R13*R40 - R14;
    double R42 = -R19*R37 + R20*R40 + R22*R38 - R25;
    double R43 = R17*(R21*R41 + R42);
    double R44 = 2*R3;
    double R45 = 2.0*R10;
    double R46 = R45*panel_straw0y;
    double R47 = R45*panel_straw0x;
    double R48 = R15*(-R2*R46 + R37*R47 + R40*R46)/pow(R16, 2);
    double R49 = R28*R48;
    double R50 = 2*R2;
    double R51 = R43*R50;
    double R52 = R17*(R15*R42 + R27*R41);
    double R53 = R49*R50;
    double R54 = R32*R48;
    double R55 = ((1.0/2.0)*R30*(-R38*R39 + R43*R44 + R44*R49) + (1.0/2.0)*R34*(R37*R39 - R47*R52 - R47*R54 - R51 - R53) + (1.0/2.0)*R35*(2*R31 - R39*R40 + R46*R52 + R46*R54 + R51*a1 + R53*a1))/sqrt(pow(R30, 2) + pow(R34, 2) + pow(R35, 2));
    double result = 16.0*((R33 > 0) ? (
   R55
)
: (
   -R55
));
    return result;
}


double CosmicTrack_DCA_Deriv_b1(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = pow(panel_straw0x, 2);
    double R4 = pow(panel_straw0y, 2);
    double R5 = R3 + R4;
    double R6 = sqrt(R5);
    double R7 = 1.0/R6;
    double R8 = 1.0/R5;
    double R9 = R7/sqrt(R3*R8 + R4*R8);
    double R10 = 1.0*R9;
    double R11 = R10*panel_straw0x;
    double R12 = R2*a1;
    double R13 = R10*panel_straw0y;
    double R14 = -R11*R2 - R12*R13;
    double R15 = 1 - pow(R14, 2);
    double R16 = 1.0/R15;
    double R17 = R7*(-R6 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R18 = -R17*panel_straw0y - panel_straw0y;
    double R19 = -R17*panel_straw0x + a0 - panel_straw0x;
    double R20 = R11*R18 - R13*R19;
    double R21 = b0 - wire_z;
    double R22 = R2*R21;
    double R23 = R22*b1;
    double R24 = R18*R2;
    double R25 = R12*R19;
    double R26 = -R23 + R24 - R25;
    double R27 = R14*R20 + R26;
    double R28 = R16*R27;
    double R29 = R2*R28;
    double R30 = R21 + R29*b1;
    double R31 = -R14*(R23 - R24 + R25) + R20;
    double R32 = R16*R31;
    double R33 = -R11*R32 + R18 - R29;
    double R34 = R12*R28 + R13*R32 + R19;
    double R35 = pow(R1, -3.0/2.0);
    double R36 = R0*R35;
    double R37 = 2*R28;
    double R38 = 2*R2;
    double R39 = R35*b1;
    double R40 = R39*a1;
    double R41 = R11*R39 + R13*R40;
    double R42 = -R18*R39 + R19*R40 + R21*R36 - R22;
    double R43 = R16*(R20*R41 + R42);
    double R44 = R38*R43;
    double R45 = 2.0*R9;
    double R46 = R45*panel_straw0x;
    double R47 = R45*panel_straw0y;
    double R48 = R14*(R39*R46 + R40*R47)/pow(R15, 2);
    double R49 = R27*R48;
    double R50 = R38*R49;
    double R51 = R16*(R14*R42 + R26*R41);
    double R52 = R31*R48;
    double R53 = 2*R12;
    double R54 = ((1.0/2.0)*R30*(2*R29 - R36*R37 + R44*b1 + R50*b1) + (1.0/2.0)*R33*(R37*R39 - R44 - R46*R51 - R46*R52 - R50) + (1.0/2.0)*R34*(-R37*R40 + R43*R53 + R47*R51 + R47*R52 + R49*R53))/sqrt(pow(R30, 2) + pow(R33, 2) + pow(R34, 2));
    double result = 16.0*((R32 > 0) ? (
   R54
)
: (
   -R54
));
    return result;
}


double CosmicTrack_DCA_Deriv_t0(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double result = 1;
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dx(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = pow(panel_straw0x, 2);
    double R3 = pow(panel_straw0y, 2);
    double R4 = R2 + R3;
    double R5 = sqrt(R4);
    double R6 = 1.0/R5;
    double R7 = 1.0/R4;
    double R8 = R6/sqrt(R2*R7 + R3*R7);
    double R9 = 1.0*R8;
    double R10 = R9*panel_straw0x;
    double R11 = R0*a1;
    double R12 = R9*panel_straw0y;
    double R13 = -R0*R10 - R11*R12;
    double R14 = 1.0/(1 - pow(R13, 2));
    double R15 = b0 - wire_z;
    double R16 = R1*R15;
    double R17 = R6*(-R5 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R18 = -R17*panel_straw0y - panel_straw0y;
    double R19 = R0*R18;
    double R20 = -R17*panel_straw0x + a0 - panel_straw0x;
    double R21 = R11*R20;
    double R22 = R10*R18 - R12*R20;
    double R23 = R14*(R13*R22 - R16 + R19 - R21);
    double R24 = R1*R23 + R15;
    double R25 = R14*(-R13*(R16 - R19 + R21) + R22);
    double R26 = -R0*R23 - R10*R25 + R18;
    double R27 = R11*R23 + R12*R25 + R20;
    double R28 = R14*(R11 + R12*R13);
    double R29 = 2*R28;
    double R30 = 2.0*R14*R8*(R11*R13 + R12);
    double R31 = (R1*R24*R28 + (1.0/2.0)*R26*(-R0*R29 - R30*panel_straw0x) + (1.0/2.0)*R27*(R11*R29 + R30*panel_straw0y - 2))/sqrt(pow(R24, 2) + pow(R26, 2) + pow(R27, 2));
    double result = 16.0*((R25 > 0) ? (
   R31
)
: (
   -R31
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dy(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = pow(panel_straw0x, 2);
    double R3 = pow(panel_straw0y, 2);
    double R4 = R2 + R3;
    double R5 = sqrt(R4);
    double R6 = 1.0/R5;
    double R7 = 1.0/R4;
    double R8 = R6/sqrt(R2*R7 + R3*R7);
    double R9 = 1.0*R8;
    double R10 = R9*panel_straw0x;
    double R11 = R0*a1;
    double R12 = R9*panel_straw0y;
    double R13 = -R0*R10 - R11*R12;
    double R14 = 1.0/(1 - pow(R13, 2));
    double R15 = b0 - wire_z;
    double R16 = R1*R15;
    double R17 = R6*(-R5 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R18 = -R17*panel_straw0y - panel_straw0y;
    double R19 = R0*R18;
    double R20 = -R17*panel_straw0x + a0 - panel_straw0x;
    double R21 = R11*R20;
    double R22 = R10*R18 - R12*R20;
    double R23 = R14*(R13*R22 - R16 + R19 - R21);
    double R24 = R1*R23 + R15;
    double R25 = R14*(-R13*(R16 - R19 + R21) + R22);
    double R26 = -R0*R23 - R10*R25 + R18;
    double R27 = R11*R23 + R12*R25 + R20;
    double R28 = R14*(-R0 - R10*R13);
    double R29 = 2*R28;
    double R30 = 2.0*R14*R8*(-R0*R13 - R10);
    double R31 = (R1*R24*R28 + (1.0/2.0)*R26*(-R0*R29 - R30*panel_straw0x - 2) + (1.0/2.0)*R27*(R11*R29 + R30*panel_straw0y))/sqrt(pow(R24, 2) + pow(R26, 2) + pow(R27, 2));
    double result = 16.0*((R25 > 0) ? (
   R31
)
: (
   -R31
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dz(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = pow(panel_straw0x, 2);
    double R5 = pow(panel_straw0y, 2);
    double R6 = R4 + R5;
    double R7 = sqrt(R6);
    double R8 = 1.0/R7;
    double R9 = 1.0/R6;
    double R10 = R8/sqrt(R4*R9 + R5*R9);
    double R11 = 1.0*R10;
    double R12 = R11*panel_straw0x;
    double R13 = R2*a1;
    double R14 = R11*panel_straw0y;
    double R15 = -R12*R2 - R13*R14;
    double R16 = 1.0/(1 - pow(R15, 2));
    double R17 = b0 - wire_z;
    double R18 = R17*R3;
    double R19 = R8*(-R7 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R20 = -R19*panel_straw0y - panel_straw0y;
    double R21 = R2*R20;
    double R22 = -R19*panel_straw0x + a0 - panel_straw0x;
    double R23 = R13*R22;
    double R24 = R12*R20 - R14*R22;
    double R25 = R16*(R15*R24 - R18 + R21 - R23);
    double R26 = R17 + R25*R3;
    double R27 = R16*(-R15*(R18 - R21 + R23) + R24);
    double R28 = -R12*R27 - R2*R25 + R20;
    double R29 = R13*R25 + R14*R27 + R22;
    double R30 = 2*R16/R1;
    double R31 = R30*b1;
    double R32 = 2.0*R10*R15*R16*R3;
    double R33 = ((1.0/2.0)*R26*(R0*R30 - 2) + (1.0/2.0)*R28*(-R31 - R32*panel_straw0x) + (1.0/2.0)*R29*(R31*a1 + R32*panel_straw0y))/sqrt(pow(R26, 2) + pow(R28, 2) + pow(R29, 2));
    double result = 16.0*((R27 > 0) ? (
   R33
)
: (
   -R33
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_a(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = pow(panel_straw0x, 2);
    double R5 = pow(panel_straw0y, 2);
    double R6 = R4 + R5;
    double R7 = sqrt(R6);
    double R8 = 1.0/R7;
    double R9 = 1.0/R6;
    double R10 = 1.0*R9;
    double R11 = R10*R4 + R10*R5;
    double R12 = R8/sqrt(R11);
    double R13 = 1.0*R12;
    double R14 = R13*panel_straw0x;
    double R15 = R2*a1;
    double R16 = R13*panel_straw0y;
    double R17 = -R14*R2 - R15*R16;
    double R18 = 1 - pow(R17, 2);
    double R19 = 1.0/R18;
    double R20 = -panel_straw0y;
    double R21 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R22 = R8*(R21 - R7);
    double R23 = R20 - R22*panel_straw0y;
    double R24 = -R22*panel_straw0x + a0 - panel_straw0x;
    double R25 = R14*R23 - R16*R24;
    double R26 = b0 - wire_z;
    double R27 = R26*R3;
    double R28 = R2*R23;
    double R29 = R15*R24;
    double R30 = -R27 + R28 - R29;
    double R31 = R17*R25 + R30;
    double R32 = R19*R31;
    double R33 = R26 + R3*R32;
    double R34 = -R17*(R27 - R28 + R29) + R25;
    double R35 = R19*R34;
    double R36 = -R14*R35 - R2*R32 + R23;
    double R37 = R15*R32 + R16*R35 + R24;
    double R38 = 2*panel_straw0y;
    double R39 = R8*(-R21 + R7);
    double R40 = 2.0*R12;
    double R41 = R40*panel_straw0x;
    double R42 = -plane_z + wire_z;
    double R43 = R14*R26 + R14*R42;
    double R44 = R14*R3;
    double R45 = R2*R42 - R3*(R20 + R39*panel_straw0y);
    double R46 = 2*R19*(R17*R43 + R25*R44 + R45);
    double R47 = 4.0*R17/pow(R18, 2);
    double R48 = R47*panel_straw0x;
    double R49 = R12*R31*R48/R1;
    double R50 = R19*(R17*R45 + R30*R44 + R43);
    double R51 = R49*b1;
    double R52 = R3*R34*R9/R11;
    double R53 = ((1.0/2.0)*R33*(R0*R49 + R3*R46 - R35*R41 + R38*R39 - R38) + (1.0/2.0)*R36*(-R2*R46 - R4*R47*R52 - R41*R50 - R51 - 2*plane_z + 2*wire_z) + (1.0/2.0)*R37*(R15*R46 + R40*R50*panel_straw0y + R48*R52*panel_straw0y + R51*a1))/sqrt(pow(R33, 2) + pow(R36, 2) + pow(R37, 2));
    double result = 16.0*((R35 > 0) ? (
   R53
)
: (
   -R53
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_b(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = pow(panel_straw0x, 2);
    double R5 = pow(panel_straw0y, 2);
    double R6 = R4 + R5;
    double R7 = sqrt(R6);
    double R8 = 1.0/R7;
    double R9 = 1.0/R6;
    double R10 = 1.0*R9;
    double R11 = R10*R4 + R10*R5;
    double R12 = R8/sqrt(R11);
    double R13 = 1.0*R12;
    double R14 = R13*panel_straw0x;
    double R15 = R2*a1;
    double R16 = R13*panel_straw0y;
    double R17 = -R14*R2 - R15*R16;
    double R18 = 1 - pow(R17, 2);
    double R19 = 1.0/R18;
    double R20 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R21 = R8*(R20 - R7);
    double R22 = -R21*panel_straw0y - panel_straw0y;
    double R23 = -R21*panel_straw0x + a0 - panel_straw0x;
    double R24 = R14*R22 - R16*R23;
    double R25 = -wire_z;
    double R26 = R25 + b0;
    double R27 = R26*R3;
    double R28 = R2*R22;
    double R29 = R15*R23;
    double R30 = -R27 + R28 - R29;
    double R31 = R17*R24 + R30;
    double R32 = R19*R31;
    double R33 = R26 + R3*R32;
    double R34 = -R17*(R27 - R28 + R29) + R24;
    double R35 = R19*R34;
    double R36 = -R14*R35 - R2*R32 + R22;
    double R37 = R15*R32 + R16*R35 + R23;
    double R38 = 2*panel_straw0x;
    double R39 = R8*(-R20 + R7);
    double R40 = 2.0*R12;
    double R41 = R40*panel_straw0y;
    double R42 = R25 + plane_z;
    double R43 = R16*R26 - R16*R42;
    double R44 = R16*R3;
    double R45 = -R15*R42 - R3*(-R39*panel_straw0x + panel_straw0x);
    double R46 = 2*R19*(R17*R43 + R24*R44 + R45);
    double R47 = 4.0*R17/pow(R18, 2);
    double R48 = R47*panel_straw0y;
    double R49 = R12*R31*R48/R1;
    double R50 = R19*(R17*R45 + R30*R44 + R43);
    double R51 = R49*b1;
    double R52 = R3*R34*R9/R11;
    double R53 = ((1.0/2.0)*R33*(R0*R49 + R3*R46 - R35*R41 - R38*R39 + R38) + (1.0/2.0)*R36*(-R2*R46 - R40*R50*panel_straw0x - R48*R52*panel_straw0x - R51) + (1.0/2.0)*R37*(R15*R46 + R41*R50 + R47*R5*R52 + R51*a1 + 2*plane_z - 2*wire_z))/sqrt(pow(R33, 2) + pow(R36, 2) + pow(R37, 2));
    double result = 16.0*((R35 > 0) ? (
   R53
)
: (
   -R53
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_g(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = pow(panel_straw0x, 2);
    double R3 = pow(panel_straw0y, 2);
    double R4 = R2 + R3;
    double R5 = sqrt(R4);
    double R6 = 1.0/R5;
    double R7 = 1.0/R4;
    double R8 = R6/sqrt(R2*R7 + R3*R7);
    double R9 = 1.0*R8;
    double R10 = R9*panel_straw0x;
    double R11 = R0*a1;
    double R12 = R9*panel_straw0y;
    double R13 = -R0*R10 - R11*R12;
    double R14 = 1 - pow(R13, 2);
    double R15 = 1.0/R14;
    double R16 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R17 = R6*(R16 - R5);
    double R18 = -R17*panel_straw0y - panel_straw0y;
    double R19 = -panel_straw0x;
    double R20 = -R17*panel_straw0x + R19 + a0;
    double R21 = R10*R18 - R12*R20;
    double R22 = b0 - wire_z;
    double R23 = R1*R22;
    double R24 = R0*R18;
    double R25 = R11*R20;
    double R26 = -R23 + R24 - R25;
    double R27 = R13*R21 + R26;
    double R28 = R15*R27;
    double R29 = R1*R28 + R22;
    double R30 = -R13*(R23 - R24 + R25) + R21;
    double R31 = R15*R30;
    double R32 = -R0*R28 - R10*R31 + R18;
    double R33 = R11*R28 + R12*R31 + R20;
    double R34 = 2.0*R8;
    double R35 = R34*panel_straw0y;
    double R36 = R34*panel_straw0x;
    double R37 = R13*(R0*R35 - R11*R36)/pow(R14, 2);
    double R38 = R27*R37;
    double R39 = 2*R38;
    double R40 = R0*R12 - R10*R11;
    double R41 = R6*(-R16 + R5);
    double R42 = -R41*panel_straw0y + panel_straw0y;
    double R43 = R19 + R41*panel_straw0x;
    double R44 = -R10*R20 + R10*R43 - R12*R18 - R12*R42;
    double R45 = R0*R43 - R11*R42;
    double R46 = R15*(R13*R44 + R21*R40 + R45);
    double R47 = 2*R46;
    double R48 = 2*panel_straw0x;
    double R49 = 2*R0;
    double R50 = R30*R37;
    double R51 = R15*(R13*R45 + R26*R40 + R44);
    double R52 = 2*panel_straw0y;
    double R53 = ((1.0/2.0)*R29*(R1*R39 + R1*R47) + (1.0/2.0)*R32*(R31*R35 - R36*R50 - R36*R51 - R38*R49 + R41*R48 - R46*R49 - R48) + (1.0/2.0)*R33*(R11*R39 + R11*R47 + R31*R36 + R35*R50 + R35*R51 - R41*R52 + R52))/sqrt(pow(R29, 2) + pow(R32, 2) + pow(R33, 2));
    double result = 16.0*((R31 > 0) ? (
   R53
)
: (
   -R53
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dx(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = pow(panel_straw0x, 2);
    double R3 = pow(panel_straw0y, 2);
    double R4 = R2 + R3;
    double R5 = sqrt(R4);
    double R6 = 1.0/R5;
    double R7 = R0*R6;
    double R8 = 1.0/R4;
    double R9 = pow(R2*R8 + R3*R8, -1.0/2.0);
    double R10 = 1.0*R9;
    double R11 = R0*a1;
    double R12 = R6*panel_straw0y;
    double R13 = R10*R12;
    double R14 = -R10*R7*panel_straw0x - R11*R13;
    double R15 = 1.0/(1 - pow(R14, 2));
    double R16 = b0 - wire_z;
    double R17 = R1*R16;
    double R18 = R6*(-R5 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R19 = -R18*panel_straw0y - panel_straw0y;
    double R20 = R0*R19;
    double R21 = -R18*panel_straw0x + a0 - panel_straw0x;
    double R22 = R11*R21;
    double R23 = R6*panel_straw0x;
    double R24 = R10*R23;
    double R25 = -R13*R21 + R19*R24;
    double R26 = R15*(R14*R25 - R17 + R20 - R22);
    double R27 = R1*R26 + R16;
    double R28 = R15*(-R14*(R17 - R20 + R22) + R25);
    double R29 = -R0*R26 + R19 - R24*R28;
    double R30 = R11*R26 + R13*R28 + R21;
    double R31 = R15*(R11*R23 - R7*panel_straw0y);
    double R32 = 2*R6;
    double R33 = 2*R31;
    double R34 = 2.0*R14*R31*R9;
    double R35 = (R1*R27*R31 + (1.0/2.0)*R29*(-R0*R33 - R23*R34 - R32*panel_straw0y) + (1.0/2.0)*R30*(R11*R33 + R12*R34 - R32*panel_straw0x))/sqrt(pow(R27, 2) + pow(R29, 2) + pow(R30, 2));
    double result = 16.0*((R28 > 0) ? (
   R35
)
: (
   -R35
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dy(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = pow(panel_straw0x, 2);
    double R3 = pow(panel_straw0y, 2);
    double R4 = R2 + R3;
    double R5 = sqrt(R4);
    double R6 = 1.0/R5;
    double R7 = R6*panel_straw0x;
    double R8 = R0*R7;
    double R9 = 1.0/R4;
    double R10 = R2*R9;
    double R11 = R3*R9;
    double R12 = pow(R10 + R11, -1.0/2.0);
    double R13 = 1.0*R12;
    double R14 = R0*a1;
    double R15 = R6*panel_straw0y;
    double R16 = R14*R15;
    double R17 = -R13*R16 - R13*R8;
    double R18 = 1.0/(1 - pow(R17, 2));
    double R19 = b0 - wire_z;
    double R20 = R1*R19;
    double R21 = R6*(-R5 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R22 = -R21*panel_straw0y - panel_straw0y;
    double R23 = R0*R22;
    double R24 = -R21*panel_straw0x + a0 - panel_straw0x;
    double R25 = R14*R24;
    double R26 = R13*R7;
    double R27 = R13*R15;
    double R28 = R22*R26 - R24*R27;
    double R29 = R18*(R17*R28 - R20 + R23 - R25);
    double R30 = R1*R29 + R19;
    double R31 = R18*(-R17*(R20 - R23 + R25) + R28);
    double R32 = -R0*R29 + R22 - R26*R31;
    double R33 = R14*R29 + R24 + R27*R31;
    double R34 = -R10*R12 - R11*R12;
    double R35 = -R16 - R8;
    double R36 = R18*(R17*R34 + R35);
    double R37 = 2*R36;
    double R38 = 2.0*R12*R18*(R17*R35 + R34);
    double R39 = (R1*R30*R36 + (1.0/2.0)*R32*(-R0*R37 - R38*R7 - 2*R7) + (1.0/2.0)*R33*(R14*R37 + R15*R38 + 2*R15))/sqrt(pow(R30, 2) + pow(R32, 2) + pow(R33, 2));
    double result = 16.0*((R31 > 0) ? (
   R39
)
: (
   -R39
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dz(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = pow(panel_straw0x, 2);
    double R5 = pow(panel_straw0y, 2);
    double R6 = R4 + R5;
    double R7 = sqrt(R6);
    double R8 = 1.0/R7;
    double R9 = 1.0/R6;
    double R10 = R8/sqrt(R4*R9 + R5*R9);
    double R11 = 1.0*R10;
    double R12 = R11*panel_straw0x;
    double R13 = R2*a1;
    double R14 = R11*panel_straw0y;
    double R15 = -R12*R2 - R13*R14;
    double R16 = 1.0/(1 - pow(R15, 2));
    double R17 = b0 - wire_z;
    double R18 = R17*R3;
    double R19 = R8*(-R7 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R20 = -R19*panel_straw0y - panel_straw0y;
    double R21 = R2*R20;
    double R22 = -R19*panel_straw0x + a0 - panel_straw0x;
    double R23 = R13*R22;
    double R24 = R12*R20 - R14*R22;
    double R25 = R16*(R15*R24 - R18 + R21 - R23);
    double R26 = R17 + R25*R3;
    double R27 = R16*(-R15*(R18 - R21 + R23) + R24);
    double R28 = -R12*R27 - R2*R25 + R20;
    double R29 = R13*R25 + R14*R27 + R22;
    double R30 = 2*R16/R1;
    double R31 = R30*b1;
    double R32 = 2.0*R10*R15*R16*R3;
    double R33 = ((1.0/2.0)*R26*(R0*R30 - 2) + (1.0/2.0)*R28*(-R31 - R32*panel_straw0x) + (1.0/2.0)*R29*(R31*a1 + R32*panel_straw0y))/sqrt(pow(R26, 2) + pow(R28, 2) + pow(R29, 2));
    double result = 16.0*((R27 > 0) ? (
   R33
)
: (
   -R33
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_a(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = pow(panel_straw0x, 2);
    double R5 = pow(panel_straw0y, 2);
    double R6 = R4 + R5;
    double R7 = sqrt(R6);
    double R8 = 1.0/R7;
    double R9 = 1.0/R6;
    double R10 = R4*R9;
    double R11 = R5*R9;
    double R12 = R10 + R11;
    double R13 = pow(R12, -1.0/2.0);
    double R14 = 1.0*R13;
    double R15 = R14*R8;
    double R16 = R15*panel_straw0x;
    double R17 = R2*a1;
    double R18 = R15*panel_straw0y;
    double R19 = -R16*R2 - R17*R18;
    double R20 = 1 - pow(R19, 2);
    double R21 = 1.0/R20;
    double R22 = R8*(-R7 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R23 = -R22*panel_straw0y - panel_straw0y;
    double R24 = -R22*panel_straw0x + a0 - panel_straw0x;
    double R25 = R16*R23 - R18*R24;
    double R26 = -wire_z;
    double R27 = R26 + b0;
    double R28 = R27*R3;
    double R29 = R2*R23;
    double R30 = R17*R24;
    double R31 = -R28 + R29 - R30;
    double R32 = R19*R25 + R31;
    double R33 = R21*R32;
    double R34 = R27 + R3*R33;
    double R35 = -R19*(R28 - R29 + R30) + R25;
    double R36 = R21*R35;
    double R37 = -R16*R36 - R2*R33 + R23;
    double R38 = R17*R33 + R18*R36 + R24;
    double R39 = 2.0*R13;
    double R40 = R26 + panel_straw0z;
    double R41 = R13*R40;
    double R42 = -R10*R41 - R11*R41 + R14*R27;
    double R43 = R14*R3;
    double R44 = R40*R8;
    double R45 = R44*panel_straw0x;
    double R46 = R44*panel_straw0y;
    double R47 = -R17*R46 - R2*R45;
    double R48 = 2*R21*(R19*R42 + R25*R43 + R47);
    double R49 = 4.0*R19/pow(R20, 2);
    double R50 = R13*R32*R49/R1;
    double R51 = R8*panel_straw0x;
    double R52 = R21*R39*(R19*R47 + R31*R43 + R42);
    double R53 = R50*b1;
    double R54 = R3*R35*R49/R12;
    double R55 = R8*panel_straw0y;
    double R56 = ((1.0/2.0)*R34*(R0*R50 + R3*R48 - R36*R39) + (1.0/2.0)*R37*(-R2*R48 - 2*R45 - R51*R52 - R51*R54 - R53) + (1.0/2.0)*R38*(R17*R48 + 2*R46 + R52*R55 + R53*a1 + R54*R55))/sqrt(pow(R34, 2) + pow(R37, 2) + pow(R38, 2));
    double result = 16.0*((R36 > 0) ? (
   R56
)
: (
   -R56
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_b(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = pow(panel_straw0x, 2);
    double R3 = pow(panel_straw0y, 2);
    double R4 = R2 + R3;
    double R5 = sqrt(R4);
    double R6 = 1.0/R5;
    double R7 = 1.0/R4;
    double R8 = R6/sqrt(R2*R7 + R3*R7);
    double R9 = 1.0*R8;
    double R10 = R9*panel_straw0x;
    double R11 = R0*a1;
    double R12 = R9*panel_straw0y;
    double R13 = -R0*R10 - R11*R12;
    double R14 = 1.0/(1 - pow(R13, 2));
    double R15 = -wire_z;
    double R16 = R15 + b0;
    double R17 = R1*R16;
    double R18 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R19 = R18 - R5;
    double R20 = R19*R6;
    double R21 = -R20*panel_straw0y - panel_straw0y;
    double R22 = R0*R21;
    double R23 = -R20*panel_straw0x + a0 - panel_straw0x;
    double R24 = R11*R23;
    double R25 = R10*R21 - R12*R23;
    double R26 = R14*(R13*R25 - R17 + R22 - R24);
    double R27 = R1*R26 + R16;
    double R28 = R14*(-R13*(R17 - R22 + R24) + R25);
    double R29 = -R0*R26 - R10*R28 + R21;
    double R30 = R11*R26 + R12*R28 + R23;
    double R31 = R6*(R15 + panel_straw0z);
    double R32 = R31*panel_straw0y;
    double R33 = R31*panel_straw0x;
    double R34 = R14*(R0*R32 - R1*R19 - R11*R33);
    double R35 = 2*R34;
    double R36 = 2.0*R13*R34*R8;
    double R37 = ((1.0/2.0)*R27*(R1*R35 + 2*R18 - 2*R5) + (1.0/2.0)*R29*(-R0*R35 + 2*R32 - R36*panel_straw0x) + (1.0/2.0)*R30*(R11*R35 + 2*R33 + R36*panel_straw0y))/sqrt(pow(R27, 2) + pow(R29, 2) + pow(R30, 2));
    double result = 16.0*((R28 > 0) ? (
   R37
)
: (
   -R37
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_g(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = pow(panel_straw0x, 2);
    double R3 = pow(panel_straw0y, 2);
    double R4 = R2 + R3;
    double R5 = sqrt(R4);
    double R6 = 1.0/R5;
    double R7 = 1.0/R4;
    double R8 = R2*R7;
    double R9 = R3*R7;
    double R10 = pow(R8 + R9, -1.0/2.0);
    double R11 = R10*R6;
    double R12 = 1.0*R11;
    double R13 = R12*panel_straw0x;
    double R14 = R0*a1;
    double R15 = R12*panel_straw0y;
    double R16 = -R0*R13 - R14*R15;
    double R17 = 1 - pow(R16, 2);
    double R18 = 1.0/R17;
    double R19 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R20 = R6*(R19 - R5);
    double R21 = -R20*panel_straw0y - panel_straw0y;
    double R22 = -R20*panel_straw0x + a0 - panel_straw0x;
    double R23 = R13*R21 - R15*R22;
    double R24 = b0 - wire_z;
    double R25 = R1*R24;
    double R26 = R0*R21;
    double R27 = R14*R22;
    double R28 = -R25 + R26 - R27;
    double R29 = R16*R23 + R28;
    double R30 = R18*R29;
    double R31 = R1*R30 + R24;
    double R32 = -R16*(R25 - R26 + R27) + R23;
    double R33 = R18*R32;
    double R34 = -R0*R30 - R13*R33 + R21;
    double R35 = R14*R30 + R15*R33 + R22;
    double R36 = 2.0*R11;
    double R37 = R36*panel_straw0y;
    double R38 = R36*panel_straw0x;
    double R39 = R16*(R0*R37 - R14*R38)/pow(R17, 2);
    double R40 = R29*R39;
    double R41 = 2*R40;
    double R42 = R0*R15 - R13*R14;
    double R43 = -R19 + R5;
    double R44 = R10*R43;
    double R45 = -R13*R22 - R15*R21 + R44*R8 + R44*R9;
    double R46 = R43*R6;
    double R47 = R46*panel_straw0x;
    double R48 = R46*panel_straw0y;
    double R49 = R0*R47 + R14*R48;
    double R50 = R18*(R16*R45 + R23*R42 + R49);
    double R51 = 2*R50;
    double R52 = 2*R0;
    double R53 = R32*R39;
    double R54 = R18*(R16*R49 + R28*R42 + R45);
    double R55 = ((1.0/2.0)*R31*(R1*R41 + R1*R51) + (1.0/2.0)*R34*(R33*R37 - R38*R53 - R38*R54 - R40*R52 + 2*R47 - R50*R52) + (1.0/2.0)*R35*(R14*R41 + R14*R51 + R33*R38 + R37*R53 + R37*R54 - 2*R48))/sqrt(pow(R31, 2) + pow(R34, 2) + pow(R35, 2));
    double result = 16.0*((R33 > 0) ? (
   R55
)
: (
   -R55
));
    return result;
}


std::vector<float> CosmicTrack_DCA_LocalDeriv(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
        std::vector<float> result = {(float)CosmicTrack_DCA_Deriv_a0(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_b0(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_a1(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_b1(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_t0(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z)};
return result;
}

std::vector<float> CosmicTrack_DCA_GlobalDeriv(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
        std::vector<float> result = {(float)CosmicTrack_DCA_Deriv_plane_dx(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_plane_dy(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_plane_dz(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_plane_a(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_plane_b(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_plane_g(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_panel_dx(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_panel_dy(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_panel_dz(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_panel_a(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_panel_b(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_panel_g(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z)};
return result;
}
std::vector<double> CosmicTrack_DCA_LocalDeriv_double(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
        std::vector<double> result = {CosmicTrack_DCA_Deriv_a0(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_b0(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_a1(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_b1(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_t0(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z)};
return result;
}

