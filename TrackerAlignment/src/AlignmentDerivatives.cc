

# include "TrackerAlignment/inc/AlignmentDerivatives.hh"
# include <math.h>
# include <vector>

double CosmicTrack_DCA(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = pow(panel_straw0x, 2);
    double R2 = pow(panel_straw0y, 2);
    double R3 = R1 + R2;
    double R4 = sqrt(R3);
    double R5 = 1.0/R4;
    double R6 = 1.0/R3;
    double R7 = 1.0*R5/sqrt(R1*R6 + R2*R6);
    double R8 = R7*panel_straw0x;
    double R9 = R0*a1;
    double R10 = R7*panel_straw0y;
    double R11 = -R0*R8 - R10*R9;
    double R12 = 1.0/(1.0 - pow(R11, 2));
    double R13 = R5*(-R4 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R14 = R13*panel_straw0y + panel_straw0y;
    double R15 = R14*R8;
    double R16 = R13*panel_straw0x - a0 + panel_straw0x;
    double R17 = R10*R16;
    double R18 = -b0 + wire_z;
    double R19 = R0*b1;
    double R20 = -R0*R14 + R16*R9 + R18*R19;
    double R21 = R12*(-R11*(R15 - R17) + R20);
    double R22 = R12*(R11*R20 - R15 + R17);
    double R23 = sqrt(pow(R18 - R19*R21, 2) + pow(R0*R21 + R14 + R22*R8, 2) + pow(-R10*R22 + R16 - R21*R9, 2));
    double result = ((R21 > 0) ? (
   R23
)
: (
   -1.0*R23
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
    double R1 = pow(panel_straw0x, 2);
    double R2 = pow(panel_straw0y, 2);
    double R3 = R1 + R2;
    double R4 = sqrt(R3);
    double R5 = 1.0/R4;
    double R6 = 1.0/R3;
    double R7 = R5/sqrt(R1*R6 + R2*R6);
    double R8 = 1.0*R7;
    double R9 = R8*panel_straw0x;
    double R10 = R0*a1;
    double R11 = R8*panel_straw0y;
    double R12 = -R0*R9 - R10*R11;
    double R13 = 1.0/(1.0 - pow(R12, 2));
    double R14 = R5*(-R4 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R15 = R14*panel_straw0y + panel_straw0y;
    double R16 = R15*R9;
    double R17 = R14*panel_straw0x - a0 + panel_straw0x;
    double R18 = R11*R17;
    double R19 = -b0 + wire_z;
    double R20 = R0*b1;
    double R21 = -R0*R15 + R10*R17 + R19*R20;
    double R22 = R13*(-R12*(R16 - R18) + R21);
    double R23 = R19 - R20*R22;
    double R24 = R13*(R12*R21 - R16 + R18);
    double R25 = R0*R22 + R15 + R24*R9;
    double R26 = -R10*R22 - R11*R24 + R17;
    double R27 = R13*(-R10 - R11*R12);
    double R28 = 2*R27;
    double R29 = 2.0*R13*R7*(-R10*R12 - R11);
    double R30 = (-R20*R23*R27 + (1.0/2.0)*R25*(R0*R28 + R29*panel_straw0x) + (1.0/2.0)*R26*(-R10*R28 - R29*panel_straw0y - 2))/sqrt(pow(R23, 2) + pow(R25, 2) + pow(R26, 2));
    double result = 16.0*((R22 > 0) ? (
   R30
)
: (
   -1.0*R30
));
    return result;
}


double CosmicTrack_DCA_Deriv_b0(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
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
    double R15 = 1.0/(1.0 - pow(R14, 2));
    double R16 = R7*(-R6 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R17 = R16*panel_straw0y + panel_straw0y;
    double R18 = R11*R17;
    double R19 = R16*panel_straw0x - a0 + panel_straw0x;
    double R20 = R13*R19;
    double R21 = -b0 + wire_z;
    double R22 = R2*b1;
    double R23 = R12*R19 - R17*R2 + R21*R22;
    double R24 = R15*(-R14*(R18 - R20) + R23);
    double R25 = R21 - R22*R24;
    double R26 = R15*(R14*R23 - R18 + R20);
    double R27 = R11*R26 + R17 + R2*R24;
    double R28 = -R12*R24 - R13*R26 + R19;
    double R29 = 2*R15/R1;
    double R30 = R29*b1;
    double R31 = 2.0*R14*R15*R22*R9;
    double R32 = ((1.0/2.0)*R25*(R0*R29 - 2) + (1.0/2.0)*R27*(-R30 - R31*panel_straw0x) + (1.0/2.0)*R28*(R30*a1 + R31*panel_straw0y))/sqrt(pow(R25, 2) + pow(R27, 2) + pow(R28, 2));
    double result = 16.0*((R24 > 0) ? (
   R32
)
: (
   -1.0*R32
));
    return result;
}


double CosmicTrack_DCA_Deriv_a1(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(a1, 2);
    double R1 = R0 + pow(b1, 2) + 1;
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
    double R12 = R10*panel_straw0y;
    double R13 = R12*R2;
    double R14 = -R11*R2 - R13*a1;
    double R15 = 1.0 - pow(R14, 2);
    double R16 = 1.0/R15;
    double R17 = R7*(-R6 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R18 = R17*panel_straw0y + panel_straw0y;
    double R19 = R11*R18;
    double R20 = R17*panel_straw0x - a0 + panel_straw0x;
    double R21 = R12*R20;
    double R22 = -b0 + wire_z;
    double R23 = R2*b1;
    double R24 = R2*R20;
    double R25 = -R18*R2 + R22*R23 + R24*a1;
    double R26 = -R14*(R19 - R21) + R25;
    double R27 = R16*R26;
    double R28 = R22 - R23*R27;
    double R29 = R2*R27;
    double R30 = -R19 + R21;
    double R31 = R14*R25 + R30;
    double R32 = R16*R31;
    double R33 = R11*R32 + R18 + R29;
    double R34 = -R12*R32 + R20 - R29*a1;
    double R35 = pow(R1, -3.0/2.0);
    double R36 = R35*a1;
    double R37 = R36*b1;
    double R38 = 2*R27;
    double R39 = R0*R35;
    double R40 = R11*R36 + R12*R39 - R13;
    double R41 = R18*R36 - R20*R39 - R22*R37 + R24;
    double R42 = R16*(R30*R40 + R41);
    double R43 = 2*R23;
    double R44 = 2.0*R9;
    double R45 = R44*panel_straw0y;
    double R46 = R44*panel_straw0x;
    double R47 = R14*(-R2*R45 + R36*R46 + R39*R45)/pow(R15, 2);
    double R48 = R26*R47;
    double R49 = 2*R2;
    double R50 = R42*R49;
    double R51 = R16*(R14*R41 + R25*R40);
    double R52 = R48*R49;
    double R53 = R31*R47;
    double R54 = ((1.0/2.0)*R28*(R37*R38 - R42*R43 - R43*R48) + (1.0/2.0)*R33*(-R36*R38 + R46*R51 + R46*R53 + R50 + R52) + (1.0/2.0)*R34*(-2*R29 + R38*R39 - R45*R51 - R45*R53 - R50*a1 - R52*a1))/sqrt(pow(R28, 2) + pow(R33, 2) + pow(R34, 2));
    double result = 16.0*((R27 > 0) ? (
   R54
)
: (
   -1.0*R54
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
    double R15 = 1.0 - pow(R14, 2);
    double R16 = 1.0/R15;
    double R17 = R7*(-R6 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R18 = R17*panel_straw0y + panel_straw0y;
    double R19 = R11*R18;
    double R20 = R17*panel_straw0x - a0 + panel_straw0x;
    double R21 = R13*R20;
    double R22 = -b0 + wire_z;
    double R23 = R2*R22;
    double R24 = R12*R20 - R18*R2 + R23*b1;
    double R25 = -R14*(R19 - R21) + R24;
    double R26 = R16*R25;
    double R27 = R2*R26;
    double R28 = R22 - R27*b1;
    double R29 = -R19 + R21;
    double R30 = R14*R24 + R29;
    double R31 = R16*R30;
    double R32 = R11*R31 + R18 + R27;
    double R33 = -R12*R26 - R13*R31 + R20;
    double R34 = pow(R1, -3.0/2.0);
    double R35 = R0*R34;
    double R36 = 2*R26;
    double R37 = 2*R2;
    double R38 = R34*b1;
    double R39 = R38*a1;
    double R40 = R11*R38 + R13*R39;
    double R41 = R18*R38 - R20*R39 - R22*R35 + R23;
    double R42 = R16*(R29*R40 + R41);
    double R43 = R37*R42;
    double R44 = 2.0*R9;
    double R45 = R44*panel_straw0x;
    double R46 = R44*panel_straw0y;
    double R47 = R14*(R38*R45 + R39*R46)/pow(R15, 2);
    double R48 = R25*R47;
    double R49 = R37*R48;
    double R50 = R16*(R14*R41 + R24*R40);
    double R51 = R30*R47;
    double R52 = 2*R12;
    double R53 = ((1.0/2.0)*R28*(-2*R27 + R35*R36 - R43*b1 - R49*b1) + (1.0/2.0)*R32*(-R36*R38 + R43 + R45*R50 + R45*R51 + R49) + (1.0/2.0)*R33*(R36*R39 - R42*R52 - R46*R50 - R46*R51 - R48*R52))/sqrt(pow(R28, 2) + pow(R32, 2) + pow(R33, 2));
    double result = 16.0*((R26 > 0) ? (
   R53
)
: (
   -1.0*R53
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
    double R1 = pow(panel_straw0x, 2);
    double R2 = pow(panel_straw0y, 2);
    double R3 = R1 + R2;
    double R4 = sqrt(R3);
    double R5 = 1.0/R4;
    double R6 = 1.0/R3;
    double R7 = R5/sqrt(R1*R6 + R2*R6);
    double R8 = 1.0*R7;
    double R9 = R8*panel_straw0x;
    double R10 = R0*a1;
    double R11 = R8*panel_straw0y;
    double R12 = -R0*R9 - R10*R11;
    double R13 = 1.0/(1.0 - pow(R12, 2));
    double R14 = R5*(-R4 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R15 = R14*panel_straw0y + panel_straw0y;
    double R16 = R15*R9;
    double R17 = R14*panel_straw0x - a0 + panel_straw0x;
    double R18 = R11*R17;
    double R19 = -b0 + wire_z;
    double R20 = R0*b1;
    double R21 = -R0*R15 + R10*R17 + R19*R20;
    double R22 = R13*(-R12*(R16 - R18) + R21);
    double R23 = R19 - R20*R22;
    double R24 = R13*(R12*R21 - R16 + R18);
    double R25 = R0*R22 + R15 + R24*R9;
    double R26 = -R10*R22 - R11*R24 + R17;
    double R27 = R13*(R10 + R11*R12);
    double R28 = 2*R27;
    double R29 = 2.0*R13*R7*(R10*R12 + R11);
    double R30 = (-R20*R23*R27 + (1.0/2.0)*R25*(R0*R28 + R29*panel_straw0x) + (1.0/2.0)*R26*(-R10*R28 - R29*panel_straw0y + 2))/sqrt(pow(R23, 2) + pow(R25, 2) + pow(R26, 2));
    double result = 16.0*((R22 > 0) ? (
   R30
)
: (
   -1.0*R30
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dy(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = pow(panel_straw0x, 2);
    double R2 = pow(panel_straw0y, 2);
    double R3 = R1 + R2;
    double R4 = sqrt(R3);
    double R5 = 1.0/R4;
    double R6 = 1.0/R3;
    double R7 = R5/sqrt(R1*R6 + R2*R6);
    double R8 = 1.0*R7;
    double R9 = R8*panel_straw0x;
    double R10 = R0*a1;
    double R11 = R8*panel_straw0y;
    double R12 = -R0*R9 - R10*R11;
    double R13 = 1.0/(1.0 - pow(R12, 2));
    double R14 = R5*(-R4 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R15 = R14*panel_straw0y + panel_straw0y;
    double R16 = R15*R9;
    double R17 = R14*panel_straw0x - a0 + panel_straw0x;
    double R18 = R11*R17;
    double R19 = -b0 + wire_z;
    double R20 = R0*b1;
    double R21 = -R0*R15 + R10*R17 + R19*R20;
    double R22 = R13*(-R12*(R16 - R18) + R21);
    double R23 = R19 - R20*R22;
    double R24 = R13*(R12*R21 - R16 + R18);
    double R25 = R0*R22 + R15 + R24*R9;
    double R26 = -R10*R22 - R11*R24 + R17;
    double R27 = R13*(-R0 - R12*R9);
    double R28 = 2*R27;
    double R29 = 2.0*R13*R7*(-R0*R12 - R9);
    double R30 = (-R20*R23*R27 + (1.0/2.0)*R25*(R0*R28 + R29*panel_straw0x + 2) + (1.0/2.0)*R26*(-R10*R28 - R29*panel_straw0y))/sqrt(pow(R23, 2) + pow(R25, 2) + pow(R26, 2));
    double result = 16.0*((R22 > 0) ? (
   R30
)
: (
   -1.0*R30
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dz(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
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
    double R15 = 1.0/(1.0 - pow(R14, 2));
    double R16 = R7*(-R6 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R17 = R16*panel_straw0y + panel_straw0y;
    double R18 = R11*R17;
    double R19 = R16*panel_straw0x - a0 + panel_straw0x;
    double R20 = R13*R19;
    double R21 = -b0 + wire_z;
    double R22 = R2*b1;
    double R23 = R12*R19 - R17*R2 + R21*R22;
    double R24 = R15*(-R14*(R18 - R20) + R23);
    double R25 = R21 - R22*R24;
    double R26 = R15*(R14*R23 - R18 + R20);
    double R27 = R11*R26 + R17 + R2*R24;
    double R28 = -R12*R24 - R13*R26 + R19;
    double R29 = 2*R15/R1;
    double R30 = R29*b1;
    double R31 = 2.0*R14*R15*R22*R9;
    double R32 = ((1.0/2.0)*R25*(-R0*R29 + 2) + (1.0/2.0)*R27*(R30 + R31*panel_straw0x) + (1.0/2.0)*R28*(-R30*a1 - R31*panel_straw0y))/sqrt(pow(R25, 2) + pow(R27, 2) + pow(R28, 2));
    double result = 16.0*((R24 > 0) ? (
   R32
)
: (
   -1.0*R32
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_a(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
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
    double R9 = 1.0*R8;
    double R10 = R3*R9 + R4*R9;
    double R11 = R7/sqrt(R10);
    double R12 = 1.0*R11;
    double R13 = R12*panel_straw0x;
    double R14 = R2*a1;
    double R15 = R12*panel_straw0y;
    double R16 = -R13*R2 - R14*R15;
    double R17 = 1.0 - pow(R16, 2);
    double R18 = 1.0/R17;
    double R19 = R7*(-R6 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R20 = R19*panel_straw0y + panel_straw0y;
    double R21 = R13*R20;
    double R22 = R19*panel_straw0x - a0 + panel_straw0x;
    double R23 = R15*R22;
    double R24 = -b0 + wire_z;
    double R25 = R2*b1;
    double R26 = R2*R20;
    double R27 = R14*R22 + R24*R25 - R26;
    double R28 = -R16*(R21 - R23) + R27;
    double R29 = R18*R28;
    double R30 = R24 - R25*R29;
    double R31 = -R21 + R23;
    double R32 = R16*R27 + R31;
    double R33 = R18*R32;
    double R34 = R13*R33 + R2*R29 + R20;
    double R35 = -R14*R29 - R15*R33 + R22;
    double R36 = 2*panel_straw0y;
    double R37 = 2.0*R11;
    double R38 = R37*panel_straw0x;
    double R39 = plane_z - wire_z;
    double R40 = -R13*R24 - R13*R39;
    double R41 = R13*R25;
    double R42 = -R2*R39 + R26*b1;
    double R43 = 2*R18*(R16*R40 + R31*R41 + R42);
    double R44 = 4.0*R16/pow(R17, 2);
    double R45 = R44*panel_straw0x;
    double R46 = R11*R28*R45/R1;
    double R47 = R18*(R16*R42 + R27*R41 + R40);
    double R48 = R46*b1;
    double R49 = R25*R32*R8/R10;
    double R50 = ((1.0/2.0)*R30*(-R0*R46 + R19*R36 - R25*R43 + R33*R38 + R36) + (1.0/2.0)*R34*(R2*R43 + R3*R44*R49 + R38*R47 + R48 + 2*plane_z - 2*wire_z) + (1.0/2.0)*R35*(-R14*R43 - R37*R47*panel_straw0y - R45*R49*panel_straw0y - R48*a1))/sqrt(pow(R30, 2) + pow(R34, 2) + pow(R35, 2));
    double result = 16.0*((R29 > 0) ? (
   R50
)
: (
   -1.0*R50
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_b(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
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
    double R9 = 1.0*R8;
    double R10 = R3*R9 + R4*R9;
    double R11 = R7/sqrt(R10);
    double R12 = 1.0*R11;
    double R13 = R12*panel_straw0x;
    double R14 = R2*a1;
    double R15 = R12*panel_straw0y;
    double R16 = -R13*R2 - R14*R15;
    double R17 = 1.0 - pow(R16, 2);
    double R18 = 1.0/R17;
    double R19 = R7*(-R6 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R20 = R19*panel_straw0y + panel_straw0y;
    double R21 = R13*R20;
    double R22 = R19*panel_straw0x;
    double R23 = R22 - a0 + panel_straw0x;
    double R24 = R15*R23;
    double R25 = -b0 + wire_z;
    double R26 = R2*b1;
    double R27 = R14*R23 - R2*R20 + R25*R26;
    double R28 = -R16*(R21 - R24) + R27;
    double R29 = R18*R28;
    double R30 = R25 - R26*R29;
    double R31 = -R21 + R24;
    double R32 = R16*R27 + R31;
    double R33 = R18*R32;
    double R34 = R13*R33 + R2*R29 + R20;
    double R35 = -R14*R29 - R15*R33 + R23;
    double R36 = 2*panel_straw0x;
    double R37 = 2.0*R11;
    double R38 = R37*panel_straw0y;
    double R39 = -plane_z + wire_z;
    double R40 = -R15*R25 + R15*R39;
    double R41 = R15*R26;
    double R42 = R14*R39 + R26*(-R22 - panel_straw0x);
    double R43 = 2*R18*(R16*R40 + R31*R41 + R42);
    double R44 = 4.0*R16/pow(R17, 2);
    double R45 = R44*panel_straw0y;
    double R46 = R11*R28*R45/R1;
    double R47 = R18*(R16*R42 + R27*R41 + R40);
    double R48 = R46*b1;
    double R49 = R26*R32*R8/R10;
    double R50 = ((1.0/2.0)*R30*(-R0*R46 - R19*R36 - R26*R43 + R33*R38 - R36) + (1.0/2.0)*R34*(R2*R43 + R37*R47*panel_straw0x + R45*R49*panel_straw0x + R48) + (1.0/2.0)*R35*(-R14*R43 - R38*R47 - R4*R44*R49 - R48*a1 - 2*plane_z + 2*wire_z))/sqrt(pow(R30, 2) + pow(R34, 2) + pow(R35, 2));
    double result = 16.0*((R29 > 0) ? (
   R50
)
: (
   -1.0*R50
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_g(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = pow(panel_straw0x, 2);
    double R2 = pow(panel_straw0y, 2);
    double R3 = R1 + R2;
    double R4 = sqrt(R3);
    double R5 = 1.0/R4;
    double R6 = 1.0/R3;
    double R7 = R5/sqrt(R1*R6 + R2*R6);
    double R8 = 1.0*R7;
    double R9 = R8*panel_straw0x;
    double R10 = R0*a1;
    double R11 = R8*panel_straw0y;
    double R12 = -R0*R9 - R10*R11;
    double R13 = 1.0 - pow(R12, 2);
    double R14 = 1.0/R13;
    double R15 = R5*(-R4 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R16 = R15*panel_straw0y;
    double R17 = R16 + panel_straw0y;
    double R18 = R17*R9;
    double R19 = R15*panel_straw0x + panel_straw0x;
    double R20 = R19 - a0;
    double R21 = R11*R20;
    double R22 = -b0 + wire_z;
    double R23 = R0*b1;
    double R24 = -R0*R17 + R10*R20 + R22*R23;
    double R25 = -R12*(R18 - R21) + R24;
    double R26 = R14*R25;
    double R27 = R22 - R23*R26;
    double R28 = -R18 + R21;
    double R29 = R12*R24 + R28;
    double R30 = R14*R29;
    double R31 = R0*R26 + R17 + R30*R9;
    double R32 = -R10*R26 - R11*R30 + R20;
    double R33 = 2.0*R7;
    double R34 = R33*panel_straw0y;
    double R35 = R33*panel_straw0x;
    double R36 = R12*(R0*R34 - R10*R35)/pow(R13, 2);
    double R37 = R25*R36;
    double R38 = 2*R37;
    double R39 = R0*R11 - R10*R9;
    double R40 = -R16 - panel_straw0y;
    double R41 = R11*R17 + R11*R40 - R19*R9 + R20*R9;
    double R42 = -R0*R19 + R10*R40;
    double R43 = R14*(R12*R41 + R28*R39 + R42);
    double R44 = 2*R43;
    double R45 = 2*panel_straw0x;
    double R46 = 2*R0;
    double R47 = R29*R36;
    double R48 = R14*(R12*R42 + R24*R39 + R41);
    double R49 = 2*panel_straw0y;
    double R50 = ((1.0/2.0)*R27*(-R23*R38 - R23*R44) + (1.0/2.0)*R31*(R15*R45 - R30*R34 + R35*R47 + R35*R48 + R37*R46 + R43*R46 + R45) + (1.0/2.0)*R32*(-R10*R38 - R10*R44 - R15*R49 - R30*R35 - R34*R47 - R34*R48 - R49))/sqrt(pow(R27, 2) + pow(R31, 2) + pow(R32, 2));
    double result = 16.0*((R26 > 0) ? (
   R50
)
: (
   -1.0*R50
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dx(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = pow(panel_straw0x, 2);
    double R2 = pow(panel_straw0y, 2);
    double R3 = R1 + R2;
    double R4 = sqrt(R3);
    double R5 = 1.0/R4;
    double R6 = R0*R5;
    double R7 = 1.0/R3;
    double R8 = pow(R1*R7 + R2*R7, -1.0/2.0);
    double R9 = 1.0*R8;
    double R10 = R0*a1;
    double R11 = R5*panel_straw0y;
    double R12 = R11*R9;
    double R13 = -R10*R12 - R6*R9*panel_straw0x;
    double R14 = 1.0/(1.0 - pow(R13, 2));
    double R15 = R5*(-R4 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R16 = R15*panel_straw0y + panel_straw0y;
    double R17 = R5*panel_straw0x;
    double R18 = R17*R9;
    double R19 = R16*R18;
    double R20 = R15*panel_straw0x - a0 + panel_straw0x;
    double R21 = R12*R20;
    double R22 = -b0 + wire_z;
    double R23 = R0*b1;
    double R24 = -R0*R16 + R10*R20 + R22*R23;
    double R25 = R14*(-R13*(R19 - R21) + R24);
    double R26 = R22 - R23*R25;
    double R27 = R14*(R13*R24 - R19 + R21);
    double R28 = R0*R25 + R16 + R18*R27;
    double R29 = -R10*R25 - R12*R27 + R20;
    double R30 = R14*(R10*R17 - R6*panel_straw0y);
    double R31 = 2*R5;
    double R32 = 2*R30;
    double R33 = 2.0*R13*R30*R8;
    double R34 = (-R23*R26*R30 + (1.0/2.0)*R28*(R0*R32 + R17*R33 + R31*panel_straw0y) + (1.0/2.0)*R29*(-R10*R32 - R11*R33 + R31*panel_straw0x))/sqrt(pow(R26, 2) + pow(R28, 2) + pow(R29, 2));
    double result = 16.0*((R25 > 0) ? (
   R34
)
: (
   -1.0*R34
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dy(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = pow(panel_straw0x, 2);
    double R2 = pow(panel_straw0y, 2);
    double R3 = R1 + R2;
    double R4 = sqrt(R3);
    double R5 = 1.0/R4;
    double R6 = R5*panel_straw0x;
    double R7 = R0*R6;
    double R8 = 1.0/R3;
    double R9 = R1*R8;
    double R10 = R2*R8;
    double R11 = pow(R10 + R9, -1.0/2.0);
    double R12 = 1.0*R11;
    double R13 = R0*a1;
    double R14 = R5*panel_straw0y;
    double R15 = R13*R14;
    double R16 = -R12*R15 - R12*R7;
    double R17 = 1.0/(1.0 - pow(R16, 2));
    double R18 = R5*(-R4 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R19 = R18*panel_straw0y + panel_straw0y;
    double R20 = R12*R6;
    double R21 = R19*R20;
    double R22 = R18*panel_straw0x - a0 + panel_straw0x;
    double R23 = R12*R14;
    double R24 = R22*R23;
    double R25 = -b0 + wire_z;
    double R26 = R0*b1;
    double R27 = -R0*R19 + R13*R22 + R25*R26;
    double R28 = R17*(-R16*(R21 - R24) + R27);
    double R29 = R25 - R26*R28;
    double R30 = R17*(R16*R27 - R21 + R24);
    double R31 = R0*R28 + R19 + R20*R30;
    double R32 = -R13*R28 + R22 - R23*R30;
    double R33 = -R10*R11 - R11*R9;
    double R34 = -R15 - R7;
    double R35 = R17*(R16*R33 + R34);
    double R36 = 2*R5;
    double R37 = 2*R35;
    double R38 = 2.0*R11*R17*(R16*R34 + R33);
    double R39 = (-R26*R29*R35 + (1.0/2.0)*R31*(R0*R37 + R36*panel_straw0x + R38*R6) + (1.0/2.0)*R32*(-R13*R37 - R14*R38 - R36*panel_straw0y))/sqrt(pow(R29, 2) + pow(R31, 2) + pow(R32, 2));
    double result = 16.0*((R28 > 0) ? (
   R39
)
: (
   -1.0*R39
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dz(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
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
    double R15 = 1.0/(1.0 - pow(R14, 2));
    double R16 = R7*(-R6 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R17 = R16*panel_straw0y + panel_straw0y;
    double R18 = R11*R17;
    double R19 = R16*panel_straw0x - a0 + panel_straw0x;
    double R20 = R13*R19;
    double R21 = -b0 + wire_z;
    double R22 = R2*b1;
    double R23 = R12*R19 - R17*R2 + R21*R22;
    double R24 = R15*(-R14*(R18 - R20) + R23);
    double R25 = R21 - R22*R24;
    double R26 = R15*(R14*R23 - R18 + R20);
    double R27 = R11*R26 + R17 + R2*R24;
    double R28 = -R12*R24 - R13*R26 + R19;
    double R29 = 2*R15/R1;
    double R30 = R29*b1;
    double R31 = 2.0*R14*R15*R22*R9;
    double R32 = ((1.0/2.0)*R25*(-R0*R29 + 2) + (1.0/2.0)*R27*(R30 + R31*panel_straw0x) + (1.0/2.0)*R28*(-R30*a1 - R31*panel_straw0y))/sqrt(pow(R25, 2) + pow(R27, 2) + pow(R28, 2));
    double result = 16.0*((R24 > 0) ? (
   R32
)
: (
   -1.0*R32
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_a(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
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
    double R9 = R3*R8;
    double R10 = R4*R8;
    double R11 = R10 + R9;
    double R12 = pow(R11, -1.0/2.0);
    double R13 = 1.0*R12;
    double R14 = R13*R7;
    double R15 = R14*panel_straw0x;
    double R16 = R2*a1;
    double R17 = R14*panel_straw0y;
    double R18 = -R15*R2 - R16*R17;
    double R19 = 1.0 - pow(R18, 2);
    double R20 = 1.0/R19;
    double R21 = R7*(-R6 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R22 = R21*panel_straw0y + panel_straw0y;
    double R23 = R15*R22;
    double R24 = R21*panel_straw0x - a0 + panel_straw0x;
    double R25 = R17*R24;
    double R26 = -b0 + wire_z;
    double R27 = R2*b1;
    double R28 = R16*R24 - R2*R22 + R26*R27;
    double R29 = -R18*(R23 - R25) + R28;
    double R30 = R20*R29;
    double R31 = R26 - R27*R30;
    double R32 = -R23 + R25;
    double R33 = R18*R28 + R32;
    double R34 = R20*R33;
    double R35 = R15*R34 + R2*R30 + R22;
    double R36 = -R16*R30 - R17*R34 + R24;
    double R37 = 2.0*R12;
    double R38 = -panel_straw0z + wire_z;
    double R39 = R12*R38;
    double R40 = R10*R39 - R13*R26 + R39*R9;
    double R41 = R13*R27;
    double R42 = R38*R7;
    double R43 = R42*panel_straw0x;
    double R44 = R42*panel_straw0y;
    double R45 = R16*R44 + R2*R43;
    double R46 = 2*R20*(R18*R40 + R32*R41 + R45);
    double R47 = 4.0*R18/pow(R19, 2);
    double R48 = R12*R29*R47/R1;
    double R49 = R7*panel_straw0x;
    double R50 = R20*R37*(R18*R45 + R28*R41 + R40);
    double R51 = R48*b1;
    double R52 = R27*R33*R47/R11;
    double R53 = R7*panel_straw0y;
    double R54 = ((1.0/2.0)*R31*(-R0*R48 - R27*R46 + R34*R37) + (1.0/2.0)*R35*(R2*R46 - 2*R43 + R49*R50 + R49*R52 + R51) + (1.0/2.0)*R36*(-R16*R46 + 2*R44 - R50*R53 - R51*a1 - R52*R53))/sqrt(pow(R31, 2) + pow(R35, 2) + pow(R36, 2));
    double result = 16.0*((R30 > 0) ? (
   R54
)
: (
   -1.0*R54
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_b(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = pow(panel_straw0x, 2);
    double R2 = pow(panel_straw0y, 2);
    double R3 = R1 + R2;
    double R4 = sqrt(R3);
    double R5 = 1.0/R4;
    double R6 = 1.0/R3;
    double R7 = R5/sqrt(R1*R6 + R2*R6);
    double R8 = 1.0*R7;
    double R9 = R8*panel_straw0x;
    double R10 = R0*a1;
    double R11 = R8*panel_straw0y;
    double R12 = -R0*R9 - R10*R11;
    double R13 = 1.0/(1.0 - pow(R12, 2));
    double R14 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R15 = R5*(R14 - R4);
    double R16 = R15*panel_straw0y + panel_straw0y;
    double R17 = R16*R9;
    double R18 = R15*panel_straw0x - a0 + panel_straw0x;
    double R19 = R11*R18;
    double R20 = -b0 + wire_z;
    double R21 = R0*b1;
    double R22 = -R0*R16 + R10*R18 + R20*R21;
    double R23 = R13*(-R12*(R17 - R19) + R22);
    double R24 = R20 - R21*R23;
    double R25 = R13*(R12*R22 - R17 + R19);
    double R26 = R0*R23 + R16 + R25*R9;
    double R27 = -R10*R23 - R11*R25 + R18;
    double R28 = R5*(-panel_straw0z + wire_z);
    double R29 = R28*panel_straw0y;
    double R30 = R28*panel_straw0x;
    double R31 = R13*(-R0*R29 + R10*R30 + R21*(-R14 + R4));
    double R32 = 2*R31;
    double R33 = 2.0*R12*R31*R7;
    double R34 = ((1.0/2.0)*R24*(-2*R14 - R21*R32 + 2*R4) + (1.0/2.0)*R26*(R0*R32 + 2*R29 + R33*panel_straw0x) + (1.0/2.0)*R27*(-R10*R32 + 2*R30 - R33*panel_straw0y))/sqrt(pow(R24, 2) + pow(R26, 2) + pow(R27, 2));
    double result = 16.0*((R23 > 0) ? (
   R34
)
: (
   -1.0*R34
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_g(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = pow(panel_straw0x, 2);
    double R2 = pow(panel_straw0y, 2);
    double R3 = R1 + R2;
    double R4 = sqrt(R3);
    double R5 = 1.0/R4;
    double R6 = 1.0/R3;
    double R7 = R1*R6;
    double R8 = R2*R6;
    double R9 = pow(R7 + R8, -1.0/2.0);
    double R10 = R5*R9;
    double R11 = 1.0*R10;
    double R12 = R11*panel_straw0x;
    double R13 = R0*a1;
    double R14 = R11*panel_straw0y;
    double R15 = -R0*R12 - R13*R14;
    double R16 = 1.0 - pow(R15, 2);
    double R17 = 1.0/R16;
    double R18 = -R4 + sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R19 = R18*R5;
    double R20 = R19*panel_straw0y;
    double R21 = R20 + panel_straw0y;
    double R22 = R12*R21;
    double R23 = R19*panel_straw0x;
    double R24 = R23 - a0 + panel_straw0x;
    double R25 = R14*R24;
    double R26 = -b0 + wire_z;
    double R27 = R0*b1;
    double R28 = -R0*R21 + R13*R24 + R26*R27;
    double R29 = -R15*(R22 - R25) + R28;
    double R30 = R17*R29;
    double R31 = R26 - R27*R30;
    double R32 = -R22 + R25;
    double R33 = R15*R28 + R32;
    double R34 = R17*R33;
    double R35 = R0*R30 + R12*R34 + R21;
    double R36 = -R13*R30 - R14*R34 + R24;
    double R37 = 2.0*R10;
    double R38 = R37*panel_straw0y;
    double R39 = R37*panel_straw0x;
    double R40 = R15*(R0*R38 - R13*R39)/pow(R16, 2);
    double R41 = R29*R40;
    double R42 = 2*R41;
    double R43 = R0*R14 - R12*R13;
    double R44 = R18*R9;
    double R45 = R12*R24 + R14*R21 - R44*R7 - R44*R8;
    double R46 = -R0*R23 - R13*R20;
    double R47 = R17*(R15*R45 + R32*R43 + R46);
    double R48 = 2*R47;
    double R49 = 2*R0;
    double R50 = R33*R40;
    double R51 = R17*(R15*R46 + R28*R43 + R45);
    double R52 = ((1.0/2.0)*R31*(-R27*R42 - R27*R48) + (1.0/2.0)*R35*(2*R23 - R34*R38 + R39*R50 + R39*R51 + R41*R49 + R47*R49) + (1.0/2.0)*R36*(-R13*R42 - R13*R48 - 2*R20 - R34*R39 - R38*R50 - R38*R51))/sqrt(pow(R31, 2) + pow(R35, 2) + pow(R36, 2));
    double result = 16.0*((R30 > 0) ? (
   R52
)
: (
   -1.0*R52
));
    return result;
}


std::vector<double> CosmicTrack_DCA_LocalDeriv(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
        std::vector<double> result = {CosmicTrack_DCA_Deriv_a0(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_b0(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_a1(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_b1(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_t0(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z)};
return result;
}

std::vector<double> CosmicTrack_DCA_GlobalDeriv(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
        std::vector<double> result = {CosmicTrack_DCA_Deriv_plane_dx(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_plane_dy(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_plane_dz(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_plane_a(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_plane_b(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_plane_g(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_panel_dx(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_panel_dy(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_panel_dz(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_panel_a(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_panel_b(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_panel_g(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z)};
return result;
}

