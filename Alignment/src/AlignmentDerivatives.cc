

# include "Alignment/inc/AlignmentDerivatives.hh"
# include <math.h>
# include <vector>

double CosmicTrack_DCA(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*a1;
    double R2 = R0*b1;
    double R3 = -R0*wdir_y + R1*wdir_x + R2*wdir_z;
    double R4 = 1.0/(1 - pow(R3, 2));
    double R5 = R0*wire_y;
    double R6 = a0 - wire_x;
    double R7 = R1*R6;
    double R8 = b0 - wire_z;
    double R9 = R2*R8;
    double R10 = R6*wdir_x + R8*wdir_z - wdir_y*wire_y;
    double R11 = R4*(R10 - R3*(R5 + R7 + R9));
    double R12 = R4*(R10*R3 - R5 - R7 - R9);
    double R13 = sqrt(pow(-R0*R12 - R11*wdir_y - wire_y, 2) + pow(R1*R12 - R11*wdir_x + R6, 2) + pow(-R11*wdir_z + R12*R2 + R8, 2));
    double result = ((R11 > 0) ? (
   R13
)
: (
   -R13
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
    double R2 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R3 = 1.0/R2;
    double R4 = 1.0*R3;
    double R5 = R4*panel_straw0x;
    double R6 = R0*R5;
    double R7 = R0*a1;
    double R8 = R4*panel_straw0y;
    double R9 = R7*R8;
    double R10 = -R6 - R9;
    double R11 = 1.0/(1 - pow(R10, 2));
    double R12 = b0 - wire_z;
    double R13 = R1*R12;
    double R14 = R3*(-R2 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R15 = -R14*panel_straw0y - panel_straw0y;
    double R16 = R0*R15;
    double R17 = -R14*panel_straw0x + a0 - panel_straw0x;
    double R18 = R17*R7;
    double R19 = R15*R5 - R17*R8;
    double R20 = R11*(R10*R19 - R13 + R16 - R18);
    double R21 = R1*R20 + R12;
    double R22 = R11*(-R10*(R13 - R16 + R18) + R19);
    double R23 = -R0*R20 + R15 - R22*R5;
    double R24 = R17 + R20*R7 + R22*R8;
    double R25 = R11*(-R10*R8 - R7);
    double R26 = 2.0*R11*R3*(R7*(R6 + R9) - R8);
    double R27 = 2*R25;
    double R28 = (R1*R21*R25 + (1.0/2.0)*R23*(-R0*R27 - R26*panel_straw0x) + (1.0/2.0)*R24*(R26*panel_straw0y + R27*R7 + 2))/sqrt(pow(R21, 2) + pow(R23, 2) + pow(R24, 2));
    double result = 16.0*((R22 > 0) ? (
   R28
)
: (
   -R28
));
    return result;
}


double CosmicTrack_DCA_Deriv_b0(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R5 = 1.0/R4;
    double R6 = 1.0*R5;
    double R7 = R6*panel_straw0x;
    double R8 = R2*R7;
    double R9 = R2*a1;
    double R10 = R6*panel_straw0y;
    double R11 = R10*R9;
    double R12 = -R11 - R8;
    double R13 = 1.0/(1 - pow(R12, 2));
    double R14 = b0 - wire_z;
    double R15 = R14*R3;
    double R16 = R5*(-R4 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R17 = -R16*panel_straw0y - panel_straw0y;
    double R18 = R17*R2;
    double R19 = -R16*panel_straw0x + a0 - panel_straw0x;
    double R20 = R19*R9;
    double R21 = -R10*R19 + R17*R7;
    double R22 = R13*(R12*R21 - R15 + R18 - R20);
    double R23 = R14 + R22*R3;
    double R24 = R13*(-R12*(R15 - R18 + R20) + R21);
    double R25 = R17 - R2*R22 - R24*R7;
    double R26 = R10*R24 + R19 + R22*R9;
    double R27 = 2*R13/R1;
    double R28 = R27*b1;
    double R29 = 2.0*R13*R3*R5*(R11 + R8);
    double R30 = ((1.0/2.0)*R23*(-R0*R27 + 2) + (1.0/2.0)*R25*(R28 - R29*panel_straw0x) + (1.0/2.0)*R26*(-R28*a1 + R29*panel_straw0y))/sqrt(pow(R23, 2) + pow(R25, 2) + pow(R26, 2));
    double result = 16.0*((R24 > 0) ? (
   R30
)
: (
   -R30
));
    return result;
}


double CosmicTrack_DCA_Deriv_a1(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(a1, 2);
    double R1 = R0 + pow(b1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R5 = 1.0/R4;
    double R6 = 1.0*R5;
    double R7 = R6*panel_straw0x;
    double R8 = R2*R7;
    double R9 = R6*panel_straw0y;
    double R10 = R2*R9;
    double R11 = R10*a1;
    double R12 = -R11 - R8;
    double R13 = 1 - pow(R12, 2);
    double R14 = 1.0/R13;
    double R15 = b0 - wire_z;
    double R16 = R15*R3;
    double R17 = R5*(-R4 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R18 = -R17*panel_straw0y - panel_straw0y;
    double R19 = R18*R2;
    double R20 = -R17*panel_straw0x + a0 - panel_straw0x;
    double R21 = R2*R20;
    double R22 = R21*a1;
    double R23 = R18*R7 - R20*R9;
    double R24 = R12*R23 - R16 + R19 - R22;
    double R25 = R14*R24;
    double R26 = R15 + R25*R3;
    double R27 = R2*R25;
    double R28 = R16 - R19 + R22;
    double R29 = -R12*R28 + R23;
    double R30 = R14*R29;
    double R31 = R18 - R27 - R30*R7;
    double R32 = R20 + R27*a1 + R30*R9;
    double R33 = pow(R1, -3.0/2.0);
    double R34 = R33*a1;
    double R35 = R34*b1;
    double R36 = 2*R25;
    double R37 = R15*R35;
    double R38 = R18*R34;
    double R39 = R0*R33;
    double R40 = R20*R39;
    double R41 = R34*R7;
    double R42 = R39*R9;
    double R43 = R14*(-R21 + R23*(-R10 + R41 + R42) + R37 - R38 + R40);
    double R44 = 2*R3;
    double R45 = 2.0*R5;
    double R46 = R45*panel_straw0y;
    double R47 = R45*panel_straw0x;
    double R48 = R12*(-R2*R46 + R34*R47 + R39*R46)/pow(R13, 2);
    double R49 = R24*R48;
    double R50 = 2*R2;
    double R51 = R43*R50;
    double R52 = R49*R50;
    double R53 = R29*R48;
    double R54 = R14*(R28*(R10 - R41 - R42) + (R11 + R8)*(R21 - R37 + R38 - R40));
    double R55 = ((1.0/2.0)*R26*(-R35*R36 + R43*R44 + R44*R49) + (1.0/2.0)*R31*(R34*R36 - R47*R53 - R47*R54 - R51 - R52) + (1.0/2.0)*R32*(2*R27 - R36*R39 + R46*R53 + R46*R54 + R51*a1 + R52*a1))/sqrt(pow(R26, 2) + pow(R31, 2) + pow(R32, 2));
    double result = 16.0*((R30 > 0) ? (
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
    double R3 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R4 = 1.0/R3;
    double R5 = 1.0*R4;
    double R6 = R5*panel_straw0x;
    double R7 = R2*R6;
    double R8 = R2*a1;
    double R9 = R5*panel_straw0y;
    double R10 = R8*R9;
    double R11 = -R10 - R7;
    double R12 = 1 - pow(R11, 2);
    double R13 = 1.0/R12;
    double R14 = b0 - wire_z;
    double R15 = R14*R2;
    double R16 = R15*b1;
    double R17 = R4*(-R3 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R18 = -R17*panel_straw0y - panel_straw0y;
    double R19 = R18*R2;
    double R20 = -R17*panel_straw0x + a0 - panel_straw0x;
    double R21 = R20*R8;
    double R22 = R18*R6 - R20*R9;
    double R23 = R11*R22 - R16 + R19 - R21;
    double R24 = R13*R23;
    double R25 = R2*R24;
    double R26 = R14 + R25*b1;
    double R27 = R16 - R19 + R21;
    double R28 = -R11*R27 + R22;
    double R29 = R13*R28;
    double R30 = R18 - R25 - R29*R6;
    double R31 = R20 + R24*R8 + R29*R9;
    double R32 = pow(R1, -3.0/2.0);
    double R33 = R0*R32;
    double R34 = 2*R24;
    double R35 = 2*R2;
    double R36 = R14*R33;
    double R37 = R32*b1;
    double R38 = R18*R37;
    double R39 = R37*a1;
    double R40 = R20*R39;
    double R41 = R37*R6;
    double R42 = R39*R9;
    double R43 = R13*(-R15 + R22*(R41 + R42) + R36 - R38 + R40);
    double R44 = R35*R43;
    double R45 = 2.0*R4;
    double R46 = R45*panel_straw0x;
    double R47 = R45*panel_straw0y;
    double R48 = R11*(R37*R46 + R39*R47)/pow(R12, 2);
    double R49 = R23*R48;
    double R50 = R35*R49;
    double R51 = R13*(R27*(-R41 - R42) + (R10 + R7)*(R15 - R36 + R38 - R40));
    double R52 = R28*R48;
    double R53 = 2*R8;
    double R54 = ((1.0/2.0)*R26*(2*R25 - R33*R34 + R44*b1 + R50*b1) + (1.0/2.0)*R30*(R34*R37 - R44 - R46*R51 - R46*R52 - R50) + (1.0/2.0)*R31*(-R34*R39 + R43*R53 + R47*R51 + R47*R52 + R49*R53))/sqrt(pow(R26, 2) + pow(R30, 2) + pow(R31, 2));
    double result = 16.0*((R29 > 0) ? (
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
    double R2 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R3 = 1.0/R2;
    double R4 = 1.0*R3;
    double R5 = R4*panel_straw0x;
    double R6 = R0*R5;
    double R7 = R0*a1;
    double R8 = R4*panel_straw0y;
    double R9 = R7*R8;
    double R10 = -R6 - R9;
    double R11 = 1.0/(1 - pow(R10, 2));
    double R12 = b0 - wire_z;
    double R13 = R1*R12;
    double R14 = R3*(-R2 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R15 = -R14*panel_straw0y - panel_straw0y;
    double R16 = R0*R15;
    double R17 = -R14*panel_straw0x + a0 - panel_straw0x;
    double R18 = R17*R7;
    double R19 = R15*R5 - R17*R8;
    double R20 = R11*(R10*R19 - R13 + R16 - R18);
    double R21 = R1*R20 + R12;
    double R22 = R11*(-R10*(R13 - R16 + R18) + R19);
    double R23 = -R0*R20 + R15 - R22*R5;
    double R24 = R17 + R20*R7 + R22*R8;
    double R25 = R11*(R10*R8 + R7);
    double R26 = 2*R25;
    double R27 = 2.0*R11*R3*(-R7*(R6 + R9) + R8);
    double R28 = (R1*R21*R25 + (1.0/2.0)*R23*(-R0*R26 - R27*panel_straw0x) + (1.0/2.0)*R24*(R26*R7 + R27*panel_straw0y - 2))/sqrt(pow(R21, 2) + pow(R23, 2) + pow(R24, 2));
    double result = 16.0*((R22 > 0) ? (
   R28
)
: (
   -R28
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dy(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R3 = 1.0/R2;
    double R4 = 1.0*R3;
    double R5 = R4*panel_straw0x;
    double R6 = R0*R5;
    double R7 = R0*a1;
    double R8 = R4*panel_straw0y;
    double R9 = R7*R8;
    double R10 = -R6 - R9;
    double R11 = 1.0/(1 - pow(R10, 2));
    double R12 = b0 - wire_z;
    double R13 = R1*R12;
    double R14 = R3*(-R2 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R15 = -R14*panel_straw0y - panel_straw0y;
    double R16 = R0*R15;
    double R17 = -R14*panel_straw0x + a0 - panel_straw0x;
    double R18 = R17*R7;
    double R19 = R15*R5 - R17*R8;
    double R20 = R11*(R10*R19 - R13 + R16 - R18);
    double R21 = R1*R20 + R12;
    double R22 = R11*(-R10*(R13 - R16 + R18) + R19);
    double R23 = -R0*R20 + R15 - R22*R5;
    double R24 = R17 + R20*R7 + R22*R8;
    double R25 = R11*(-R0 - R10*R5);
    double R26 = 2.0*R11*R3*(R0*(R6 + R9) - R5);
    double R27 = 2*R25;
    double R28 = (R1*R21*R25 + (1.0/2.0)*R23*(-R0*R27 - R26*panel_straw0x - 2) + (1.0/2.0)*R24*(R26*panel_straw0y + R27*R7))/sqrt(pow(R21, 2) + pow(R23, 2) + pow(R24, 2));
    double result = 16.0*((R22 > 0) ? (
   R28
)
: (
   -R28
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dz(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R5 = 1.0/R4;
    double R6 = 1.0*R5;
    double R7 = R6*panel_straw0x;
    double R8 = R2*R7;
    double R9 = R2*a1;
    double R10 = R6*panel_straw0y;
    double R11 = R10*R9;
    double R12 = -R11 - R8;
    double R13 = 1.0/(1 - pow(R12, 2));
    double R14 = b0 - wire_z;
    double R15 = R14*R3;
    double R16 = R5*(-R4 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R17 = -R16*panel_straw0y - panel_straw0y;
    double R18 = R17*R2;
    double R19 = -R16*panel_straw0x + a0 - panel_straw0x;
    double R20 = R19*R9;
    double R21 = -R10*R19 + R17*R7;
    double R22 = R13*(R12*R21 - R15 + R18 - R20);
    double R23 = R14 + R22*R3;
    double R24 = R13*(-R12*(R15 - R18 + R20) + R21);
    double R25 = R17 - R2*R22 - R24*R7;
    double R26 = R10*R24 + R19 + R22*R9;
    double R27 = 2*R13/R1;
    double R28 = R27*b1;
    double R29 = 2.0*R13*R3*R5*(R11 + R8);
    double R30 = ((1.0/2.0)*R23*(R0*R27 - 2) + (1.0/2.0)*R25*(-R28 + R29*panel_straw0x) + (1.0/2.0)*R26*(R28*a1 - R29*panel_straw0y))/sqrt(pow(R23, 2) + pow(R25, 2) + pow(R26, 2));
    double result = 16.0*((R24 > 0) ? (
   R30
)
: (
   -R30
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
    double R5 = R4 + pow(panel_straw0y, 2);
    double R6 = sqrt(R5);
    double R7 = 1.0/R6;
    double R8 = R7*panel_straw0x;
    double R9 = 1.0*R8;
    double R10 = R2*R9;
    double R11 = R2*a1;
    double R12 = R7*panel_straw0y;
    double R13 = 1.0*R12;
    double R14 = R11*R13;
    double R15 = -R10 - R14;
    double R16 = 1 - pow(R15, 2);
    double R17 = 1.0/R16;
    double R18 = b0 - wire_z;
    double R19 = R18*R3;
    double R20 = -panel_straw0y;
    double R21 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R22 = R7*(R21 - R6);
    double R23 = R20 - R22*panel_straw0y;
    double R24 = R2*R23;
    double R25 = -R22*panel_straw0x + a0 - panel_straw0x;
    double R26 = R11*R25;
    double R27 = -R13*R25 + R23*R9;
    double R28 = R15*R27 - R19 + R24 - R26;
    double R29 = R17*R28;
    double R30 = R18 + R29*R3;
    double R31 = R19 - R24 + R26;
    double R32 = -R15*R31 + R27;
    double R33 = R17*R32;
    double R34 = -R2*R29 + R23 - R33*R9;
    double R35 = R11*R29 + R13*R33 + R25;
    double R36 = 2*panel_straw0y;
    double R37 = -R21 + R6;
    double R38 = 2.0*R8;
    double R39 = -plane_z + wire_z;
    double R40 = R2*R39;
    double R41 = R3*(R12*R37 + R20);
    double R42 = R18*R9 + R39*R9;
    double R43 = R3*R9;
    double R44 = 2*R17*(R15*R42 + R27*R43 + R40 - R41);
    double R45 = 4.0*R15/pow(R16, 2);
    double R46 = R28*R45*R8/R1;
    double R47 = R17*(-R31*R43 + R42 + (R10 + R14)*(-R40 + R41));
    double R48 = R46*b1;
    double R49 = R3*R32*R45/R5;
    double R50 = ((1.0/2.0)*R30*(R0*R46 + R3*R44 - R33*R38 + R36*R37*R7 - R36) + (1.0/2.0)*R34*(-R2*R44 - R38*R47 - R4*R49 - R48 - 2*plane_z + 2*wire_z) + (1.0/2.0)*R35*(R11*R44 + 2.0*R12*R47 + R48*a1 + R49*panel_straw0x*panel_straw0y))/sqrt(pow(R30, 2) + pow(R34, 2) + pow(R35, 2));
    double result = 16.0*((R33 > 0) ? (
   R50
)
: (
   -R50
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_b(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = pow(panel_straw0y, 2);
    double R5 = R4 + pow(panel_straw0x, 2);
    double R6 = sqrt(R5);
    double R7 = 1.0/R6;
    double R8 = R7*panel_straw0x;
    double R9 = 1.0*R8;
    double R10 = R2*R9;
    double R11 = R2*a1;
    double R12 = R7*panel_straw0y;
    double R13 = 1.0*R12;
    double R14 = R11*R13;
    double R15 = -R10 - R14;
    double R16 = 1 - pow(R15, 2);
    double R17 = 1.0/R16;
    double R18 = -wire_z;
    double R19 = R18 + b0;
    double R20 = R19*R3;
    double R21 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R22 = R7*(R21 - R6);
    double R23 = -R22*panel_straw0y - panel_straw0y;
    double R24 = R2*R23;
    double R25 = -R22*panel_straw0x + a0 - panel_straw0x;
    double R26 = R11*R25;
    double R27 = -R13*R25 + R23*R9;
    double R28 = R15*R27 - R20 + R24 - R26;
    double R29 = R17*R28;
    double R30 = R19 + R29*R3;
    double R31 = R20 - R24 + R26;
    double R32 = -R15*R31 + R27;
    double R33 = R17*R32;
    double R34 = -R2*R29 + R23 - R33*R9;
    double R35 = R11*R29 + R13*R33 + R25;
    double R36 = 2*panel_straw0x;
    double R37 = -R21 + R6;
    double R38 = 2.0*R12;
    double R39 = R18 + plane_z;
    double R40 = R11*R39;
    double R41 = R3*(-R37*R8 + panel_straw0x);
    double R42 = R13*R19 - R13*R39;
    double R43 = R13*R3;
    double R44 = 2*R17*(R15*R42 + R27*R43 - R40 - R41);
    double R45 = 4.0*R15/pow(R16, 2);
    double R46 = R12*R28*R45/R1;
    double R47 = R17*(-R31*R43 + R42 + (R10 + R14)*(R40 + R41));
    double R48 = R46*b1;
    double R49 = R3*R32*R45/R5;
    double R50 = ((1.0/2.0)*R30*(R0*R46 + R3*R44 - R33*R38 - R36*R37*R7 + R36) + (1.0/2.0)*R34*(-R2*R44 - 2.0*R47*R8 - R48 - R49*panel_straw0x*panel_straw0y) + (1.0/2.0)*R35*(R11*R44 + R38*R47 + R4*R49 + R48*a1 + 2*plane_z - 2*wire_z))/sqrt(pow(R30, 2) + pow(R34, 2) + pow(R35, 2));
    double result = 16.0*((R33 > 0) ? (
   R50
)
: (
   -R50
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_g(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R3 = 1.0/R2;
    double R4 = 1.0*R3;
    double R5 = R4*panel_straw0x;
    double R6 = R0*R5;
    double R7 = R0*a1;
    double R8 = R4*panel_straw0y;
    double R9 = R7*R8;
    double R10 = -R6 - R9;
    double R11 = 1 - pow(R10, 2);
    double R12 = 1.0/R11;
    double R13 = b0 - wire_z;
    double R14 = R1*R13;
    double R15 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R16 = R3*(R15 - R2);
    double R17 = -R16*panel_straw0y - panel_straw0y;
    double R18 = R0*R17;
    double R19 = -panel_straw0x;
    double R20 = -R16*panel_straw0x + R19 + a0;
    double R21 = R20*R7;
    double R22 = R17*R5 - R20*R8;
    double R23 = R10*R22 - R14 + R18 - R21;
    double R24 = R12*R23;
    double R25 = R1*R24 + R13;
    double R26 = R14 - R18 + R21;
    double R27 = -R10*R26 + R22;
    double R28 = R12*R27;
    double R29 = -R0*R24 + R17 - R28*R5;
    double R30 = R20 + R24*R7 + R28*R8;
    double R31 = 2.0*R3;
    double R32 = R31*panel_straw0y;
    double R33 = R31*panel_straw0x;
    double R34 = R10*(R0*R32 - R33*R7)/pow(R11, 2);
    double R35 = R23*R34;
    double R36 = 2*R35;
    double R37 = R3*(-R15 + R2);
    double R38 = R19 + R37*panel_straw0x;
    double R39 = R0*R38;
    double R40 = -R37*panel_straw0y + panel_straw0y;
    double R41 = R40*R7;
    double R42 = R0*R8;
    double R43 = R5*R7;
    double R44 = -R17*R8 - R20*R5 + R38*R5 - R40*R8;
    double R45 = R12*(R10*R44 + R22*(R42 - R43) + R39 - R41);
    double R46 = 2*R45;
    double R47 = 2*panel_straw0x;
    double R48 = 2*R0;
    double R49 = R27*R34;
    double R50 = R12*(R26*(-R42 + R43) + R44 + (-R39 + R41)*(R6 + R9));
    double R51 = 2*panel_straw0y;
    double R52 = ((1.0/2.0)*R25*(R1*R36 + R1*R46) + (1.0/2.0)*R29*(R28*R32 - R33*R49 - R33*R50 - R35*R48 + R37*R47 - R45*R48 - R47) + (1.0/2.0)*R30*(R28*R33 + R32*R49 + R32*R50 + R36*R7 - R37*R51 + R46*R7 + R51))/sqrt(pow(R25, 2) + pow(R29, 2) + pow(R30, 2));
    double result = 16.0*((R28 > 0) ? (
   R52
)
: (
   -R52
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dx(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R3 = 1.0/R2;
    double R4 = R0*R3;
    double R5 = 1.0*panel_straw0x;
    double R6 = R4*R5;
    double R7 = R0*a1;
    double R8 = R3*panel_straw0y;
    double R9 = 1.0*R8;
    double R10 = R7*R9;
    double R11 = -R10 - R6;
    double R12 = 1.0/(1 - pow(R11, 2));
    double R13 = b0 - wire_z;
    double R14 = R1*R13;
    double R15 = R3*(-R2 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R16 = -R15*panel_straw0y - panel_straw0y;
    double R17 = R0*R16;
    double R18 = -R15*panel_straw0x + a0 - panel_straw0x;
    double R19 = R18*R7;
    double R20 = R3*R5;
    double R21 = R16*R20 - R18*R9;
    double R22 = R12*(R11*R21 - R14 + R17 - R19);
    double R23 = R1*R22 + R13;
    double R24 = R12*(-R11*(R14 - R17 + R19) + R21);
    double R25 = -R0*R22 + R16 - R20*R24;
    double R26 = R18 + R22*R7 + R24*R9;
    double R27 = R4*panel_straw0y;
    double R28 = R3*panel_straw0x;
    double R29 = R28*R7;
    double R30 = R12*(-R27 + R29);
    double R31 = 2*R3;
    double R32 = 2*R30;
    double R33 = 2.0*R12*(R10 + R6)*(R27 - R29);
    double R34 = (R1*R23*R30 + (1.0/2.0)*R25*(-R0*R32 - R28*R33 - R31*panel_straw0y) + (1.0/2.0)*R26*(-R31*panel_straw0x + R32*R7 + R33*R8))/sqrt(pow(R23, 2) + pow(R25, 2) + pow(R26, 2));
    double result = 16.0*((R24 > 0) ? (
   R34
)
: (
   -R34
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
    double R9 = 1.0*R8;
    double R10 = R0*a1;
    double R11 = R6*panel_straw0y;
    double R12 = 1.0*R11;
    double R13 = R10*R12;
    double R14 = -R13 - R9;
    double R15 = 1.0/(1 - pow(R14, 2));
    double R16 = b0 - wire_z;
    double R17 = R1*R16;
    double R18 = R6*(-R5 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R19 = -R18*panel_straw0y - panel_straw0y;
    double R20 = R0*R19;
    double R21 = -R18*panel_straw0x + a0 - panel_straw0x;
    double R22 = R10*R21;
    double R23 = 1.0*R7;
    double R24 = -R12*R21 + R19*R23;
    double R25 = R15*(R14*R24 - R17 + R20 - R22);
    double R26 = R1*R25 + R16;
    double R27 = R15*(-R14*(R17 - R20 + R22) + R24);
    double R28 = -R0*R25 + R19 - R23*R27;
    double R29 = R10*R25 + R12*R27 + R21;
    double R30 = R10*R11;
    double R31 = 1.0/R4;
    double R32 = -R2*R31 - R3*R31;
    double R33 = R15*(R14*R32 - R30 - R8);
    double R34 = 2.0*R15*(R32 + (R13 + R9)*(R30 + R8));
    double R35 = 2*R33;
    double R36 = (R1*R26*R33 + (1.0/2.0)*R28*(-R0*R35 - R34*R7 - 2*R7) + (1.0/2.0)*R29*(R10*R35 + R11*R34 + 2*R11))/sqrt(pow(R26, 2) + pow(R28, 2) + pow(R29, 2));
    double result = 16.0*((R27 > 0) ? (
   R36
)
: (
   -R36
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dz(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R5 = 1.0/R4;
    double R6 = 1.0*R5;
    double R7 = R6*panel_straw0x;
    double R8 = R2*R7;
    double R9 = R2*a1;
    double R10 = R6*panel_straw0y;
    double R11 = R10*R9;
    double R12 = -R11 - R8;
    double R13 = 1.0/(1 - pow(R12, 2));
    double R14 = b0 - wire_z;
    double R15 = R14*R3;
    double R16 = R5*(-R4 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R17 = -R16*panel_straw0y - panel_straw0y;
    double R18 = R17*R2;
    double R19 = -R16*panel_straw0x + a0 - panel_straw0x;
    double R20 = R19*R9;
    double R21 = -R10*R19 + R17*R7;
    double R22 = R13*(R12*R21 - R15 + R18 - R20);
    double R23 = R14 + R22*R3;
    double R24 = R13*(-R12*(R15 - R18 + R20) + R21);
    double R25 = R17 - R2*R22 - R24*R7;
    double R26 = R10*R24 + R19 + R22*R9;
    double R27 = 2*R13/R1;
    double R28 = R27*b1;
    double R29 = 2.0*R13*R3*R5*(R11 + R8);
    double R30 = ((1.0/2.0)*R23*(R0*R27 - 2) + (1.0/2.0)*R25*(-R28 + R29*panel_straw0x) + (1.0/2.0)*R26*(R28*a1 - R29*panel_straw0y))/sqrt(pow(R23, 2) + pow(R25, 2) + pow(R26, 2));
    double result = 16.0*((R24 > 0) ? (
   R30
)
: (
   -R30
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
    double R9 = 1.0*R8;
    double R10 = R9*panel_straw0x;
    double R11 = R10*R2;
    double R12 = R2*a1;
    double R13 = R9*panel_straw0y;
    double R14 = R12*R13;
    double R15 = -R11 - R14;
    double R16 = 1 - pow(R15, 2);
    double R17 = 1.0/R16;
    double R18 = -wire_z;
    double R19 = R18 + b0;
    double R20 = R19*R3;
    double R21 = R8*(-R7 + sqrt(pow(wire_x, 2) + pow(wire_y, 2)));
    double R22 = -R21*panel_straw0y - panel_straw0y;
    double R23 = R2*R22;
    double R24 = -R21*panel_straw0x + a0 - panel_straw0x;
    double R25 = R12*R24;
    double R26 = R10*R22 - R13*R24;
    double R27 = R15*R26 - R20 + R23 - R25;
    double R28 = R17*R27;
    double R29 = R19 + R28*R3;
    double R30 = R20 - R23 + R25;
    double R31 = -R15*R30 + R26;
    double R32 = R17*R31;
    double R33 = -R10*R32 - R2*R28 + R22;
    double R34 = R12*R28 + R13*R32 + R24;
    double R35 = R18 + panel_straw0z;
    double R36 = R35*R8;
    double R37 = R36*panel_straw0x;
    double R38 = R2*R37;
    double R39 = R36*panel_straw0y;
    double R40 = R12*R39;
    double R41 = 1.0*R35/R6;
    double R42 = -R4*R41 - R41*R5 + 1.0*b0 - 1.0*wire_z;
    double R43 = 1.0*R3;
    double R44 = 2*R17*(R15*R42 + R26*R43 - R38 - R40);
    double R45 = 4.0*R15/pow(R16, 2);
    double R46 = R27*R45/R1;
    double R47 = R8*panel_straw0x;
    double R48 = 2.0*R17*(-R30*R43 + R42 + (R11 + R14)*(R38 + R40));
    double R49 = R46*b1;
    double R50 = R3*R31*R45;
    double R51 = R8*panel_straw0y;
    double R52 = ((1.0/2.0)*R29*(R0*R46 + R3*R44 - 2.0*R32) + (1.0/2.0)*R33*(-R2*R44 - 2*R37 - R47*R48 - R47*R50 - R49) + (1.0/2.0)*R34*(R12*R44 + 2*R39 + R48*R51 + R49*a1 + R50*R51))/sqrt(pow(R29, 2) + pow(R33, 2) + pow(R34, 2));
    double result = 16.0*((R32 > 0) ? (
   R52
)
: (
   -R52
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_b(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R3 = 1.0/R2;
    double R4 = 1.0*R3;
    double R5 = R4*panel_straw0x;
    double R6 = R0*R5;
    double R7 = R0*a1;
    double R8 = R4*panel_straw0y;
    double R9 = R7*R8;
    double R10 = -R6 - R9;
    double R11 = 1.0/(1 - pow(R10, 2));
    double R12 = -wire_z;
    double R13 = R12 + b0;
    double R14 = R1*R13;
    double R15 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R16 = R15 - R2;
    double R17 = R16*R3;
    double R18 = -R17*panel_straw0y - panel_straw0y;
    double R19 = R0*R18;
    double R20 = -R17*panel_straw0x + a0 - panel_straw0x;
    double R21 = R20*R7;
    double R22 = R18*R5 - R20*R8;
    double R23 = R11*(R10*R22 - R14 + R19 - R21);
    double R24 = R1*R23 + R13;
    double R25 = R11*(-R10*(R14 - R19 + R21) + R22);
    double R26 = -R0*R23 + R18 - R25*R5;
    double R27 = R20 + R23*R7 + R25*R8;
    double R28 = R3*(R12 + panel_straw0z);
    double R29 = R28*panel_straw0y;
    double R30 = R0*R29;
    double R31 = R28*panel_straw0x;
    double R32 = R31*R7;
    double R33 = R1*R16;
    double R34 = 2*R11*(R30 - R32 - R33);
    double R35 = 2.0*R11*R3*(R6 + R9)*(-R30 + R32 + R33);
    double R36 = ((1.0/2.0)*R24*(R1*R34 + 2*R15 - 2*R2) + (1.0/2.0)*R26*(-R0*R34 + 2*R29 - R35*panel_straw0x) + (1.0/2.0)*R27*(2*R31 + R34*R7 + R35*panel_straw0y))/sqrt(pow(R24, 2) + pow(R26, 2) + pow(R27, 2));
    double result = 16.0*((R25 > 0) ? (
   R36
)
: (
   -R36
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
    double R7 = 1.0*R6;
    double R8 = R7*panel_straw0x;
    double R9 = R0*R8;
    double R10 = R0*a1;
    double R11 = R7*panel_straw0y;
    double R12 = R10*R11;
    double R13 = -R12 - R9;
    double R14 = 1 - pow(R13, 2);
    double R15 = 1.0/R14;
    double R16 = b0 - wire_z;
    double R17 = R1*R16;
    double R18 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R19 = R6*(R18 - R5);
    double R20 = -R19*panel_straw0y - panel_straw0y;
    double R21 = R0*R20;
    double R22 = -R19*panel_straw0x + a0 - panel_straw0x;
    double R23 = R10*R22;
    double R24 = -R11*R22 + R20*R8;
    double R25 = R13*R24 - R17 + R21 - R23;
    double R26 = R15*R25;
    double R27 = R1*R26 + R16;
    double R28 = R17 - R21 + R23;
    double R29 = -R13*R28 + R24;
    double R30 = R15*R29;
    double R31 = -R0*R26 + R20 - R30*R8;
    double R32 = R10*R26 + R11*R30 + R22;
    double R33 = R0*R6;
    double R34 = 2.0*panel_straw0y;
    double R35 = R10*R6;
    double R36 = 2.0*panel_straw0x;
    double R37 = R13*(R33*R34 - R35*R36)/pow(R14, 2);
    double R38 = R25*R37;
    double R39 = 2*R38;
    double R40 = -R18 + R5;
    double R41 = R40*panel_straw0x;
    double R42 = R33*R41;
    double R43 = R40*panel_straw0y;
    double R44 = R35*R43;
    double R45 = R0*R11;
    double R46 = R10*R8;
    double R47 = 1.0*R40/R4;
    double R48 = -R11*R20 + R2*R47 - R22*R8 + R3*R47;
    double R49 = R15*(R13*R48 + R24*(R45 - R46) + R42 + R44);
    double R50 = 2*R49;
    double R51 = 2*R6;
    double R52 = R30*R6;
    double R53 = 2*R0;
    double R54 = R36*R6;
    double R55 = R29*R37;
    double R56 = R15*(R28*(-R45 + R46) + R48 + (R12 + R9)*(-R42 - R44));
    double R57 = R34*R6;
    double R58 = ((1.0/2.0)*R27*(R1*R39 + R1*R50) + (1.0/2.0)*R31*(R34*R52 - R38*R53 + R41*R51 - R49*R53 - R54*R55 - R54*R56) + (1.0/2.0)*R32*(R10*R39 + R10*R50 + R36*R52 - R43*R51 + R55*R57 + R56*R57))/sqrt(pow(R27, 2) + pow(R31, 2) + pow(R32, 2));
    double result = 16.0*((R30 > 0) ? (
   R58
)
: (
   -R58
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

std::vector<double> CosmicTrack_DCA_LocalDeriv_double(double a0, double b0, double a1, double b1, double t0, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
        std::vector<double> result = {CosmicTrack_DCA_Deriv_a0(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_b0(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_a1(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_b1(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
CosmicTrack_DCA_Deriv_t0(a0,b0,a1,b1,t0,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z)};
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

