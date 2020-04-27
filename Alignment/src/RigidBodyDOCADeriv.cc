

# include "Alignment/inc/RigidBodyDOCADeriv.hh"
# include <math.h>
# include <vector>

double CosmicTrack_DCA(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
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
    double R0 = atan2(panel_straw0y, panel_straw0x);
    double result = R0*panel_dy + panel_dx + panel_straw0x + plane_b*(panel_dz + panel_straw0z - plane_z) + plane_dx + plane_g*(-R0*panel_dx + panel_dy + panel_straw0y) + (-panel_straw0z + wire_z)*(R0*panel_a + panel_b + plane_b + plane_g*(-R0*panel_b + panel_a)) + (-sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2)) + sqrt(pow(wire_x, 2) + pow(wire_y, 2)))*(-R0*panel_g + panel_b*plane_b + plane_g*(-R0 - panel_g) + 1);
    return result;
}


double CosmicTrack_DCAalignpos_y(double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
{
    double R0 = atan2(panel_straw0y, panel_straw0x);
    double result = -R0*panel_dx + panel_dy + panel_straw0y + plane_a*(panel_dz + panel_straw0z - plane_z) + plane_dy - plane_g*(R0*panel_dy + panel_dx + panel_straw0x) + (-panel_straw0z + wire_z)*(-R0*panel_b + panel_a + plane_a - plane_g*(R0*panel_a + panel_b)) + (-sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2)) + sqrt(pow(wire_x, 2) + pow(wire_y, 2)))*(-R0 + panel_b*plane_a - panel_g - plane_g*(-R0*panel_g + 1));
    return result;
}


double CosmicTrack_DCAalignpos_z(double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
{
    double R0 = atan2(panel_straw0y, panel_straw0x);
    double result = panel_dz + panel_straw0z - plane_a*(-R0*panel_dx + panel_dy + panel_straw0y) + plane_b*(R0*panel_dy + panel_dx + panel_straw0x) + plane_dz + (-panel_straw0z + wire_z)*(-plane_a*(-R0*panel_b + panel_a) + plane_b*(R0*panel_a + panel_b) + 1) + (-sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2)) + sqrt(pow(wire_x, 2) + pow(wire_y, 2)))*(panel_b - plane_a*(-R0 - panel_g) + plane_b*(-R0*panel_g + 1));
    return result;
}


double CosmicTrack_DCAaligndir_x(double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
{
    double R0 = atan2(panel_straw0y, panel_straw0x);
    double result = 1.0*R0 - 1.0*panel_a*plane_b + 1.0*panel_g + 1.0*plane_g*(-R0*panel_g + 1);
    return result;
}


double CosmicTrack_DCAaligndir_y(double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
{
    double R0 = atan2(panel_straw0y, panel_straw0x);
    double result = -1.0*R0*panel_g - 1.0*panel_a*plane_a - 1.0*plane_g*(R0 + panel_g) + 1.0;
    return result;
}


double CosmicTrack_DCAaligndir_z(double plane_dx, double plane_dy, double plane_dz, double plane_a, double plane_b, double plane_g, double panel_dx, double panel_dy, double panel_dz, double panel_a, double panel_b, double panel_g, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z)
{
    double R0 = atan2(panel_straw0y, panel_straw0x);
    double result = -1.0*panel_a - 1.0*plane_a*(-R0*panel_g + 1) + 1.0*plane_b*(R0 + panel_g);
    return result;
}


double CosmicTrack_DCA_Deriv_a0(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = 1.0*R0;
    double R3 = atan2(panel_straw0y, panel_straw0x);
    double R4 = R2*R3*a1;
    double R5 = -R2 + R4;
    double R6 = 1.0/(1 - pow(R5, 2));
    double R7 = b0 - wire_z;
    double R8 = R1*R7;
    double R9 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R10 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R11 = -R10 + R9 + a0 - panel_straw0x;
    double R12 = R0*a1;
    double R13 = R11*R12;
    double R14 = R3*(R10 - R9);
    double R15 = R14 - panel_straw0y;
    double R16 = R0*R15;
    double R17 = 1.0*R3;
    double R18 = R11*R17 + 1.0*R14 - 1.0*panel_straw0y;
    double R19 = R6*(-R13 + R16 + R18*R5 - R8);
    double R20 = R1*R19 + R7;
    double R21 = R6*(R18 - R5*(R13 - R16 + R8));
    double R22 = R11 + R12*R19 - R17*R21;
    double R23 = -R0*R19 + R15 - 1.0*R21;
    double R24 = R6*(-R12 + R17*R5);
    double R25 = 2.0*R6*(R12*(R2 - R4) + R17);
    double R26 = 2*R24;
    double R27 = (R1*R20*R24 + (1.0/2.0)*R22*(R12*R26 - R25*R3 + 2) + (1.0/2.0)*R23*(-R0*R26 - R25))/sqrt(pow(R20, 2) + pow(R22, 2) + pow(R23, 2));
    double result = ((R21 > 0) ? (
   R27
)
: (
   -R27
));
    return result;
}


double CosmicTrack_DCA_Deriv_b0(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = 1.0*R2;
    double R5 = atan2(panel_straw0y, panel_straw0x);
    double R6 = R4*R5*a1;
    double R7 = -R4 + R6;
    double R8 = 1.0/(1 - pow(R7, 2));
    double R9 = b0 - wire_z;
    double R10 = R3*R9;
    double R11 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R12 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R13 = R11 - R12 + a0 - panel_straw0x;
    double R14 = R2*a1;
    double R15 = R13*R14;
    double R16 = R5*(-R11 + R12);
    double R17 = R16 - panel_straw0y;
    double R18 = R17*R2;
    double R19 = 1.0*R5;
    double R20 = R13*R19 + 1.0*R16 - 1.0*panel_straw0y;
    double R21 = R8*(-R10 - R15 + R18 + R20*R7);
    double R22 = R21*R3 + R9;
    double R23 = R8*(R20 - R7*(R10 + R15 - R18));
    double R24 = R13 + R14*R21 - R19*R23;
    double R25 = R17 - R2*R21 - 1.0*R23;
    double R26 = 2*R8/R1;
    double R27 = R26*b1;
    double R28 = 2.0*R3*R8*(R4 - R6);
    double R29 = ((1.0/2.0)*R22*(-R0*R26 + 2) + (1.0/2.0)*R24*(-R27*a1 - R28*R5) + (1.0/2.0)*R25*(R27 - R28))/sqrt(pow(R22, 2) + pow(R24, 2) + pow(R25, 2));
    double result = ((R23 > 0) ? (
   R29
)
: (
   -R29
));
    return result;
}


double CosmicTrack_DCA_Deriv_a1(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(a1, 2);
    double R1 = R0 + pow(b1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = 1.0*R2;
    double R5 = atan2(panel_straw0y, panel_straw0x);
    double R6 = R4*R5;
    double R7 = R6*a1;
    double R8 = -R4 + R7;
    double R9 = 1 - pow(R8, 2);
    double R10 = 1.0/R9;
    double R11 = b0 - wire_z;
    double R12 = R11*R3;
    double R13 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R14 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R15 = R13 - R14 + a0 - panel_straw0x;
    double R16 = R15*R2;
    double R17 = R16*a1;
    double R18 = R5*(-R13 + R14);
    double R19 = R18 - panel_straw0y;
    double R20 = R19*R2;
    double R21 = 1.0*R5;
    double R22 = R15*R21 + 1.0*R18 - 1.0*panel_straw0y;
    double R23 = -R12 - R17 + R20 + R22*R8;
    double R24 = R10*R23;
    double R25 = R11 + R24*R3;
    double R26 = R12 + R17 - R20;
    double R27 = R22 - R26*R8;
    double R28 = R10*R27;
    double R29 = R2*R24;
    double R30 = R15 - R21*R28 + R29*a1;
    double R31 = R19 - 1.0*R28 - R29;
    double R32 = pow(R1, -3.0/2.0);
    double R33 = R32*a1;
    double R34 = R33*b1;
    double R35 = 2*R24;
    double R36 = R11*R34;
    double R37 = R0*R32;
    double R38 = R15*R37;
    double R39 = R19*R33;
    double R40 = 1.0*R33;
    double R41 = R21*R37;
    double R42 = R10*(-R16 + R22*(R40 - R41 + R6) + R36 + R38 - R39);
    double R43 = 2*R3;
    double R44 = 2.0*R5;
    double R45 = R8*(R2*R44 + 2.0*R33 - R37*R44)/pow(R9, 2);
    double R46 = R23*R45;
    double R47 = 2*R2;
    double R48 = R42*R47;
    double R49 = R27*R45;
    double R50 = R46*R47;
    double R51 = R10*(R26*(-R40 + R41 - R6) + (R4 - R7)*(R16 - R36 - R38 + R39));
    double R52 = ((1.0/2.0)*R25*(-R34*R35 + R42*R43 + R43*R46) + (1.0/2.0)*R30*(2*R29 - R35*R37 - R44*R49 - R44*R51 + R48*a1 + R50*a1) + (1.0/2.0)*R31*(R33*R35 - R48 - 2.0*R49 - R50 - 2.0*R51))/sqrt(pow(R25, 2) + pow(R30, 2) + pow(R31, 2));
    double result = ((R28 > 0) ? (
   R52
)
: (
   -R52
));
    return result;
}


double CosmicTrack_DCA_Deriv_b1(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = 1.0*R2;
    double R4 = atan2(panel_straw0y, panel_straw0x);
    double R5 = R4*a1;
    double R6 = R3*R5;
    double R7 = -R3 + R6;
    double R8 = 1 - pow(R7, 2);
    double R9 = 1.0/R8;
    double R10 = b0 - wire_z;
    double R11 = R10*R2;
    double R12 = R11*b1;
    double R13 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R14 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R15 = R13 - R14 + a0 - panel_straw0x;
    double R16 = R2*a1;
    double R17 = R15*R16;
    double R18 = R4*(-R13 + R14);
    double R19 = R18 - panel_straw0y;
    double R20 = R19*R2;
    double R21 = 1.0*R4;
    double R22 = R15*R21 + 1.0*R18 - 1.0*panel_straw0y;
    double R23 = -R12 - R17 + R20 + R22*R7;
    double R24 = R23*R9;
    double R25 = R2*R24;
    double R26 = R10 + R25*b1;
    double R27 = R12 + R17 - R20;
    double R28 = R22 - R27*R7;
    double R29 = R28*R9;
    double R30 = R15 + R16*R24 - R21*R29;
    double R31 = R19 - R25 - 1.0*R29;
    double R32 = pow(R1, -3.0/2.0);
    double R33 = R0*R32;
    double R34 = 2*R24;
    double R35 = 2*R2;
    double R36 = R10*R33;
    double R37 = R32*b1;
    double R38 = R37*a1;
    double R39 = R15*R38;
    double R40 = R19*R37;
    double R41 = 1.0*R37;
    double R42 = R21*R38;
    double R43 = R9*(-R11 + R22*(R41 - R42) + R36 + R39 - R40);
    double R44 = R35*R43;
    double R45 = 2.0*R37;
    double R46 = R7*(-R45*R5 + R45)/pow(R8, 2);
    double R47 = R23*R46;
    double R48 = R35*R47;
    double R49 = 2.0*R28*R46;
    double R50 = 2.0*R9*(R27*(-R41 + R42) + (R3 - R6)*(R11 - R36 - R39 + R40));
    double R51 = 2*R16;
    double R52 = ((1.0/2.0)*R26*(2*R25 - R33*R34 + R44*b1 + R48*b1) + (1.0/2.0)*R30*(-R34*R38 - R4*R49 - R4*R50 + R43*R51 + R47*R51) + (1.0/2.0)*R31*(R34*R37 - R44 - R48 - R49 - R50))/sqrt(pow(R26, 2) + pow(R30, 2) + pow(R31, 2));
    double result = ((R29 > 0) ? (
   R52
)
: (
   -R52
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dx(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = 1.0*R0;
    double R3 = atan2(panel_straw0y, panel_straw0x);
    double R4 = R2*R3*a1;
    double R5 = -R2 + R4;
    double R6 = 1.0/(1 - pow(R5, 2));
    double R7 = b0 - wire_z;
    double R8 = R1*R7;
    double R9 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R10 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R11 = -R10 + R9 + a0 - panel_straw0x;
    double R12 = R0*a1;
    double R13 = R11*R12;
    double R14 = R3*(R10 - R9);
    double R15 = R14 - panel_straw0y;
    double R16 = R0*R15;
    double R17 = 1.0*R3;
    double R18 = R11*R17 + 1.0*R14 - 1.0*panel_straw0y;
    double R19 = R6*(-R13 + R16 + R18*R5 - R8);
    double R20 = R1*R19 + R7;
    double R21 = R6*(R18 - R5*(R13 - R16 + R8));
    double R22 = R11 + R12*R19 - R17*R21;
    double R23 = -R0*R19 + R15 - 1.0*R21;
    double R24 = R6*(R12 - R17*R5);
    double R25 = 2.0*R6*(-R12*(R2 - R4) - R17);
    double R26 = 2*R24;
    double R27 = (R1*R20*R24 + (1.0/2.0)*R22*(R12*R26 - R25*R3 - 2) + (1.0/2.0)*R23*(-R0*R26 - R25))/sqrt(pow(R20, 2) + pow(R22, 2) + pow(R23, 2));
    double result = ((R21 > 0) ? (
   R27
)
: (
   -R27
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dy(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(a1, 2);
    double R1 = R0 + pow(b1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = 1.0*R2;
    double R5 = atan2(panel_straw0y, panel_straw0x);
    double R6 = R5*a1;
    double R7 = R4*R6;
    double R8 = -R4 + R7;
    double R9 = 1.0/(1 - pow(R8, 2));
    double R10 = b0 - wire_z;
    double R11 = R10*R3;
    double R12 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R13 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R14 = R12 - R13 + a0 - panel_straw0x;
    double R15 = R2*a1;
    double R16 = R14*R15;
    double R17 = R5*(-R12 + R13);
    double R18 = R17 - panel_straw0y;
    double R19 = R18*R2;
    double R20 = 1.0*R5;
    double R21 = R14*R20 + 1.0*R17 - 1.0*panel_straw0y;
    double R22 = R9*(-R11 - R16 + R19 + R21*R8);
    double R23 = R10 + R22*R3;
    double R24 = R9*(R21 - R8*(R11 + R16 - R19));
    double R25 = R14 + R15*R22 - R20*R24;
    double R26 = R18 - R2*R22 - 1.0*R24;
    double R27 = 1.0/R1;
    double R28 = 2.0*R9;
    double R29 = R27*R28;
    double R30 = R28*(R2*(R4 - R7) - 1.0);
    double R31 = (-R20*R23*R27*R9*a1*b1 + (1.0/2.0)*R25*(-R0*R29*R5 - R30*R5) + (1.0/2.0)*R26*(R29*R6 - R30 - 2))/sqrt(pow(R23, 2) + pow(R25, 2) + pow(R26, 2));
    double result = ((R24 > 0) ? (
   R31
)
: (
   -R31
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_dz(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = 1.0*R2;
    double R5 = atan2(panel_straw0y, panel_straw0x);
    double R6 = R4*R5*a1;
    double R7 = -R4 + R6;
    double R8 = 1.0/(1 - pow(R7, 2));
    double R9 = b0 - wire_z;
    double R10 = R3*R9;
    double R11 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R12 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R13 = R11 - R12 + a0 - panel_straw0x;
    double R14 = R2*a1;
    double R15 = R13*R14;
    double R16 = R5*(-R11 + R12);
    double R17 = R16 - panel_straw0y;
    double R18 = R17*R2;
    double R19 = 1.0*R5;
    double R20 = R13*R19 + 1.0*R16 - 1.0*panel_straw0y;
    double R21 = R8*(-R10 - R15 + R18 + R20*R7);
    double R22 = R21*R3 + R9;
    double R23 = R8*(R20 - R7*(R10 + R15 - R18));
    double R24 = R13 + R14*R21 - R19*R23;
    double R25 = R17 - R2*R21 - 1.0*R23;
    double R26 = 2*R8/R1;
    double R27 = R26*b1;
    double R28 = 2.0*R3*R8*(R4 - R6);
    double R29 = ((1.0/2.0)*R22*(R0*R26 - 2) + (1.0/2.0)*R24*(R27*a1 + R28*R5) + (1.0/2.0)*R25*(-R27 + R28))/sqrt(pow(R22, 2) + pow(R24, 2) + pow(R25, 2));
    double result = ((R23 > 0) ? (
   R29
)
: (
   -R29
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_a(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = 1.0*R2;
    double R5 = atan2(panel_straw0y, panel_straw0x);
    double R6 = R4*R5*a1;
    double R7 = -R4 + R6;
    double R8 = 1 - pow(R7, 2);
    double R9 = 1.0/R8;
    double R10 = -wire_z;
    double R11 = R10 + b0;
    double R12 = R11*R3;
    double R13 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R14 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R15 = R13 - R14;
    double R16 = R15 + a0 - panel_straw0x;
    double R17 = R2*a1;
    double R18 = R16*R17;
    double R19 = R5*(-R13 + R14);
    double R20 = R19 - panel_straw0y;
    double R21 = R2*R20;
    double R22 = 1.0*R5;
    double R23 = R16*R22 + 1.0*R19 - 1.0*panel_straw0y;
    double R24 = -R12 - R18 + R21 + R23*R7;
    double R25 = R24*R9;
    double R26 = R11 + R25*R3;
    double R27 = R12 + R18 - R21;
    double R28 = R23 - R27*R7;
    double R29 = R28*R9;
    double R30 = R16 + R17*R25 - R22*R29;
    double R31 = -R2*R25 + R20 - 1.0*R29;
    double R32 = R15*R5;
    double R33 = R2*(R10 + plane_z);
    double R34 = -1.0*b0 + 1.0*plane_z;
    double R35 = R3*(R32 + panel_straw0y);
    double R36 = R4*b1;
    double R37 = 2*R9*(-R23*R36 + R33 + R34*R7 - R35);
    double R38 = 4.0*R7/pow(R8, 2);
    double R39 = R24*R38/R1;
    double R40 = 2.0*R9*(R27*R36 + R34 + (-R33 + R35)*(R4 - R6));
    double R41 = R39*b1;
    double R42 = R28*R3*R38;
    double R43 = ((1.0/2.0)*R26*(-R0*R39 + 2.0*R29 + R3*R37 + 2*R32 + 2*panel_straw0y) + (1.0/2.0)*R30*(R17*R37 - R40*R5 - R41*a1 + R42*R5) + (1.0/2.0)*R31*(-R2*R37 - R40 + R41 + R42 + 2*plane_z - 2*wire_z))/sqrt(pow(R26, 2) + pow(R30, 2) + pow(R31, 2));
    double result = ((R29 > 0) ? (
   R43
)
: (
   -R43
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_b(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = 1.0*R2;
    double R5 = atan2(panel_straw0y, panel_straw0x);
    double R6 = R4*R5;
    double R7 = R6*a1;
    double R8 = -R4 + R7;
    double R9 = 1 - pow(R8, 2);
    double R10 = 1.0/R9;
    double R11 = -wire_z;
    double R12 = R11 + b0;
    double R13 = R12*R3;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R16 = R14 - R15 - panel_straw0x;
    double R17 = R16 + a0;
    double R18 = R2*a1;
    double R19 = R17*R18;
    double R20 = R5*(-R14 + R15);
    double R21 = R20 - panel_straw0y;
    double R22 = R2*R21;
    double R23 = 1.0*R5;
    double R24 = R17*R23 + 1.0*R20 - 1.0*panel_straw0y;
    double R25 = -R13 - R19 + R22 + R24*R8;
    double R26 = R10*R25;
    double R27 = R12 + R26*R3;
    double R28 = R13 + R19 - R22;
    double R29 = R24 - R28*R8;
    double R30 = R10*R29;
    double R31 = R17 + R18*R26 - R23*R30;
    double R32 = -R2*R26 + R21 - 1.0*R30;
    double R33 = R11 + plane_z;
    double R34 = R18*R33;
    double R35 = R16*R3;
    double R36 = R12*R23 + R23*R33;
    double R37 = R6*b1;
    double R38 = 2*R10*(R24*R37 - R34 - R35 + R36*R8);
    double R39 = 4.0*R8/pow(R9, 2);
    double R40 = R39*R5;
    double R41 = R25*R40/R1;
    double R42 = 2.0*R10*(-R28*R37 + R36 + (R34 + R35)*(R4 - R7));
    double R43 = R41*b1;
    double R44 = R29*R3;
    double R45 = ((1.0/2.0)*R27*(R0*R41 + 2*R14 - 2*R15 + R3*R38 - 2.0*R30*R5 - 2*panel_straw0x) + (1.0/2.0)*R31*(R18*R38 - R39*R44*pow(R5, 2) - R42*R5 + R43*a1 + 2*plane_z - 2*wire_z) + (1.0/2.0)*R32*(-R2*R38 - R40*R44 - R42 - R43))/sqrt(pow(R27, 2) + pow(R31, 2) + pow(R32, 2));
    double result = ((R30 > 0) ? (
   R45
)
: (
   -R45
));
    return result;
}


double CosmicTrack_DCA_Deriv_plane_g(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = 1.0*R0;
    double R3 = atan2(panel_straw0y, panel_straw0x);
    double R4 = R2*a1;
    double R5 = R3*R4;
    double R6 = -R2 + R5;
    double R7 = 1 - pow(R6, 2);
    double R8 = 1.0/R7;
    double R9 = b0 - wire_z;
    double R10 = R1*R9;
    double R11 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R12 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R13 = R11 - R12;
    double R14 = R13 + a0 - panel_straw0x;
    double R15 = R0*a1;
    double R16 = R14*R15;
    double R17 = -panel_straw0y;
    double R18 = -R11 + R12;
    double R19 = R18*R3;
    double R20 = R17 + R19;
    double R21 = R0*R20;
    double R22 = 1.0*R3;
    double R23 = R14*R22 + 1.0*R19 - 1.0*panel_straw0y;
    double R24 = -R10 - R16 + R21 + R23*R6;
    double R25 = R24*R8;
    double R26 = R1*R25 + R9;
    double R27 = R10 + R16 - R21;
    double R28 = R23 - R27*R6;
    double R29 = R28*R8;
    double R30 = R14 + R15*R25 - R22*R29;
    double R31 = -R0*R25 + R20 - 1.0*R29;
    double R32 = 2.0*R3;
    double R33 = R6*(R0*R32 + 2.0*R15)/pow(R7, 2);
    double R34 = R24*R33;
    double R35 = 2*R34;
    double R36 = R0*(R18 + panel_straw0x);
    double R37 = R13*R3;
    double R38 = R17 - R37;
    double R39 = R15*R38;
    double R40 = R2*R3;
    double R41 = -R20*R22 + R22*R38 + 1.0*a0;
    double R42 = R8*(R23*(R4 + R40) + R36 - R39 + R41*R6);
    double R43 = 2*R42;
    double R44 = R28*R33;
    double R45 = 2*R0;
    double R46 = R8*(R27*(-R4 - R40) + R41 + (R2 - R5)*(-R36 + R39));
    double R47 = ((1.0/2.0)*R26*(R1*R35 + R1*R43) + (1.0/2.0)*R30*(R15*R35 + R15*R43 - 2.0*R29 - R32*R44 - R32*R46 - 2*R37 - 2*panel_straw0y) + (1.0/2.0)*R31*(-2*R11 + 2*R12 + R29*R32 - R34*R45 - R42*R45 - 2.0*R44 - 2.0*R46 + 2*panel_straw0x))/sqrt(pow(R26, 2) + pow(R30, 2) + pow(R31, 2));
    double result = ((R29 > 0) ? (
   R47
)
: (
   -R47
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dx(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = 1.0*R0;
    double R3 = atan2(panel_straw0y, panel_straw0x);
    double R4 = R2*R3*a1;
    double R5 = -R2 + R4;
    double R6 = 1.0/(1 - pow(R5, 2));
    double R7 = b0 - wire_z;
    double R8 = R1*R7;
    double R9 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R10 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R11 = -R10 + R9 + a0 - panel_straw0x;
    double R12 = R0*a1;
    double R13 = R11*R12;
    double R14 = R3*(R10 - R9);
    double R15 = R14 - panel_straw0y;
    double R16 = R0*R15;
    double R17 = 1.0*R3;
    double R18 = R11*R17 + 1.0*R14 - 1.0*panel_straw0y;
    double R19 = R6*(-R13 + R16 + R18*R5 - R8);
    double R20 = R1*R19 + R7;
    double R21 = R6*(R18 - R5*(R13 - R16 + R8));
    double R22 = R11 + R12*R19 - R17*R21;
    double R23 = -R0*R19 + R15 - 1.0*R21;
    double R24 = R0*R3;
    double R25 = R6*(R12 + R24);
    double R26 = 2*R25;
    double R27 = 2.0*R6*(-R12 - R24)*(R2 - R4);
    double R28 = (R1*R20*R25 + (1.0/2.0)*R22*(R12*R26 - R27*R3 - 2) + (1.0/2.0)*R23*(-R0*R26 - R27 + 2*R3))/sqrt(pow(R20, 2) + pow(R22, 2) + pow(R23, 2));
    double result = ((R21 > 0) ? (
   R28
)
: (
   -R28
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dy(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = 1.0*R0;
    double R3 = atan2(panel_straw0y, panel_straw0x);
    double R4 = R2*R3*a1;
    double R5 = -R2 + R4;
    double R6 = 1.0/(1 - pow(R5, 2));
    double R7 = b0 - wire_z;
    double R8 = R1*R7;
    double R9 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R10 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R11 = -R10 + R9 + a0 - panel_straw0x;
    double R12 = R0*a1;
    double R13 = R11*R12;
    double R14 = R3*(R10 - R9);
    double R15 = R14 - panel_straw0y;
    double R16 = R0*R15;
    double R17 = 1.0*R3;
    double R18 = R11*R17 + 1.0*R14 - 1.0*panel_straw0y;
    double R19 = R6*(-R13 + R16 + R18*R5 - R8);
    double R20 = R1*R19 + R7;
    double R21 = R6*(R18 - R5*(R13 - R16 + R8));
    double R22 = R11 + R12*R19 - R17*R21;
    double R23 = -R0*R19 + R15 - 1.0*R21;
    double R24 = R12*R3;
    double R25 = -1.0*pow(R3, 2) - 1.0;
    double R26 = R6*(-R0 + R24 + R25*R5);
    double R27 = 2.0*R6*(R25 + (R0 - R24)*(R2 - R4));
    double R28 = 2*R26;
    double R29 = (R1*R20*R26 + (1.0/2.0)*R22*(R12*R28 - R27*R3 - 2*R3) + (1.0/2.0)*R23*(-R0*R28 - R27 - 2))/sqrt(pow(R20, 2) + pow(R22, 2) + pow(R23, 2));
    double result = ((R21 > 0) ? (
   R29
)
: (
   -R29
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_dz(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = 1.0*R2;
    double R5 = atan2(panel_straw0y, panel_straw0x);
    double R6 = R4*R5*a1;
    double R7 = -R4 + R6;
    double R8 = 1.0/(1 - pow(R7, 2));
    double R9 = b0 - wire_z;
    double R10 = R3*R9;
    double R11 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R12 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R13 = R11 - R12 + a0 - panel_straw0x;
    double R14 = R2*a1;
    double R15 = R13*R14;
    double R16 = R5*(-R11 + R12);
    double R17 = R16 - panel_straw0y;
    double R18 = R17*R2;
    double R19 = 1.0*R5;
    double R20 = R13*R19 + 1.0*R16 - 1.0*panel_straw0y;
    double R21 = R8*(-R10 - R15 + R18 + R20*R7);
    double R22 = R21*R3 + R9;
    double R23 = R8*(R20 - R7*(R10 + R15 - R18));
    double R24 = R13 + R14*R21 - R19*R23;
    double R25 = R17 - R2*R21 - 1.0*R23;
    double R26 = 2*R8/R1;
    double R27 = R26*b1;
    double R28 = 2.0*R3*R8*(R4 - R6);
    double R29 = ((1.0/2.0)*R22*(R0*R26 - 2) + (1.0/2.0)*R24*(R27*a1 + R28*R5) + (1.0/2.0)*R25*(-R27 + R28))/sqrt(pow(R22, 2) + pow(R24, 2) + pow(R25, 2));
    double result = ((R23 > 0) ? (
   R29
)
: (
   -R29
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_a(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0/2.0);
    double R3 = R2*b1;
    double R4 = 1.0*R2;
    double R5 = atan2(panel_straw0y, panel_straw0x);
    double R6 = R5*a1;
    double R7 = R4*R6;
    double R8 = -R4 + R7;
    double R9 = 1 - pow(R8, 2);
    double R10 = 1.0/R9;
    double R11 = -wire_z;
    double R12 = R11 + b0;
    double R13 = R12*R3;
    double R14 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R15 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R16 = R14 - R15 + a0 - panel_straw0x;
    double R17 = R2*a1;
    double R18 = R16*R17;
    double R19 = R5*(-R14 + R15);
    double R20 = R19 - panel_straw0y;
    double R21 = R2*R20;
    double R22 = 1.0*R5;
    double R23 = R16*R22 + 1.0*R19 - 1.0*panel_straw0y;
    double R24 = -R13 - R18 + R21 + R23*R8;
    double R25 = R10*R24;
    double R26 = R12 + R25*R3;
    double R27 = R13 + R18 - R21;
    double R28 = R23 - R27*R8;
    double R29 = R10*R28;
    double R30 = R16 + R17*R25 - R22*R29;
    double R31 = -R2*R25 + R20 - 1.0*R29;
    double R32 = R11 + panel_straw0z;
    double R33 = R2*R32;
    double R34 = R33*R6;
    double R35 = 1.0*R32*pow(R5, 2) - 1.0*b0 + 1.0*panel_straw0z;
    double R36 = R4*b1;
    double R37 = 2*R10*(-R23*R36 + R33 - R34 + R35*R8);
    double R38 = 4.0*R8/pow(R9, 2);
    double R39 = R24*R38/R1;
    double R40 = 2.0*R10*(R27*R36 + R35 + (-R33 + R34)*(R4 - R7));
    double R41 = R39*b1;
    double R42 = R28*R3*R38;
    double R43 = ((1.0/2.0)*R26*(-R0*R39 + 2.0*R29 + R3*R37) + (1.0/2.0)*R30*(R17*R37 + 2*R32*R5 - R40*R5 - R41*a1 + R42*R5) + (1.0/2.0)*R31*(-R2*R37 - R40 + R41 + R42 + 2*panel_straw0z - 2*wire_z))/sqrt(pow(R26, 2) + pow(R30, 2) + pow(R31, 2));
    double result = ((R29 > 0) ? (
   R43
)
: (
   -R43
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_b(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = 1.0*R0;
    double R3 = atan2(panel_straw0y, panel_straw0x);
    double R4 = R2*R3*a1;
    double R5 = -R2 + R4;
    double R6 = 1.0/(1 - pow(R5, 2));
    double R7 = -wire_z;
    double R8 = R7 + b0;
    double R9 = R1*R8;
    double R10 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R11 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R12 = R10 - R11;
    double R13 = R12 + a0 - panel_straw0x;
    double R14 = R0*a1;
    double R15 = R13*R14;
    double R16 = R3*(-R10 + R11);
    double R17 = R16 - panel_straw0y;
    double R18 = R0*R17;
    double R19 = 1.0*R3;
    double R20 = R13*R19 + 1.0*R16 - 1.0*panel_straw0y;
    double R21 = R6*(-R15 + R18 + R20*R5 - R9);
    double R22 = R1*R21 + R8;
    double R23 = R6*(R20 - R5*(R15 - R18 + R9));
    double R24 = R13 + R14*R21 - R19*R23;
    double R25 = -R0*R21 + R17 - 1.0*R23;
    double R26 = R7 + panel_straw0z;
    double R27 = R14*R26;
    double R28 = R26*R3;
    double R29 = R0*R28;
    double R30 = R1*R12;
    double R31 = 2*R6*(-R27 - R29 - R30);
    double R32 = 2.0*R6*(R2 - R4)*(R27 + R29 + R30);
    double R33 = ((1.0/2.0)*R22*(R1*R31 + 2*R10 - 2*R11) + (1.0/2.0)*R24*(R14*R31 - R3*R32 + 2*panel_straw0z - 2*wire_z) + (1.0/2.0)*R25*(-R0*R31 - 2*R28 - R32))/sqrt(pow(R22, 2) + pow(R24, 2) + pow(R25, 2));
    double result = ((R23 > 0) ? (
   R33
)
: (
   -R33
));
    return result;
}


double CosmicTrack_DCA_Deriv_panel_g(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0/2.0);
    double R1 = R0*b1;
    double R2 = 1.0*R0;
    double R3 = atan2(panel_straw0y, panel_straw0x);
    double R4 = R2*a1;
    double R5 = R3*R4;
    double R6 = -R2 + R5;
    double R7 = 1 - pow(R6, 2);
    double R8 = 1.0/R7;
    double R9 = b0 - wire_z;
    double R10 = R1*R9;
    double R11 = sqrt(pow(panel_straw0x, 2) + pow(panel_straw0y, 2));
    double R12 = sqrt(pow(wire_x, 2) + pow(wire_y, 2));
    double R13 = R11 - R12;
    double R14 = R13 + a0 - panel_straw0x;
    double R15 = R0*a1;
    double R16 = R14*R15;
    double R17 = R3*(-R11 + R12);
    double R18 = R17 - panel_straw0y;
    double R19 = R0*R18;
    double R20 = 1.0*R3;
    double R21 = R14*R20 + 1.0*R17 - 1.0*panel_straw0y;
    double R22 = -R10 - R16 + R19 + R21*R6;
    double R23 = R22*R8;
    double R24 = R1*R23 + R9;
    double R25 = R10 + R16 - R19;
    double R26 = R21 - R25*R6;
    double R27 = R26*R8;
    double R28 = R14 + R15*R23 - R20*R27;
    double R29 = -R0*R23 + R18 - 1.0*R27;
    double R30 = 2.0*R3;
    double R31 = R6*(R0*R30 + 2.0*R15)/pow(R7, 2);
    double R32 = R22*R31;
    double R33 = 2*R32;
    double R34 = R0*R13;
    double R35 = R3*R34*a1;
    double R36 = R2*R3;
    double R37 = -1.0*R13*pow(R3, 2) - R18*R20 + 1.0*a0 - 1.0*panel_straw0x;
    double R38 = R8*(R21*(R36 + R4) - R34 + R35 + R37*R6);
    double R39 = 2*R38;
    double R40 = R26*R31;
    double R41 = 2*R0;
    double R42 = R8*(R25*(-R36 - R4) + R37 + (R2 - R5)*(R34 - R35));
    double R43 = ((1.0/2.0)*R24*(R1*R33 + R1*R39) + (1.0/2.0)*R28*(-2*R13*R3 + R15*R33 + R15*R39 - 2.0*R27 - R30*R40 - R30*R42) + (1.0/2.0)*R29*(-2*R11 + 2*R12 + R27*R30 - R32*R41 - R38*R41 - 2.0*R40 - 2.0*R42))/sqrt(pow(R24, 2) + pow(R28, 2) + pow(R29, 2));
    double result = ((R27 > 0) ? (
   R43
)
: (
   -R43
));
    return result;
}


std::vector<float> CosmicTrack_DCA_LocalDeriv(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
        std::vector<float> result = {(float)CosmicTrack_DCA_Deriv_a0(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_b0(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_a1(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_b1(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z)};
return result;
}

std::vector<float> CosmicTrack_DCA_GlobalDeriv(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z)
{
        std::vector<float> result = {(float)CosmicTrack_DCA_Deriv_plane_dx(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_plane_dy(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_plane_dz(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_plane_a(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_plane_b(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_plane_g(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_panel_dx(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_panel_dy(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_panel_dz(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_panel_a(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_panel_b(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z),
(float)CosmicTrack_DCA_Deriv_panel_g(a0,b0,a1,b1,wire_x,wire_y,wire_z,wdir_x,wdir_y,wdir_z,plane_x,plane_y,plane_z,panel_straw0x,panel_straw0y,panel_straw0z)};
return result;
}

