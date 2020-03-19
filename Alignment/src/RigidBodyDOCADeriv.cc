

#include "CosmicAlignment/inc/RigidBodyDOCADeriv.hh"
#include <math.h>
#include <vector>

double CosmicTrack_DCA(double a0, double b0, double a1, double b1, double wire_x, double wire_y,
                       double wire_z, double wdir_x, double wdir_y, double wdir_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R1 = R0 * a1;
    double R2 = R0 * b1;
    double R3 = R0 * wdir_z + R1 * wdir_x + R2 * wdir_y;
    double R4 = 1.0 / (1 - pow(R3, 2));
    double R5 = R0 * wire_z;
    double R6 = a0 - wire_x;
    double R7 = R1 * R6;
    double R8 = b0 - wire_y;
    double R9 = R2 * R8;
    double R10 = R6 * wdir_x + R8 * wdir_y - wdir_z * wire_z;
    double R11 = R4 * (R10 - R3 * (-R5 + R7 + R9));
    double R12 = R4 * (R10 * R3 + R5 - R7 - R9);
    double R13 = sqrt(pow(R0 * R12 - R11 * wdir_z - wire_z, 2) +
                      pow(R1 * R12 - R11 * wdir_x + R6, 2) + pow(-R11 * wdir_y + R12 * R2 + R8, 2));
    double result = ((R11 > 0) ? (R13) : (-R13));
    return result;
}

double CosmicTrack_DCA_Deriv_a0(double a0, double b0, double a1, double b1, double wire_x,
                                double wire_y, double wire_z, double wdir_x, double wdir_y,
                                double wdir_z, double plane_x, double plane_y, double plane_z,
                                double panel_x, double panel_y, double panel_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R1 = R0 * wdir_z;
    double R2 = R0 * a1;
    double R3 = R2 * wdir_x;
    double R4 = R0 * b1;
    double R5 = R4 * wdir_y;
    double R6 = R1 + R3 + R5;
    double R7 = 1.0 / (1 - pow(R6, 2));
    double R8 = R0 * wire_z;
    double R9 = a0 - wire_x;
    double R10 = R2 * R9;
    double R11 = b0 - wire_y;
    double R12 = R11 * R4;
    double R13 = R11 * wdir_y + R9 * wdir_x - wdir_z * wire_z;
    double R14 = R7 * (R13 - R6 * (R10 + R12 - R8));
    double R15 = R7 * (-R10 - R12 + R13 * R6 + R8);
    double R16 = R0 * R15 - R14 * wdir_z - wire_z;
    double R17 = -R14 * wdir_x + R15 * R2 + R9;
    double R18 = R11 - R14 * wdir_y + R15 * R4;
    double R19 = 2 * R7;
    double R20 = R19 * (R2 * (-R1 - R3 - R5) + wdir_x);
    double R21 = R19 * (-R2 + R6 * wdir_x);
    double R22 = ((1.0 / 2.0) * R16 * (R0 * R21 - R20 * wdir_z) +
                  (1.0 / 2.0) * R17 * (R2 * R21 - R20 * wdir_x + 2) +
                  (1.0 / 2.0) * R18 * (-R20 * wdir_y + R21 * R4)) /
                 sqrt(pow(R16, 2) + pow(R17, 2) + pow(R18, 2));
    double result = ((R14 > 0) ? (R22) : (-R22));
    return result;
}

double CosmicTrack_DCA_Deriv_b0(double a0, double b0, double a1, double b1, double wire_x,
                                double wire_y, double wire_z, double wdir_x, double wdir_y,
                                double wdir_z, double plane_x, double plane_y, double plane_z,
                                double panel_x, double panel_y, double panel_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R1 = R0 * wdir_z;
    double R2 = R0 * a1;
    double R3 = R2 * wdir_x;
    double R4 = R0 * b1;
    double R5 = R4 * wdir_y;
    double R6 = R1 + R3 + R5;
    double R7 = 1.0 / (1 - pow(R6, 2));
    double R8 = R0 * wire_z;
    double R9 = a0 - wire_x;
    double R10 = R2 * R9;
    double R11 = b0 - wire_y;
    double R12 = R11 * R4;
    double R13 = R11 * wdir_y + R9 * wdir_x - wdir_z * wire_z;
    double R14 = R7 * (R13 - R6 * (R10 + R12 - R8));
    double R15 = R7 * (-R10 - R12 + R13 * R6 + R8);
    double R16 = R0 * R15 - R14 * wdir_z - wire_z;
    double R17 = -R14 * wdir_x + R15 * R2 + R9;
    double R18 = R11 - R14 * wdir_y + R15 * R4;
    double R19 = 2 * R7;
    double R20 = R19 * (R4 * (-R1 - R3 - R5) + wdir_y);
    double R21 = R19 * (-R4 + R6 * wdir_y);
    double R22 = ((1.0 / 2.0) * R16 * (R0 * R21 - R20 * wdir_z) +
                  (1.0 / 2.0) * R17 * (R2 * R21 - R20 * wdir_x) +
                  (1.0 / 2.0) * R18 * (-R20 * wdir_y + R21 * R4 + 2)) /
                 sqrt(pow(R16, 2) + pow(R17, 2) + pow(R18, 2));
    double result = ((R14 > 0) ? (R22) : (-R22));
    return result;
}

double CosmicTrack_DCA_Deriv_a1(double a0, double b0, double a1, double b1, double wire_x,
                                double wire_y, double wire_z, double wdir_x, double wdir_y,
                                double wdir_z, double plane_x, double plane_y, double plane_z,
                                double panel_x, double panel_y, double panel_z)
{
    double R0 = pow(a1, 2);
    double R1 = R0 + pow(b1, 2) + 1;
    double R2 = pow(R1, -1.0 / 2.0);
    double R3 = R2 * wdir_z;
    double R4 = R2 * wdir_x;
    double R5 = R4 * a1;
    double R6 = R2 * b1;
    double R7 = R6 * wdir_y;
    double R8 = R3 + R5 + R7;
    double R9 = 1 - pow(R8, 2);
    double R10 = 1.0 / R9;
    double R11 = R2 * wire_z;
    double R12 = a0 - wire_x;
    double R13 = R12 * R2;
    double R14 = R13 * a1;
    double R15 = b0 - wire_y;
    double R16 = R15 * R6;
    double R17 = -R11 + R14 + R16;
    double R18 = R12 * wdir_x + R15 * wdir_y - wdir_z * wire_z;
    double R19 = -R17 * R8 + R18;
    double R20 = R10 * R19;
    double R21 = R11 - R14 - R16 + R18 * R8;
    double R22 = R10 * R21;
    double R23 = R2 * R22;
    double R24 = -R20 * wdir_z + R23 - wire_z;
    double R25 = R12 - R20 * wdir_x + R23 * a1;
    double R26 = R15 - R20 * wdir_y + R22 * R6;
    double R27 = pow(R1, -3.0 / 2.0);
    double R28 = R27 * a1;
    double R29 = 2 * R22;
    double R30 = R28 * wire_z;
    double R31 = R28 * b1;
    double R32 = R15 * R31;
    double R33 = R0 * R27;
    double R34 = R12 * R33;
    double R35 = R28 * wdir_z;
    double R36 = R31 * wdir_y;
    double R37 = R33 * wdir_x;
    double R38 = 2 * R10;
    double R39 = R38 * (-R13 + R18 * (-R35 - R36 - R37 + R4) - R30 + R32 + R34);
    double R40 = R2 * R39;
    double R41 = 2 * R8 * (-2 * R35 - 2 * R36 - 2 * R37 + 2 * R4) / pow(R9, 2);
    double R42 = R19 * R41;
    double R43 = R21 * R41;
    double R44 = R2 * R43;
    double R45 = R38 * (R17 * (R35 + R36 + R37 - R4) + (-R3 - R5 - R7) * (R13 + R30 - R32 - R34));
    double R46 =
        ((1.0 / 2.0) * R24 * (-R28 * R29 + R40 - R42 * wdir_z + R44 - R45 * wdir_z) +
         (1.0 / 2.0) * R25 *
             (2 * R23 - R29 * R33 + R40 * a1 - R42 * wdir_x + R44 * a1 - R45 * wdir_x) +
         (1.0 / 2.0) * R26 * (-R29 * R31 + R39 * R6 - R42 * wdir_y + R43 * R6 - R45 * wdir_y)) /
        sqrt(pow(R24, 2) + pow(R25, 2) + pow(R26, 2));
    double result = ((R20 > 0) ? (R46) : (-R46));
    return result;
}

double CosmicTrack_DCA_Deriv_b1(double a0, double b0, double a1, double b1, double wire_x,
                                double wire_y, double wire_z, double wdir_x, double wdir_y,
                                double wdir_z, double plane_x, double plane_y, double plane_z,
                                double panel_x, double panel_y, double panel_z)
{
    double R0 = pow(b1, 2);
    double R1 = R0 + pow(a1, 2) + 1;
    double R2 = pow(R1, -1.0 / 2.0);
    double R3 = R2 * wdir_z;
    double R4 = R2 * a1;
    double R5 = R4 * wdir_x;
    double R6 = R2 * wdir_y;
    double R7 = R6 * b1;
    double R8 = R3 + R5 + R7;
    double R9 = 1 - pow(R8, 2);
    double R10 = 1.0 / R9;
    double R11 = R2 * wire_z;
    double R12 = a0 - wire_x;
    double R13 = R12 * R4;
    double R14 = b0 - wire_y;
    double R15 = R14 * R2;
    double R16 = R15 * b1;
    double R17 = -R11 + R13 + R16;
    double R18 = R12 * wdir_x + R14 * wdir_y - wdir_z * wire_z;
    double R19 = -R17 * R8 + R18;
    double R20 = R10 * R19;
    double R21 = R11 - R13 - R16 + R18 * R8;
    double R22 = R10 * R21;
    double R23 = R2 * R22;
    double R24 = -R20 * wdir_z + R23 - wire_z;
    double R25 = R12 - R20 * wdir_x + R22 * R4;
    double R26 = R14 - R20 * wdir_y + R23 * b1;
    double R27 = pow(R1, -3.0 / 2.0);
    double R28 = R27 * b1;
    double R29 = 2 * R22;
    double R30 = R28 * wire_z;
    double R31 = R28 * a1;
    double R32 = R12 * R31;
    double R33 = R0 * R27;
    double R34 = R14 * R33;
    double R35 = R28 * wdir_z;
    double R36 = R31 * wdir_x;
    double R37 = R33 * wdir_y;
    double R38 = 2 * R10;
    double R39 = R38 * (-R15 + R18 * (-R35 - R36 - R37 + R6) - R30 + R32 + R34);
    double R40 = R2 * R39;
    double R41 = 2 * R8 * (-2 * R35 - 2 * R36 - 2 * R37 + 2 * R6) / pow(R9, 2);
    double R42 = R19 * R41;
    double R43 = R21 * R41;
    double R44 = R2 * R43;
    double R45 = R38 * (R17 * (R35 + R36 + R37 - R6) + (-R3 - R5 - R7) * (R15 + R30 - R32 - R34));
    double R46 =
        ((1.0 / 2.0) * R24 * (-R28 * R29 + R40 - R42 * wdir_z + R44 - R45 * wdir_z) +
         (1.0 / 2.0) * R25 * (-R29 * R31 + R39 * R4 + R4 * R43 - R42 * wdir_x - R45 * wdir_x) +
         (1.0 / 2.0) * R26 *
             (2 * R23 - R29 * R33 + R40 * b1 - R42 * wdir_y + R44 * b1 - R45 * wdir_y)) /
        sqrt(pow(R24, 2) + pow(R25, 2) + pow(R26, 2));
    double result = ((R20 > 0) ? (R46) : (-R46));
    return result;
}

double CosmicTrack_DCA_Deriv_plane_dx(double a0, double b0, double a1, double b1, double wire_x,
                                      double wire_y, double wire_z, double wdir_x, double wdir_y,
                                      double wdir_z, double plane_x, double plane_y, double plane_z,
                                      double panel_x, double panel_y, double panel_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R1 = R0 * wdir_z;
    double R2 = R0 * a1;
    double R3 = R2 * wdir_x;
    double R4 = R0 * b1;
    double R5 = R4 * wdir_y;
    double R6 = R1 + R3 + R5;
    double R7 = 1.0 / (1 - pow(R6, 2));
    double R8 = R0 * wire_z;
    double R9 = a0 - wire_x;
    double R10 = R2 * R9;
    double R11 = b0 - wire_y;
    double R12 = R11 * R4;
    double R13 = R11 * wdir_y + R9 * wdir_x - wdir_z * wire_z;
    double R14 = R7 * (R13 - R6 * (R10 + R12 - R8));
    double R15 = R7 * (-R10 - R12 + R13 * R6 + R8);
    double R16 = R0 * R15 - R14 * wdir_z - wire_z;
    double R17 = -R14 * wdir_x + R15 * R2 + R9;
    double R18 = R11 - R14 * wdir_y + R15 * R4;
    double R19 = 2 * R7;
    double R20 = R19 * (-R2 * (-R1 - R3 - R5) - wdir_x);
    double R21 = R19 * (R2 - R6 * wdir_x);
    double R22 = ((1.0 / 2.0) * R16 * (R0 * R21 - R20 * wdir_z) +
                  (1.0 / 2.0) * R17 * (R2 * R21 - R20 * wdir_x - 2) +
                  (1.0 / 2.0) * R18 * (-R20 * wdir_y + R21 * R4)) /
                 sqrt(pow(R16, 2) + pow(R17, 2) + pow(R18, 2));
    double result = ((R14 > 0) ? (R22) : (-R22));
    return result;
}

double CosmicTrack_DCA_Deriv_plane_dy(double a0, double b0, double a1, double b1, double wire_x,
                                      double wire_y, double wire_z, double wdir_x, double wdir_y,
                                      double wdir_z, double plane_x, double plane_y, double plane_z,
                                      double panel_x, double panel_y, double panel_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R1 = R0 * wdir_z;
    double R2 = R0 * a1;
    double R3 = R2 * wdir_x;
    double R4 = R0 * b1;
    double R5 = R4 * wdir_y;
    double R6 = R1 + R3 + R5;
    double R7 = 1.0 / (1 - pow(R6, 2));
    double R8 = R0 * wire_z;
    double R9 = a0 - wire_x;
    double R10 = R2 * R9;
    double R11 = b0 - wire_y;
    double R12 = R11 * R4;
    double R13 = R11 * wdir_y + R9 * wdir_x - wdir_z * wire_z;
    double R14 = R7 * (R13 - R6 * (R10 + R12 - R8));
    double R15 = R7 * (-R10 - R12 + R13 * R6 + R8);
    double R16 = R0 * R15 - R14 * wdir_z - wire_z;
    double R17 = -R14 * wdir_x + R15 * R2 + R9;
    double R18 = R11 - R14 * wdir_y + R15 * R4;
    double R19 = 2 * R7;
    double R20 = R19 * (-R4 * (-R1 - R3 - R5) - wdir_y);
    double R21 = R19 * (R4 - R6 * wdir_y);
    double R22 = ((1.0 / 2.0) * R16 * (R0 * R21 - R20 * wdir_z) +
                  (1.0 / 2.0) * R17 * (R2 * R21 - R20 * wdir_x) +
                  (1.0 / 2.0) * R18 * (-R20 * wdir_y + R21 * R4 - 2)) /
                 sqrt(pow(R16, 2) + pow(R17, 2) + pow(R18, 2));
    double result = ((R14 > 0) ? (R22) : (-R22));
    return result;
}

double CosmicTrack_DCA_Deriv_plane_dz(double a0, double b0, double a1, double b1, double wire_x,
                                      double wire_y, double wire_z, double wdir_x, double wdir_y,
                                      double wdir_z, double plane_x, double plane_y, double plane_z,
                                      double panel_x, double panel_y, double panel_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R1 = R0 * wdir_z;
    double R2 = R0 * a1;
    double R3 = R2 * wdir_x;
    double R4 = R0 * b1;
    double R5 = R4 * wdir_y;
    double R6 = R1 + R3 + R5;
    double R7 = 1.0 / (1 - pow(R6, 2));
    double R8 = R0 * wire_z;
    double R9 = a0 - wire_x;
    double R10 = R2 * R9;
    double R11 = b0 - wire_y;
    double R12 = R11 * R4;
    double R13 = R11 * wdir_y + R9 * wdir_x - wdir_z * wire_z;
    double R14 = R7 * (R13 - R6 * (R10 + R12 - R8));
    double R15 = R7 * (-R10 - R12 + R13 * R6 + R8);
    double R16 = R0 * R15 - R14 * wdir_z - wire_z;
    double R17 = -R14 * wdir_x + R15 * R2 + R9;
    double R18 = R11 - R14 * wdir_y + R15 * R4;
    double R19 = 2 * R7;
    double R20 = R19 * (-R0 * (-R1 - R3 - R5) - wdir_z);
    double R21 = R19 * (R0 - R6 * wdir_z);
    double R22 = ((1.0 / 2.0) * R16 * (R0 * R21 - R20 * wdir_z - 2) +
                  (1.0 / 2.0) * R17 * (R2 * R21 - R20 * wdir_x) +
                  (1.0 / 2.0) * R18 * (-R20 * wdir_y + R21 * R4)) /
                 sqrt(pow(R16, 2) + pow(R17, 2) + pow(R18, 2));
    double result = ((R14 > 0) ? (R22) : (-R22));
    return result;
}

double CosmicTrack_DCA_Deriv_plane_a(double a0, double b0, double a1, double b1, double wire_x,
                                     double wire_y, double wire_z, double wdir_x, double wdir_y,
                                     double wdir_z, double plane_x, double plane_y, double plane_z,
                                     double panel_x, double panel_y, double panel_z)
{
    double R0 = -wire_z;
    double R1 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R2 = R1 * wdir_z;
    double R3 = R1 * a1;
    double R4 = R3 * wdir_x;
    double R5 = R1 * wdir_y;
    double R6 = R5 * b1;
    double R7 = R2 + R4 + R6;
    double R8 = 1 - pow(R7, 2);
    double R9 = 1.0 / R8;
    double R10 = R1 * wire_z;
    double R11 = a0 - wire_x;
    double R12 = R11 * R3;
    double R13 = b0 - wire_y;
    double R14 = R1 * b1;
    double R15 = R13 * R14;
    double R16 = -R10 + R12 + R15;
    double R17 = R11 * wdir_x + R13 * wdir_y - wdir_z * wire_z;
    double R18 = -R16 * R7 + R17;
    double R19 = R18 * R9;
    double R20 = R19 * wdir_z;
    double R21 = R10 - R12 - R15 + R17 * R7;
    double R22 = R21 * R9;
    double R23 = R0 + R1 * R22 - R20;
    double R24 = R11 - R19 * wdir_x + R22 * R3;
    double R25 = R19 * wdir_y;
    double R26 = R13 + R14 * R22 - R25;
    double R27 = -plane_y + wire_y;
    double R28 = R1 * R27;
    double R29 = R0 + plane_z;
    double R30 = R14 * R29;
    double R31 = R2 * b1;
    double R32 = R13 * wdir_z + R27 * wdir_z + R29 * wdir_y + wdir_y * wire_z;
    double R33 = 2 * R9;
    double R34 = R33 * (R17 * (R31 - R5) - R28 - R30 + R32 * R7);
    double R35 = 2 * R7 * (2 * R31 - 2 * R5) / pow(R8, 2);
    double R36 = R18 * R35;
    double R37 = R33 * (R16 * (-R31 + R5) + R32 + (R28 + R30) * (-R2 - R4 - R6));
    double R38 = R21 * R35;
    double R39 = ((1.0 / 2.0) * R23 *
                      (R1 * R34 + R1 * R38 + 2 * R25 - R36 * wdir_z - R37 * wdir_z - 2 * plane_y +
                       2 * wire_y) +
                  (1.0 / 2.0) * R24 * (R3 * R34 + R3 * R38 - R36 * wdir_x - R37 * wdir_x) +
                  (1.0 / 2.0) * R26 *
                      (R14 * R34 + R14 * R38 - 2 * R20 - R36 * wdir_y - R37 * wdir_y + 2 * plane_z -
                       2 * wire_z)) /
                 sqrt(pow(R23, 2) + pow(R24, 2) + pow(R26, 2));
    double result = ((R19 > 0) ? (R39) : (-R39));
    return result;
}

double CosmicTrack_DCA_Deriv_plane_b(double a0, double b0, double a1, double b1, double wire_x,
                                     double wire_y, double wire_z, double wdir_x, double wdir_y,
                                     double wdir_z, double plane_x, double plane_y, double plane_z,
                                     double panel_x, double panel_y, double panel_z)
{
    double R0 = -wire_z;
    double R1 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R2 = R1 * wdir_z;
    double R3 = R1 * wdir_x;
    double R4 = R3 * a1;
    double R5 = R1 * b1;
    double R6 = R5 * wdir_y;
    double R7 = R2 + R4 + R6;
    double R8 = 1 - pow(R7, 2);
    double R9 = 1.0 / R8;
    double R10 = R1 * wire_z;
    double R11 = -wire_x;
    double R12 = R11 + a0;
    double R13 = R1 * a1;
    double R14 = R12 * R13;
    double R15 = b0 - wire_y;
    double R16 = R15 * R5;
    double R17 = -R10 + R14 + R16;
    double R18 = R12 * wdir_x + R15 * wdir_y - wdir_z * wire_z;
    double R19 = -R17 * R7 + R18;
    double R20 = R19 * R9;
    double R21 = R20 * wdir_z;
    double R22 = R10 - R14 - R16 + R18 * R7;
    double R23 = R22 * R9;
    double R24 = R0 + R1 * R23 - R21;
    double R25 = R20 * wdir_x;
    double R26 = R12 + R13 * R23 - R25;
    double R27 = R15 - R20 * wdir_y + R23 * R5;
    double R28 = R11 + plane_x;
    double R29 = R1 * R28;
    double R30 = R0 + plane_z;
    double R31 = R13 * R30;
    double R32 = R2 * a1;
    double R33 = R12 * wdir_z + R28 * wdir_z + R30 * wdir_x - wdir_x * wire_z;
    double R34 = 2 * R9;
    double R35 = R34 * (R18 * (R3 + R32) - R29 - R31 + R33 * R7);
    double R36 = 2 * R7 * (2 * R3 + 2 * R32) / pow(R8, 2);
    double R37 = R19 * R36;
    double R38 = R34 * (R17 * (-R3 - R32) + R33 + (R29 + R31) * (-R2 - R4 - R6));
    double R39 = R22 * R36;
    double R40 = ((1.0 / 2.0) * R24 *
                      (R1 * R35 + R1 * R39 - 2 * R25 - R37 * wdir_z - R38 * wdir_z + 2 * plane_x -
                       2 * wire_x) +
                  (1.0 / 2.0) * R26 *
                      (R13 * R35 + R13 * R39 - 2 * R21 - R37 * wdir_x - R38 * wdir_x + 2 * plane_z -
                       2 * wire_z) +
                  (1.0 / 2.0) * R27 * (R35 * R5 - R37 * wdir_y - R38 * wdir_y + R39 * R5)) /
                 sqrt(pow(R24, 2) + pow(R26, 2) + pow(R27, 2));
    double result = ((R20 > 0) ? (R40) : (-R40));
    return result;
}

double CosmicTrack_DCA_Deriv_plane_g(double a0, double b0, double a1, double b1, double wire_x,
                                     double wire_y, double wire_z, double wdir_x, double wdir_y,
                                     double wdir_z, double plane_x, double plane_y, double plane_z,
                                     double panel_x, double panel_y, double panel_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R1 = R0 * wdir_z;
    double R2 = R0 * a1;
    double R3 = R2 * wdir_x;
    double R4 = R0 * b1;
    double R5 = R4 * wdir_y;
    double R6 = R1 + R3 + R5;
    double R7 = 1 - pow(R6, 2);
    double R8 = 1.0 / R7;
    double R9 = R0 * wire_z;
    double R10 = a0 - wire_x;
    double R11 = R10 * R2;
    double R12 = -wire_y;
    double R13 = R12 + b0;
    double R14 = R13 * R4;
    double R15 = R11 + R14 - R9;
    double R16 = R10 * wdir_x + R13 * wdir_y - wdir_z * wire_z;
    double R17 = -R15 * R6 + R16;
    double R18 = R17 * R8;
    double R19 = -R11 - R14 + R16 * R6 + R9;
    double R20 = R19 * R8;
    double R21 = R0 * R20 - R18 * wdir_z - wire_z;
    double R22 = R18 * wdir_x;
    double R23 = R10 + R2 * R20 - R22;
    double R24 = R18 * wdir_y;
    double R25 = R13 + R20 * R4 - R24;
    double R26 = R12 + plane_y;
    double R27 = R2 * R26;
    double R28 = -plane_x + wire_x;
    double R29 = R28 * R4;
    double R30 = R2 * wdir_y;
    double R31 = R4 * wdir_x;
    double R32 = R10 * wdir_y - R13 * wdir_x + R26 * wdir_x + R28 * wdir_y;
    double R33 = 2 * R8;
    double R34 = R33 * (R16 * (R30 - R31) - R27 - R29 + R32 * R6);
    double R35 = 2 * R6 * (2 * R30 - 2 * R31) / pow(R7, 2);
    double R36 = R17 * R35;
    double R37 = R33 * (R15 * (-R30 + R31) + R32 + (R27 + R29) * (-R1 - R3 - R5));
    double R38 = R19 * R35;
    double R39 = ((1.0 / 2.0) * R21 * (R0 * R34 + R0 * R38 - R36 * wdir_z - R37 * wdir_z) +
                  (1.0 / 2.0) * R23 *
                      (R2 * R34 + R2 * R38 - 2 * R24 - R36 * wdir_x - R37 * wdir_x + 2 * plane_y -
                       2 * wire_y) +
                  (1.0 / 2.0) * R25 *
                      (2 * R22 + R34 * R4 - R36 * wdir_y - R37 * wdir_y + R38 * R4 - 2 * plane_x +
                       2 * wire_x)) /
                 sqrt(pow(R21, 2) + pow(R23, 2) + pow(R25, 2));
    double result = ((R18 > 0) ? (R39) : (-R39));
    return result;
}

double CosmicTrack_DCA_Deriv_panel_dx(double a0, double b0, double a1, double b1, double wire_x,
                                      double wire_y, double wire_z, double wdir_x, double wdir_y,
                                      double wdir_z, double plane_x, double plane_y, double plane_z,
                                      double panel_x, double panel_y, double panel_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R1 = R0 * wdir_z;
    double R2 = R0 * a1;
    double R3 = R2 * wdir_x;
    double R4 = R0 * b1;
    double R5 = R4 * wdir_y;
    double R6 = R1 + R3 + R5;
    double R7 = 1.0 / (1 - pow(R6, 2));
    double R8 = R0 * wire_z;
    double R9 = a0 - wire_x;
    double R10 = R2 * R9;
    double R11 = b0 - wire_y;
    double R12 = R11 * R4;
    double R13 = R11 * wdir_y + R9 * wdir_x - wdir_z * wire_z;
    double R14 = R7 * (R13 - R6 * (R10 + R12 - R8));
    double R15 = R7 * (-R10 - R12 + R13 * R6 + R8);
    double R16 = R0 * R15 - R14 * wdir_z - wire_z;
    double R17 = -R14 * wdir_x + R15 * R2 + R9;
    double R18 = R11 - R14 * wdir_y + R15 * R4;
    double R19 = 2 * R7;
    double R20 = R19 * (-R2 * (-R1 - R3 - R5) - wdir_x);
    double R21 = R19 * (R2 - R6 * wdir_x);
    double R22 = ((1.0 / 2.0) * R16 * (R0 * R21 - R20 * wdir_z) +
                  (1.0 / 2.0) * R17 * (R2 * R21 - R20 * wdir_x - 2) +
                  (1.0 / 2.0) * R18 * (-R20 * wdir_y + R21 * R4)) /
                 sqrt(pow(R16, 2) + pow(R17, 2) + pow(R18, 2));
    double result = ((R14 > 0) ? (R22) : (-R22));
    return result;
}

double CosmicTrack_DCA_Deriv_panel_dy(double a0, double b0, double a1, double b1, double wire_x,
                                      double wire_y, double wire_z, double wdir_x, double wdir_y,
                                      double wdir_z, double plane_x, double plane_y, double plane_z,
                                      double panel_x, double panel_y, double panel_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R1 = R0 * wdir_z;
    double R2 = R0 * a1;
    double R3 = R2 * wdir_x;
    double R4 = R0 * b1;
    double R5 = R4 * wdir_y;
    double R6 = R1 + R3 + R5;
    double R7 = 1.0 / (1 - pow(R6, 2));
    double R8 = R0 * wire_z;
    double R9 = a0 - wire_x;
    double R10 = R2 * R9;
    double R11 = b0 - wire_y;
    double R12 = R11 * R4;
    double R13 = R11 * wdir_y + R9 * wdir_x - wdir_z * wire_z;
    double R14 = R7 * (R13 - R6 * (R10 + R12 - R8));
    double R15 = R7 * (-R10 - R12 + R13 * R6 + R8);
    double R16 = R0 * R15 - R14 * wdir_z - wire_z;
    double R17 = -R14 * wdir_x + R15 * R2 + R9;
    double R18 = R11 - R14 * wdir_y + R15 * R4;
    double R19 = 2 * R7;
    double R20 = R19 * (-R4 * (-R1 - R3 - R5) - wdir_y);
    double R21 = R19 * (R4 - R6 * wdir_y);
    double R22 = ((1.0 / 2.0) * R16 * (R0 * R21 - R20 * wdir_z) +
                  (1.0 / 2.0) * R17 * (R2 * R21 - R20 * wdir_x) +
                  (1.0 / 2.0) * R18 * (-R20 * wdir_y + R21 * R4 - 2)) /
                 sqrt(pow(R16, 2) + pow(R17, 2) + pow(R18, 2));
    double result = ((R14 > 0) ? (R22) : (-R22));
    return result;
}

double CosmicTrack_DCA_Deriv_panel_dz(double a0, double b0, double a1, double b1, double wire_x,
                                      double wire_y, double wire_z, double wdir_x, double wdir_y,
                                      double wdir_z, double plane_x, double plane_y, double plane_z,
                                      double panel_x, double panel_y, double panel_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R1 = R0 * wdir_z;
    double R2 = R0 * a1;
    double R3 = R2 * wdir_x;
    double R4 = R0 * b1;
    double R5 = R4 * wdir_y;
    double R6 = R1 + R3 + R5;
    double R7 = 1.0 / (1 - pow(R6, 2));
    double R8 = R0 * wire_z;
    double R9 = a0 - wire_x;
    double R10 = R2 * R9;
    double R11 = b0 - wire_y;
    double R12 = R11 * R4;
    double R13 = R11 * wdir_y + R9 * wdir_x - wdir_z * wire_z;
    double R14 = R7 * (R13 - R6 * (R10 + R12 - R8));
    double R15 = R7 * (-R10 - R12 + R13 * R6 + R8);
    double R16 = R0 * R15 - R14 * wdir_z - wire_z;
    double R17 = -R14 * wdir_x + R15 * R2 + R9;
    double R18 = R11 - R14 * wdir_y + R15 * R4;
    double R19 = 2 * R7;
    double R20 = R19 * (-R0 * (-R1 - R3 - R5) - wdir_z);
    double R21 = R19 * (R0 - R6 * wdir_z);
    double R22 = ((1.0 / 2.0) * R16 * (R0 * R21 - R20 * wdir_z - 2) +
                  (1.0 / 2.0) * R17 * (R2 * R21 - R20 * wdir_x) +
                  (1.0 / 2.0) * R18 * (-R20 * wdir_y + R21 * R4)) /
                 sqrt(pow(R16, 2) + pow(R17, 2) + pow(R18, 2));
    double result = ((R14 > 0) ? (R22) : (-R22));
    return result;
}

double CosmicTrack_DCA_Deriv_panel_a(double a0, double b0, double a1, double b1, double wire_x,
                                     double wire_y, double wire_z, double wdir_x, double wdir_y,
                                     double wdir_z, double plane_x, double plane_y, double plane_z,
                                     double panel_x, double panel_y, double panel_z)
{
    double R0 = -wire_z;
    double R1 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R2 = R1 * wdir_z;
    double R3 = R1 * a1;
    double R4 = R3 * wdir_x;
    double R5 = R1 * wdir_y;
    double R6 = R5 * b1;
    double R7 = R2 + R4 + R6;
    double R8 = 1 - pow(R7, 2);
    double R9 = 1.0 / R8;
    double R10 = R1 * wire_z;
    double R11 = a0 - wire_x;
    double R12 = R11 * R3;
    double R13 = b0 - wire_y;
    double R14 = R1 * b1;
    double R15 = R13 * R14;
    double R16 = -R10 + R12 + R15;
    double R17 = R11 * wdir_x + R13 * wdir_y - wdir_z * wire_z;
    double R18 = -R16 * R7 + R17;
    double R19 = R18 * R9;
    double R20 = R19 * wdir_z;
    double R21 = R10 - R12 - R15 + R17 * R7;
    double R22 = R21 * R9;
    double R23 = R0 + R1 * R22 - R20;
    double R24 = R11 - R19 * wdir_x + R22 * R3;
    double R25 = R19 * wdir_y;
    double R26 = R13 + R14 * R22 - R25;
    double R27 = -panel_y + wire_y;
    double R28 = R1 * R27;
    double R29 = R0 + panel_z;
    double R30 = R14 * R29;
    double R31 = R2 * b1;
    double R32 = R13 * wdir_z + R27 * wdir_z + R29 * wdir_y + wdir_y * wire_z;
    double R33 = 2 * R9;
    double R34 = R33 * (R17 * (R31 - R5) - R28 - R30 + R32 * R7);
    double R35 = 2 * R7 * (2 * R31 - 2 * R5) / pow(R8, 2);
    double R36 = R18 * R35;
    double R37 = R33 * (R16 * (-R31 + R5) + R32 + (R28 + R30) * (-R2 - R4 - R6));
    double R38 = R21 * R35;
    double R39 = ((1.0 / 2.0) * R23 *
                      (R1 * R34 + R1 * R38 + 2 * R25 - R36 * wdir_z - R37 * wdir_z - 2 * panel_y +
                       2 * wire_y) +
                  (1.0 / 2.0) * R24 * (R3 * R34 + R3 * R38 - R36 * wdir_x - R37 * wdir_x) +
                  (1.0 / 2.0) * R26 *
                      (R14 * R34 + R14 * R38 - 2 * R20 - R36 * wdir_y - R37 * wdir_y + 2 * panel_z -
                       2 * wire_z)) /
                 sqrt(pow(R23, 2) + pow(R24, 2) + pow(R26, 2));
    double result = ((R19 > 0) ? (R39) : (-R39));
    return result;
}

double CosmicTrack_DCA_Deriv_panel_b(double a0, double b0, double a1, double b1, double wire_x,
                                     double wire_y, double wire_z, double wdir_x, double wdir_y,
                                     double wdir_z, double plane_x, double plane_y, double plane_z,
                                     double panel_x, double panel_y, double panel_z)
{
    double R0 = -wire_z;
    double R1 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R2 = R1 * wdir_z;
    double R3 = R1 * wdir_x;
    double R4 = R3 * a1;
    double R5 = R1 * b1;
    double R6 = R5 * wdir_y;
    double R7 = R2 + R4 + R6;
    double R8 = 1 - pow(R7, 2);
    double R9 = 1.0 / R8;
    double R10 = R1 * wire_z;
    double R11 = -wire_x;
    double R12 = R11 + a0;
    double R13 = R1 * a1;
    double R14 = R12 * R13;
    double R15 = b0 - wire_y;
    double R16 = R15 * R5;
    double R17 = -R10 + R14 + R16;
    double R18 = R12 * wdir_x + R15 * wdir_y - wdir_z * wire_z;
    double R19 = -R17 * R7 + R18;
    double R20 = R19 * R9;
    double R21 = R20 * wdir_z;
    double R22 = R10 - R14 - R16 + R18 * R7;
    double R23 = R22 * R9;
    double R24 = R0 + R1 * R23 - R21;
    double R25 = R20 * wdir_x;
    double R26 = R12 + R13 * R23 - R25;
    double R27 = R15 - R20 * wdir_y + R23 * R5;
    double R28 = R11 + panel_x;
    double R29 = R1 * R28;
    double R30 = R0 + panel_z;
    double R31 = R13 * R30;
    double R32 = R2 * a1;
    double R33 = R12 * wdir_z + R28 * wdir_z + R30 * wdir_x - wdir_x * wire_z;
    double R34 = 2 * R9;
    double R35 = R34 * (R18 * (R3 + R32) - R29 - R31 + R33 * R7);
    double R36 = 2 * R7 * (2 * R3 + 2 * R32) / pow(R8, 2);
    double R37 = R19 * R36;
    double R38 = R34 * (R17 * (-R3 - R32) + R33 + (R29 + R31) * (-R2 - R4 - R6));
    double R39 = R22 * R36;
    double R40 = ((1.0 / 2.0) * R24 *
                      (R1 * R35 + R1 * R39 - 2 * R25 - R37 * wdir_z - R38 * wdir_z + 2 * panel_x -
                       2 * wire_x) +
                  (1.0 / 2.0) * R26 *
                      (R13 * R35 + R13 * R39 - 2 * R21 - R37 * wdir_x - R38 * wdir_x + 2 * panel_z -
                       2 * wire_z) +
                  (1.0 / 2.0) * R27 * (R35 * R5 - R37 * wdir_y - R38 * wdir_y + R39 * R5)) /
                 sqrt(pow(R24, 2) + pow(R26, 2) + pow(R27, 2));
    double result = ((R20 > 0) ? (R40) : (-R40));
    return result;
}

double CosmicTrack_DCA_Deriv_panel_g(double a0, double b0, double a1, double b1, double wire_x,
                                     double wire_y, double wire_z, double wdir_x, double wdir_y,
                                     double wdir_z, double plane_x, double plane_y, double plane_z,
                                     double panel_x, double panel_y, double panel_z)
{
    double R0 = pow(pow(a1, 2) + pow(b1, 2) + 1, -1.0 / 2.0);
    double R1 = R0 * wdir_z;
    double R2 = R0 * a1;
    double R3 = R2 * wdir_x;
    double R4 = R0 * b1;
    double R5 = R4 * wdir_y;
    double R6 = R1 + R3 + R5;
    double R7 = 1 - pow(R6, 2);
    double R8 = 1.0 / R7;
    double R9 = R0 * wire_z;
    double R10 = a0 - wire_x;
    double R11 = R10 * R2;
    double R12 = -wire_y;
    double R13 = R12 + b0;
    double R14 = R13 * R4;
    double R15 = R11 + R14 - R9;
    double R16 = R10 * wdir_x + R13 * wdir_y - wdir_z * wire_z;
    double R17 = -R15 * R6 + R16;
    double R18 = R17 * R8;
    double R19 = -R11 - R14 + R16 * R6 + R9;
    double R20 = R19 * R8;
    double R21 = R0 * R20 - R18 * wdir_z - wire_z;
    double R22 = R18 * wdir_x;
    double R23 = R10 + R2 * R20 - R22;
    double R24 = R18 * wdir_y;
    double R25 = R13 + R20 * R4 - R24;
    double R26 = R12 + panel_y;
    double R27 = R2 * R26;
    double R28 = -panel_x + wire_x;
    double R29 = R28 * R4;
    double R30 = R2 * wdir_y;
    double R31 = R4 * wdir_x;
    double R32 = R10 * wdir_y - R13 * wdir_x + R26 * wdir_x + R28 * wdir_y;
    double R33 = 2 * R8;
    double R34 = R33 * (R16 * (R30 - R31) - R27 - R29 + R32 * R6);
    double R35 = 2 * R6 * (2 * R30 - 2 * R31) / pow(R7, 2);
    double R36 = R17 * R35;
    double R37 = R33 * (R15 * (-R30 + R31) + R32 + (R27 + R29) * (-R1 - R3 - R5));
    double R38 = R19 * R35;
    double R39 = ((1.0 / 2.0) * R21 * (R0 * R34 + R0 * R38 - R36 * wdir_z - R37 * wdir_z) +
                  (1.0 / 2.0) * R23 *
                      (R2 * R34 + R2 * R38 - 2 * R24 - R36 * wdir_x - R37 * wdir_x + 2 * panel_y -
                       2 * wire_y) +
                  (1.0 / 2.0) * R25 *
                      (2 * R22 + R34 * R4 - R36 * wdir_y - R37 * wdir_y + R38 * R4 - 2 * panel_x +
                       2 * wire_x)) /
                 sqrt(pow(R21, 2) + pow(R23, 2) + pow(R25, 2));
    double result = ((R18 > 0) ? (R39) : (-R39));
    return result;
}

std::vector<float> CosmicTrack_DCA_LocalDeriv(double a0, double b0, double a1, double b1,
                                              double wire_x, double wire_y, double wire_z,
                                              double wdir_x, double wdir_y, double wdir_z,
                                              double plane_x, double plane_y, double plane_z,
                                              double panel_x, double panel_y, double panel_z)
{
    std::vector<float> result = {
        (float)CosmicTrack_DCA_Deriv_a0(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x, wdir_y,
                                        wdir_z, plane_x, plane_y, plane_z, panel_x, panel_y,
                                        panel_z),
        (float)CosmicTrack_DCA_Deriv_b0(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x, wdir_y,
                                        wdir_z, plane_x, plane_y, plane_z, panel_x, panel_y,
                                        panel_z),
        (float)CosmicTrack_DCA_Deriv_a1(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x, wdir_y,
                                        wdir_z, plane_x, plane_y, plane_z, panel_x, panel_y,
                                        panel_z),
        (float)CosmicTrack_DCA_Deriv_b1(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x, wdir_y,
                                        wdir_z, plane_x, plane_y, plane_z, panel_x, panel_y,
                                        panel_z)};
    return result;
}

std::vector<float> CosmicTrack_DCA_GlobalDeriv(double a0, double b0, double a1, double b1,
                                               double wire_x, double wire_y, double wire_z,
                                               double wdir_x, double wdir_y, double wdir_z,
                                               double plane_x, double plane_y, double plane_z,
                                               double panel_x, double panel_y, double panel_z)
{
    std::vector<float> result = {
        (float)CosmicTrack_DCA_Deriv_plane_dx(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x,
                                              wdir_y, wdir_z, plane_x, plane_y, plane_z, panel_x,
                                              panel_y, panel_z),
        (float)CosmicTrack_DCA_Deriv_plane_dy(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x,
                                              wdir_y, wdir_z, plane_x, plane_y, plane_z, panel_x,
                                              panel_y, panel_z),
        (float)CosmicTrack_DCA_Deriv_plane_dz(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x,
                                              wdir_y, wdir_z, plane_x, plane_y, plane_z, panel_x,
                                              panel_y, panel_z),
        (float)CosmicTrack_DCA_Deriv_plane_a(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x, wdir_y,
                                             wdir_z, plane_x, plane_y, plane_z, panel_x, panel_y,
                                             panel_z),
        (float)CosmicTrack_DCA_Deriv_plane_b(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x, wdir_y,
                                             wdir_z, plane_x, plane_y, plane_z, panel_x, panel_y,
                                             panel_z),
        (float)CosmicTrack_DCA_Deriv_plane_g(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x, wdir_y,
                                             wdir_z, plane_x, plane_y, plane_z, panel_x, panel_y,
                                             panel_z),
        (float)CosmicTrack_DCA_Deriv_panel_dx(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x,
                                              wdir_y, wdir_z, plane_x, plane_y, plane_z, panel_x,
                                              panel_y, panel_z),
        (float)CosmicTrack_DCA_Deriv_panel_dy(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x,
                                              wdir_y, wdir_z, plane_x, plane_y, plane_z, panel_x,
                                              panel_y, panel_z),
        (float)CosmicTrack_DCA_Deriv_panel_dz(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x,
                                              wdir_y, wdir_z, plane_x, plane_y, plane_z, panel_x,
                                              panel_y, panel_z),
        (float)CosmicTrack_DCA_Deriv_panel_a(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x, wdir_y,
                                             wdir_z, plane_x, plane_y, plane_z, panel_x, panel_y,
                                             panel_z),
        (float)CosmicTrack_DCA_Deriv_panel_b(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x, wdir_y,
                                             wdir_z, plane_x, plane_y, plane_z, panel_x, panel_y,
                                             panel_z),
        (float)CosmicTrack_DCA_Deriv_panel_g(a0, b0, a1, b1, wire_x, wire_y, wire_z, wdir_x, wdir_y,
                                             wdir_z, plane_x, plane_y, plane_z, panel_x, panel_y,
                                             panel_z)};
    return result;
}
