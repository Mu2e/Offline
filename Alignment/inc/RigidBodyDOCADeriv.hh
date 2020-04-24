
# ifndef RIGIDBODYDOCADERIV_H
# define RIGIDBODYDOCADERIV_H
# include <vector>
double CosmicTrack_DCA(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z);

double CosmicTrack_DCA_Deriv_a0(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_b0(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_a1(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_b1(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_plane_dx(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_plane_dy(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_plane_dz(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_plane_a(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_plane_b(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_plane_g(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_panel_dx(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_panel_dy(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_panel_dz(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_panel_a(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_panel_b(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

double CosmicTrack_DCA_Deriv_panel_g(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

std::vector<float> CosmicTrack_DCA_LocalDeriv(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

std::vector<float> CosmicTrack_DCA_GlobalDeriv(double a0, double b0, double a1, double b1, double wire_x, double wire_y, double wire_z, double wdir_x, double wdir_y, double wdir_z, double plane_x, double plane_y, double plane_z, double panel_straw0x, double panel_straw0y, double panel_straw0z);

# endif

