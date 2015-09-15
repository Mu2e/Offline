// the following approximation is from Czarnecki etal, 'Muon decay in orbit:spectrum of high-energy electrons',
// for E>85 MeV
#include "Rtypes.h"
#include "math.h"
Double_t DIOCZ(Double_t *x, Double_t *par) {
  double ee = x[0];
  double norm = par[0];
  double mal(25133);
  //    double mmu(105.654);
  double emu(105.194);
  //    double emue(104.973);
  //    double me(0.511);
  double a5(8.6434e-17);
  double a6(1.16874e-17);
  double a7(-1.87828e-19);
  double a8(9.16327e-20);
  double delta = emu - ee - ee*ee/(2*mal);
  if (delta > 0.0)
    return norm*(a5*pow(delta,5) + a6*pow(delta,6) + a7*pow(delta,7) + a8*pow(delta,8));
  else
    return 0.0;
}


