#include "Rtypes.h"
#include "TMath.h"
#include "math.h"
Double_t DIOSC(Double_t *x, Double_t *par) {
  double ee = x[0];
  double norm = par[0];
  double expon(0.023);
  double efactor(1.44);
  double offset(-0.22);
  double emue(104.973);
  double mmu(105.6583715);
  double delta = (emue - ee)/emue;
  double GF(1.1663787e-5);
  static double gamma0 = GF*GF*pow(mmu,5)/(192*pow(TMath::Pi(),3));
  return norm * (1.0/mmu) * (efactor * pow(delta,expon) + offset)* 1e-4 * pow(delta,5);
}

