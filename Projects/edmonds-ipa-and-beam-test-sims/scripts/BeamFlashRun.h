#ifndef BeamFlashRun_h_
#define BeamFlashRun_h_


#include "BaseRun.h"

class BeamFlashRun : public BaseRun {

 public:
  BeamFlashRun(std::string filename); 

};

BeamFlashRun::BeamFlashRun(std::string filename) : BaseRun(filename, "Flash") {

  std::cout << "Mean EDep = " << GetSHEDepPlot()->GetMean() << std::endl;
  std::cout << "Flash Charge Deposit per Hit = " << fChargeDepositPerHit << std::endl;
  double n_flashes_per_POT = 0.0022923; // from Doc-DB 3774: Table 1
  fNParticlesPerMicrobunch = n_POT_per_microbunch * n_flashes_per_POT;
}

#endif
