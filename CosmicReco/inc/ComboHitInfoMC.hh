#ifndef CosmicReco_ComboHitInfoMC_HH
namespace mu2e {
  struct ComboHitInfoMC {
    Int_t _rel; // relation to the 1st hit
    XYZVectorF _mcpos; // position of this hit
  };
}
#endif
