#include "Offline/GeneralUtilities/inc/compressPdgId.hh"

namespace mu2e {

  int compressPdgIdCosmic(int pdgId) {

    int ind;

    switch (pdgId) {
    case PDGCode::mu_minus: // mu-
      ind = 0; break;
    case PDGCode::mu_plus: // mu+
      ind = 0; break;
    case PDGCode::gamma: // photon
      ind = 1; break;
    case PDGCode::e_plus: // e+
      ind = 2; break;
    case PDGCode::e_minus: // e-
      ind = 2; break;
    case PDGCode::n0: // neutron
      ind = 3; break;
    case PDGCode::anti_n0: // neutron
      ind = 3; break;
    case PDGCode::proton: // proton
      ind = 4; break;
    case PDGCode::anti_proton: // proton
      ind = 4; break;
    case PDGCode::pi0: // pi0
      ind = 5; break;
    case PDGCode::pi_plus: // pi+
      ind = 5; break;
    case PDGCode::pi_minus: // pi-
      ind = 5; break;
    case PDGCode::K_L0: // k0 L
      ind = 6; break;
    case PDGCode::K_S0: // k0 S
      ind = 6; break;
    case PDGCode::K0: // k0
      ind = 6; break;
    case PDGCode::anti_K0: // anti k0
      ind = 6; break;
    case PDGCode::K_plus: // k+
      ind = 6; break;
    case PDGCode::K_minus: // k-
      ind = 6; break;
    default: // others
      ind = 7; break;
    }

    return ind;
  }


}
