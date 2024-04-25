#include "Offline/Mu2eUtilities/inc/compressPdgId.hh"

namespace mu2e {

  CompressedPDGCode::enum_type compressPDGCode(PDGCode::enum_type pdgId) {

    int q3 = (pdgId / 10) % 10;
    int q2 = (pdgId / 100) % 10;
    int q1 = (pdgId / 1000) % 10;
    int rr = pdgId / 10000;

    CompressedPDGCode::enum_type code = CompressedPDGCode::unknown;

    PDGCode::enum_type apdgId = PDGCode::enum_type( abs(pdgId) );

    if (pdgId == PDGCode::gamma) {
      code = CompressedPDGCode::gamma;
    } else if (pdgId == PDGCode::e_plus) {
      code = CompressedPDGCode::e_plus;
    } else if (pdgId == PDGCode::e_minus) {
      code = CompressedPDGCode::e_minus;
    } else if (pdgId == PDGCode::mu_plus) {
      code = CompressedPDGCode::mu_plus;
    } else if (pdgId == PDGCode::mu_minus) {
      code = CompressedPDGCode::mu_minus;
    } else if (apdgId >= PDGCode::e_minus && apdgId <= PDGCode::nu_tau) {
      code = (pdgId>0 ?  CompressedPDGCode::lepton :
              CompressedPDGCode::anti_lepton);
    } else if (pdgId == PDGCode::pi_plus) {
      code = CompressedPDGCode::pi_plus;
    } else if (pdgId == PDGCode::pi_minus) {
      code = CompressedPDGCode::pi_minus;
    } else if (pdgId == PDGCode::K_plus) {
      code = CompressedPDGCode::K_plus;
    } else if (pdgId == PDGCode::K_minus) {
      code = CompressedPDGCode::K_minus;
    } else if (pdgId == PDGCode::proton) {
      code = CompressedPDGCode::proton;
    } else if (pdgId == PDGCode::anti_proton) {
      code = CompressedPDGCode::anti_proton;
    } else if (pdgId == PDGCode::n0) {
      code = CompressedPDGCode::n0;
    } else if (pdgId == PDGCode::anti_n0) {
      code = CompressedPDGCode::anti_n0;
    } else if (pdgId == PDGCode::pi0) {
      code = CompressedPDGCode::pi0;
    } else if (pdgId == PDGCode::K_S0) {
      code = CompressedPDGCode::K_S0;
    } else if (pdgId == PDGCode::K_L0) {
      code = CompressedPDGCode::K_L0;
    } else if (q1 == 0 && q2 > 0 && q3 > 0 && rr < 1000) {  // meson
      if (q2 <= 2) {
        code = CompressedPDGCode::ud_meson;
      } else if (q2 == 3) {
        code = CompressedPDGCode::s_meson;
      } else {
        code = CompressedPDGCode::cb_meson;
      }
    } else if (q1 > 0 && q2 > 0 && q3 > 0 && rr < 1000) {  // baryon
      if (q1 <= 2) {
        code = CompressedPDGCode::ud_baryon;
      } else if (q1 == 3) {
        code = CompressedPDGCode::s_baryon;
      } else {
        code = CompressedPDGCode::cb_baryon;
      }
    } else if (pdgId > PDGCode::G4Threshold) {  // nuclei
      code = CompressedPDGCode::nuclei;
    } else if (pdgId == PDGCode::unknown) {
      code = CompressedPDGCode::unknown;
    } else {
      code = CompressedPDGCode::other;
    }
    return code;
  }

  TH1D* compressPDGCodeHisto(art::ServiceHandle<art::TFileService> tfs,
                             std::string name, std::string title) {

    float low = float(CompressedPDGCode::minBin) - 0.5;
    float hgh = float(CompressedPDGCode::maxBin) + 0.5;
    int nbin = CompressedPDGCode::maxBin - CompressedPDGCode::minBin +1;

    auto hp = tfs->make<TH1D>(name.c_str(),title.c_str(),nbin,low,hgh);

    int rootBin = 0;
    for(int ibin=CompressedPDGCode::minBin;
        ibin<=CompressedPDGCode::maxBin; ibin++) {
      rootBin++;
      std::string name = CompressedPDGCode(ibin,false).name();
      hp->GetXaxis()->SetBinLabel(rootBin,name.c_str());
    }

    return hp;
  }


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
