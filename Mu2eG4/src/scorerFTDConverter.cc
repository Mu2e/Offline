#include "Offline/Mu2eG4/inc/scorerFTDConverter.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "cetlib_except/exception.h"

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>


namespace mu2e{

  scorerFTDConverter::scorerFTDConverter(const scorerDoseType& type, const std::string& method)
  {
    if (type == scorerDoseType::Effective){
       photon_   = scorerFTDTable("Offline/Mu2eG4/data/Photon_effective_dose.dat",  method);
       electron_ = scorerFTDTable("Offline/Mu2eG4/data/Electron_effective_dose.dat",method);
       positron_ = scorerFTDTable("Offline/Mu2eG4/data/Positron_effective_dose.dat",method);
       muminus_  = scorerFTDTable("Offline/Mu2eG4/data/Muminus_effective_dose.dat", method);
       muplus_   = scorerFTDTable("Offline/Mu2eG4/data/Muplus_effective_dose.dat",  method);
       piminus_  = scorerFTDTable("Offline/Mu2eG4/data/Piminus_effective_dose.dat", method);
       piplus_   = scorerFTDTable("Offline/Mu2eG4/data/Piplus_effective_dose.dat",  method);
       proton_   = scorerFTDTable("Offline/Mu2eG4/data/Proton_effective_dose.dat",  method);
       neutron_  = scorerFTDTable("Offline/Mu2eG4/data/Neutron_effective_dose.dat", method);
    } else if (type == scorerDoseType::Ambient){
       photon_   = scorerFTDTable("Offline/Mu2eG4/data/Photon_ambient_dose.dat",  "ISO");
       electron_ = scorerFTDTable("Offline/Mu2eG4/data/Electron_ambient_dose.dat","ISO");
       positron_ = scorerFTDTable("Offline/Mu2eG4/data/Positron_ambient_dose.dat","ISO");
       muminus_  = scorerFTDTable("Offline/Mu2eG4/data/Muon_ambient_dose.dat",    "ISO");
       muplus_   = scorerFTDTable("Offline/Mu2eG4/data/Muon_ambient_dose.dat",    "ISO");
       piminus_  = scorerFTDTable("Offline/Mu2eG4/data/Piminus_ambient_dose.dat", "ISO");
       piplus_   = scorerFTDTable("Offline/Mu2eG4/data/Piplus_ambient_dose.dat",  "ISO");
       proton_   = scorerFTDTable("Offline/Mu2eG4/data/Proton_ambient_dose.dat",  "ISO");
       neutron_  = scorerFTDTable("Offline/Mu2eG4/data/Neutron_ambient_dose.dat", "ISO");
    } else {
       throw cet::exception("INIT")<<"scorerFTDConverter unknown type\n";
    }
  }



  void scorerFTDConverter::print()
  {
    photon_.print();
    electron_.print();
    positron_.print();
    muplus_.print();
    muminus_.print();
    piplus_.print();
    piminus_.print();
    proton_.print();
    neutron_.print();
  }

  double scorerFTDConverter::evaluate(int pdgCode, double energy)
  {
    switch(pdgCode) {
      case PDGCode::gamma:
        return photon_.evaluate(energy);
      case PDGCode::e_minus:
        return electron_.evaluate(energy);
      case PDGCode::e_plus:
        return positron_.evaluate(energy);
      case PDGCode::mu_minus:
        return muminus_.evaluate(energy);
      case PDGCode::mu_plus:
        return muplus_.evaluate(energy);
      case PDGCode::pi_minus:
        return piminus_.evaluate(energy);
      case PDGCode::pi_plus:
        return piplus_.evaluate(energy);
      case PDGCode::n0:
        return neutron_.evaluate(energy);
      case PDGCode::proton:
        return proton_.evaluate(energy);
      default:
        //use proton as proxy for helium and heavier nuclei
        return proton_.evaluate(energy);
    }
 }

}
