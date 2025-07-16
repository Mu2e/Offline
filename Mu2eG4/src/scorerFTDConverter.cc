#include "Offline/Mu2eG4/inc/scorerFTDConverter.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "cetlib_except/exception.h"

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>


namespace mu2e{

  scorerFTDConverter::scorerFTDConverter(const std::string& method) :
    photon  ("Offline/Mu2eG4/data/Photon_effective_dose.dat",  method),
    electron("Offline/Mu2eG4/data/Electron_effective_dose.dat",method),
    positron("Offline/Mu2eG4/data/Positron_effective_dose.dat",method),
    muminus ("Offline/Mu2eG4/data/Muminus_effective_dose.dat", method),
    muplus  ("Offline/Mu2eG4/data/Muplus_effective_dose.dat",  method),
    piminus ("Offline/Mu2eG4/data/Piminus_effective_dose.dat", method),
    piplus  ("Offline/Mu2eG4/data/Piplus_effective_dose.dat",  method),
    proton  ("Offline/Mu2eG4/data/Proton_effective_dose.dat",  method),
    neutron ("Offline/Mu2eG4/data/Neutron_effective_dose.dat", method),
    helium  ("Offline/Mu2eG4/data/Helium_effective_dose.dat",  method)
  {}


  void scorerFTDConverter::print()
  {
    photon.print();
    electron.print();
    positron.print();
    muplus.print();
    muminus.print();
    piplus.print();
    piminus.print();
    proton.print();
    neutron.print();
    helium.print();
  }

  double scorerFTDConverter::evaluate(int pdgCode, double energy)
  {
    switch(pdgCode) {
      case 22:
        return photon.evaluate(energy);
      case 11:
        return positron.evaluate(energy);
      case -11:
        return electron.evaluate(energy);
      case -13:
        return muminus.evaluate(energy);
      case 13:
        return muplus.evaluate(energy);
      case -211:
        return piminus.evaluate(energy);
      case 211:
        return piplus.evaluate(energy);
      case 2112:
        return neutron.evaluate(energy);
      case 2212:
        return proton.evaluate(energy);
      default:
        //use proton as proxy for heavier nuclei and other particles
        return proton.evaluate(energy);
    }
 }

}
