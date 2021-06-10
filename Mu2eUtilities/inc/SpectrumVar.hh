#ifndef Mu2eUtilities_inc_SpectrumVar_hh
#define Mu2eUtilities_inc_SpectrumVar_hh

#include <string>

namespace mu2e {

  enum SpectrumVar  { TOTAL_ENERGY, KINETIC_ENERY, MOMENTUM };

  inline SpectrumVar    parseSpectrumVar(const std::string& name) {
    if (name == "totalEnergy"  )  return TOTAL_ENERGY;
    if (name == "kineticEnergy")  return KINETIC_ENERY;
    if (name == "momentum"     )  return MOMENTUM;
    throw cet::exception("BADCONFIG")<<"parseSpectrumVar(): unknown spectrum variable "<<name<<"\n";
  }
}

#endif/*Mu2eUtilities_inc_SpectrumVar_hh*/
