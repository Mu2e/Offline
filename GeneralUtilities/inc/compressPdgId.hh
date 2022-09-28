#ifndef Mu2eUtilities_compressPdgId_hh
#define Mu2eUtilities_compressPdgId_hh

#include "Offline/DataProducts/inc/PDGCode.hh"


namespace mu2e {

  // compress PDG e,mu,gamma, pions and kaons produced
  // in cosmic generators into a index 0-7
  int compressPdgIdCosmic(int pdgId);

}

#endif
