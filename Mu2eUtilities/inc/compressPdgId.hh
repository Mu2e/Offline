#ifndef Mu2eUtilities_compressPdgId_hh
#define Mu2eUtilities_compressPdgId_hh

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/CompressedPDGCode.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

  // compress all PDGId's into 27 categories,
  // useful for histogramming (see also next)
  CompressedPDGCode::enum_type compressPDGCode(PDGCode::enum_type pdgId);

  // a histogram with text labels, ready for use with CompressedPDGCode
  TH1D* compressPDGCodeHisto(art::ServiceHandle<art::TFileService> tfs,
                             std::string name="compPdgId",
                             std::string title="Compressed PDG ID");

  // compress PDG e,mu,gamma, pions and kaons produced
  // in cosmic generators into a index 0-7
  int compressPdgIdCosmic(int pdgId);

}

#endif
