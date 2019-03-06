//
//  Collection of tools useful for dealing with various tracking functions
//
#ifndef TrkReco_TrkPrintUtils_HH
#define TrkReco_TrkPrintUtils_HH

#include "fhiclcpp/ParameterSet.h"

#include "BTrk/KalmanTrack/KalRep.hh"
// #include "Mu2eUtilities/inc/McUtilsToolBase.hh"

#include <string>

namespace art {
  class Event;
}

namespace mu2e {

  class McUtilsToolBase;

  class TrkPrintUtils {
  protected:
    int                              _mcTruth;
    std::unique_ptr<McUtilsToolBase> _mcUtils;
    std::string                      _strawHitCollTag;
    
  public:
    TrkPrintUtils(const fhicl::ParameterSet& PSet);
    ~TrkPrintUtils();
//-----------------------------------------------------------------------------
// Option = ""              : print banner + track parameters (but not hits)
//        includes "banner" : print banner
//        includes "data"   : print track parameters
//        includes "hits    : print hits
// if Message is not empty, it is also printed
//-----------------------------------------------------------------------------
    void       printTrack(const art::Event* event, const KalRep* Krep, const char* Option = "", const char* Message = "");
  };
}

#endif
