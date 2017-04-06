#ifndef __CalPatRec_CalTrkFit_hist_hh__
#define __CalPatRec_CalTrkFit_hist_hh__

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "TH1.h"

namespace mu2e {

class CalTrkFit_hist {
public:

  std::string                           _shDigiLabel;
  const PtrStepPointMCVectorCollection* _listOfMCStrawHits;
  static CalTrkFit_hist*                _instance;

  int (*_NGenHits) (const art::Event*, fhicl::ParameterSet*, const char*, const StrawHitCollection*);
private:
  CalTrkFit_hist();
public:
  ~CalTrkFit_hist();
  
  static int book(art::ServiceHandle<art::TFileService> & Tfs, TObject* Hist);
  static int fill(const TObject* Data, TObject* Hist);

  static int nGenHits(const art::Event*, fhicl::ParameterSet*, const char*, const StrawHitCollection*) {
    return -1;
  }

  static CalTrkFit_hist* instance();

  friend class Cleaner;

  class Cleaner {
  public:
    Cleaner() {};
    ~Cleaner() { if (_instance) delete _instance; }
  };

};

}
#endif
