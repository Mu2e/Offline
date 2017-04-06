#ifndef __CalPatRec_CalTrkFit_types_hh__
#define __CalPatRec_CalTrkFit_types_hh__

#include "TObject.h"
#include "TH1.h"

#include <vector>

namespace art {
  class Event;
};

namespace mu2e {
  
  class KalFitResultNew;
  class TTracker;
  class Doublet;

  struct CalTrkFit_Data_t : public TObject {
    const art::Event*     event;
    KalFitResultNew*      result;
    fhicl::ParameterSet*  timeOffsets;
    const TTracker*       tracker;
    std::vector<Doublet>* listOfDoublets;

    int                   ntracks[2];
  };

  struct CalTrkFit_Hist_t : public TObject {
    TH1F*  nhits;
    TH1F*  chi2   [2];
    TH1F*  p      [2];
    TH1F*  kaldoca[2];
    TH1F*  ntracks[2];

    TH1F*  dSlopeOS[2];
    TH1F*  dSlopeSS[2];
    TH1F*  chi2bOS [2];
    TH1F*  chi2bSS [2];
    TH1F*  chi2rOS [2];
    TH1F*  chi2rSS [2];
  };

}
#endif
