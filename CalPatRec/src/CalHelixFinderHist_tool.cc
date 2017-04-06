#ifndef __CalPatRec_CalHelixFinderHist_hh__
#define __CalPatRec_CalHelixFinderHist_hh__

#include "TH1.h"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/CalHelixFinder_types.hh"
#include "CalPatRec/inc/CalHelixFinderHist.hh"
#include "CalPatRec/inc/CalHelixFinder_types.hh"
#include "CalPatRec/inc/CprModuleHistBase.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

namespace mu2e {

class CalHelixFinderHist : public mu2e::CprModuleHistBase {
protected:
  std::string  _shDigiLabel;
  const        mu2e::PtrStepPointMCVectorCollection* _listOfMCStrawHits;
  
public:

  CalHelixFinderHist(const fhicl::ParameterSet& PSet);
  ~CalHelixFinderHist();

private:

  virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs, TObject* Hist) override ;
  virtual int fillHistograms(int Mode, const TObject* Data, TObject* Hist) override ;
};


//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------

CalHelixFinderHist::CalHelixFinderHist(const fhicl::ParameterSet& PSet) {
  printf(" CalHelixFinderHist::CalHelixFinderHist : HOORAY! \n");
}

CalHelixFinderHist::~CalHelixFinderHist() {
}
  
//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
int CalHelixFinderHist::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs, TObject* Hist) {

  mu2e::CalHelixFinder_Hist_t* hist = ( mu2e::CalHelixFinder_Hist_t*) Hist;

  hist->nseeds[0]   = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events", 21, -0.5, 20.5);
  hist->nseeds[1]   = Tfs->make<TH1F>("nseeds1"  , "number of track candidates: nhits > 15", 21, -0.5, 20.5);
  hist->nhits       = Tfs->make<TH1F>("nhits"    , "number of hits within a track candidate; nHits", 101, -0.5, 100.5);
  hist->radius[0]   = Tfs->make<TH1F>("radius0"  , "helix radius; r [mm]"                  , 401, -0.5, 400.5);
  hist->radius[1]   = Tfs->make<TH1F>("radius1"  , "helix radius nhits > 15; r [mm]"       , 401, -0.5, 400.5);
  hist->pT [0]      = Tfs->make<TH1F>("pT0"      , "transverse momentum; pT [MeV/c]"       , 400, -0.5, 200.5);
  hist->p  [0]      = Tfs->make<TH1F>("p0"       , "momentum; p [MeV/c]"                   , 400, -0.5, 200.5);
  hist->pT [1]      = Tfs->make<TH1F>("pT1"      , "transverse momentum nhits > 15; pT [MeV/c]"       , 400, -0.5, 200.5);
  hist->p  [1]      = Tfs->make<TH1F>("p1"       , "momentum nhits > 15; p [MeV/c]"                   , 400, -0.5, 200.5);
  hist->chi2XY[0]   = Tfs->make<TH1F>("chi2XY0"  , "normalized chi2-XY"                   , 200, 0., 20.);
  hist->chi2XY[1]   = Tfs->make<TH1F>("chi2XY1"  , "normalized chi2-XY: nhits>15"         , 200, 0., 20.);
  hist->chi2ZPhi[0] = Tfs->make<TH1F>("chi2ZPhi0", "normalized chi2-ZPhi"             , 200, 0., 20.);
  hist->chi2ZPhi[1] = Tfs->make<TH1F>("chi2ZPhi1", "normalized chi2-ZPhi: nhits>15"   , 200, 0., 20.);
  hist->nhitsvspT   = Tfs->make<TH2F>("nhitsvspT", "nhits vs pT", 100, 0, 100, 400, 0, 200);
  hist->nhitsvsp    = Tfs->make<TH2F>("nhitsvsp" , "nhits vs p" , 100, 0, 100, 400, 0, 200);

  return 0;
}


//-----------------------------------------------------------------------------
  int CalHelixFinderHist::fillHistograms(int Mode, const TObject* Data, TObject* Hist) {

  static int                              first_call(1);
  static int                              event_number(-1);
  static  mu2e::SimParticleTimeOffset*    timeOffsets(NULL);
  //    static PtrStepPointMCVectorCollection*  listOfMCStrawHits;

  mu2e::CalHelixFinder_Hist_t* hist = ( mu2e::CalHelixFinder_Hist_t*) Hist;
  mu2e::CalHelixFinder_Data_t* data = ( mu2e::CalHelixFinder_Data_t*) Data;

  if (first_call == 1) {
    timeOffsets = new  mu2e::SimParticleTimeOffset(*data->timeOffsets);
    first_call  = 0;
  }
  //-----------------------------------------------------------------------------
  // update if new event
  //-----------------------------------------------------------------------------
  int en = data->event->event();
  if ( en != event_number) {
    event_number = en;

    timeOffsets->updateMap(*data->event);

    // art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandle;
    // Data->event->getByLabel(Data->shLabel,"",mcptrHandle);
    // if (mcptrHandle.isValid()) {
    // 	listOfMCStrawHits = (mu2e::PtrStepPointMCVectorCollection*) mcptrHandle.product();
    // }
    // else {
    // 	listOfMCStrawHits = NULL;
    // }
  }

  hist->nseeds[0]->Fill(data->nseeds[0]);
  hist->nseeds[1]->Fill(data->nseeds[1]);

  for (int i=0; i<data->nseeds[0]; i++) {
    hist->nhits->Fill(data->nhits[i]);
      
    hist->p[0]->Fill(data->p[i]);
    hist->pT[0]->Fill(data->pT[i]);
    hist->radius[0]->Fill(data->radius[i]);
    hist->chi2XY[0]->Fill(data->chi2XY[i]);
    hist->chi2ZPhi[0]->Fill(data->chi2ZPhi[i]);

    if (data->good[i] != 0) {
      hist->p[1]->Fill(data->p[i]);
      hist->pT[1]->Fill(data->pT[i]);
      hist->radius[1]->Fill(data->radius[i]);
      hist->chi2XY[1]->Fill(data->chi2XY[i]);
      hist->chi2ZPhi[1]->Fill(data->chi2ZPhi[i]);
    }

    hist->nhitsvspT ->Fill(data->nhits[i],data->pT[i]);
    hist->nhitsvsp  ->Fill(data->nhits[i],data->p [i]);
  }

  return 0;
}

}

DEFINE_ART_CLASS_TOOL(mu2e::CalHelixFinderHist)

#endif
