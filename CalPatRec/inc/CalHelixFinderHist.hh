#ifndef __CalPatRec_CalHelixFinderHist_hh__
#define __CalPatRec_CalHelixFinderHist_hh__

#include "TH1.h"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/CalPatRec_utils.hh"
#include "CalPatRec/inc/CalHelixFinder_types.hh"
#include "CalPatRec/inc/CalHelixFinderHist.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "CalPatRec/inc/CalHelixFinderTypes.hh"

namespace mu2e {
  
  class CalHelixFinderHist : public ModuleHist {
  public:

    std::string                           _shDigiLabel;
    const PtrStepPointMCVectorCollection* _listOfMCStrawHits;
    
    CalHelixFinderHist();
    ~CalHelixFinderHist();
    
    int book(art::ServiceHandle<art::TFileService> & Tfs, CalHelixFinder_Hist_t* Hist);
    int fill(const CalHelixFinder_Data_t* Data          , CalHelixFinder_Hist_t* Hist);
  };




  // static XXX_hist XXX_CalHelixFinder("CalHelixFinder",CalHelixFinder_hist::book,CalHelixFinder_hist::fill);
//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
  int CalHelixFinderHist::book(art::ServiceHandle<art::TFileService>& Tfs, CalHelixFinder_Hist_t* Hist) {

    Hist->nseeds[0]   = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events", 21, -0.5, 20.5);
    Hist->nseeds[1]   = Tfs->make<TH1F>("nseeds1"  , "number of track candidates: nhits > 15", 21, -0.5, 20.5);
    Hist->nhits       = Tfs->make<TH1F>("nhits"    , "number of hits within a track candidate; nHits", 101, -0.5, 100.5);
    Hist->radius[0]   = Tfs->make<TH1F>("radius0"  , "helix radius; r [mm]"                  , 401, -0.5, 400.5);
    Hist->radius[1]   = Tfs->make<TH1F>("radius1"  , "helix radius nhits > 15; r [mm]"       , 401, -0.5, 400.5);
    Hist->pT [0]      = Tfs->make<TH1F>("pT0"      , "transverse momentum; pT [MeV/c]"       , 400, -0.5, 200.5);
    Hist->p  [0]      = Tfs->make<TH1F>("p0"       , "momentum; p [MeV/c]"                   , 400, -0.5, 200.5);
    Hist->pT [1]      = Tfs->make<TH1F>("pT1"      , "transverse momentum nhits > 15; pT [MeV/c]"       , 400, -0.5, 200.5);
    Hist->p  [1]      = Tfs->make<TH1F>("p1"       , "momentum nhits > 15; p [MeV/c]"                   , 400, -0.5, 200.5);
    Hist->chi2XY[0]   = Tfs->make<TH1F>("chi2XY0"  , "normalized chi2-XY"                   , 200, 0., 20.);
    Hist->chi2XY[1]   = Tfs->make<TH1F>("chi2XY1"  , "normalized chi2-XY: nhits>15"         , 200, 0., 20.);
    Hist->chi2ZPhi[0] = Tfs->make<TH1F>("chi2ZPhi0", "normalized chi2-ZPhi"             , 200, 0., 20.);
    Hist->chi2ZPhi[1] = Tfs->make<TH1F>("chi2ZPhi1", "normalized chi2-ZPhi: nhits>15"   , 200, 0., 20.);
    Hist->nhitsvspT   = Tfs->make<TH2F>("nhitsvspT", "nhits vs pT", 100, 0, 100, 400, 0, 200);
    Hist->nhitsvsp    = Tfs->make<TH2F>("nhitsvsp" , "nhits vs p" , 100, 0, 100, 400, 0, 200);

    return 0;
  }


  //-----------------------------------------------------------------------------
  int CalHelixFinder_hist::fill(const CalHelixFinder_Data_t* Data, CalHelixFinder_Hist_t* Hist) {

    static int                              first_call(1);
    static int                              event_number(-1);
    static SimParticleTimeOffset*           timeOffsets(NULL);
    //    static PtrStepPointMCVectorCollection*  listOfMCStrawHits;

    if (first_call == 1) {
      timeOffsets = new SimParticleTimeOffset(*Data->timeOffsets);
      first_call  = 0;
    }
//-----------------------------------------------------------------------------
// update if new event
//-----------------------------------------------------------------------------
    int en = Data->event->event();
    if ( en != event_number) {
      event_number = en;

      timeOffsets->updateMap(*Data->event);

      // art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandle;
      // Data->event->getByLabel(Data->shLabel,"",mcptrHandle);
      // if (mcptrHandle.isValid()) {
      // 	listOfMCStrawHits = (mu2e::PtrStepPointMCVectorCollection*) mcptrHandle.product();
      // }
      // else {
      // 	listOfMCStrawHits = NULL;
      // }
    }

    Hist->nseeds[0]->Fill(Data->nseeds[0]);
    Hist->nseeds[1]->Fill(Data->nseeds[1]);

    for (int i=0; i<Data->nseeds[0]; i++) {
      Hist->nhits->Fill(Data->nhits[i]);
      
      Hist->p[0]->Fill(Data->p[i]);
      Hist->pT[0]->Fill(Data->pT[i]);
      Hist->radius[0]->Fill(Data->radius[i]);
      Hist->chi2XY[0]->Fill(Data->chi2XY[i]);
      Hist->chi2ZPhi[0]->Fill(Data->chi2ZPhi[i]);

      if (Data->good[i] != 0) {
	Hist->p[1]->Fill(Data->p[i]);
	Hist->pT[1]->Fill(Data->pT[i]);
	Hist->radius[1]->Fill(Data->radius[i]);
	Hist->chi2XY[1]->Fill(Data->chi2XY[i]);
	Hist->chi2ZPhi[1]->Fill(Data->chi2ZPhi[i]);
      }

      Hist->nhitsvspT ->Fill(Data->nhits[i],Data->pT[i]);
      Hist->nhitsvsp  ->Fill(Data->nhits[i],Data->p [i]);
    }

    return 0;
  }

}
#endif
