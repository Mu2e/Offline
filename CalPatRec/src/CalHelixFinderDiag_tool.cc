#ifndef __CalPatRec_CalHelixFinderDiag_hh__
#define __CalPatRec_CalHelixFinderDiag_hh__

#include "TH1.h"
#include "TH2.h"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/CalHelixFinder_types.hh"
#include "CalPatRec/inc/ModuleHistToolBase.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

namespace mu2e {

  using namespace CalHelixFinderTypes;
  
  class CalHelixFinderDiag : public mu2e::ModuleHistToolBase {

  protected:
    const PtrStepPointMCVectorCollection* _listOfMCStrawHits;
    Hist_t                     _hist;
    Data_t*                    _data;
    int                        _first_call;
    int                        _event_number;
    SimParticleTimeOffset*     _timeOffsets;
    
  public:
    
    CalHelixFinderDiag(const fhicl::ParameterSet& PSet);
    ~CalHelixFinderDiag();
    
  private:
    
    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1) override ;
  };

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
CalHelixFinderDiag::CalHelixFinderDiag(const fhicl::ParameterSet& PSet) {
  printf(" CalHelixFinderDiag::CalHelixFinderDiag : HOORAY! \n");
  _first_call   = 1;
  _event_number = -1;
  _timeOffsets  = NULL;
}


CalHelixFinderDiag::~CalHelixFinderDiag() {
}
  
//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
int CalHelixFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {

  _hist.nseeds[0]     = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events"    , 21, -0.5, 20.5);
  _hist.nseeds[1]     = Tfs->make<TH1F>("nseeds1"  , "number of track candidates: nhits > 15"    , 21, -0.5, 20.5);
  _hist.nhits         = Tfs->make<TH1F>("nhits"    , "number of hits on a track candidate"       , 101, -0.5, 100.5);
  _hist.radius[0]     = Tfs->make<TH1F>("radius0"  , "helix radius; r [mm]"                      , 401, -0.5, 400.5);
  _hist.radius[1]     = Tfs->make<TH1F>("radius1"  , "helix radius nhits > 15; r [mm]"           , 401, -0.5, 400.5);
  _hist.pT [0]        = Tfs->make<TH1F>("pT0"      , "transverse momentum; pT [MeV/c]"           , 400, -0.5, 200.5);
  _hist.p  [0]        = Tfs->make<TH1F>("p0"       , "momentum; p [MeV/c]"                       , 400, -0.5, 200.5);
  _hist.pT [1]        = Tfs->make<TH1F>("pT1"      , "transverse momentum nhits > 15; pT [MeV/c]", 400, -0.5, 200.5);
  _hist.p  [1]        = Tfs->make<TH1F>("p1"       , "momentum nhits > 15; p [MeV/c]"            , 400, -0.5, 200.5);
  _hist.chi2XY[0]     = Tfs->make<TH1F>("chi2XY0"  , "normalized chi2-XY"                        , 200, 0., 20.);
  _hist.chi2XY[1]     = Tfs->make<TH1F>("chi2XY1"  , "normalized chi2-XY: nhits>15"              , 200, 0., 20.);
  _hist.chi2ZPhi[0]   = Tfs->make<TH1F>("chi2ZPhi0", "normalized chi2-ZPhi"                      , 200, 0., 20.);
  _hist.chi2ZPhi[1]   = Tfs->make<TH1F>("chi2ZPhi1", "normalized chi2-ZPhi: nhits>15"            , 200, 0., 20.);
  _hist.nhitsvspT     = Tfs->make<TH2F>("nhitsvspT", "nhits vs pT"                               , 100, 0, 100, 400, 0, 200);
  _hist.nhitsvsp      = Tfs->make<TH2F>("nhitsvsp" , "nhits vs p"                                , 100, 0, 100, 400, 0, 200);
  _hist.nStationPairs = Tfs->make<TH1F>("nstpairs" , "n-station pairs used in findDfDz"          , 200, 0, 200);
  return 0;
}


//-----------------------------------------------------------------------------
// called once per event, Mode is not used..
//-----------------------------------------------------------------------------
  int CalHelixFinderDiag::fillHistograms(void * Data, int Mode) {

    _data = (Data_t*) Data;

    int en = _data->event->event();
    if ( en == _event_number) return -1;

    if (_first_call == 1) {
      _timeOffsets = new SimParticleTimeOffset(*_data->timeOffsets);
      _first_call = 0;
    }
    
    _event_number = en;
    _timeOffsets->updateMap(*_data->event);
//-----------------------------------------------------------------------------
// this part is so far simple
//-----------------------------------------------------------------------------
    _hist.nseeds[0]->Fill(_data->nseeds[0]);
    _hist.nseeds[1]->Fill(_data->nseeds[1]);

    for (int i=0; i<_data->nseeds[0]; i++) {
      _hist.nhits->Fill(_data->nhits[i]);
      
      _hist.p[0]->Fill(_data->p[i]);
      _hist.pT[0]->Fill(_data->pT[i]);
      _hist.radius[0]->Fill(_data->radius[i]);
      _hist.chi2XY[0]->Fill(_data->chi2XY[i]);
      _hist.chi2ZPhi[0]->Fill(_data->chi2ZPhi[i]);

      if (_data->good[i] != 0) {
	_hist.p[1]->Fill(_data->p[i]);
	_hist.pT[1]->Fill(_data->pT[i]);
	_hist.radius[1]->Fill(_data->radius[i]);
	_hist.chi2XY[1]->Fill(_data->chi2XY[i]);
	_hist.chi2ZPhi[1]->Fill(_data->chi2ZPhi[i]);
      }

      _hist.nhitsvspT ->Fill(_data->nhits[i],_data->pT[i]);
      _hist.nhitsvsp  ->Fill(_data->nhits[i],_data->p [i]);

      _hist.nStationPairs->Fill(_data->nStationPairs[i]);
    }

    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::CalHelixFinderDiag)

#endif
