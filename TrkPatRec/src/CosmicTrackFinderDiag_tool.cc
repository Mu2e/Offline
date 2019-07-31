#ifndef __TrkPatRec_CosmicTrackFinderDiag_hh__
#define __TrkPatRec_CosmicTrackFinderDiag_hh__

#include "TrkPatRec/inc/CosmicTrackFinder_types.hh"
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "TrkReco/inc/CosmicTrackFit.hh"

#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrkData/inc/TrkStrawHit.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

#include "TH1F.h"

namespace mu2e {

  using namespace CosmicTrackFinderTypes;
  
  class CosmicTrackFinderDiag : public mu2e::ModuleHistToolBase {
  public:

    enum {
      kNEventHistSets = 10,
      kNHelixHistSets = 10,
      kNHitHistSets   = 10
    };

    
    struct Hist_t {
      TH1F*  nTimePeaks;
      TH1F*  nChPerPanel;
      TH1F*  nChHits;
      TH1F*  ntclhits;
      TH1F*  nhits   ;           // number of hits on a helix  
      TH1F*  nseeds;
      TH1F*  niters;
      TH1F*  nShFit;
      TH1F*  nChFit;  
      TH1F*  nXYSh;
      TH1F*  chi2;
      TH1F*  pT;
      TH1F*  p; 
      TH1F*  chi2d_track;
      TH1F*  residualX;
      TH1F*  residualY;
      TH1F*  pullX;
      TH1F*  pullY;
      
    };

  protected:
    int                              _mcTruth;
    //    std::unique_ptr<McUtilsToolBase> _mcUtils;
    std::string                      _shDigiLabel;
    Hist_t                           _hist;            // owned
    Data_t*                          _data;            // cached

  public:

    CosmicTrackFinderDiag(const fhicl::ParameterSet& PSet);
    ~CosmicTrackFinderDiag();

  private:
    
    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1) override ;
  };


//-----------------------------------------------------------------------------
  CosmicTrackFinderDiag::CosmicTrackFinderDiag(const fhicl::ParameterSet& PSet) {
    printf(" CosmicTrackFinderDiag::CosmicTrackFinderDiag : HOORAY! \n");
  }

//-----------------------------------------------------------------------------
  CosmicTrackFinderDiag::~CosmicTrackFinderDiag() {
  }

    
//-----------------------------------------------------------------------------
  int CosmicTrackFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
  
    _hist.nTimePeaks    = Tfs->make<TH1F>("ntpeaks"  , "number of time peaks"                      , 51, -0.5, 5.5);
    _hist.nChPerPanel   = Tfs->make<TH1F>("nchppanel", "number of ComboHits per panel"             , 101, -0.5, 50.5);
    _hist.nChHits       = Tfs->make<TH1F>("nchhits"  , "number of ComboHits processed"             , 101, -0.5, 50.5);
    _hist.nseeds     = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events"    , 21, -0.5, 1000);
    
    _hist.ntclhits   = Tfs->make<TH1F>("ntclhits0" , "number of hits on a time peak - no delta"  , 201, -0.5, 50.5);
   
    
    _hist.nhits       = Tfs->make<TH1F>("nhits"    , "number of hits on a track candidate"       , 401, -0.5, 50.5);
    
    _hist.niters    = Tfs->make<TH1F>("niter"     , "number of iterations  fit"   , 401, -0.5, 800.5);

    _hist.nShFit    = Tfs->make<TH1F>("nShFit"     , "number of strawhits after fit "   , 101, -0.5, 100.5);
   
    _hist.nChFit    = Tfs->make<TH1F>("nChFit"     , "number of combo hits after fit,  "   , 101, -0.5, 50.5);

    _hist.nXYSh     = Tfs->make<TH1F>("nXYSh"     , "number of strawHits from the circle fit,  "   , 201, -0.5, 200.5);
    
    _hist.pT         = Tfs->make<TH1F>("pT0"      , "transverse momentum; pT [MeV/c]"           , 400, -0.5, 200.5);
    _hist.p          = Tfs->make<TH1F>("p0"       , "momentum; p [MeV/c]"                       , 400, -0.5, 200.5);
    
    _hist.chi2d_track= Tfs->make<TH1F>("chi2dtrack" , "global chi2d; #chi^{2}/ndof"                   , 50, 0, 50); 
   
    _hist.residualX     = Tfs->make<TH1F>("resid_X" , "track hit X dr, "               , 50,   -500., 500.);
    _hist.residualY     = Tfs->make<TH1F>("resid_Y" , "helix hit dr, "               , 50,   -500., 500.);
    _hist.pullX     = Tfs->make<TH1F>("pull_X" , "track hit X dr/error, "               , 50,   -10., 10.);
    _hist.pullY     = Tfs->make<TH1F>("pull_Y" , "track hit dr/error, "               , 50,   -10., 10.);
    
    return 0;
  }


//-----------------------------------------------------------------------------
// Mode is not used
//-----------------------------------------------------------------------------
  int CosmicTrackFinderDiag::fillHistograms(void* Data, int Mode) {
    _data = (Data_t*) Data;
//-----------------------------------------------------------------------------
// fill track-level histograms
//-----------------------------------------------------------------------------    
    _hist.nTimePeaks ->Fill(_data->nTimePeaks);
    _hist.nseeds->Fill(_data->nseeds);
    for (int i=0; i<_data->nseeds; i++) {
	_hist.nChPerPanel->Fill(_data->nChPPanel[i]);
	_hist.nChHits    ->Fill(_data->nChHits[i]);
	_hist.ntclhits   ->Fill(_data->ntclhits[i]);
	_hist.nhits      ->Fill(_data->nhits[i]);	
	_hist.nShFit    ->Fill(_data->nShFit[i]);
	_hist.nChFit    ->Fill(_data->nChFit[i]);
	_hist.nXYSh     ->Fill(_data->nsh[i]);       	
	//_hist.p       ->Fill(_data->p[i]);
	//_hist.pT      ->Fill(_data->pT[i]);	
	_hist.chi2d_track->Fill(_data->chi2d_track[i]);
	for (int j=0; j<_data->nChFit[i];++j){
	  _hist.residualX->Fill(_data->hit_residualX[i][j]);
	  _hist.residualY->Fill(_data->hit_residualX[i][j]);
	  _hist.pullX->Fill(_data->hit_residualX[i][j]/_data->hit_errorX[i][j]);
	  _hist.pullY->Fill(_data->hit_residualY[i][j]/_data->hit_errorY[i][j]);
	}//end CH loop
	
      }//end seed loop
    

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(mu2e::CosmicTrackFinderDiag)

#endif
