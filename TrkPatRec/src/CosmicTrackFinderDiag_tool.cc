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
    
    struct Hist_t {
      //General Track info:
      TH1F*  nTimePeaks;
      //TH1F*  nChPerPanel;
      //TH1F*  nChHits;
      TH1F*  ntclhits;
      TH1F*  nhits   ;           // number of hits on a helix  
      TH1F*  nseeds;
      TH1F*  niters;
      TH1F*  nShFit;
      TH1F*  nChFit;  
      //TH1F*  nXYSh;
      TH1F*  chi2;
      TH1F*  pT;
      TH1F*  p; 
      //Chi2:
      TH1F*  Final_chi2dX_track;
      TH1F*  Initial_chi2dX_track;
      TH1F*  Final_chi2dY_track;
      TH1F*  Initial_chi2dY_track;
      TH1F*  Final_chi2d_track;
      TH1F*  Initial_chi2d_track;
      //Residuals and Pulls:
      TH1F*  Final_residualX;
      TH1F*  Final_residualY;
      TH1F*  Final_pullX;
      TH1F*  Final_pullY;
      TH1F*  Initial_residualX;
      TH1F*  Initial_residualY;
      TH1F*  Initial_pullX;
      TH1F*  Initial_pullY;
     
      
    };

  protected:
    //int                              _mcTruth;
    //    std::unique_ptr<McUtilsToolBase> _mcUtils;
    //std::string                      _shDigiLabel;
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
  
    _hist.nTimePeaks    = Tfs->make<TH1F>("ntpeaks"  , "number of time peaks"                      , 1000,0,1000);
   
    _hist.nseeds     = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events"    , 1000,0,1000);
    
    _hist.ntclhits   = Tfs->make<TH1F>("ntclhits0" , "number of hits on a time peak "  , 201, -0.5, 50.5);
   
    
    _hist.nhits       = Tfs->make<TH1F>("nhits"    , "number of hits on a track candidate"       , 100,0,100);
    
    _hist.niters    = Tfs->make<TH1F>("niter"     , "number of iterations  fit"   , 1000, 0, 1000);

    _hist.nShFit    = Tfs->make<TH1F>("nShFit"     , "number of strawhits after fit "   , 101, 0, 100);
   
    _hist.nChFit    = Tfs->make<TH1F>("nChFit"     , "number of combo hits after fit,  "   , 101, 0,100);

    _hist.pT         = Tfs->make<TH1F>("pT0"      , "transverse momentum; pT [MeV/c]"           , 400, -0.5, 200.5);
    _hist.p          = Tfs->make<TH1F>("p0"       , "momentum; p [MeV/c]"                       , 400, -0.5, 200.5);
    
    _hist.Final_chi2d_track = Tfs->make<TH1F>("Final_chi2dtrack" , "End total chi2", 50, 0, 50); 
    _hist.Final_chi2dX_track = Tfs->make<TH1F>("Final_chi2dtrackX" , "End X chi2"  , 50, 0, 50);
    _hist.Final_chi2dY_track = Tfs->make<TH1F>("Final_chi2dtrackY" , "End Y chi2"   , 50, 0, 50);
   
    _hist.Final_residualX  = Tfs->make<TH1F>("final_resid_X" , "Final Track Residual X, ", 100,   -100., 100.);
    _hist.Final_residualY   = Tfs->make<TH1F>("final_resid_Y" , "Final Track Residual Y, " , 100,   -100., 100.);
    _hist.Final_pullX  = Tfs->make<TH1F>("final_pull_X" , "Final Track Pull X"               , 50,   -10., 10.);
    _hist.Final_pullY  = Tfs->make<TH1F>("final_pull_Y" , "Final Track Pull Y"               , 50,   -10., 10.);
    
    _hist.Initial_chi2d_track = Tfs->make<TH1F>("Initial_chi2dtrack" , "End total chi2", 50, 0, 50); 
    _hist.Initial_chi2dX_track = Tfs->make<TH1F>("Initial_chi2dtrackX" , "End X chi2", 50, 0, 50); 
    _hist.Initial_chi2dY_track = Tfs->make<TH1F>("Initial_chi2dtrackY" , "End Y chi2"  , 50, 0, 50); 
    
    _hist.Initial_residualX  = Tfs->make<TH1F>("Initial_resid_X" , "Initial Track Residual X, " , 50,   -500., 500.);
    _hist.Initial_residualY  = Tfs->make<TH1F>("Initial_resid_Y" , "Initial Track Residual Y, "  , 50,   -500., 500.);
    _hist.Initial_pullX     = Tfs->make<TH1F>("Initial_pull_X" , "Initial Track Pull X" , 100,   -50., 50.);
    _hist.Initial_pullY     = Tfs->make<TH1F>("Initial_pull_Y" , "Initial Track Pull Y"  , 100,   -50., 50.);
    
    return 0;
  }

int CosmicTrackFinderDiag::fillHistograms(void* Data,  int Mode) {
    _data = (Data_t*) Data;
    _hist.nTimePeaks ->Fill(_data->nTimePeaks);
    _hist.nseeds->Fill(_data->nseeds);
   //for (int i=0; i<_data->nseeds; i++) {
	
	_hist.ntclhits   ->Fill(_data->ntclhits);
	_hist.nhits      ->Fill(_data->nhits);	
	_hist.nShFit    ->Fill(_data->nShFit);
	_hist.nChFit    ->Fill(_data->nChFit);
	  	
	//Diagnostics	
	_hist.Final_chi2d_track->Fill(_data->Final_chi2d_track);
	_hist.Initial_chi2d_track->Fill(_data->Final_chi2d_track);
	_hist.Final_chi2dX_track->Fill(_data->Final_chi2dX_track);
	_hist.Initial_chi2dX_track->Fill(_data->Final_chi2dX_track);
	_hist.Final_chi2d_track->Fill(_data->Final_chi2dY_track);
	_hist.Initial_chi2d_track->Fill(_data->Final_chi2dY_track);
	
	for (int j=0; j<_data->nChFit;++j){
	  _hist.Final_residualX->Fill(_data->Final_hit_residualX[j]);
	  _hist.Final_residualY->Fill(_data->Final_hit_residualY[j]);
	  _hist.Final_pullX->Fill(_data->Final_hit_residualX[j]/_data->Final_hit_errorX[j]);
	  _hist.Final_pullY->Fill(_data->Final_hit_residualY[j]/_data->Final_hit_errorY[j]);
	  _hist.Initial_residualX->Fill(_data->Initial_hit_residualX[j]);
	  _hist.Initial_residualY->Fill(_data->Initial_hit_residualX[j]);
	  _hist.Initial_pullX->Fill(_data->Initial_hit_residualX[j]/_data->Initial_hit_errorX[j]);
	  _hist.Initial_pullY->Fill(_data->Initial_hit_residualY[j]/_data->Initial_hit_errorY[j]);
	//}//end CH loop
	
      }//end seed loop
   
    return 0;
  }
  
}

DEFINE_ART_CLASS_TOOL(mu2e::CosmicTrackFinderDiag)

#endif
