//

#include "TrkPatRec/inc/RobustHelixFinder_types.hh"
#include "CalPatRec/inc/ModuleHistToolBase.hh"
#include "TrkReco/inc/RobustHelixFit.hh"

#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrkData/inc/TrkStrawHit.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

#include "TH1F.h"

namespace mu2e {

  using namespace RobustHelixFinderTypes;
  
  class RobustHelixFinderDiag : public mu2e::ModuleHistToolBase {
  public:

    enum {
      kNEventHistSets = 10,
      kNHelixHistSets = 10,
      kNHitHistSets   = 10
    };

    
    struct Hist_t {
      TH1F*  nTimePeaks;
      TH1F*  ntclhits[2];
      TH1F*  nhits;           // number of hits on a helix  
      TH1F*  nseeds  [2];
      TH1F*  radius  [2];   
      TH1F*  chi2XY  [2];
      TH1F*  chi2ZPhi[2];
      TH1F*  pT      [2];
      TH1F*  p       [2];
      TH1F*  dr[2];
      TH1F*  chi2d_helix[2];
    };

  protected:
    int                              _mcTruth;
    std::unique_ptr<McUtilsToolBase> _mcUtils;
    std::string                      _shDigiLabel;
    Hist_t                           _hist;            // owned
    Data_t*                          _data;            // cached

    //    const        mu2e::PtrStepPointMCVectorCollection* _listOfMCStrawHits;

  public:

    RobustHelixFinderDiag(const fhicl::ParameterSet& PSet);
    ~RobustHelixFinderDiag();

  private:
    
    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1) override ;
  };


//-----------------------------------------------------------------------------
  RobustHelixFinderDiag::RobustHelixFinderDiag(const fhicl::ParameterSet& PSet) {
    printf(" RobustHelixFinderDiag::RobustHelixFinderDiag : HOORAY! \n");
  }

//-----------------------------------------------------------------------------
  RobustHelixFinderDiag::~RobustHelixFinderDiag() {
  }

    
//-----------------------------------------------------------------------------
  int RobustHelixFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
  
    _hist.nTimePeaks    = Tfs->make<TH1F>("ntpeaks"  , "number of time peaks"                      , 11, -0.5, 10.5);
    _hist.nseeds[0]     = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events"    , 21, -0.5, 20.5);
    _hist.nseeds[1]     = Tfs->make<TH1F>("nseeds1"  , "number of track candidates: nhits > 15;"    , 21, -0.5, 20.5);
    _hist.ntclhits[0]   = Tfs->make<TH1F>("ntclhits0" , "number of hits on a time peak - no delta"  , 101, -0.5, 100.5);
    _hist.ntclhits[1]   = Tfs->make<TH1F>("ntclhits1" , "number of hits on a time peak - no delta: nhits > 15"  , 101, -0.5, 100.5);
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
    _hist.dr  [0]       = Tfs->make<TH1F>("dr0"      , "dr; r - r_{no-target} [mm]"                , 800, -200, 200);
    _hist.dr  [1]       = Tfs->make<TH1F>("dr1"      , "dr: nhits>15; r - r_{no-target} [mm]"      , 800, -200, 200);
    _hist.shmeanr  [0]  = Tfs->make<TH1F>("shmeanr0" , "straw hit mean radius; r_{sh} [mm]"          , 1800, 0, 900);
    _hist.shmeanr  [1]  = Tfs->make<TH1F>("shmeanr1" , "straw hit mean radius: nhits>15; r_{sh} [mm]", 1800, 0, 900);
    _hist.chi2d_helix[0]= Tfs->make<TH1F>("chi2dhel0" , "global chi2d; #chi^{2}/ndof"                   , 100, 0, 10); 
    _hist.chi2d_helix[1]= Tfs->make<TH1F>("chi2dhel1" , "global chi2d: nhits>15; #chi^{2}/ndof"         , 100, 0, 10); 

    return 0;
  }

xs
//-----------------------------------------------------------------------------
// Mode is not used
//-----------------------------------------------------------------------------
  int RobustHelixFinderDiag::fillHistograms(void* Data, int Mode) {

    _data = (Data_t*) Data;

//-----------------------------------------------------------------------------
// fill helix-level histograms
//-----------------------------------------------------------------------------
    
    _hist.nTimePeaks->Fill(_data->nTimePeaks);

    for (int k=0; k<nhelicities; ++k){
      _hist.nseeds[k]->Fill(_data->nseeds[k]);
      _hist.nseeds[k]->Fill(_data->nseeds[k]);

      for (int i=0; i<_data->nseeds[k]; i++) {
	_hist.ntclhits   [k]->Fill(_data->ntclhits[k][i]   );
	_hist.nhits      [k]->Fill(_data->nhits[k][i]      );
	_hist.p          [k]->Fill(_data->p[k][i]          );
	_hist.pT         [k]->Fill(_data->pT[k][i]         );
	_hist.radius     [k]->Fill(_data->radius[k][i]     );
	_hist.chi2XY     [k]->Fill(_data->chi2XY[k][i]     );
	_hist.chi2ZPhi   [k]->Fill(_data->chi2ZPhi[k][i]   );
	_hist.dr         [k]->Fill(_data->dr[k][i]         );
	_hist.chi2d_helix[k]->Fill(_data->chi2d_helix[k][i]);
      
      }
    }

    return 0;
  }

  DEFINE_ART_CLASS_TOOL(RobustHelixFinderDiag)

}
