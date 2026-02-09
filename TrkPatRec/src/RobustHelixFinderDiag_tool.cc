#include "Offline/TrkPatRec/inc/RobustHelixFinder_types.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/TrkReco/inc/RobustHelixFit.hh"


#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "TH1F.h"

namespace mu2e {

  using namespace RobustHelixFinderTypes;
  using RobustHelixFinderTypes::Config;

  class RobustHelixFinderDiag : public mu2e::ModuleHistToolBase {
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
        TH1F*  ntclhits[2];
        TH1F*  nhits   [2];           // number of hits on a helix
        TH1F*  nseeds  [2];
        TH1F*  ntripl0 ;
        TH1F*  ntripl1 [2];
        TH1F*  xyniter [2];
        TH1F*  fzniter [2];
        TH1F*  niter   [2];

        TH1F*  nShFitXY    [2];
        TH1F*  nChFitXY    [2];
        TH1F*  nShFitCircle[2];
        TH1F*  nChFitCircle[2];

        TH1F*  nrescuedhits[2];

        TH1F*  nXYSh   [2];
        TH1F*  nZPhiSh [2];

        TH1F*  nfz0counter[2];
        TH1F*  nfz0counter_0[2];

        TH1F*  nshsxy_0    [2];
        TH1F*  rsxy_0      [2];
        TH1F*  chi2dsxy_0  [2];

        TH1F*   nshsxy_1    [2];
        TH1F*   rsxy_1      [2];
        TH1F*   chi2dsxy_1  [2];

        TH1F*   nshszphi_0   [2];
        TH1F*   lambdaszphi_0[2];
        TH1F*   chi2dszphi_0 [2];

        TH1F*   nshszphi_1   [2];
        TH1F*   lambdaszphi_1[2];
        TH1F*   chi2dszphi_1 [2];


        TH1F*  lambda0 [2];
        TH1F*  lambda1 [2];

        TH1F*  rinit   [2];
        TH1F*  radius  [2];
        TH1F*  chi2XY  [2];
        TH1F*  chi2ZPhi[2];

        TH1F*  pT      [2];
        TH1F*  p       [2];
        TH1F*  dr[2];
        TH1F*  chi2d_helix[2];

        TH1F*  resid   [2];
        TH1F*  rwdot   [2];

        TH1F*  nLoops  [2];
        TH1F*  nHitsLoopFailed[2];
        TH1F*  eDepAvg[2];
        TH1F*  chRadialDist[2];
      };

    protected:
      int                              _mcTruth;
      //    std::unique_ptr<McUtilsToolBase> _mcUtils;
      std::string                      _shDigiLabel;
      Hist_t                           _hist;            // owned
      Data_t*                          _data;            // cached

      //    const        mu2e::PtrStepPointMCVectorCollection* _listOfMCStrawHits;

    public:

      explicit RobustHelixFinderDiag(const Config& config);
      explicit RobustHelixFinderDiag(const fhicl::ParameterSet& PSet);
      ~RobustHelixFinderDiag();

    private:

      virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
      virtual int fillHistograms(void* Data, int Mode = -1) override ;
  };


  //-----------------------------------------------------------------------------
  RobustHelixFinderDiag::RobustHelixFinderDiag(const Config& config) {
  }

  //-----------------------------------------------------------------------------
  RobustHelixFinderDiag::RobustHelixFinderDiag(const fhicl::ParameterSet& PSet) {
  }

  //-----------------------------------------------------------------------------
  RobustHelixFinderDiag::~RobustHelixFinderDiag() {
  }


  //-----------------------------------------------------------------------------
  int RobustHelixFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {

    _hist.nTimePeaks    = Tfs->make<TH1F>("ntpeaks"  , "number of time peaks"                      , 51, -0.5, 50.5);
    _hist.nChPerPanel   = Tfs->make<TH1F>("nchppanel", "number of ComboHits per panel"             , 101, -0.5, 100.5);
    _hist.nChHits       = Tfs->make<TH1F>("nchhits"  , "number of ComboHits processed"             , 101, -0.5, 100.5);
    _hist.nseeds[0]     = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events"    , 21, -0.5, 20.5);
    _hist.nseeds[1]     = Tfs->make<TH1F>("nseeds1"  , "number of track candidates: nhits > 15;"    , 21, -0.5, 20.5);
    _hist.ntclhits[0]   = Tfs->make<TH1F>("ntclhits0" , "number of hits on a time peak - no delta"  , 201, -0.5, 200.5);
    _hist.ntclhits[1]   = Tfs->make<TH1F>("ntclhits1" , "number of hits on a time peak - no delta: nhits > 15"  , 201, -0.5, 200.5);
    _hist.nhits[0]       = Tfs->make<TH1F>("nhitsNeg"    , "number of hits on a track candidate"       , 401, -0.5, 800.5);
    _hist.nhits[1]       = Tfs->make<TH1F>("nhitsPos"    , "number of hits on a track candidate"       , 401, -0.5, 800.5);
    _hist.ntripl0       = Tfs->make<TH1F>("ntripl0"    , "number of triplets fit Circle"           , 401, -0.5, 800.5);

    _hist.lambda0[0]    = Tfs->make<TH1F>("lambda0Neg"    , "initFZ, Neg;#lambda [mm/rad] "   , 401, -0.5, 400.5);
    _hist.lambda0[1]    = Tfs->make<TH1F>("lambda0Pos"    , "initFZ, Pos;#lambda [mm/rad] "   , 401, -0.5, 400.5);
    _hist.lambda1[0]    = Tfs->make<TH1F>("lambda1Neg"    , "fitFZ, Neg;#lambda [mm/rad] "    , 401, -0.5, 400.5);
    _hist.lambda1[1]    = Tfs->make<TH1F>("lambda1Pos"    , "fitFZ, Pos;#lambda [mm/rad]"     , 401, -0.5, 400.5);

    _hist.ntripl1[0]    = Tfs->make<TH1F>("ntripl1Neg"    , "number of triplets fit Helix, Neg "   , 401, -0.5, 800.5);
    _hist.ntripl1[1]    = Tfs->make<TH1F>("ntripl1Pos"    , "number of triplets fit Helix, Pos "   , 401, -0.5, 800.5);

    _hist.xyniter[0]    = Tfs->make<TH1F>("xyniterNeg"   , "number of iterations XY fit, Neg "   , 401, -0.5, 800.5);
    _hist.xyniter[1]    = Tfs->make<TH1F>("xyniterPos"   , "number of iterations XY fit, Pos "   , 401, -0.5, 800.5);

    _hist.fzniter[0]    = Tfs->make<TH1F>("fzniterNeg"   , "number of iterations z#phi fit, Neg "   , 401, -0.5, 800.5);
    _hist.fzniter[1]    = Tfs->make<TH1F>("fzniterPos"   , "number of iterations z#phi fit, Pos "   , 401, -0.5, 800.5);

    _hist.niter  [0]    = Tfs->make<TH1F>("niterNeg"     , "number of iterations helix fit, Neg "   , 401, -0.5, 800.5);
    _hist.niter  [1]    = Tfs->make<TH1F>("niterPos"     , "number of iterations helix fit, Pos "   , 401, -0.5, 800.5);

    _hist.nShFitXY  [0]    = Tfs->make<TH1F>("nShFitXYNeg"     , "number of strawhits after the XY fit, Neg "   , 101, -0.5, 100.5);
    _hist.nShFitXY  [1]    = Tfs->make<TH1F>("nShFitXYPos"     , "number of strawhits after the XY fit, Pos "   , 101, -0.5, 100.5);

    _hist.nChFitXY  [0]    = Tfs->make<TH1F>("nChFitXYNeg"     , "number of combo hits after the XY fit, Neg "   , 101, -0.5, 100.5);
    _hist.nChFitXY  [1]    = Tfs->make<TH1F>("nChFitXYPos"     , "number of combo hits after the XY fit, Pos "   , 101, -0.5, 100.5);

    _hist.nShFitCircle  [0]    = Tfs->make<TH1F>("nShFitCircleNeg"     , "number of strawhits after the Circle fit, Neg "   , 101, -0.5, 100.5);
    _hist.nShFitCircle  [1]    = Tfs->make<TH1F>("nShFitCirclePos"     , "number of strawhits after the Circle fit, Pos "   , 101, -0.5, 100.5);

    _hist.nChFitCircle  [0]    = Tfs->make<TH1F>("nChFitCircleNeg"     , "number of combo hits after the Circle fit, Neg "   , 101, -0.5, 100.5);
    _hist.nChFitCircle  [1]    = Tfs->make<TH1F>("nChFitCirclePos"     , "number of combo hits after the Circle fit, Pos "   , 101, -0.5, 100.5);

    _hist.nrescuedhits  [0]    = Tfs->make<TH1F>("nrescuedhitsNeg"     , "number of rescued hits, Neg "   , 801, -0.5, 800.5);
    _hist.nrescuedhits  [1]    = Tfs->make<TH1F>("nrescuedhitsyPos"    , "number of rescued hits, Pos "   , 801, -0.5, 800.5);

    _hist.nXYSh  [0]    = Tfs->make<TH1F>("nXYShNeg"     , "number of strawHits from the circle fit, Neg "   , 201, -0.5, 200.5);
    _hist.nZPhiSh[0]    = Tfs->make<TH1F>("nZPhiShNeg"   , "number of strawHits from the z#phi fit, Neg "   , 201, -0.5, 200.5);
    _hist.nXYSh  [1]    = Tfs->make<TH1F>("nXYShPos"     , "number of strawHits from the circle fit, Pos "   , 201, -0.5, 200.5);
    _hist.nZPhiSh[1]    = Tfs->make<TH1F>("nZPhiShPos"   , "number of strawHits from the z#phi fit, Pos "   , 201, -0.5, 200.5);

    _hist.nfz0counter  [0]    = Tfs->make<TH1F>("nfz0counter_Neg"   , "counts used in initFZ for fz0 eval, Neg "   , 201, -0.5, 200.5);
    _hist.nfz0counter  [1]    = Tfs->make<TH1F>("nfz0counter_Pos"   , "counts used in initFZ for fz0 eval, Pos "   , 201, -0.5, 200.5);

    _hist.nfz0counter_0  [0]    = Tfs->make<TH1F>("nfz0counter0_Neg"   , "counts used in initFZ for fz0 eval, nsh #geq 15, Neg "   , 201, -0.5, 200.5);
    _hist.nfz0counter_0  [1]    = Tfs->make<TH1F>("nfz0counter0_Pos"   , "counts used in initFZ for fz0 eval, nsh #geq 15, Pos "   , 201, -0.5, 200.5);
    _hist.nshsxy_0     [0]   = Tfs->make<TH1F>("nshsxy_0Neg"   , "number of strawHits from the xy fit, Neg "   , 201, -0.5, 200.5);
    _hist.nshsxy_0     [1]   = Tfs->make<TH1F>("nshsxy_0Pos"   , "number of strawHits from the xy fit, Pos "   , 201, -0.5, 200.5);

    _hist.rsxy_0       [0]   = Tfs->make<TH1F>("rsxy_0Neg"     , "radius from the xy fit, Neg "   , 301, 99.5, 400.5);
    _hist.rsxy_0       [1]   = Tfs->make<TH1F>("rsxy_0Pos"     , "radius  from the xy fit, Pos "   , 301, 99.5, 400.5);
    _hist.chi2dsxy_0   [0]   = Tfs->make<TH1F>("chi2dsxy_0Neg", "#chi^{2}/ndof from the xy fit, Neg "   , 2001, -0.05, 200.05);
    _hist.chi2dsxy_0   [1]   = Tfs->make<TH1F>("chi2dsxy_0Pos", "#chi^{2}/ndof from the xy fit, Pos "   , 2001, -0.05, 200.05);

    _hist.nshsxy_1     [0]   = Tfs->make<TH1F>("nshsxy_1Neg"   , "number of strawHits from the xy fit, Neg "   , 201, -0.5, 200.5);
    _hist.nshsxy_1     [1]   = Tfs->make<TH1F>("nshsxy_1Pos"   , "number of strawHits from the xy fit, Pos "   , 201, -0.5, 200.5);

    _hist.rsxy_1       [0]   = Tfs->make<TH1F>("rsxy_1Neg"     , "number of strawHits from the xy fit, Neg "   , 301, 99.5, 400.5);
    _hist.rsxy_1       [1]   = Tfs->make<TH1F>("rsxy_1Pos"     , "number of strawHits from the xy fit, Pos "   , 301, 99.5, 400.5);
    _hist.chi2dsxy_1   [0]   = Tfs->make<TH1F>("chi2dsxy_1Neg", "#chi^{2}/ndof from the xy fit, Neg "   , 201, -0.05, 20.05);
    _hist.chi2dsxy_1   [1]   = Tfs->make<TH1F>("chi2dsxy_1Pos", "#chi^{2}/ndof from the xy fit, Pos "   , 201, -0.5, 20.5);

    _hist.nshszphi_0   [0] = Tfs->make<TH1F>("nshszphi_0Neg"   , "number of strawHits from the zphi fit, Neg "   , 201, -0.5, 200.5);
    _hist.nshszphi_0   [1] = Tfs->make<TH1F>("nshszphi_0Pos"   , "number of strawHits from the zphi fit, Pos "   , 201, -0.5, 200.5);

    _hist.lambdaszphi_0[0] = Tfs->make<TH1F>("lambdaszphi_0Neg", "#lambda from the zphi fit, Neg "   , 301, 99.5, 400.5);
    _hist.lambdaszphi_0[1] = Tfs->make<TH1F>("lambdaszphi_0Pos", "#lambda from the zphi fit, Pos "   , 301, 99.5, 400.5);
    _hist.chi2dszphi_0 [0] = Tfs->make<TH1F>("chi2dszphi_0Neg", "#chi^{2}/ndof from the zphi fit, Neg "   , 201, -0.5, 20.5);
    _hist.chi2dszphi_0 [1] = Tfs->make<TH1F>("chi2dszphi_0Pos", "#chi^{2}/ndof from the zphi fit, Pos "   , 201, -0.5, 20.5);

    _hist.nshszphi_1   [0] = Tfs->make<TH1F>("nshszphi_1Neg"   , "number of strawHits from the zphi fit, Neg "   , 201, -0.5, 200.5);
    _hist.nshszphi_1   [1] = Tfs->make<TH1F>("nshszphi_1Pos"   , "number of strawHits from the zphi fit, Pos "   , 201, -0.5, 200.5);

    _hist.lambdaszphi_1[0] = Tfs->make<TH1F>("lambdaszphi_1Neg", "#lambda from the zphi fit, Neg "   , 301, 99.5, 400.5);
    _hist.lambdaszphi_1[1] = Tfs->make<TH1F>("lambdaszphi_1Pos", "#lambda from the zphi fit, Pos "   , 301, 99.5, 400.5);
    _hist.chi2dszphi_1 [0] = Tfs->make<TH1F>("chi2dszphi_1Neg", "#chi^{2}/ndof from the zphi fit, Neg "   , 201, -0.5, 20.5);
    _hist.chi2dszphi_1 [1] = Tfs->make<TH1F>("chi2dszphi_1Pos", "#chi^{2}/ndof from the zphi fit, Pos "   , 201, -0.5, 20.5);


    _hist.rinit  [0]    = Tfs->make<TH1F>("rinitNeg"       , "helix radius fitCircle, Neg; r [mm]"      , 401, -0.5, 400.5);

    _hist.rinit  [1]    = Tfs->make<TH1F>("rinitPos"       , "helix radius fitCircle, Pos; r [mm]"      , 401, -0.5, 400.5);

    _hist.radius[0]     = Tfs->make<TH1F>("radius0"  , "helix radius; r [mm]"                      , 401, -0.5, 400.5);
    _hist.radius[1]     = Tfs->make<TH1F>("radius1"  , "helix radius nhits > 15; r [mm]"           , 401, -0.5, 400.5);
    _hist.pT [0]        = Tfs->make<TH1F>("pT0"      , "transverse momentum; pT [MeV/c]"           , 400, -0.5, 200.5);
    _hist.p  [0]        = Tfs->make<TH1F>("p0"       , "momentum; p [MeV/c]"                       , 400, -0.5, 200.5);
    _hist.pT [1]        = Tfs->make<TH1F>("pT1"      , "transverse momentum nhits > 15; pT [MeV/c]", 400, -0.5, 200.5);
    _hist.p  [1]        = Tfs->make<TH1F>("p1"       , "momentum nhits > 15; p [MeV/c]"            , 400, -0.5, 200.5);

    _hist.chi2XY[0]     = Tfs->make<TH1F>("chi2XY0_Neg" , "normalized chi2-XY, Neg"               , 500,   0., 50.);
    _hist.chi2XY[1]     = Tfs->make<TH1F>("chi2XY0_Pos" , "normalized chi2-XY, Pos"               , 500,   0., 50.);
    _hist.chi2ZPhi[0]   = Tfs->make<TH1F>("chi2ZPhi_Neg", "normalized chi2-ZPhi, Neg"                 , 500,   0., 50.);
    _hist.chi2ZPhi[1]   = Tfs->make<TH1F>("chi2ZPhi_Pos", "normalized chi2-ZPhi, Pos"                 , 500,   0., 50.);

    _hist.dr  [0]       = Tfs->make<TH1F>("dr0"      , "dr; r - r_{no-target} [mm]"                , 800, -200, 200);
    _hist.dr  [1]       = Tfs->make<TH1F>("dr1"      , "dr: nhits>15; r - r_{no-target} [mm]"      , 800, -200, 200);
    // _hist.shmeanr  [0]  = Tfs->make<TH1F>("shmeanr0" , "straw hit mean radius; r_{sh} [mm]"          , 1800, 0, 900);
    // _hist.shmeanr  [1]  = Tfs->make<TH1F>("shmeanr1" , "straw hit mean radius: nhits>15; r_{sh} [mm]", 1800, 0, 900);
    _hist.chi2d_helix[0]= Tfs->make<TH1F>("chi2dhel0" , "global chi2d; #chi^{2}/ndof"                   , 100, 0, 10);
    _hist.chi2d_helix[1]= Tfs->make<TH1F>("chi2dhel1" , "global chi2d: nhits>15; #chi^{2}/ndof"         , 100, 0, 10);

    _hist.resid[0]     = Tfs->make<TH1F>("resid_Neg" , "helix hit dr, Neg"               , 500,   -250., 250.);
    _hist.resid[1]     = Tfs->make<TH1F>("resid_Pos" , "helix hit dr, Pos"               , 500,   -250., 250.);

    _hist.rwdot[0]     = Tfs->make<TH1F>("rwdot_Neg" , "helix hit rwdot, Neg"               , 200,   -2., 2.);
    _hist.rwdot[1]     = Tfs->make<TH1F>("rwdot_Pos" , "helix hit rwdot, Pos"               , 200,   -2., 2.);

    _hist.nLoops[0]     = Tfs->make<TH1F>("nLoops_Neg" , "helix  nLoops, Neg"               , 7,   -0.5, 6.5);
    _hist.nLoops[1]     = Tfs->make<TH1F>("nLoops_Pos" , "helix  nLoops, Pos"               , 7,   -0.5, 6.5);

    _hist.nHitsLoopFailed[0]     = Tfs->make<TH1F>("nHitsLoopFailed_Neg" , "helix  nHitsLoopFailed, Neg", 21,   -0.5, 20.5);
    _hist.nHitsLoopFailed[1]     = Tfs->make<TH1F>("nHitsLoopFailed_Pos" , "helix  nHitsLoopFailed, Pos", 21,   -0.5, 20.5);

    _hist.eDepAvg[0]    = Tfs->make<TH1F>("eDepAvg_Neg" , "helix  eDepAvg, Neg", 1000,   0, 0.01);
    _hist.eDepAvg[1]    = Tfs->make<TH1F>("eDepAvg_Pos" , "helix  eDepAvg, Pos", 1000,   0, 0.01);

    _hist.chRadialDist[0]     = Tfs->make<TH1F>("meanHitRadialDist_Neg" , "helix hit meanHitRadialDist, Neg"               , 60,   100, 700.);
    _hist.chRadialDist[1]     = Tfs->make<TH1F>("meanHitRadialDist_Pos" , "helix hit meanHitRadialDist, Pos"               , 60,   100, 700.);

    return 0;
  }


  //-----------------------------------------------------------------------------
  // Mode is not used
  //-----------------------------------------------------------------------------
  int RobustHelixFinderDiag::fillHistograms(void* Data, int Mode) {

    _data = (Data_t*) Data;

    //-----------------------------------------------------------------------------
    // fill helix-level histograms
    //-----------------------------------------------------------------------------

    _hist.nTimePeaks ->Fill(_data->nTimePeaks);

    int   nhelicities(2);

    for (int k=0; k<nhelicities; ++k){
      _hist.nseeds[k]->Fill(_data->nseeds[k]);
      _hist.nseeds[k]->Fill(_data->nseeds[k]);

      for (int i=0; i<_data->nseeds[k]; i++) {
        _hist.nChPerPanel->Fill(_data->nChPPanel[k][i] );
        _hist.nChHits    ->Fill(_data->nChHits  [k][i] );

        _hist.ntclhits   [k]->Fill(_data->ntclhits[k][i]   );
        if (k==1){
          _hist.ntripl0 ->Fill(_data->ntriplet0[k][i]      );
        }
        _hist.nhits      [k]->Fill(_data->nhits[k][i]      );
        _hist.rinit      [k]->Fill(_data->rinit[k][i]      );

        _hist.nrescuedhits[k]->Fill(_data->nrescuedhits[k][i]);

        _hist.lambda0    [k]->Fill(fabs(_data->lambda0[k][i]));
        _hist.lambda1    [k]->Fill(fabs(_data->lambda1[k][i]));

        _hist.ntripl1    [k]->Fill(_data->ntriplet1[k][i]  );
        _hist.xyniter    [k]->Fill(_data->xyniter[k][i]    );
        _hist.fzniter    [k]->Fill(_data->fzniter[k][i]    );

        _hist.nShFitXY    [k]->Fill(_data->nShFitXY[k][i]    );
        _hist.nChFitXY    [k]->Fill(_data->nChFitXY[k][i]    );

        _hist.nShFitCircle    [k]->Fill(_data->nShFitCircle[k][i]    );
        _hist.nChFitCircle    [k]->Fill(_data->nChFitCircle[k][i]    );

        _hist.nXYSh      [k]->Fill(_data->nXYSh[k][i]      );
        _hist.nZPhiSh    [k]->Fill(_data->nZPhiSh[k][i]   );

        if (_data->nXYSh[k][i] >= 15){
          _hist.nfz0counter_0  [k]->Fill(_data->nfz0counter[k][i]);
        }

        _hist.p          [k]->Fill(_data->p[k][i]          );
        _hist.pT         [k]->Fill(_data->pT[k][i]         );
        _hist.radius     [k]->Fill(_data->radius[k][i]     );
        _hist.chi2XY     [k]->Fill(_data->chi2XY[k][i]     );
        _hist.chi2ZPhi   [k]->Fill(_data->chi2ZPhi[k][i]   );
        _hist.dr         [k]->Fill(_data->dr[k][i]         );
        _hist.chi2d_helix[k]->Fill(_data->chi2d_helix[k][i]);

        _hist.nfz0counter  [k]->Fill(_data->nfz0counter[k][i]);

        _hist.nshsxy_0     [k] ->Fill(_data->nshsxy_0[k][i]     );

        _hist.rsxy_0       [k] ->Fill(_data->rsxy_0[k][i]       );
        _hist.chi2dsxy_0   [k] ->Fill(_data->chi2dsxy_0[k][i]   );

        _hist.nshsxy_1     [k] ->Fill(_data->nshsxy_1[k][i]     );

        _hist.rsxy_1       [k] ->Fill(_data->rsxy_1[k][i]       );
        _hist.chi2dsxy_1   [k] ->Fill(_data->chi2dsxy_1[k][i]   );

        _hist.nshszphi_0   [k] ->Fill(_data->nshszphi_0[k][i]   );

        _hist.lambdaszphi_0[k] ->Fill(_data->lambdaszphi_0[k][i]);
        _hist.chi2dszphi_0 [k] ->Fill(_data->chi2dszphi_0[k][i] );


        _hist.nshszphi_1   [k] ->Fill(_data->nshszphi_1[k][i]   );

        _hist.lambdaszphi_1[k] ->Fill(_data->lambdaszphi_1[k][i]);
        _hist.chi2dszphi_1 [k] ->Fill(_data->chi2dszphi_1[k][i] );

        for (int j=0; j<_data->nXYCh[k][i];++j){
          _hist.rwdot[k]->Fill(_data->hitRWDot[k][i][j]);
          _hist.resid[k]->Fill(_data->hitDr   [k][i][j]);
        }

        if (_data->nZPhiSh[k][i] >= 15){
          _hist.nLoops       [k] ->Fill(_data->nLoops[k][i]);
          _hist.nHitsLoopFailed       [k] ->Fill(_data->nHitsLoopFailed[k][i]);
          _hist.chRadialDist [k] ->Fill(_data->meanHitRadialDist[k][i]);
        }

        _hist.eDepAvg      [k] ->Fill(_data->eDepAvg[k][i]);

      }
    }

    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::RobustHelixFinderDiag)
