#include "TH1.h"
#include "TH2.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "Offline/CalPatRec/inc/CalHelixFinder_types.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {

  using namespace CalHelixFinderTypes;

  class CalHelixFinderDiag : public mu2e::ModuleHistToolBase {

  protected:
    Hist_t                     _hist;
    Data_t*                    _data;
    int                        _event_number;

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
  _event_number = -1;
}

CalHelixFinderDiag::~CalHelixFinderDiag() {
}

//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
int CalHelixFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {

  _hist.nTimePeaks    = Tfs->make<TH1F>("ntpeaks"  , "number of time peaks"                      , 11, -0.5, 10.5);
  _hist.nseeds[0]     = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events"    , 21, -0.5, 20.5);
  _hist.nseeds[1]     = Tfs->make<TH1F>("nseeds1"  , "number of track candidates: nhits > 15;"    , 21, -0.5, 20.5);
  _hist.drVsDzSeed[0] = Tfs->make<TH2F>("drVsDzSeed0"  , "Dist from prediction vs dz from seed-hit: all events; #Delta z from seed-hit [mm]; dr from prediction [mm]", 100, -1e3, 3e3, 60, 0, 300);
  _hist.drVsDzSeed[1] = Tfs->make<TH2F>("drVsDzSeed1"  , "Dist from prediction vs dz from seed-hit: nhits > 15; #Delta z from seed-hit [mm]; dr from prediction [mm]", 100, -1e3, 3e3, 60, 0, 300);
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
  _hist.nhitsvspT     = Tfs->make<TH2F>("nhitsvspT", "nhits vs pT"                               , 100, 0, 100, 400, 0, 200);
  _hist.nhitsvsp      = Tfs->make<TH2F>("nhitsvsp" , "nhits vs p"                                , 100, 0, 100, 400, 0, 200);
  _hist.nStationPairs = Tfs->make<TH1F>("nstpairs" , "n-station pairs used in findDfDz"          , 200, 0, 200);
  _hist.dr  [0]       = Tfs->make<TH1F>("dr0"      , "dr; r - r_{no-target} [mm]"                , 800, -200, 200);
  _hist.dr  [1]       = Tfs->make<TH1F>("dr1"      , "dr: nhits>15; r - r_{no-target} [mm]"      , 800, -200, 200);
  _hist.shmeanr  [0]  = Tfs->make<TH1F>("shmeanr0" , "straw hit mean radius; r_{sh} [mm]"          , 1800, 0, 900);
  _hist.shmeanr  [1]  = Tfs->make<TH1F>("shmeanr1" , "straw hit mean radius: nhits>15; r_{sh} [mm]", 1800, 0, 900);
  _hist.chi2d_helix[0]= Tfs->make<TH1F>("chi2dhel0" , "global chi2d; #chi^{2}/ndof"                   , 100, 0, 10);
  _hist.chi2d_helix[1]= Tfs->make<TH1F>("chi2dhel1" , "global chi2d: nhits>15; #chi^{2}/ndof"         , 100, 0, 10);
  _hist.chi2d_loop0[0] = Tfs->make<TH1F>("chi2dloop0", "XY chi2d: loop 0; #chi^{2}/ndof"              , 100, 0, 10);
  _hist.chi2d_loop0[1] = Tfs->make<TH1F>("chi2dloop0Good", "XY chi2d: loop 0: nhits>15; #chi^{2}/ndof", 100, 0, 10);
  _hist.chi2d_loop1[0] = Tfs->make<TH1F>("chi2dloop1", "XY chi2d: loop 1; #chi^{2}/ndof"              , 100, 0, 10);
  _hist.chi2d_loop1[1] = Tfs->make<TH1F>("chi2dloop1Good", "XY chi2d: loop 1: nhits>15; #chi^{2}/ndof", 100, 0, 10);
  _hist.chi2d_line_loop0[0] = Tfs->make<TH1F>("chi2dZPhiloop0", "Z#phi chi2d: loop 0; {#chi^{2}/ndof}_{z#phi}"              , 1000, 0, 10);
  _hist.chi2d_line_loop0[1] = Tfs->make<TH1F>("chi2dZPhiloop0Good", "Z#phi chi2d: loop 0: nhits>15; {#chi^{2}/ndof}_{z#phi}", 1000, 0, 10);
  _hist.chi2d_line_loop1[0] = Tfs->make<TH1F>("chi2dZPhiloop1", "Z#phi chi2d: loop 1; {#chi^{2}/ndof}_{z#phi}"              , 1000, 0, 10);
  _hist.chi2d_line_loop1[1] = Tfs->make<TH1F>("chi2dZPhiloop1Good", "Z#phi chi2d: loop 1: nhits>15; {#chi^{2}/ndof}_{z#phi}", 1000, 0, 10);
  _hist.npoints_loop0 = Tfs->make<TH1F>("npointsloop0", "XY npoints: loop 0; nhits"           , 101, -0.5, 100.5);
  _hist.npoints_loop1 = Tfs->make<TH1F>("npointsloop1", "XY npoints: loop 1; nhits"           , 101, -0.5, 100.5);
  _hist.loopId[0]     = Tfs->make<TH1F>("loopAll"   , "loopId; loopId"                           , 10, 0, 10);
  _hist.loopId[1]     = Tfs->make<TH1F>("loopGood"  , "loopId: nhits>15: loopId"                 , 10, 0, 10);
  _hist.nHitsRatio[0] = Tfs->make<TH1F>("nHitsRatio0","circleHits/phiHits"                       ,1000,0, 5);
  _hist.nHitsRatio[1] = Tfs->make<TH1F>("nHitsRatio1","circleHits/phiHits:nhits>15"              ,1000,0, 5);
  _hist.eDepAvg[0]    = Tfs->make<TH1F>("eDepAvg0","avg comboHit eDep"                          ,1000,0,0.01);
  _hist.eDepAvg[1]    = Tfs->make<TH1F>("eDepAvg1","avg comboHit eDep:nhits>15"                 ,1000,0,0.01);
  return 0;
}


//-----------------------------------------------------------------------------
// called once per event, Mode is not used..
//-----------------------------------------------------------------------------
  int CalHelixFinderDiag::fillHistograms(void * Data, int Mode) {

    _data = (Data_t*) Data;

    int en = _data->event->event();
    if ( en == _event_number) return -1;

    _event_number = en;
//-----------------------------------------------------------------------------
// this part is so far simple
//-----------------------------------------------------------------------------
    _hist.nTimePeaks->Fill(_data->nTimePeaks);
    _hist.nseeds[0]->Fill(_data->nseeds[0]);
    _hist.nseeds[1]->Fill(_data->nseeds[1]);

    for (int i=0; i<_data->nseeds[0]; i++) {
      _hist.ntclhits[0]->Fill(_data->ntclhits[i]);
      _hist.nhits   ->Fill(_data->nhits   [i]);
      for (int j=0; j<_data->nhits   [i]; ++j){
        _hist.drVsDzSeed[0]->Fill(_data->hitDzSeed[i][j], _data->hitDrPred[i][j]);
      }
      _hist.p[0]->Fill(_data->p[i]);
      _hist.pT[0]->Fill(_data->pT[i]);
      _hist.radius[0]->Fill(_data->radius[i]);
      _hist.chi2XY[0]->Fill(_data->chi2XY[i]);
      _hist.chi2ZPhi[0]->Fill(_data->chi2ZPhi[i]);
      _hist.dr[0]->Fill(_data->dr[i]);
      _hist.shmeanr[0]->Fill(_data->shmeanr[i]);
      _hist.chi2d_helix[0]->Fill(_data->chi2d_helix[i]);
      _hist.loopId     [0]->Fill(_data->loopId[i]);
      _hist.nHitsRatio[0]->Fill(_data->nHitsRatio[i]);
      _hist.eDepAvg[0]->Fill(_data->eDepAvg[i]);
      if (_data->loopId[i] == 1){
        _hist.chi2d_loop0 [0]->Fill(_data->chi2d_loop0[i]);
        _hist.chi2d_line_loop0 [0]->Fill(_data->chi2d_line_loop0[i]);
        _hist.npoints_loop0  ->Fill(_data->npoints_loop0[i]);
      }
      if (_data->loopId[i] == 2) {
        _hist.chi2d_loop1 [0]->Fill(_data->chi2d_loop1[i]);
        _hist.chi2d_line_loop1 [0]->Fill(_data->chi2d_line_loop1[i]);
        _hist.npoints_loop1  ->Fill(_data->npoints_loop1[i]);
      }

      if (_data->good[i] != 0) {
        _hist.ntclhits[1]->Fill(_data->ntclhits[i]);
        _hist.loopId  [1]->Fill(_data->loopId[i]);
        _hist.p       [1]->Fill(_data->p[i]);
        _hist.pT      [1]->Fill(_data->pT[i]);
        _hist.radius  [1]->Fill(_data->radius[i]);
        _hist.chi2XY  [1]->Fill(_data->chi2XY[i]);
        _hist.chi2ZPhi[1]->Fill(_data->chi2ZPhi[i]);
        _hist.dr      [1]->Fill(_data->dr[i]);
        _hist.shmeanr [1]->Fill(_data->shmeanr[i]);
        _hist.chi2d_helix[1]->Fill(_data->chi2d_helix[i]);
        _hist.nHitsRatio[1]->Fill(_data->nHitsRatio[i]);
        _hist.eDepAvg[1]->Fill(_data->eDepAvg[i]);
        if (_data->loopId[i] == 1) {
          _hist.chi2d_loop0   [1]->Fill(_data->chi2d_loop0[i]);
        }
        if (_data->loopId[i] == 2) {
          _hist.chi2d_loop1   [1]->Fill(_data->chi2d_loop1[i]);
        }

        for (int j=0; j<_data->nhits   [i]; ++j){
        _hist.drVsDzSeed[1]->Fill(_data->hitDzSeed[i][j], _data->hitDrPred[i][j]);
      }
      }

      _hist.nhitsvspT ->Fill(_data->nhits[i],_data->pT[i]);
      _hist.nhitsvsp  ->Fill(_data->nhits[i],_data->p [i]);

      _hist.nStationPairs->Fill(_data->nStationPairs[i]);
    }

    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::CalHelixFinderDiag)
