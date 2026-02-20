#include <iostream>
#include <string>
#include <cmath>

// Cosmic Tracks:
#include "Offline/CosmicReco/inc/CosmicTrackFit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/CosmicReco/inc/CosmicTrackMCInfo.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/CosmicReco/inc/PDFFit.hh"
#include "Offline/CosmicReco/inc/MinuitDriftFitter.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

//Mu2e Data Prods:
//
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

//Utilities
#include "Offline/Mu2eUtilities/inc/ParametricFit.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"

// Mu2e diagnostics
#include "Offline/GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
//ROOT
#include "TStyle.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TTree.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TProfile.h"
using namespace std;

namespace mu2e
{
  class CosmicTrackDiag : public art::EDAnalyzer {
    public:
      struct Config{
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diag{Name("diagLevel"), Comment("set to 1 for info"),2};
        fhicl::Atom<int> printlevel{Name("printLevel"), Comment("print level"), 0};
        fhicl::Atom<bool> mcdiag{Name("mcdiag"), Comment("set on for MC info"),false};
        fhicl::Atom<bool> shdiag{Name("shdiag"), Comment("shdiag"),false};
        fhicl::Atom<bool> ubdiag{Name("ubresids"), Comment("Calculate unbiased residuals"),false};
        fhicl::Atom<art::InputTag> phtag{Name("PanelHitCollection"),Comment("tag for panel hits")};
        fhicl::Atom<art::InputTag> chtag{Name("ComboHitCollection"),Comment("tag for combo hit collection (single straw hit level)")};
        fhicl::Atom<art::InputTag> shtag{Name("StrawHitCollection"),Comment("tag for straw hit objects")};
        fhicl::Atom<art::InputTag> tctag{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
        fhicl::Atom<art::InputTag> tstag{Name("CosmicTrackSeedCollection"),Comment("CosmicTrackSeed collection tag")};
        fhicl::Atom<art::InputTag> mcdigistag{Name("StrawDigiMCCollection"),Comment("StrawDigi collection tag")};
        fhicl::Atom<art::InputTag> pbtmcTag{Name("ProtonBunchTimeMC"),Comment("ProtonBunchTimeMC tag")};
      };
      typedef art::EDAnalyzer::Table<Config> Parameters;

      explicit CosmicTrackDiag(const Parameters& conf);
      virtual ~CosmicTrackDiag();
      virtual void beginJob() override;
      virtual void analyze(const art::Event& e) override;
    private:

      Config _conf;

      int  _diag;
      int _printlevel;
      bool _mcdiag;
      bool _shdiag;
      bool _ubdiag;
      art::InputTag   _phtag;
      art::InputTag   _chtag;//combo
      art::InputTag   _shtag;
      art::InputTag   _tctag;//timeclusters
      art::InputTag   _tstag;//Striaght tracks
      art::InputTag   _mcdigistag; //MC digis
      art::InputTag _pbtmcTag;
      const ComboHitCollection* _phcol;
      const ComboHitCollection* _chcol;
      const StrawHitCollection* _shcol;
      const TimeClusterCollection* _tccol;
      const CosmicTrackSeedCollection* _tscol;
      const StrawDigiMCCollection* _mcdigis;
      CLHEP::Hep3Vector _mcpos, _mcdir;

      Float_t _pbtmc;
      ProditionsHandle<Tracker> _alignedTracker_h;
      ProditionsHandle<StrawResponse> _strawResponse_h;


      //TTree Info:
      TTree* _trackT;
      TTree* _hitT;

      // track tree
      Int_t _evt;
      Int_t _run;
      Int_t _subrun;
      Int_t _ntrack;
      Int_t _nsh; // # associated straw hits / event
      Int_t _ntc; // # clusters/event
      Int_t _n_panels; // # panels
      Int_t _n_stations; // # stations
      Int_t _n_planes; // # stations
      Int_t _mcnsh;
      Int_t _outsidehits;
      Float_t _minuitdoca, _minuitangle;
      Int_t _tcnhits, _ontrackhits;
      Int_t _converged, _minuitconverged;
      Float_t _maxllike, _maxtresid, _llike;
      Float_t _a0,_b0,_a1,_b1,_t0;
      Float_t _a0err,_b0err,_a1err,_b1err,_t0err;
      Float_t _mca0,_mcb0,_mca1,_mcb1, _mct0;

      // hit tree
      Float_t _hitminuitdoca, _hitminuitlong, _hitminuitdpocat;
      Float_t _hitublong, _hitubdoca, _hitubtresid, _hitublongrms, _hitubtresidrms;
      Float_t _hitubt0;
      Float_t _hitfiterr,_hitfiterr2;
      Float_t _hitdz;
      Float_t _hittruedoca, _hitmcdoca, _hitrecodoca, _hitdeltat;
      Float_t _hittot[2];

      Float_t _hitmcdpocat;
      Float_t _hittruelong, _hitmclong, _hitrecolong;
      Float_t _hitmctresid, _hitmctresidrms, _hitmcdrifttime, _hitmcdrifttimeoffset, _hitmctrajtime;
      Float_t _hitmcmom;
      Int_t _hitmcpdg, _hitmccc;
      Int_t _hitbackground;
      Int_t _hitused, _hitintclust;
      Int_t _hitambig, _hitmcambig;
      Float_t _hitrecolongrms;
      Float_t _hittresid, _hittresidrms, _hitdrifttime, _hitdrifttimeoffset, _hittrajtime, _hitllike;
      Float_t _hitdresid, _hitdresiderr;
      Float_t _hitedep, _hitproptime, _hittime;

      Int_t _hitplane, _hitpanel, _hitstraw;
      Float_t _strawlen;

      Float_t _ewMarkerOffset;

      void GetMCTrack(const art::Event& event, const StrawDigiMCCollection& mccol, const Tracker* tracker);
      void hitDiag(const art::Event& event, StrawResponse const& srep, int its, const Tracker* tracker);
      bool findData(const art::Event& evt);
  };

  CosmicTrackDiag::CosmicTrackDiag(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _diag (conf().diag()),
    _printlevel (conf().printlevel()),
    _mcdiag (conf().mcdiag()),
    _shdiag (conf().shdiag()),
    _ubdiag (conf().ubdiag()),
    _phtag (conf().phtag()),
    _chtag (conf().chtag()),
    _shtag (conf().shtag()),
    _tctag (conf().tctag()),
    _tstag (conf().tstag()),
    _mcdigistag (conf().mcdigistag()),
    _pbtmcTag (conf().pbtmcTag())
  {
    consumes<CosmicTrackSeedCollection>(_tstag);
  }

  CosmicTrackDiag::~CosmicTrackDiag(){}


  void CosmicTrackDiag::beginJob() {

    // create diagnostics if requested...
    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      //Tree for detailed diagnostics
      _trackT=tfs->make<TTree>("trackT","Track level diagnostics");

      //Create branches:
      _trackT->Branch("evt",&_evt,"evt/I");  // add event id
      _trackT->Branch("run",&_run,"run/I");
      _trackT->Branch("subrun",&_subrun,"subrun/I");
      _trackT->Branch("ntrack",&_ntrack,"ntrack/I");
      _trackT->Branch("outsidehits",&_outsidehits,"outsidehits/I");
      _trackT->Branch("StrawHitsInEvent", &_nsh, "StrawHitsInEvent/I");
      _trackT->Branch("PanelsCrossedInEvent", &_n_panels, "PanelsCrossedInEvent/I");
      _trackT->Branch("PlanesCrossedInEvent", &_n_planes, "PlanesCrossedInEvent/I");
      _trackT->Branch("StationsCrossedInEvent", &_n_stations, "StationsCrossedInEvent/I");
      _trackT->Branch("TimeClustersInEvent", &_ntc, "TimeClusterInEvent/I");
      _trackT->Branch("a0",&_a0,"a0/F");
      _trackT->Branch("b0",&_b0,"b0/F");
      _trackT->Branch("a1",&_a1,"a1/F");
      _trackT->Branch("b1",&_b1,"b1/F");
      _trackT->Branch("t0",&_t0,"t0/F");
      _trackT->Branch("a0err",&_a0err,"a0err/F");
      _trackT->Branch("b0err",&_b0err,"b0err/F");
      _trackT->Branch("a1err",&_a1err,"a1err/F");
      _trackT->Branch("b1err",&_b1err,"b1err/F");
      _trackT->Branch("t0err",&_t0err,"t0err/F");
      _trackT->Branch("tcnhits",&_tcnhits,"tcnhits/I");
      _trackT->Branch("ontrackhits",&_ontrackhits,"ontrackhits/I");
      _trackT->Branch("ewmoffset",&_ewMarkerOffset,"ewmoffset/F");
      _trackT->Branch("converged",&_converged,"converged/I");
      _trackT->Branch("minuitconverged",&_minuitconverged,"minuitconverged/I");
      _trackT->Branch("llike",&_llike,"llike/F");
      _trackT->Branch("maxtresid",&_maxtresid,"maxtresid/F");
      _trackT->Branch("maxllike",&_maxllike,"maxllike/F");
      if (_mcdiag){
        _trackT->Branch("mcnsh",&_mcnsh,"mcnsh/I");
        _trackT->Branch("mca0",&_mca0,"mca0/F");
        _trackT->Branch("mcb0",&_mcb0,"mcb0/F");
        _trackT->Branch("mca1",&_mca1,"mca1/F");
        _trackT->Branch("mcb1",&_mcb1,"mcb1/F");
        _trackT->Branch("mct0",&_mct0,"mct0/F");
        _trackT->Branch("minuitdoca",&_minuitdoca,"minuitdoca/F");
        _trackT->Branch("minuitangle",&_minuitangle,"minuitangle/F");
      }

      _hitT=tfs->make<TTree>("hitT","Hit tree");
      _hitT->Branch("evt",&_evt,"evt/I");  // add event id
      _hitT->Branch("run",&_run,"run/I");
      _hitT->Branch("subrun",&_subrun,"subrun/I");
      _hitT->Branch("ntrack",&_ntrack,"ntrack/I");
      _hitT->Branch("outsidehits",&_outsidehits,"outsidehits/I");
      _hitT->Branch("doca",&_hitminuitdoca,"doca/F");
      _hitT->Branch("hitused",&_hitused,"hitused/I");
      _hitT->Branch("hitintclust",&_hitintclust,"hitintclust/I");
      _hitT->Branch("long",&_hitminuitlong,"long/F");
      _hitT->Branch("time",&_hittime,"time/F");
      _hitT->Branch("proptime",&_hitproptime,"proptime/F");
      _hitT->Branch("trajtime",&_hittrajtime,"trajtime/F");
      _hitT->Branch("drifttime",&_hitdrifttime,"drifttime/F");
      _hitT->Branch("drifttimeoffset",&_hitdrifttimeoffset,"drifttimeoffset/F");
      _hitT->Branch("recolong",&_hitrecolong,"recolong/F");
      _hitT->Branch("recodoca",&_hitrecodoca,"recodoca/F");
      _hitT->Branch("recolongrms",&_hitrecolongrms,"recolongrms/F");
      _hitT->Branch("tresid",&_hittresid,"tresid/F");
      _hitT->Branch("tresidrms",&_hittresidrms,"tresidrms/F");
      _hitT->Branch("dresid",&_hitdresid,"dresid/F");
      _hitT->Branch("dresiderr",&_hitdresiderr,"dresiderr/F");
      _hitT->Branch("hitllike",&_hitllike,"hitllike/F");
      _hitT->Branch("tcnhits",&_tcnhits,"tcnhits/I");
      _hitT->Branch("ontrackhits",&_ontrackhits,"ontrackhits/I");
      _hitT->Branch("t0",&_t0,"t0/F");
      _hitT->Branch("ewmoffset",&_ewMarkerOffset,"ewmoffset/F");
      _hitT->Branch("converged",&_converged,"converged/I");
      _hitT->Branch("minuitconverged",&_minuitconverged,"minuitconverged/I");
      _hitT->Branch("llike",&_llike,"llike/F");
      _hitT->Branch("maxtresid",&_maxtresid,"maxtresid/F");
      _hitT->Branch("maxllike",&_maxllike,"maxllike/F");
      _hitT->Branch("PanelsCrossedInEvent", &_n_panels, "PanelsCrossedInEvent/I");
      _hitT->Branch("plane",&_hitplane,"plane/I");
      _hitT->Branch("panel",&_hitpanel,"panel/I");
      _hitT->Branch("straw",&_hitstraw,"straw/I");
      _hitT->Branch("strawlen",&_strawlen,"strawlen/F");
      _hitT->Branch("edep",&_hitedep,"edep/F");
      _hitT->Branch("ambig",&_hitambig,"ambig/I");
      _hitT->Branch("dz",&_hitdz,"dz/F");
      _hitT->Branch("edep",&_hitedep,"edep/F");

      if (_ubdiag){
        _hitT->Branch("ublong",&_hitublong,"ublong/F");
        _hitT->Branch("ublongrms",&_hitublongrms,"ublongrms/F");
        _hitT->Branch("ubdoca",&_hitubdoca,"ubdoca/F");
        _hitT->Branch("ubt0",&_hitubt0,"ubt0/F");
        _hitT->Branch("ubtresid",&_hitubtresid,"ubtresid/F");
        _hitT->Branch("ubtresidrms",&_hitubtresidrms,"ubtresidrms/F");
      }
      if (_shdiag){
        _hitT->Branch("deltat",&_hitdeltat,"deltat/F");
        _hitT->Branch("tot",&_hittot,"totcal/F:tothv/F");
      }


      if (_mcdiag){
        _hitT->Branch("background",&_hitbackground,"background/I");
        _hitT->Branch("mct0",&_mct0,"mct0/F");
        _hitT->Branch("truedoca",&_hittruedoca,"truedoca/F");
        _hitT->Branch("mcdoca",&_hitmcdoca,"mcdoca/F");
        _hitT->Branch("mctresid",&_hitmctresid,"mctresid/F");
        _hitT->Branch("mctresidrms",&_hitmctresidrms,"mctresidrms/F");
        _hitT->Branch("mcdrifttime",&_hitmcdrifttime,"mcdrifttime/F");
        _hitT->Branch("mcdrifttimeoffset",&_hitmcdrifttimeoffset,"mcdrifttimeoffset/F");
        _hitT->Branch("mctrajtime",&_hitmctrajtime,"mctrajtime/F");
        _hitT->Branch("mcdpocat",&_hitmcdpocat,"mcdpocat/F");
        _hitT->Branch("minuitdpocat",&_hitminuitdpocat,"minuitdpocat/F");
        _hitT->Branch("truelong",&_hittruelong,"truelong/F");
        _hitT->Branch("mclong",&_hitmclong,"mclong/F");
        _hitT->Branch("mcnsh",&_mcnsh,"mcnsh/I");
        _hitT->Branch("mct0",&_mct0,"mct0/F");
        _hitT->Branch("minuitdoca",&_minuitdoca,"minuitdoca/F");
        _hitT->Branch("minuitangle",&_minuitangle,"minuitangle/F");
        _hitT->Branch("mcambig",&_hitmcambig,"mcambig/I");
        _hitT->Branch("mcpdg",&_hitmcpdg,"mcpdg/I");
        _hitT->Branch("mccc",&_hitmccc,"mccc/I");
        _hitT->Branch("mcmom",&_hitmcmom,"mcmom/F");
      }
    }
  }

  void CosmicTrackDiag::analyze(const art::Event& event) {

    _evt = event.id().event();  // add event id
    _run = event.id().run();
    _subrun = event.id().subRun();
    if(!findData(event)) // find data
      throw cet::exception("RECO")<<"No Time Clusters in event"<< endl;

    Tracker const& tracker = _alignedTracker_h.get(event.id());

    StrawResponse const& srep = _strawResponse_h.get(event.id());

    //find time clusters:
    _ntc = _tccol->size();


    if (_mcdiag){
      GetMCTrack(event, *_mcdigis,&tracker);
    }

    std::vector<uint16_t> panels, planes, stations;
    _nsh = 0;
    for(size_t ich = 0;ich < _chcol->size(); ++ich){
      ComboHit const& chit =(*_chcol)[ich];
      _nsh += chit.nStrawHits();
      uint16_t panelid = chit.strawId().uniquePanel();
      if (std::find(panels.begin(),panels.end(),panelid) == panels.end())
        panels.push_back(panelid);
      uint16_t planeid = chit.strawId().plane();
      if (std::find(planes.begin(),planes.end(),planeid) == planes.end())
        planes.push_back(planeid);
      uint16_t stationid = chit.strawId().station();
      if (std::find(stations.begin(),stations.end(),stationid) == stations.end())
        stations.push_back(stationid);
    }
    _n_panels = panels.size();
    _n_planes = planes.size();
    _n_stations = stations.size();
    _ntrack = -1;

    _minuitconverged = 0;
    _converged = 0;
    for (size_t its=0;its<_tscol->size();its++){
      auto tseed = _tscol->at(its);

      auto tclust = tseed._timeCluster;
      _tcnhits = 0;
      const std::vector<StrawHitIndex>& shIndices = tclust->hits();
      for (size_t i=0; i<shIndices.size(); ++i) {
        int loc = shIndices[i];
        const ComboHit& ch  = _phcol->at(loc);
        _tcnhits += ch.nStrawHits();
      }

      _a0 = tseed._track.MinuitParams.A0;
      _b0 = tseed._track.MinuitParams.B0;
      _a1 = tseed._track.MinuitParams.A1;
      _b1 = tseed._track.MinuitParams.B1;
      _a0err = tseed._track.MinuitParams.deltaA0;
      _b0err = tseed._track.MinuitParams.deltaB0;
      _a1err = tseed._track.MinuitParams.deltaA1;
      _b1err = tseed._track.MinuitParams.deltaB1;
      _t0 = tseed._t0._t0;
      _t0err = tseed._t0._t0err;
      _converged = tseed._track.converged;
      _minuitconverged = tseed._track.minuit_converged;
      _ntrack = 0;
      auto minuitpos = GenVector::Hep3Vec(tseed._track.MinuitEquation.Pos);
      auto minuitdir = GenVector::Hep3Vec(tseed._track.MinuitEquation.Dir.unit());
      // convert to y orientation
      if (minuitdir.y() > 0)
        minuitdir *= -1;
      if (minuitdir.y() != 0)
        minuitpos -= minuitdir*minuitpos.y()/minuitdir.y();

//      const std::vector<StrawHitIndex>& lseedShIndices = lseed._strawHitIdxs;
      _maxtresid = 0;
      _maxllike = 0;
      _llike = 0;
      _ontrackhits = 0;
      _outsidehits = 0;

      for (size_t i=0; i<tseed._straw_chits.size(); ++i) {
        auto sh = tseed._straw_chits[i];

        double llike = 0;
        Straw const& straw = tracker.getStraw(sh.strawId());
        TwoLinePCA pca(straw.getMidPoint(), straw.getDirection(), minuitpos, minuitdir);

        double longdist = (pca.point1() - straw.getMidPoint()).dot(straw.getDirection());
        if (fabs(longdist) > straw.halfLength()){
          llike += pow((fabs(longdist)-straw.halfLength()),2);
          longdist = std::copysign(straw.halfLength(),longdist);
        }
        //double longres = srep.wpRes(sh.energyDep() * 1000., fabs(longdist));
        double longres = srep.wpRes(sh.energyDep()*1000., sh.wireDist());

        //FIXME ignoring likelihood contribution from outside straw
        //llike += pow(longdist - sh.wireDist(), 2) / pow(longres, 2);
        double trunc_wireDist = std::copysign(std::min(straw.halfLength(),fabs(sh.wireDist())),sh.wireDist());
        llike += pow(longdist - trunc_wireDist, 2) / pow(longres, 2);


        double even_z = tracker.getStraw(StrawId(sh.strawId().plane(),sh.strawId().panel(),0)).getMidPoint().z();
        double odd_z = tracker.getStraw(StrawId(sh.strawId().plane(),sh.strawId().panel(),1)).getMidPoint().z();
        if ((pca.point2().z() > even_z && pca.point2().z() > odd_z) || (pca.point2().z() < even_z && pca.point2().z() < odd_z)){
          _outsidehits++;
        }

        double drift_time = srep.driftDistanceToTime(sh.strawId(), pca.dca(), 0);
        drift_time += srep.driftTimeOffset(sh.strawId(), pca.dca(), 0);

        double drift_res = srep.driftTimeError(sh.strawId(), pca.dca(), 0);

        double traj_time = ((pca.point2() - minuitpos).dot(minuitdir))/299.9;
        double hit_t0 = sh.time() - sh.propTime() - traj_time - drift_time;
        double tresid = hit_t0-_t0;
        llike += pow(tresid,2)/pow(drift_res,2);

        if (pca.dca() < 2.55 && fabs(tresid) < 100)
          _ontrackhits += 1;

        if (fabs(tresid) > _maxtresid)
          _maxtresid = fabs(tresid);
        if (llike > _maxllike)
          _maxllike = llike;
        _llike += llike;
      }
      _llike /= tseed._straw_chits.size();


      TwoLinePCA minuitpca( _mcpos, _mcdir, minuitpos, minuitdir);
      _minuitdoca = minuitpca.dca();
      _minuitangle = _mcdir.dot(minuitdir);

      if (_printlevel > 0){
        std::cout << "Track found! angle: " << _minuitangle << " hits on track: " << _ontrackhits << " / " << _nsh << std::endl;
      }
      _trackT->Fill();

      if (_diag > 1)
        hitDiag(event, srep, its, &tracker);
    }

    if (_tscol->size() == 0){
      if (_diag > 1)
        hitDiag(event, srep, -1, &tracker);
      _trackT->Fill();
    }
  }

  bool CosmicTrackDiag::findData(const art::Event& evt){
    _phcol = 0;
    _chcol = 0;
    _shcol = 0;
    _tccol = 0;
    _tscol = 0;
    _ewMarkerOffset = 0;
    _pbtmc = 0;

    auto phH = evt.getValidHandle<ComboHitCollection>(_phtag);
    _phcol = phH.product();
    auto chH = evt.getValidHandle<ComboHitCollection>(_chtag);
    _chcol = chH.product();
    auto tcH = evt.getValidHandle<TimeClusterCollection>(_tctag);
    _tccol =tcH.product();
    auto tsH = evt.getValidHandle<CosmicTrackSeedCollection>(_tstag);
    _tscol =tsH.product();
    auto pbtmcHandle = evt.getValidHandle<ProtonBunchTimeMC>(_pbtmcTag);
    _pbtmc = pbtmcHandle.product()->pbtime_;
    _ewMarkerOffset = -1*_pbtmc;
    if(_mcdiag){
      _mcdigis=0;
      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigistag);
      _mcdigis = mcdH.product();
    }

    if (_shdiag){
      auto shH = evt.getValidHandle<StrawHitCollection>(_shtag);
      _shcol = shH.product();
    }
    return _phcol != 0 && _chcol != 0 && _tccol!=0 && _tscol !=0 && (_mcdigis != 0 || !_mcdiag) && (_shcol != 0 || !_shdiag);
  }


  void CosmicTrackDiag::GetMCTrack(const art::Event& event, const StrawDigiMCCollection& mccol, const Tracker* tracker) {
    // get all possible directions
    std::vector<CLHEP::Hep3Vector> pppos;
    std::vector<CLHEP::Hep3Vector> ppdir;
    for (size_t i=0;i<mccol.size();i++){
      StrawDigiMC mcdigi = mccol[i];
      auto const& sgsptr = mcdigi.earlyStrawGasStep();
      auto const& sgs = *sgsptr;
      auto const& sp = *sgs.simParticle();
      auto posi = GenVector::Hep3Vec(sgs.startPosition());
      if ((sp.pdgId() == PDGCode::mu_minus || sp.pdgId() == PDGCode::mu_plus) && sp.creationCode() == 56){
        for (size_t j=i+1;j<mccol.size();j++){
          StrawDigiMC jmcdigi = mccol[j];
          auto const& jsgsptr = jmcdigi.earlyStrawGasStep();
          auto const& jsgs = *jsgsptr;
          auto const& jsp = *jsgs.simParticle();
          auto posj = GenVector::Hep3Vec(jsgs.startPosition());
          if ((jsp.pdgId() == PDGCode::mu_minus || jsp.pdgId() == PDGCode::mu_plus) && jsp.creationCode() == 56){
            pppos.push_back(posi);
            ppdir.push_back((posi-posj).unit());
          }
        }
      }
    }

    // get the one that has the most hits within 500 microns
    int max = 0;
    for (size_t j=0;j<pppos.size();j++){
      int count = 0;
      double avg_t0 = 0;
      CLHEP::Hep3Vector ppintercept(0,0,0);
      CLHEP::Hep3Vector ppdirection(0,1,0);
      for (size_t i=0;i<mccol.size();i++){
        StrawDigiMC mcdigi = mccol[i];

        const Straw& straw = tracker->getStraw( mcdigi.strawId() );
        auto const& sgsptr = mcdigi.earlyStrawGasStep();
        auto const& sgs = *sgsptr;
        auto const& sp = *sgs.simParticle();

        if ((sp.pdgId() == PDGCode::mu_minus || sp.pdgId() == PDGCode::mu_plus) && sp.creationCode() == 56){
          TwoLinePCA pca( straw.getMidPoint(), straw.getDirection(),
              GenVector::Hep3Vec(sgs.startPosition()), GenVector::Hep3Vec(sgs.endPosition()-sgs.startPosition()) );
          double true_doca = pca.dca();

          TwoLinePCA pca2( straw.getMidPoint(), straw.getDirection(),
              pppos[j], ppdir[j]);

          double mctrack_doca = pca2.dca();
          if (fabs(true_doca - mctrack_doca) < 0.5){
            count++;
            if (ppdir[j].y() != 0)
              ppintercept = pppos[j] - ppdir[j]*pppos[j].y()/ppdir[j].y();
            ppdirection = ppdir[j];
            if (ppdirection.y() > 0)
              ppdirection *= -1;
            double mctime = sgs.time();// - _ewMarkerOffset;
            double trajtime = (pca2.point2()-ppintercept).dot(ppdirection.unit())/299.9;
            mctime -= trajtime;
            avg_t0 += mctime;
          }
        }
      }
      if (count > max){
        max = count;
        _mcpos = ppintercept;
        _mcdir = ppdirection;
        if (ppdirection.y() != 0){
          ppdirection /= -1*ppdirection.y();
          ppintercept -= ppdirection*ppintercept.y()/ppdirection.y();
        }
        _mca0 = ppintercept.x();
        _mcb0 = ppintercept.z();
        _mca1 = ppdirection.x();
        _mcb1 = ppdirection.z();
        _mct0 = avg_t0/count;
      }
    }
    if (_mcdir.y() > 0)
      _mcdir *= -1;

    _mcnsh = 0;
    // now lets count the number of straw hits in this mc event
    // loop over combohits, for each combo hit get all the sub strawhits
    //
    for (size_t ich=0;ich<_chcol->size();ich++){
        const ComboHit &sh = _chcol->at(ich);
        const Straw& straw = tracker->getStraw( sh.strawId() );

        // shids.size() should be 1 if this is really all single StrawHits
        if (sh.nStrawHits() != 1){
          std::cout << "INCORRECT NUMBER OF StrawDigis " << sh.nStrawHits() << std::endl;
          continue;
        }
        const StrawDigiMC &mcdigi = _mcdigis->at(ich);

        auto const& sgsptr = mcdigi.earlyStrawGasStep();
        auto const& sgs = *sgsptr;
        auto const& sp = *sgs.simParticle();
        if ((sp.pdgId() == PDGCode::mu_minus || sp.pdgId() == PDGCode::mu_plus) && sp.creationCode() == 56){
          TwoLinePCA pca( straw.getMidPoint(), straw.getDirection(),
              GenVector::Hep3Vec(sgs.startPosition()), GenVector::Hep3Vec(sgs.endPosition()-sgs.startPosition()) );
          double true_doca = pca.dca();

          TwoLinePCA pca2( straw.getMidPoint(), straw.getDirection(),
              _mcpos, _mcdir);

          double mctrack_doca = pca2.dca();
          if (fabs(true_doca - mctrack_doca) < 0.5){
            _mcnsh += 1;
          }
        }
    }
  }

  void CosmicTrackDiag::hitDiag(const art::Event& event, StrawResponse const& srep, int its, const Tracker* tracker) {
    GaussianDriftFit gdf(*_chcol, srep, tracker);
    // loop over combohits
    for (size_t ich=0;ich<_chcol->size();ich++){
      const ComboHit &sh = _chcol->at(ich);
      const Straw& straw = tracker->getStraw( sh.strawId() );


      _hitminuitdoca = -999;
      _hitbackground = 1;
      _hittruedoca = -999;
      _hitmcdoca = -999;
      _hitmcdpocat = -999;
      _hitminuitdpocat = -999;
      _hitused = 0;
      _hittruelong = -999;
      _hitmclong = -999;
      _hitminuitlong = -999;
      _hitrecolong = -999;
      _hitrecolongrms = -999;
      _hittresid = -999;
      _hittresidrms = -999;
      _hitdresid = -999;
      _hitdresiderr = -999;
      _hitllike = 0;
      _hitrecodoca = -999;
      _hitambig = -999;
      _hitmcambig = -999;
      _hitedep = -999;
      _hitublong = -999;
      _hitublongrms = -999;
      _hitubdoca = -999;
      _hitubt0 = -999;
      _hitdz = -999;
      _hitubtresid = -999;
      _hitubtresidrms = -999;
      _hitdeltat = -999;
      _hittot[0] = -999;
      _hittot[1] = -999;
      _hitintclust = 0;


      _hitplane = sh.strawId().getPlane();
      _hitpanel = sh.strawId().getPanel();
      _hitstraw = sh.strawId().getStraw();
      _strawlen = straw.halfLength();

      _hitedep = sh.energyDep()*1000.;
      _hitrecolong = sh.wireDist();
      _hitrecolongrms = sh.wireRes();
      _hitrecodoca = sh.driftTime();

      _hitproptime = sh.propTime();
      _hittime = sh.time();

      if( _shcol){
        const StrawHit &strawh = _shcol->at(ich);
        _hitdeltat = strawh.time(StrawEnd::hv)-strawh.time(StrawEnd::cal);
        _hittot[0] = strawh.TOT(StrawEnd::cal);
        _hittot[1] = strawh.TOT(StrawEnd::hv);
      }


      CLHEP::Hep3Vector pcapoint2(0,0,0);
      if (_mcdiag){
        if (sh.nStrawHits() != 1){
          std::cout << "INCORRECT NUMBER OF StrawDigis " << sh.nStrawHits() << std::endl;
          continue;
        }

        const StrawDigiMC &mcdigi = _mcdigis->at(ich);

        auto const& sgsptr = mcdigi.earlyStrawGasStep();
        auto const& sgs = *sgsptr;
        auto const& sp = *sgs.simParticle();

        // get true ambiguity
        CLHEP::Hep3Vector mcsep = sgs.position()-straw.getMidPoint();
        CLHEP::Hep3Vector sgsdir = GenVector::Hep3Vec(sgs.momentum()).unit();
        CLHEP::Hep3Vector mcperp = (sgsdir.cross(straw.getDirection())).unit();
        double mcdperp = mcperp.dot(mcsep);
        _hitmcambig = mcdperp > 0 ? -1 : 1;

        if (sp.creationCode() == 56){
          _hitbackground = 0;
        }

        _hitmcpdg = sp.pdgId();
        _hitmccc = sp.creationCode();
        _hitmcmom = sgs.momvec().mag();

        // Get the actual DOCA of the MC step
        TwoLinePCA pca( straw.getMidPoint(), straw.getDirection(),
            GenVector::Hep3Vec(sgs.startPosition()), GenVector::Hep3Vec(sgs.endPosition()-sgs.startPosition()) );
        _hittruedoca = pca.dca();
        _hittruelong = (GenVector::Hep3Vec(sgs.startPosition())-straw.getMidPoint()).dot(straw.getDirection());

        // Get the DOCA from the MC straight line track
        TwoLinePCA pca1( straw.getMidPoint(), straw.getDirection(),
            _mcpos, _mcdir);
        pcapoint2 = pca.point2();
        _hitmcdoca = pca1.dca();
        // get the delta transverse distance between POCA of step and POCA of track
        _hitmcdpocat = sqrt((pca.point2()-pca1.point2()).mag2()-pow((pca.point2()-pca1.point2()).dot(straw.getDirection()),2));
        _hitmclong = (pca1.point1()-straw.getMidPoint()).dot(straw.getDirection());

        TwoLinePCA mcpca( _mcpos, _mcdir,
            GenVector::Hep3Vec(sgs.startPosition()), GenVector::Hep3Vec(sgs.endPosition()-sgs.startPosition()) );
        _hitmctrajtime = (mcpca.point1()-_mcpos).dot(_mcdir.unit())/299.9;
        if (mcpca.closeToParallel()){
          _hitmctrajtime = (GenVector::Hep3Vec(sgs.startPosition())-_mcpos).dot(_mcdir.unit())/299.9;
        }
        _hitmcdrifttime = srep.driftDistanceToTime(sh.strawId(), pca1.dca(), 0);
        _hitmcdrifttimeoffset = srep.driftTimeOffset(sh.strawId(), pca1.dca(), 0);
        _hitmctresid = sh.time() - (_mct0 + _hitmctrajtime + sh.propTime() + _hitmcdrifttime + _hitmcdrifttimeoffset);

      }


      if (its >= 0){
        auto tseed = _tscol->at(its);
        auto minuitpos = GenVector::Hep3Vec(tseed._track.MinuitEquation.Pos);
        auto minuitdir = GenVector::Hep3Vec(tseed._track.MinuitEquation.Dir.unit());
        // convert to y orientation
        if (minuitdir.y() > 0)
          minuitdir *= -1;
        if (minuitdir.y() != 0)
          minuitpos -= minuitdir*minuitpos.y()/minuitdir.y();


        size_t found_k = 0;
        bool found = false;
        for (size_t k=0;k<tseed._straw_chits.size();k++){
          if (tseed._straw_chits[k].strawId() == sh.strawId() && tseed._straw_chits[k].time() == sh.time()){
            found = true;
            found_k = k;
            break;
          }
        }
        if (found)
          _hitused = 1;

        auto tclust = tseed._timeCluster;
        ComboHitCollection tchits;
        std::vector<ComboHitCollection::const_iterator> chids;
        _phcol->fillComboHits(tclust->hits(), chids);
        for (auto const& it : chids){
          if (it[0].strawId() == sh.strawId() && it[0].time() == sh.time()){
            _hitintclust = 1;
            break;
          }
        }


        TwoLinePCA pca3( straw.getMidPoint(), straw.getDirection(),
            minuitpos, minuitdir);
        _hitminuitdoca = pca3.dca();
        _hitminuitlong = (pca3.point1()-straw.getMidPoint()).dot(straw.getDirection());

        CLHEP::Hep3Vector sep = (pca3.point2()-straw.getMidPoint()); // vector from PCA to straw
        CLHEP::Hep3Vector perp = (minuitdir.cross(straw.getDirection())).unit();
        double dperp = perp.dot(sep);
        _hitambig = dperp > 0 ? -1 : 1;
        _hitdz = sep.z();

        double trunc_hitminuitlong = _hitminuitlong;
        if (fabs(_hitminuitlong) > straw.halfLength()){
          _hitllike += pow((fabs(_hitminuitlong)-straw.halfLength()),2);
          trunc_hitminuitlong = std::copysign(straw.halfLength(),_hitminuitlong);
        }
//        _hitrecolongrms = srep.wpRes(sh.energyDep() * 1000., fabs(trunc_hitminuitlong));
        _hitrecolongrms = srep.wpRes(sh.energyDep() * 1000., fabs(sh.wireDist()));

        //FIXME ignoring likelihood contribution from beyond wire
        //_hitllike += pow(trunc_hitminuitlong-_hitrecolong,2)/pow(_hitrecolongrms,2);
        double trunc_wireDist = std::copysign(std::min(straw.halfLength(),fabs(_hitrecolong)),_hitrecolong);
        _hitllike += pow(trunc_hitminuitlong-trunc_wireDist,2)/pow(_hitrecolongrms,2);


        _hitdrifttime = srep.driftDistanceToTime(sh.strawId(), pca3.dca(), 0);
        _hitdrifttimeoffset = srep.driftTimeOffset(sh.strawId(), pca3.dca(), 0);
        _hittrajtime = ((pca3.point2() - minuitpos).dot(minuitdir))/299.9;
        _hittresidrms = srep.driftTimeError(sh.strawId(), pca3.dca(), 0);

        //double hit_t0 = sh.time() - sh.propTime() - _hittrajtime - _hitdrifttime - _hitdrifttimeoffset;

        //_hittresid = hit_t0-tseed._t0._t0;
        _hittresid = -1*gdf.TimeResidual(sh,tseed);
        _hitllike += pow(_hittresid,2)/pow(_hittresidrms,2);

        _hitdresid = gdf.DOCAresidual(sh,tseed);
        _hitdresiderr = gdf.DOCAresidualError(sh,tseed);

        if (_printlevel > 1){
          if (_hitused){
            std::cout << "   " << _hitminuitdoca << " " << _hitminuitlong << " " << _hitmcdoca << " " << _hitmclong << " " << sh.strawId() << std::endl;
          }else{
            std::cout << " X " << _hitminuitdoca << " " << _hitminuitlong << " " << _hitmcdoca << " " << _hitmclong << " " << sh.strawId() << std::endl;
          }
        }


        if (_mcdiag){
          _hitminuitdpocat = sqrt((pcapoint2-pca3.point2()).mag2()-pow((pcapoint2-pca3.point2()).dot(straw.getDirection()),2));
        }

        // redo the full fit without this hit to get unbiased residual
        if (_ubdiag && found){
          auto dir = tseed._track.FitEquation.Dir;
          auto intercept = tseed._track.FitEquation.Pos;
          dir /= -1 * dir.y();
          intercept -= dir * intercept.y() / dir.y();

          std::vector<double> errors(5, 0);
          std::vector<double> pars(5, 0);

          pars[0] = intercept.x();
          pars[1] = intercept.z();
          pars[2] = dir.x();
          pars[3] = dir.z();
          pars[4] = tseed._track.FitParams.T0;
          errors[0] = tseed._track.FitParams.Covarience.sigA0;
          errors[1] = tseed._track.FitParams.Covarience.sigB0;
          errors[2] = tseed._track.FitParams.Covarience.sigA1;
          errors[3] = tseed._track.FitParams.Covarience.sigB1;
//          errors[4] = tseed._t0.t0Err();
          errors[4] = 1.0;

          std::vector<double> temp_cov;
          bool temp_converged;
          // Define the PDF used by Minuit:
          GaussianDriftFit fit2(tseed._straw_chits, srep, tracker);
          fit2.setExcludeHit(found_k);
          MinuitDriftFitter::DoDriftTimeFit(pars, errors, temp_cov, temp_converged, fit2, 10, 0, 0.5);
          if (temp_converged){
            CLHEP::Hep3Vector temppos(pars[0],0,pars[1]);
            CLHEP::Hep3Vector tempdir(pars[2],-1,pars[3]);
            tempdir = tempdir.unit();
            TwoLinePCA ubpca( straw.getMidPoint(), straw.getDirection(),
                temppos, tempdir);
            _hitublong = (ubpca.point1()-straw.getMidPoint()).dot(straw.getDirection());
            _hitubdoca = ubpca.dca();
            _hitubt0 = pars[4];
            double uhitdrifttime = srep.driftDistanceToTime(sh.strawId(), ubpca.dca(), 0);
            double uhitdrifttimeoffset = srep.driftTimeOffset(sh.strawId(), ubpca.dca(), 0);
            double uhittrajtime = ((ubpca.point2() - temppos).dot(tempdir))/299.9;
            _hitubtresidrms = srep.driftTimeError(sh.strawId(), ubpca.dca(), 0);
            double uhit_t0 = sh.time() - sh.propTime() - uhittrajtime - uhitdrifttime - uhitdrifttimeoffset;
            _hitubtresid = uhit_t0-_hitubt0;
          }
        }
      }
      _hitT->Fill();
    }
  }

}

using mu2e::CosmicTrackDiag;
DEFINE_ART_MODULE(CosmicTrackDiag)
