#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <cmath>

// Cosmic Tracks:
#include "CosmicReco/inc/CosmicTrackFit.hh"
#include "CosmicReco/inc/CosmicTrackFinderData.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "CosmicReco/inc/CosmicTrackMCInfo.hh"
#include "DataProducts/inc/EventWindowMarker.hh"

//Mu2e Data Prods:
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "DataProducts/inc/XYZVec.hh"

//Utilities
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
#include "CosmicReco/inc/DriftFitUtils.hh"
#include "Mu2eUtilities/inc/ParametricFit.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

// Mu2e diagnostics
#include "TrkDiag/inc/ComboHitInfo.hh"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
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
        fhicl::Atom<bool> mcdiag{Name("mcdiag"), Comment("set on for MC info"),true};
        fhicl::Atom<art::InputTag> chtag{Name("ComboHitCollection"),Comment("tag for combo hit collection")};
        fhicl::Atom<art::InputTag> tctag{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
        fhicl::Atom<art::InputTag> costag{Name("CosmicTrackSeedCollection"),Comment("tag for cosmci track seed collection")};
        fhicl::Atom<art::InputTag> mcdigistag{Name("StrawDigiMCCollection"),Comment("StrawDigi collection tag"),"makeSD"};
        fhicl::Atom<art::InputTag> ewMarkerTag{Name("EventWindowMarkerLabel"),Comment("Event window marker tag"),"EWMProducer"};
        fhicl::Table<SimParticleTimeOffset::Config> toff{Name("TimeOffsets"), Comment("Sim particle time offset ")};
      };
      typedef art::EDAnalyzer::Table<Config> Parameters;

      explicit CosmicTrackDiag(const Parameters& conf);
      virtual ~CosmicTrackDiag();
      virtual void beginJob() override;
      virtual void beginRun(const art::Run& r) override;
      virtual void analyze(const art::Event& e) override;
    private: 

      Config _conf;

      int  _diag;
      bool _mcdiag;
      std::ofstream outputfile;
      art::InputTag   _chtag;//combo
      art::InputTag   _tctag;//timeclusters
      art::InputTag   _costag;//Striaght tracks
      art::InputTag   _mcdigistag; //MC digis
      art::InputTag _ewMarkerTag;
      SimParticleTimeOffset _toff;
      const ComboHitCollection* _chcol;
      const CosmicTrackSeedCollection* _coscol;
      const TimeClusterCollection* _tccol;
      const StrawDigiMCCollection* _mcdigis;
      CLHEP::Hep3Vector _mcpos, _mcdir;

      Float_t _ewMarkerOffset;
      const Tracker* tracker;

      //TTree Info:
      TTree* _trackT;
      TTree* _hitT;

      // track tree 
      Int_t _evt; 
      Int_t _ntrack;
      Int_t _nsh, _nch; // # associated straw hits / event
      Int_t _ntc; // # clusters/event
      Int_t _n_panels; // # panels
      Int_t _n_stations; // # stations
      Int_t _n_planes; // # stations
      Float_t _mct0;
      Int_t _mcnsh;
      Int_t _hitsok, _chifitok, _chifitconverged, _minuitok, _minuitconverged;

      // hit tree 
      Float_t _chidoca, _chiangle, _seeddoca, _seedangle, _minuitdoca, _minuitangle;
      Float_t _hittruedoca, _hitmcdoca, _hitseeddoca, _hitminuitdoca;
      Float_t _hitmcdpocat, _hitseeddpocat, _hitminuitdpocat;
      Int_t _hitbackground;


      void GetMCTrack(const StrawDigiMCCollection& mccol);
      void hitDiag(const art::Event& event);
      bool findData(const art::Event& evt);
  };

  CosmicTrackDiag::CosmicTrackDiag(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _diag (conf().diag()),
    _mcdiag (conf().mcdiag()),
    _chtag (conf().chtag()),
    _tctag (conf().tctag()),
    _costag (conf().costag()),
    _mcdigistag (conf().mcdigistag()),
    _ewMarkerTag (conf().ewMarkerTag()),
    _toff (conf().toff())
  {
    if(_mcdiag){
      for (auto const& tag : conf().toff().inputs()) {
        consumes<SimParticleTimeMap>(tag);
      }
    }
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
      _trackT->Branch("ntrack",&_ntrack,"ntrack/I");
      _trackT->Branch("StrawHitsInEvent", &_nsh, "StrawHitsInEvent/I");
      _trackT->Branch("ComboHitsInEvent", &_nch, "ComboHitsInEvent/I");
      _trackT->Branch("PanelsCrossedInEvent", &_n_panels, "PanelsCrossedInEvent/I");
      _trackT->Branch("PlanesCrossedInEvent", &_n_planes, "PlanesCrossedInEvent/I");
      _trackT->Branch("StatonsCrossedInEvent", &_n_stations, "StationsCrossedInEvent/I");
      _trackT->Branch("TimeClustersInEvent", &_ntc, "TimeClusterInEvent/I"); 
      _trackT->Branch("hitsok",&_hitsok,"hitsok/I");
      _trackT->Branch("chifitok",&_chifitok,"chifitok/I");
      _trackT->Branch("chifitconverged",&_chifitconverged,"chifitconverged/I");
      _trackT->Branch("minuitok",&_minuitok,"minuitok/I");
      _trackT->Branch("minuitconverged",&_minuitconverged,"minuitconverged/I");
      if (_mcdiag){
        _trackT->Branch("mcnsh",&_mcnsh,"mcnsh/I");
        _trackT->Branch("mct0",&_mct0,"mct0/F");
        _trackT->Branch("chidoca",&_chidoca,"chidoca/F");
        _trackT->Branch("chiangle",&_chiangle,"chiangle/F");
        _trackT->Branch("seeddoca",&_seeddoca,"seeddoca/F");
        _trackT->Branch("seedangle",&_seedangle,"seedangle/F");
        _trackT->Branch("minuitdoca",&_minuitdoca,"minuitdoca/F");
        _trackT->Branch("minuitangle",&_minuitangle,"minuitangle/F");
      }

      _hitT=tfs->make<TTree>("hitT","Hit tree");
      _hitT->Branch("seeddoca",&_hitseeddoca,"seeddoca/F");
      _hitT->Branch("minuitdoca",&_hitminuitdoca,"minuitdoca/F");
      if (_mcdiag){
        _hitT->Branch("background",&_hitbackground,"background/I");
        _hitT->Branch("truedoca",&_hittruedoca,"truedoca/F");
        _hitT->Branch("mcdoca",&_hitmcdoca,"mcdoca/F");
        _hitT->Branch("mcdpocat",&_hitmcdpocat,"mcdpocat/F");
        _hitT->Branch("seeddpocat",&_hitseeddpocat,"seeddpocat/F");
        _hitT->Branch("minuitdpocat",&_hitminuitdpocat,"minuitdpocat/F");
      }
    }
  } 

  void CosmicTrackDiag::beginRun(const art::Run& run){
    mu2e::GeomHandle<mu2e::Tracker> th;
    tracker = th.get();
  }

  void CosmicTrackDiag::analyze(const art::Event& event) {

    _evt = event.id().event();  // add event id
    if(!findData(event)) // find data
      throw cet::exception("RECO")<<"No Time Clusters in event"<< endl; 

    //find time clusters:
    _ntc = _tccol->size();
    _nch = _chcol->size();


    if (_mcdiag){
      GetMCTrack(*_mcdigis);
      _mct0 = 0;
      _mcnsh = 0;
      for (size_t i=0;i<_mcdigis->size();i++){
        StrawDigiMC mcdigi = _mcdigis->at(i); 
        auto const& sgsptr = mcdigi.earlyStrawGasStep();
        auto const& sgs = *sgsptr;
        auto const& sp = *sgs.simParticle();
        if (sp.creationCode() == 56){
          double mctime = sgs.time() + _toff.totalTimeOffset(sgs.simParticle()) -_ewMarkerOffset; 
          TwoLinePCA mcpca( _mcpos, _mcdir,
              Geom::Hep3Vec(sgs.startPosition()), Geom::Hep3Vec(sgs.endPosition()-sgs.startPosition()) );
          double trajtime = (mcpca.point1()-_mcpos).dot(_mcdir.unit())/299.9;
          if (mcpca.closeToParallel()){
            trajtime = (Geom::Hep3Vec(sgs.startPosition())-_mcpos).dot(_mcdir.unit())/299.9;
          }
          mctime -= trajtime;
          _mct0 += mctime;
          _mcnsh += 1;
        }
      }
      _mct0 /= _mcnsh;
    }

    std::vector<int> panels, planes, stations;
    _nsh = 0;
    for(size_t ich = 0;ich < _chcol->size(); ++ich){
      ComboHit const& chit =(*_chcol)[ich];
      _nsh += chit.nStrawHits(); 
      panels.push_back(chit.strawId().panel());
      planes.push_back(chit.strawId().plane());
      stations.push_back(chit.strawId().station());
    }
    _n_panels = std::set<float>( panels.begin(), panels.end() ).size();
    _n_planes = std::set<float>( planes.begin(), planes.end() ).size();
    _n_stations = std::set<float>( stations.begin(), stations.end() ).size();


    //loop over tracks
    for(size_t ist = 0;ist < _coscol->size(); ++ist){
      _chidoca = -999;
      _chiangle = -999;
      _seeddoca = -999;
      _seedangle = -999;
      _minuitdoca = -999;
      _minuitangle = -999;
      _hitsok = 0;
      _chifitok = 0;
      _chifitconverged = 0;
      _minuitok = 0;
      _minuitconverged = 0;

      _ntrack = ist;

      CosmicTrackSeed sts =(*_coscol)[ist];
      CosmicTrack st = sts._track;

      TrkFitFlag const& status = sts._status;
      if (status.hasAllProperties(TrkFitFlag::hitsOK))
        _hitsok = 1;
      if (status.hasAllProperties(TrkFitFlag::helixOK) && status.hasAllProperties(TrkFitFlag::helixConverged))
        _chifitok = 1;
      if (status.hasAllProperties(TrkFitFlag::seedOK) && status.hasAllProperties(TrkFitFlag::seedConverged))
        _chifitconverged = 1;
      if (status.hasAllProperties(TrkFitFlag::kalmanOK))
        _minuitok = 1;
      if (status.hasAllProperties(TrkFitFlag::kalmanConverged))
        _minuitconverged = 1;
      //      if (!status.hasAllProperties(TrkFitFlag::helixOK) ){ continue;}
      //      if(st.converged == false or st.minuit_converged  == false) { continue;}

      if(_mcdiag){
        auto chipos = CLHEP::Hep3Vector(st.FitEquationXYZ.Pos.X(),st.FitEquationXYZ.Pos.Y(),0);
        auto chidir = CLHEP::Hep3Vector(st.FitEquationXYZ.Dir.X(),st.FitEquationXYZ.Dir.Y(),1).unit();

        auto seedpos = CLHEP::Hep3Vector(st.FitParams.A0,st.FitParams.B0,0);
        auto seeddir = CLHEP::Hep3Vector(st.FitParams.A1,st.FitParams.B1,1).unit();

        auto minuitpos = CLHEP::Hep3Vector(st.MinuitFitParams.A0,st.MinuitFitParams.B0,0);
        auto minuitdir = CLHEP::Hep3Vector(st.MinuitFitParams.A1,st.MinuitFitParams.B1,1).unit();

        // get angle between reco and mc track, and doca between them
        TwoLinePCA chipca( _mcpos, _mcdir, chipos, chidir);
        _chidoca = chipca.dca(); 
        _chiangle = _mcdir.dot(chidir);
        TwoLinePCA seedpca( _mcpos, _mcdir, seedpos, seeddir);
        _seeddoca = seedpca.dca(); 
        _seedangle = _mcdir.dot(seeddir);
        TwoLinePCA minuitpca( _mcpos, _mcdir, minuitpos, minuitdir);
        _minuitdoca = minuitpca.dca(); 
        _minuitangle = _mcdir.dot(minuitdir);
      }

      //      for(auto const& tseed : *_coscol) {   
      //        TrkFitFlag const& status = tseed._status;
      //        _hitsOK = status.hasAllProperties(TrkFitFlag::hitsOK);
      //        _StraightTrackOK = status.hasAllProperties(TrkFitFlag::helixOK);
      //        _StraightTrackConverged = status.hasAllProperties(TrkFitFlag::helixConverged);
      //        _StraightTrackInit = status.hasAllProperties(TrkFitFlag::circleInit);
      //      }


      _trackT->Fill();
    }//end analyze

    if (_diag > 1)
      hitDiag(event);

  }

  bool CosmicTrackDiag::findData(const art::Event& evt){
    _chcol = 0; 
    _tccol = 0;
    _coscol = 0; 
    _ewMarkerOffset = 0;
    auto chH = evt.getValidHandle<ComboHitCollection>(_chtag);
    _chcol = chH.product();
    auto tcH = evt.getValidHandle<TimeClusterCollection>(_tctag);
    _tccol =tcH.product();
    auto stH = evt.getValidHandle<CosmicTrackSeedCollection>(_costag);
    _coscol =stH.product();
    auto ewMarkerHandle = evt.getValidHandle<EventWindowMarker>(_ewMarkerTag);
    auto ewMarker = ewMarkerHandle.product();
    _ewMarkerOffset = ewMarker->timeOffset();
    if(_mcdiag){
      _mcdigis=0;
      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigistag);
      _mcdigis = mcdH.product();
      _toff.updateMap(evt);
    }
    return _chcol != 0 && _tccol!=0 && _coscol !=0 && (_mcdigis != 0 || !_mcdiag);
  }


  void CosmicTrackDiag::GetMCTrack(const StrawDigiMCCollection& mccol) {
    // get all possible directions
    std::vector<CLHEP::Hep3Vector> pppos;
    std::vector<CLHEP::Hep3Vector> ppdir;
    for (size_t i=0;i<mccol.size();i++){
      StrawDigiMC mcdigi = mccol[i]; 
      auto const& sgsptr = mcdigi.earlyStrawGasStep();
      auto const& sgs = *sgsptr;
      auto const& sp = *sgs.simParticle();
      auto posi = Geom::Hep3Vec(sgs.startPosition());
      if ((sp.pdgId() == 13 || sp.pdgId() == -13) && sp.creationCode() == 56){
        for (size_t j=i+1;j<mccol.size();j++){
          StrawDigiMC jmcdigi = mccol[j]; 
          auto const& jsgsptr = jmcdigi.earlyStrawGasStep();
          auto const& jsgs = *jsgsptr;
          auto const& jsp = *jsgs.simParticle();
          auto posj = Geom::Hep3Vec(jsgs.startPosition());
          if ((jsp.pdgId() == 13 || jsp.pdgId() == -13) && jsp.creationCode() == 56){
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
      for (size_t i=0;i<mccol.size();i++){
        StrawDigiMC mcdigi = mccol[i]; 

        const Straw& straw = tracker->getStraw( mcdigi.strawId() );
        auto const& sgsptr = mcdigi.earlyStrawGasStep();
        auto const& sgs = *sgsptr;
        auto const& sp = *sgs.simParticle();

        if ((sp.pdgId() == 13 || sp.pdgId() == -13) && sp.creationCode() == 56){
          TwoLinePCA pca( straw.getMidPoint(), straw.getDirection(),
              Geom::Hep3Vec(sgs.startPosition()), Geom::Hep3Vec(sgs.endPosition()-sgs.startPosition()) );
          double true_doca = pca.dca(); 

          TwoLinePCA pca2( straw.getMidPoint(), straw.getDirection(),
              pppos[j], ppdir[j]);

          double mctrack_doca = pca2.dca(); 
          if (fabs(true_doca - mctrack_doca) < 0.5)
            count++;
        }
      }
      if (count > max){
        max = count;
        _mcpos = pppos[j];
        _mcdir = ppdir[j];
      }
    }
    if (_mcdir.y() > 0)
      _mcdir *= -1;
  }

  void CosmicTrackDiag::hitDiag(const art::Event& event) {
    // loop over combohits, for each combo hit get all the sub strawhits
    for (size_t ich=0;ich<_chcol->size();ich++){
      // get the straw hits for this combohit
      std::vector<ComboHitCollection::const_iterator> chids;  
      ComboHitCollection strawhits;
      _chcol->fillComboHits(event, ich, chids); 
      for (auto const& it : chids){
        const mu2e::ComboHit chit = it[0];
        strawhits.push_back(chit);
      }
      // get the StrawDigi indices for this combohit
      std::vector<StrawDigiIndex> shids;
      _chcol->fillStrawDigiIndices(event,ich,shids);


      for (size_t k=0;k<shids.size();k++){
        _hitseeddoca = -999;
        _hitminuitdoca = -999;
        _hitbackground = 1;
        _hittruedoca = -999;
        _hitmcdoca = -999;
        _hitmcdpocat = -999;
        _hitminuitdpocat = -999;


        const ComboHit &sh = strawhits[k];
        const Straw& straw = tracker->getStraw( sh.strawId() );

        for(size_t ist = 0;ist < _coscol->size(); ++ist){
          CosmicTrackSeed sts =(*_coscol)[ist];
          CosmicTrack st = sts._track;

          auto seedpos = CLHEP::Hep3Vector(st.FitParams.A0,st.FitParams.B0,0);
          auto seeddir = CLHEP::Hep3Vector(st.FitParams.A1,st.FitParams.B1,1).unit();

          auto minuitpos = CLHEP::Hep3Vector(st.MinuitFitParams.A0,st.MinuitFitParams.B0,0);
          auto minuitdir = CLHEP::Hep3Vector(st.MinuitFitParams.A1,st.MinuitFitParams.B1,1).unit();

          TwoLinePCA pca4( straw.getMidPoint(), straw.getDirection(),
              seedpos, seeddir);
          _hitseeddoca = pca4.dca();
//          if (_mcdiag)
//            _hitseeddpocat = sqrt((pca.point2()-pca4.point2()).mag2()-pow((pca.point2()-pca4.point2()).dot(straw.getDirection()),2));

          TwoLinePCA pca5( straw.getMidPoint(), straw.getDirection(),
              minuitpos, minuitdir);
          _hitminuitdoca = pca5.dca();
//          if (_mcdiag)
//            _hitminuitdpocat = sqrt((pca.point2()-pca5.point2()).mag2()-pow((pca.point2()-pca5.point2()).dot(straw.getDirection()),2));
        }


        if (_mcdiag){
          const StrawDigiMC &mcdigi = _mcdigis->at(shids[k]);

          auto const& sgsptr = mcdigi.earlyStrawGasStep();
          auto const& sgs = *sgsptr;
          auto const& sp = *sgs.simParticle();

          if (sp.creationCode() == 56){
            _hitbackground = 0;
          }

          // Get the actual DOCA of the MC step
          TwoLinePCA pca( straw.getMidPoint(), straw.getDirection(),
              Geom::Hep3Vec(sgs.startPosition()), Geom::Hep3Vec(sgs.endPosition()-sgs.startPosition()) );
          _hittruedoca = pca.dca(); 
          //double true_long = (Geom::Hep3Vec(sgs.startPosition())-straw.getMidPoint()).dot(straw.getDirection());

          // Get the DOCA from the MC straight line track
          TwoLinePCA pca1( straw.getMidPoint(), straw.getDirection(),
              _mcpos, _mcdir);
          _hitmcdoca = pca1.dca(); 
          // get the delta transverse distance between POCA of step and POCA of track
          _hitmcdpocat = sqrt((pca.point2()-pca1.point2()).mag2()-pow((pca.point2()-pca1.point2()).dot(straw.getDirection()),2));
        }
        _hitT->Fill();
      }
    }
  }

}

using mu2e::CosmicTrackDiag;
DEFINE_ART_MODULE(CosmicTrackDiag);
