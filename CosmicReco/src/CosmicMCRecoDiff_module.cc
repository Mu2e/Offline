//Author: S Middleton
//Date:  2021
//Purpose: For getting seed cov estimates

#include <iostream>
#include <string>
#include <cmath>

// Cosmic Tracks:
#include "Offline/CosmicReco/inc/CosmicTrackFit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/CosmicReco/inc/CosmicTrackMCInfo.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/CosmicReco/inc/ComboHitInfoMC.hh"

//Mu2e Data Prods:
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"
//Utilities
#include "Offline/CosmicReco/inc/DriftFitUtils.hh"
#include "Offline/Mu2eUtilities/inc/ParametricFit.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/Mu2eUtilities/inc/BuildLinearFitMatrixSums.hh"
#include "Offline/Mu2eUtilities/inc/CosmicTrackUtils.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA_XYZ.hh"

// Mu2e diagnostics
#include "Offline/GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

// ROOT incldues
#include "TTree.h"

//Geom
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"

using namespace std;

namespace mu2e
{
        class CosmicMCRecoDiff : public art::EDAnalyzer {
    public:
      struct Config{
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> chtag{Name("ComboHitCollection"),Comment("tag for combo hit collection")};
        fhicl::Atom<art::InputTag> tctag{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
        fhicl::Atom<art::InputTag> costag{Name("CosmicTrackSeedCollection"),Comment("tag for cosmci track seed collection")};
        fhicl::Atom<art::InputTag> mcdigistag{Name("StrawDigiMCCollection"),Comment("StrawDigi collection tag"),"makeSD"};
      };
      typedef art::EDAnalyzer::Table<Config> Parameters;

      explicit CosmicMCRecoDiff(const Parameters& conf);
      virtual ~CosmicMCRecoDiff();
      virtual void beginJob() override;
      virtual void beginRun(const art::Run& r) override;
      virtual void analyze(const art::Event& e) override;
      virtual void endJob() override;
    private:

      Config _conf;

      std::ofstream outputfile;
      art::InputTag   _chtag; //ComboHits
      art::InputTag   _tctag; //Timeclusters
      art::InputTag   _costag; //Straight tracks
      art::InputTag   _mcdigistag; //MC Digis

      const ComboHitCollection* _chcol;
      const CosmicTrackSeedCollection* _coscol;
      const TimeClusterCollection* _tccol;
      const StrawDigiMCCollection* _mcdigis;
      CosmicTrackMCInfo trueinfo;

      TTree* _cosmic_tree;

      Float_t _RecoA0;
      Float_t _RecoA1;
      Float_t _RecoB1;
      Float_t _RecoB0;

      Float_t _Recod0;
      Float_t _Recoz0;
      Float_t _RecoPhi0;
      Float_t _RecoCosT;

      Float_t _RecoErrorA0;
      Float_t _RecoErrorA1;
      Float_t _RecoErrorB1;
      Float_t _RecoErrorB0;

      Float_t _TrueA0;
      Float_t _TrueA1;
      Float_t _TrueB1;
      Float_t _TrueB0;

      Float_t _TrueErrorA0;
      Float_t _TrueErrorA1;
      Float_t _TrueErrorB1;
      Float_t _TrueErrorB0;

      Float_t _Trued0;
      Float_t _Truez0;
      Float_t _TruePhi0;
      Float_t _TrueCosT;

      Float_t _Diffd0;
      Float_t _Diffz0;
      Float_t _DiffPhi0;
      Float_t _DiffCosT;

      Int_t _evt;

      ProditionsHandle<Tracker> _alignedTracker_h;
      ProditionsHandle<StrawResponse> _strawResponse_h;
      Int_t _strawid;
      vector<ComboHitInfoMC> _chinfomc;
      bool findData(const art::Event& evt);
      std::tuple <double,double, double, double, double> GetMCTrack(const art::Event& event, const StrawDigiMCCollection& mccol);
    };

    CosmicMCRecoDiff::CosmicMCRecoDiff(const Parameters& conf) :
      art::EDAnalyzer(conf),
      _chtag (conf().chtag()),
      _tctag (conf().tctag()),
      _costag (conf().costag()),
      _mcdigistag (conf().mcdigistag())
    {}

    CosmicMCRecoDiff::~CosmicMCRecoDiff(){}

    void CosmicMCRecoDiff::beginJob() {
      art::ServiceHandle<art::TFileService> tfs;
      _cosmic_tree=tfs->make<TTree>("cosmic_tree"," Diagnostics for Cosmic Track Fitting");

      _cosmic_tree->Branch("RecoA0",&_RecoA0,"RecoA0/F");
      _cosmic_tree->Branch("RecoA1",&_RecoA1,"RecoA1/F");
      _cosmic_tree->Branch("RecoB0",&_RecoA0,"RecoB0/F");
      _cosmic_tree->Branch("RecoB1",&_RecoA1,"RecoB1/F");

      _cosmic_tree->Branch("Recod0",&_Recod0,"Recod0/F");
      _cosmic_tree->Branch("Recoz0",&_Recoz0,"Recoz0/F");
      _cosmic_tree->Branch("RecoPhi0",&_RecoPhi0,"RecoPhi0/F");
      _cosmic_tree->Branch("RecoCosT",&_RecoCosT,"RecoCosT/F");

      _cosmic_tree->Branch("RecoErrorA0",&_RecoErrorA0,"RecoErrorA0/F");
      _cosmic_tree->Branch("RecoErrorA1",&_RecoErrorA1,"RecoErrorA1/F");
      _cosmic_tree->Branch("RecoErrorB0",&_RecoErrorA0,"RecoErrorB0/F");
      _cosmic_tree->Branch("RecoErrorB1",&_RecoErrorA1,"RecoErrorB1/F");

      _cosmic_tree->Branch("TrueA0",&_TrueA0,"TrueA0/F");
      _cosmic_tree->Branch("TrueA1",&_TrueA1,"TrueA1/F");
      _cosmic_tree->Branch("TrueB0",&_TrueA0,"TrueB0/F");
      _cosmic_tree->Branch("TrueB1",&_TrueA1,"TrueB1/F");

      _cosmic_tree->Branch("Trued0",&_Trued0,"Trued0/F");
      _cosmic_tree->Branch("Truez0",&_Truez0,"Truez0/F");
      _cosmic_tree->Branch("TruePhi0",&_TruePhi0,"TruePhi0/F");
      _cosmic_tree->Branch("TrueCosT",&_TrueCosT,"TrueCosT/F");

      _cosmic_tree->Branch("TrueErrorA0",&_RecoErrorA0,"TrueErrorA0/F");
      _cosmic_tree->Branch("TrueErrorA1",&_RecoErrorA1,"TrueErrorA1/F");
      _cosmic_tree->Branch("TrueErrorB0",&_RecoErrorA0,"TrueErrorB0/F");
      _cosmic_tree->Branch("TrueErrorB1",&_RecoErrorA1,"TrueErrorB1/F");

      _cosmic_tree->Branch("Diffd0",&_Diffd0,"Diffd0/F");
      _cosmic_tree->Branch("Diffz0",&_Diffz0,"Diffz0/F");
      _cosmic_tree->Branch("DiffPhi0",&_DiffPhi0,"DiffPhi0/F");
      _cosmic_tree->Branch("DiffCosT",&_DiffCosT,"DiffCosT/F");
    }

    void CosmicMCRecoDiff::beginRun(const art::Run& run){}

      void CosmicMCRecoDiff::analyze(const art::Event& event) {


      _evt = event.id().event();

      if(!findData(event))
      //throw cet::exception("RECO")<<"No Time Clusters in event"<< endl;
      std::cout<<"Size "<<_coscol->size()<<std::endl;
      for(unsigned ist = 0;ist < _coscol->size(); ++ist){
        std::cout<<"Looping "<<std::endl;
        CosmicTrackSeed sts =(*_coscol)[ist];
        CosmicTrack st = sts._track;
        //double t0 = sts._t0.t0();
        //TrkFitFlag const& status = sts._status;

        //if (!status.hasAllProperties(TrkFitFlag::helixOK) ){ continue; }
        //if(st.converged == false or st.minuit_converged  == false) { continue; }

        _RecoA0=(st.MinuitParams.A0);
        _RecoA1=(st.MinuitParams.A1);
        _RecoB1=(st.MinuitParams.B1);
        _RecoB0=(st.MinuitParams.B0);

        // Get KinKal:
        std::tuple <double, double, double, double, double, double> KinKalParams = KinKalTrackParams(st);
        _Recod0 = get<0>(KinKalParams);
        _Recoz0 = get<2>(KinKalParams);
        _RecoCosT = get<3>(KinKalParams);
        _RecoPhi0 = get<1>(KinKalParams);
        std::cout<<_Recod0<<" "<<_Recoz0<<" "<<_RecoCosT<<" "<<_RecoPhi0<<std::endl;
        _RecoErrorA0=(st.MinuitParams.deltaA0);
        _RecoErrorA1=(st.MinuitParams.deltaA1);
        _RecoErrorB1=(st.MinuitParams.deltaB1);
        _RecoErrorB0=(st.MinuitParams.deltaB0);

              //StrawDigiMC hitP1;
              //StrawDigiMC first = (*_mcdigis)[0];

              XYZVectorF zpos(0.,0.,0);
        XYZVectorF  zdir(0.,0.,1.);
        std::tuple <double,double, double, double, double> info = GetMCTrack(event, *_mcdigis);
              XYZVectorF pos0(get<0>(info),0, get<2>(info));//a0,0,b0
              XYZVectorF dir(get<1>(info), -1, get<3>(info));//a1,-1,b1

        TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(pos0, dir, zpos, zdir);
        XYZVectorF POCA = PCA.point1()-PCA.point2();
        double DOCA = PCA.dca();
        double amsign = copysign(1.0, -(zdir.Cross(POCA)).Dot(dir));

        _TrueA0 = get<0>(info);
        _TrueA1 = get<1>(info);
        _TrueB0 = get<2>(info);
        _TrueB1 = get<3>(info);

        _Trued0 = amsign*DOCA;
        _Truez0 = PCA.point1().Z();
        _TrueCosT = dir.Z();
        _TruePhi0 = dir.Phi();

        _Diffd0 = _Recod0 - _Trued0;
        _Diffz0 = _Recoz0 - _Truez0;
        _DiffCosT =  _RecoCosT - _TrueCosT;
        _DiffPhi0 =  _RecoPhi0 - _TruePhi0;

        _cosmic_tree->Fill();

      }

    }

 std::tuple <double,double, double, double, double> CosmicMCRecoDiff::GetMCTrack(const art::Event& event, const StrawDigiMCCollection& mccol) {
    // get all possible directions
    double _mca0 = 0;
    double _mca1 = 0;
    double _mcb0 = 0;
    double _mcb1 = 0;
    double _mct0 = 0;
    CLHEP::Hep3Vector _mcpos, _mcdir;
    std::vector<CLHEP::Hep3Vector> pppos;
    std::vector<CLHEP::Hep3Vector> ppdir;
    const Tracker *tracker = _alignedTracker_h.getPtr(event.id()).get();
    //StrawResponse const& srep = _strawResponse_h.get(event.id());

    for (size_t i=0;i<mccol.size();i++){
      StrawDigiMC mcdigi = mccol[i];
      auto const& sgsptr = mcdigi.earlyStrawGasStep();
      auto const& sgs = *sgsptr;
      auto const& sp = *sgs.simParticle();
      auto posi = GenVector::Hep3Vec(sgs.startPosition());
      if ((sp.pdgId() == PDGCode::mu_minus || sp.pdgId() == PDGCode::mu_plus) && sp.creationCode() == 56){
        for (size_t j=i+1; j<mccol.size();j++){
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
    if (_mcdir.y() > 0){
      _mcdir *= -1;
    }
    std::tuple <double, double, double, double, double> info;
    info = make_tuple(_mca0, _mca1, _mcb0, _mcb1, _mct0);
    return info;
    }

    void CosmicMCRecoDiff::endJob() {}

    bool CosmicMCRecoDiff::findData(const art::Event& evt){

      _coscol = 0;
      auto stH = evt.getValidHandle<CosmicTrackSeedCollection>(_costag);
      _coscol =stH.product();

      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigistag);
      _mcdigis = mcdH.product();

      return  _coscol !=0 && _mcdigis != 0 ;
    }

}

using mu2e::CosmicMCRecoDiff;
DEFINE_ART_MODULE(CosmicMCRecoDiff)
