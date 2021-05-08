//Author: S Middleton
//Date: Oct 2021
//Purpose: For getting seed cov estimates

#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <cmath>

// Cosmic Tracks:
#include "CosmicReco/inc/CosmicTrackFit.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "CosmicReco/inc/CosmicTrackMCInfo.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"

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
#include "TrackerConditions/inc/StrawResponse.hh"
#include "Mu2eUtilities/inc/BuildLinearFitMatrixSums.hh"

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

// ROOT incldues
#include "TTree.h"

//Geom
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"

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
        fhicl::Table<SimParticleTimeOffset::Config> toff{Name("TimeOffsets"), Comment("Sim particle time offset ")};
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
      SimParticleTimeOffset _toff;

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

    };

    CosmicMCRecoDiff::CosmicMCRecoDiff(const Parameters& conf) :
      art::EDAnalyzer(conf),
      _chtag (conf().chtag()),
      _tctag (conf().tctag()),
      _costag (conf().costag()),
      _mcdigistag (conf().mcdigistag()),
      _toff (conf().toff())
    {
        for (auto const& tag : conf().toff().inputs()) {
          consumes<SimParticleTimeMap>(tag);
      }
    }

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

      //const Tracker *tracker = _alignedTracker_h.getPtr(event.id()).get();
      //StrawResponse const& srep = _strawResponse_h.get(event.id());

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
        std::tuple <double, double, double, double, double, double> KinKalParams = st.KinKalTrackParams();
        _Recod0 = get<0>(KinKalParams);
        _Recoz0 = get<2>(KinKalParams);
        _RecoCosT = get<3>(KinKalParams);
        _RecoPhi0 = get<1>(KinKalParams);
        std::cout<<_Recod0<<" "<<_Recoz0<<" "<<_RecoCosT<<" "<<_RecoPhi0<<std::endl;
        _RecoErrorA0=(st.MinuitParams.deltaA0);
        _RecoErrorA1=(st.MinuitParams.deltaA1);
        _RecoErrorB1=(st.MinuitParams.deltaB1);
        _RecoErrorB0=(st.MinuitParams.deltaB0);
        
        //const StrawDigiMCCollection*& _mcdigis;

	      StrawDigiMC hitP1;
	      StrawDigiMC first = (*_mcdigis)[0];

	      auto const& spmcp0= first.earlyStrawGasStep();
	      XYZVec zpos(0.,0.,0);
        XYZVec  zdir(0.,0.,1.);
	      XYZVec pos0(spmcp0->position().x(), spmcp0->position().y(), spmcp0->position().z());
	      XYZVec dir(spmcp0->momentum().x(), spmcp0->momentum().y(), spmcp0->momentum().z());
        std::tuple <double,double, double, double, double, double> info;
        TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(pos0, dir, zpos, zdir);
        XYZVec POCA = PCA.point1()-PCA.point2();
        double DOCA = PCA.dca();
        double amsign = copysign(1.0, -(zdir.Cross(POCA)).Dot(dir));

        _Trued0 = amsign*DOCA;
        _Truez0 = dir.Z();
        _TrueCosT = PCA.point1().Z();
        _TruePhi0 = dir.Phi();
        
        _Diffd0 = _Recod0 - _Trued0;
        _Diffz0 = _Recoz0 - _Truez0;
        _DiffCosT =  _RecoCosT - _TrueCosT;
        _DiffPhi0 =  _RecoPhi0 - _TruePhi0;


        _cosmic_tree->Fill();
        
        
      
      }

    }

    void CosmicMCRecoDiff::endJob() {}

    bool CosmicMCRecoDiff::findData(const art::Event& evt){
   
      _coscol = 0;
      auto stH = evt.getValidHandle<CosmicTrackSeedCollection>(_costag);
      _coscol =stH.product();
      
      auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigistag);
      _mcdigis = mcdH.product();
      _toff.updateMap(evt);

      return  _coscol !=0 && _mcdigis != 0 ;
    }

}

using mu2e::CosmicMCRecoDiff;
DEFINE_ART_MODULE(CosmicMCRecoDiff);
