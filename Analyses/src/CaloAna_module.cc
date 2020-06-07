// Author : S. middleton
// Date: March 2020
// Purpose: IPA Analysis
#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/DetectorSystem.hh"

#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrackClusterMatch.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"

// BaBar Kalman filter includes
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/BaBar/BaBar.hh"

// ROOT incldues
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TH2D.h"
#include "TF1.h"

#include "Rtypes.h"
#include "TApplication.h"
#include "TArc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TNtuple.h"

#include "TStyle.h"
#include "TText.h"
#include "TRotMatrix.h"
#include "TColor.h"
#include "TLorentzVector.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <algorithm> 
using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::keV;

#define me 0.511 //MeV

namespace mu2e {

  class CaloAna : public art::EDAnalyzer {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("diag level"),0};
        fhicl::Atom<int> mcdiag{Name("mcdiag"),Comment("mc diag level"),0};
        fhicl::Atom<art::InputTag> kalrepTag{Name("KalRepPtrCollection"),Comment("outcome of Kalman filter (for tracker momentum info)")};
        fhicl::Atom<art::InputTag> tcmatchTag{Name("TrackClusterMatchCollection"), Comment("track calo match"), "TrackCaloMatching"};

      };
    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit CaloAna(const Parameters& conf);
    virtual ~CaloAna() {}


    virtual void beginJob();
    virtual void endJob();
    virtual void analyze(const art::Event& e) override;

  private:
		std::ofstream Newfile, CrystalMap, DescisionFile;
		Config _conf;
		int _diagLevel;
		int _mcdiag;
    art::InputTag _kalrepTag;
    art::InputTag _tcmatchTag;
    const KalRepPtrCollection* _kalrepcol;
//    const TrackClusterMatchCollection* _tcmatchcol;

    TTree* _Ntup;

    Int_t   _nEvents = 0;
    Int_t _evt, _run, _nTracks, _nMatches, _nTrackMatched, _nHits;
    Float_t _TrackT0, _TrackT0Err, _TrackMom, _MaxEoP, _TrackBackTime ,
    _TrackBackOmega ,_TrackBackD0 , _TrackBackZ0, _TrackBackPhi0,
    _TrackBackTanDip, _TrackChi2, _TrackChi2DOF, _TrackCosTheta, _TrackEnergy,
    _TrackEoP,  _TrackMomErr;

    Float_t _matchChi2, _matchEDep, _matchPosXCl, _matchPosYCl, _matchPosZCl,
    _matchPathLen, _matchR, _matchDt, _matchPosXtrk, _matchPosYtrk,
    _matchPosZtrk,_matchTtrk, _matchClusterSize, _CaloEoP, _matchcrySum, _matchErrEdep, _matchSecondMoment, _matchAngle;

    Int_t   _cryId[16384], _crySectionId[16384], _crySimIdx[16384], _crySimLen[16384];
    Float_t _cryTime[16384], _cryEdep[16384],_cryDose[16384], _cryPosX[16384],
    _cryPosY[16384], _cryPosZ[16384], _cryLeak[16384], _cryTotE[16384],
    _cryTotSum[16384], _cryTotEErr[16384], _cryRadius[16384],  _cryMaxR[16384],
    _cryEdepErr[16384];

    Float_t _EoP_diff, _E_diff;

    bool findData(const art::Event& evt);

	};

  CaloAna::CaloAna(const Parameters& conf):
    art::EDAnalyzer(conf),
    _diagLevel(conf().diagLevel()),
    _mcdiag(conf().mcdiag()),
    _kalrepTag(conf().kalrepTag()),
    _tcmatchTag(conf().tcmatchTag())
    {}

  void CaloAna::beginJob(){
    //std::cout<<"[In BeginJob()] Beginning ..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;
    _Ntup  = tfs->make<TTree>("CaloMatchAna", "CaloMatchAna");
    _Ntup->Branch("evt",          	&_evt ,        "evt/I");
    _Ntup->Branch("run",          	&_run ,        "run/I");

    _Ntup->Branch("nTracks",	&_nTracks,		"nTracks/I");
    _Ntup->Branch("TrackT0", 	&_TrackT0, 		"TrackT0/F");
    _Ntup->Branch("TrackT0Err", 	&_TrackT0Err,		"TrackT0Err/F");
    _Ntup->Branch("TrackBackTime", 	&_TrackBackTime,   	"TrackBackTime/F");
    _Ntup->Branch("TrackBackOmega", &_TrackBackOmega,	"TrackBackOmega/F");
    _Ntup->Branch("TrackBackD0",  	&_TrackBackD0, 		"TrackBackD0/F");
    _Ntup->Branch("TrackBackZ0", 	&_TrackBackZ0, 		"TrackBackZ0/F");
    _Ntup->Branch("TrackBackPhi0",	&_TrackBackPhi0,	"TrackBackPhi0/F");
    _Ntup->Branch("TrackBackTanDip",&_TrackBackTanDip,	"TrackBackTanDip/F");
    _Ntup->Branch("TrackChi2",	&_TrackChi2,		"TrackChi2/F");
    _Ntup->Branch("TrackChi2DOF", 	&_TrackChi2DOF, 	"TrackCho2DOF/F");
    _Ntup->Branch("TrackCosTheta",  &_TrackCosTheta,	"TrackCosTheta/F");
    _Ntup->Branch("TrackEnergy",  	&_TrackEnergy,		"TrackEnergy/F");
    _Ntup->Branch("TrackMom",	&_TrackMom, 		"TrackMom/F");
    _Ntup->Branch("TrackEoP",		&_TrackEoP,			"TrackEoP/F");
    _Ntup->Branch("TrackMomErr",		&_TrackMomErr,			"TrackMomErr/F");

    _Ntup->Branch("nMatches",	&_nMatches,		"nMatches/I");
    _Ntup->Branch("matchPosXtrk",	&_matchPosXtrk,		"matchPosXtrk/F");
    _Ntup->Branch("matchPosYtrk",	&_matchPosYtrk,		"matchPosYtrk/F");
    _Ntup->Branch("matchPosZtrk",	&_matchPosZtrk,		"matchPosZtrk/F");
    _Ntup->Branch("matchTtrk",	&_matchTtrk,		"matchTtrk/F");
    _Ntup->Branch("matchChi2", 	&_matchChi2, 		"matchChi/F");
    _Ntup->Branch("matchEDep", 	&_matchEDep, 		"matchEDep/F");
    _Ntup->Branch("matchPosXcl", 	&_matchPosXCl, 		"matchPosXCl/F");
    _Ntup->Branch("matchPosYCl", 	&_matchPosYCl, 		"matchPosYCl/F");
    _Ntup->Branch("matchPosZCl", 	&_matchPosZCl, 		"matchPosZCl/F");
    _Ntup->Branch("matchPathLen", 	&_matchPathLen,	        "matchPathLen/F");
    _Ntup->Branch("matchR",		&_matchR, 		"matchR/F");
    _Ntup->Branch("matchDt",	&_matchDt,		"matchDt/F");
    _Ntup->Branch("matchClusterSize",	&_matchClusterSize,		"matchClusterSize/F");
    _Ntup->Branch("CaloEoP", 	&_CaloEoP, 		"CaloEoP/F");
    _Ntup->Branch("matchcrySum", 	&_matchcrySum, 		"_matchcrySum/F");
    _Ntup->Branch("matchErrEdep", 	&_matchErrEdep, 		"_matchErrEdep/F");
    _Ntup->Branch("matchAngle", 	&_matchAngle, 		"_matchAngle/F");
    _Ntup->Branch("matchSecondMoment", 	&_matchSecondMoment, 		"_matchSecondMoment/F");

    _Ntup->Branch("EoP_diff",	&_EoP_diff,		"_EoP_diff/F");
    _Ntup->Branch("E_diff",	&_E_diff,		"_E_diff/F");

    _Ntup->Branch("nCry",         	&_nHits ,       "nCry/I");
    _Ntup->Branch("cryId",        	&_cryId ,       "cryId[nCry]/I");
    _Ntup->Branch("crySectionId", 	&_crySectionId, "crySectionId[nCry]/I");
    _Ntup->Branch("cryPosX",      	&_cryPosX ,     "cryPosX[nCry]/F");
    _Ntup->Branch("cryPosY",      	&_cryPosY ,     "cryPosY[nCry]/F");
    _Ntup->Branch("cryPosZ",      	&_cryPosZ ,     "cryPosZ[nCry]/F");
    _Ntup->Branch("cryEdep",      	&_cryEdep ,     "cryEdep[nCry]/F");
    _Ntup->Branch("cryEdepErr",      	&_cryEdepErr ,     "cryEdepErr[nCry]/F");
    _Ntup->Branch("cryTime",      	&_cryTime ,     "cryTime[nCry]/F");
    _Ntup->Branch("cryDose",      	&_cryDose ,     "cryDose[nCry]/F");
    _Ntup->Branch("cryRadius",	&_cryRadius,	"cryRadius[nCry]/F");
    Newfile.open("Combined.csv.root");
    DescisionFile.open("List.csv.root");
  }


  void CaloAna::analyze(const art::Event& event) {
    _evt = event.id().event();
    _run = event.run();

    if(!findData(event))
    throw cet::exception("RECO")<<"No data in  event"<< endl;
   
    art::ServiceHandle<mu2e::GeometryService>   geom;
    /* const mu2e::Calorimeter* bc(nullptr);
    if (geom->hasElement<mu2e::DiskCalorimeter>() ) {
      mu2e::GeomHandle<mu2e::DiskCalorimeter> h;
      const mu2e::Calorimeter*bc = (const mu2e::Calorimeter*) h.get();
    }*/
    mu2e::GeomHandle<mu2e::DetectorSystem>      ds;
    mu2e::GeomHandle<mu2e::VirtualDetector>     vdet;
    _nTracks = 0;
    _nMatches =0;

    Hep3Vector vd_tt_back = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_Back));
    double     Z      = vd_tt_back.z();

    for(unsigned int i=0;i<_kalrepcol->size();i++){
     // int iv = 0;
      art::Ptr<KalRep> const& ptr = _kalrepcol->at(i);
      const KalRep* TrackKrep = ptr.get();
   //   const CaloCluster* ClosestCluster(nullptr);
  //    double best_chi2_match(1.e6); //high number to start
      double  ds(10.), s0, s1, s2, z0, z1, z2, dzds, sz, sz1, z01;
      const TrkHitVector* hots = &TrackKrep->hitVector();
      int nh = hots->size();
      const TrkHit *first(nullptr), *last(nullptr);

      for (int ih=0; ih<nh; ++ih) {
        const TrkHit* hit = hots->at(ih);
        if (hit  != nullptr) {
          if (first == nullptr) first = hit;
            last = hit;
          }
        }

      s1 = first->fltLen();
      s2 = last ->fltLen();

      z1     = TrackKrep->position(s1).z();
      z2     = TrackKrep->position(s2).z();

      dzds   = (z2-z1)/(s2-s1);

      if (fabs(Z-z1) > fabs(Z-z2)) {
        z0 = z2;
        s0 = s2;
      }
      else {
        z0 = z1;
        s0 = s1;
      }

    sz    = s0+(Z-z0)/dzds;

    z0     = TrackKrep->position(sz).z();     // z0 has to be close to Z(TT_FrontPA)
    z01    = TrackKrep->position(sz+ds).z();

    dzds   = (z01-z0)/ds;
    sz1    = sz+(Z-z0)/dzds;	          // should be good enough

   // double EndMom= TrackKrep->momentum(sz1).mag();//TODO
    //=========================== Add Tracks =============================//
    _TrackT0 = TrackKrep->t0().t0();
    _TrackT0Err = TrackKrep->t0().t0Err();
    _TrackMom = TrackKrep->momentum(sz1).mag();
    _TrackBackTime =   TrackKrep->arrivalTime(sz1);
    HelixParams helx  = TrackKrep->helix(sz1);
    _TrackBackOmega       = helx.omega();
    _TrackBackD0       = helx.d0();
    _TrackBackZ0       = helx.z0();
    _TrackBackPhi0     = helx.phi0();
    _TrackBackTanDip   = helx.tanDip();
    BbrVectorErr TrackMomErrVec    = TrackKrep->momentumErr(sz1);
    _TrackMomErr = TrackMomErrVec.mag();
    _TrackEnergy = sqrt(TrackKrep->momentum(sz1).mag()*TrackKrep->momentum(sz1).mag() + me*me);
    double entlen         = std::min(s1,s2);
    CLHEP::Hep3Vector fitmom = TrackKrep->momentum(entlen);
    TLorentzVector  Momentum(fitmom.x(),fitmom.y(),fitmom.z(),0.511);
    _TrackCosTheta = Momentum.CosTheta();
    _TrackChi2 =TrackKrep->chisq();
    _TrackChi2DOF= TrackKrep->chisq()/TrackKrep->nActive();
    _nTracks ++;
    //============================= Add Matches ===============================//
   /* if(_tcmatchcol->size() ==0) continue;
    for(unsigned int c=0;c<_tcmatchcol->size();c++){

      TrackClusterMatch const& tcm = (*_tcmatchcol)[c];
      const TrkCaloIntersect* extrk = tcm.textrapol();
      const KalRep* Krep  = extrk->trk().get();
      if (Krep == TrackKrep) {
          const mu2e::CaloCluster* cl = tcm.caloCluster();

          _matchcrySum = 0;
          
          iv   = cl->diskId();
          CLHEP::Hep3Vector x1   = bc->geomUtil().mu2eToDisk(iv,cl->cog3Vector());

          if ((ClosestCluster == nullptr) || (tcm.chi2() < best_chi2_match )) {
            ClosestCluster = cl;
            best_chi2_match    = tcm.chi2();
          }
          _matchClusterSize = cl->size();
          _matchEDep = cl->energyDep();
          _matchErrEdep = cl->energyDepErr();
          _matchSecondMoment = cl->secondMoment();
          _matchAngle = cl->angle();
          _matchPosXtrk = tcm.xtrk();
          _matchPosYtrk = tcm.ytrk();
          _matchPosZtrk = tcm.ztrk();
          _matchTtrk	 = tcm.ttrk();
          _matchChi2 = tcm.chi2();
          _matchPosXCl =x1.x();
          _matchPosYCl =x1.y();
          _matchPosZCl = x1.z();
          _matchPathLen = tcm.ds();
          _matchDt = tcm.dt();
          _matchR = sqrt(x1.x()*x1.x() + x1.y()*x1.y());

          art::ServiceHandle<GeometryService> geom;
      	  if( ! geom->hasElement<Calorimeter>() ) return;
      	  Calorimeter const & cal = *(GeomHandle<Calorimeter>());
          _nHits = cl->caloCrystalHitsPtrVector().size();

          
          unsigned int i =0;
          double Emax = 0;
          double Emin =1000;
          double Xmax = -100000; double Xmin=100000; double Ymin=100000; double Ymax=-100000; 
		      for(unsigned int ic=0 ; ic< cl->caloCrystalHitsPtrVector().size();ic++){
            art::Ptr<CaloCrystalHit>  hit=cl->caloCrystalHitsPtrVector()[ic] ;
            if( hit->energyDep() > Emax) Emax = hit->energyDep();
            int diskId                     = cal.crystal(hit->id()).diskId();
            CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit->id()).position());
            if(crystalPos.x()>Xmax) Xmax = crystalPos.x();
            if(crystalPos.x()<Xmin) Xmin = crystalPos.x();
            if(crystalPos.y()>Ymax) Ymax = crystalPos.y();
            if(crystalPos.y()<Ymin) Ymin = crystalPos.y();
            if(hit->energyDep()<Emin) Emin = hit->energyDep();
            i++;
          }
          double EnergySpread = Emax - Emin;
          double ShowerArea = sqrt((Xmax-Xmin)*(Xmax-Xmin)+(Ymax-Ymin)*(Ymax-Ymin));
          double EnergyDensity = EnergySpread/ShowerArea;
          double EnergyPerCrystal = ClosestCluster->energyDep()/_matchClusterSize;
          if(_matchClusterSize>2){
            DescisionFile<<_evt<<","<<_run<<","<<abs(ClosestCluster->energyDep()-_TrackEnergy)<<","<<_TrackEnergy<<","<<_matchPathLen<<","<<_matchAngle<<","
  <<_matchSecondMoment<<","<<_matchClusterSize<<","<<_TrackBackD0<<","<<_TrackBackPhi0     <<","<<_TrackChi2DOF<<","<<Emax<<","<<EnergySpread<<","<<ShowerArea<<","
  <<EnergyDensity<<","<<EnergyPerCrystal<<endl;
        }
         
          _nMatches++;
        
          //================ Add in Crystal Hits ==========================//
          for(unsigned int ic=0 ; ic< cl->caloCrystalHitsPtrVector().size();ic++){
            art::Ptr<CaloCrystalHit>  cry=cl->caloCrystalHitsPtrVector()[ic] ;
            _cryId[ic] = cry->id();
            int diskId                     = cal.crystal(cry->id()).diskId();
            CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(cry->id()).position());
            _cryTime[ic]       	= cry->time();
            _cryEdep[ic]       	= cry->energyDep();
            _cryEdepErr[ic]     = cry->energyDepErr();
            _cryTotE[ic]  		  = cry->energyDepTot();
            _cryTotEErr[ic]  		= cry->energyDepTotErr();
            _cryPosX[ic]      	= crystalPos.x();
            _cryPosY[ic]       	= crystalPos.y();
            _cryPosZ[ic]       	= crystalPos.z();
            _cryRadius[ic]  	 	= sqrt(crystalPos.x()*crystalPos.x() + crystalPos.y()*crystalPos.y());
            _matchcrySum += cry->energyDep(); 
          if(_matchClusterSize>2){
            Newfile<<_evt<<","<<_run<<","<<_nHits<<","<<cry->id()<<","
            <<cry->energyDep()<<","<<cry->energyDepErr()<<","<<ClosestCluster->energyDep()<<","<<_TrackEnergy<<","<<_TrackMom<<","<<_TrackMomErr<<std::endl;
}
          }
        

        }
      }
      if(EndMom!=0 and ClosestCluster != nullptr) {
        _CaloEoP = ClosestCluster->energyDep()/EndMom;
        _TrackEoP = _TrackEnergy/EndMom;
        _EoP_diff = _CaloEoP - _TrackEoP;
        _E_diff = ClosestCluster->energyDep() - _TrackEnergy;
      }
      
      _nTrackMatched++;
*/
    }
    _Ntup->Fill();
    _nEvents++;
}


       

bool CaloAna::findData(const art::Event& evt){
 // _tcmatchcol=0;
  _kalrepcol = 0;
  auto kalrep = evt.getValidHandle<KalRepPtrCollection>(_kalrepTag);
  _kalrepcol =kalrep.product();
//  auto tcmatch = evt.getValidHandle<TrackClusterMatchCollection>(_tcmatchTag);
//  _tcmatchcol = tcmatch.product();
  return _kalrepcol!=0;// and _tcmatchcol!=0;
}

void CaloAna::endJob(){}

}

DEFINE_ART_MODULE(mu2e::CaloAna);


