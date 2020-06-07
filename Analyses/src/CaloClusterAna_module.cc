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
//#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

//#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"

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

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <fstream>

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

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  class CaloClusterAna : public art::EDAnalyzer {
     public:
	struct Config {
		     using Name=fhicl::Name;
		     using Comment=fhicl::Comment;
	 	     fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("diag level"),0};
		     fhicl::Atom<int> mcdiag{Name("mcdiag"),Comment("mc diag level"),0};
		     fhicl::Atom<art::InputTag> calocrysTag{Name("CaloCrystalHitCollection"),Comment("cal reco crystal hit info")};
		     fhicl::Atom<art::InputTag> caloclusterTag{Name("CaloClusterCollection"),Comment("cal reco cluster info")};
		     fhicl::Atom<art::InputTag> genTag{Name("GenParticleCollection"), Comment("gen particle info")};
		     
	    	};
		typedef art::EDAnalyzer::Table<Config> Parameters;

		explicit CaloClusterAna(const Parameters& conf);
		virtual ~CaloClusterAna() {}

		virtual void beginJob();
		virtual void endJob();
		virtual void analyze(const art::Event& e) override;

     	private:
		std::ofstream outputfile;
		Config _conf;
		int _diagLevel;
		int _mcdiag;

		art::InputTag _calocrysTag;
		art::InputTag _caloclusterTag;
		art::InputTag _genTag;

		const CaloCrystalHitCollection*  _calcryhitcol;	
		const CaloClusterCollection* _calclustercol;
		const GenParticleCollection *_gencol;

		TTree* _Ntup;

		Int_t   _nEvents = 0;
		Int_t _evt, _run, _nHits, _nClusters, _nTracks, _nMatches, _nTrackMatched, _nGen;
		
		Int_t _genPdgId, _genCrCode;

		Float_t _genmomX,_genmomY, _genmomZ, _genStartX,  _genStartY, _genStartZ, _genStartT;

		Int_t   _cryId[8129], _crySectionId[8129], _crySimIdx[8129], _crySimLen[8129];
		Float_t _cryTime[8129], _cryEdep[8129],_cryDose[8129], _cryPosX[8129], _cryPosY[8129], _cryPosZ[8129], _cryLeak[8129], _cryTotE[8129], _cryTotSum[8129], _cryTotEErr[8129], _cryRadius[8129],  _cryMaxR[8129];

		Int_t   _clusterId, _clusterdiskId;
		Float_t _clustertime, _clustertimeErr, _clusterEdep, _clusterEDepErr, _clusterangle, _clustercog3VectorX, _clustercog3VectorY, _clustercog3VectorZ,_clusterR,_clusterNHits, _clustermaxECrystal, _clusterindexMaxECrystal, _clusterERatio;

		bool findData(const art::Event& evt);

	};

   CaloClusterAna::CaloClusterAna(const Parameters& conf):
		art::EDAnalyzer(conf),
		_diagLevel(conf().diagLevel()),
		_mcdiag(conf().mcdiag()),
		_calocrysTag(conf().calocrysTag()),
		_caloclusterTag(conf().caloclusterTag()),
		_genTag(conf().genTag())
	{}

  void CaloClusterAna::beginJob(){
	//std:cout<<"[In BeginJob()] Beginning ..."<<std::endl;
	art::ServiceHandle<art::TFileService> tfs;
	_Ntup  = tfs->make<TTree>("CaloCalibAna", "CaloCalibAna");
	_Ntup->Branch("evt",          	&_evt ,        "evt/I");
	_Ntup->Branch("run",          	&_run ,        "run/I");

	_Ntup->Branch("nGen",         	&_nGen ,        "nGen/I");
	_Ntup->Branch("genId",        	&_genPdgId,     "genId/I");
	_Ntup->Branch("genCrCode",    	&_genCrCode,    "genCrCode/I");
	_Ntup->Branch("genMomX",      	&_genmomX,      "genMomX/F");
	_Ntup->Branch("genMomY",      	&_genmomY,      "genMomY/F");
	_Ntup->Branch("genMomZ",      	&_genmomZ,      "genMomZ/F");
	_Ntup->Branch("genStartX",    	&_genStartX,    "genStartX/F");
	_Ntup->Branch("genStartY",    	&_genStartY,    "genStartY/F");
	_Ntup->Branch("genStartZ",    	&_genStartZ,    "genStartZ/F");
	_Ntup->Branch("genStartT",    	&_genStartT,    "genStartT/F");

	_Ntup->Branch("nCry",         	&_nHits ,       "nCry/I");
	_Ntup->Branch("cryId",        	&_cryId ,       "cryId[nCry]/I");
	_Ntup->Branch("crySectionId", 	&_crySectionId, "crySectionId[nCry]/I");
	_Ntup->Branch("cryPosX",      	&_cryPosX ,     "cryPosX[nCry]/F");
	_Ntup->Branch("cryPosY",      	&_cryPosY ,     "cryPosY[nCry]/F");
	_Ntup->Branch("cryPosZ",      	&_cryPosZ ,     "cryPosZ[nCry]/F");
	_Ntup->Branch("cryEdep",      	&_cryEdep ,     "cryEdep[nCry]/F");
	_Ntup->Branch("cryTime",      	&_cryTime ,     "cryTime[nCry]/F");
	_Ntup->Branch("cryDose",      	&_cryDose ,     "cryDose[nCry]/F");
	_Ntup->Branch("cryRadius",	&_cryRadius,	"cryRadius[nCry]/F");
	
	_Ntup->Branch("nClu",         	&_nClusters ,       	"nClu/I");
	_Ntup->Branch("clustertime", 	&_clustertime,		"clustertime/F");
	_Ntup->Branch("clustertimeErr", &_clustertimeErr,	"clustertimeErr/F");
	_Ntup->Branch("clusterEdep", 	&_clusterEdep,		"clusterEdep/F");
	_Ntup->Branch("clusterEDepErr", &_clusterEDepErr, 	"clusterEdepErr/F");
	_Ntup->Branch("clusterangle", 	&_clusterangle, 	"clusterangle/F");
	_Ntup->Branch("clusterPosX",	&_clustercog3VectorX, 	"clustercog3VectorX/F");
	_Ntup->Branch("clusterPosY",	&_clustercog3VectorY, 	"clustercog3VectorY/F");
	_Ntup->Branch("clusterPosZ",	&_clustercog3VectorZ, 	"clustercog3VectorZ/F");
	_Ntup->Branch("clusterR", 	&_clusterR,		"clusterR/F");
	_Ntup->Branch("clusterNHits", 	&_clusterNHits,		"clusterNHits/F");
	outputfile.open("IPAAnaCaloClusters.csv");
        outputfile<<"event,run,cluster_size"<<std::endl;
  }


  void CaloClusterAna::analyze(const art::Event& event) {
	//std:cout<<"[In Analyze()] Beginning ..."<<std::endl;
	_evt = event.id().event();
	_run = event.run();

	if(!findData(event)) 
		throw cet::exception("RECO")<<"No data in  event"<< endl; 

	//std:cout<<"[In Analyze()] Found Data ..."<<std::endl;

      	art::ServiceHandle<GeometryService> geom;
      	if( ! geom->hasElement<Calorimeter>() ) return;
      	Calorimeter const & cal = *(GeomHandle<Calorimeter>());

//================== GenPartile Info======================//
	//std:cout<<"[In Analyze()] Getting GenInfo ..."<<std::endl;
	_nGen = _gencol->size();
	for (unsigned int i=0; i <_gencol->size(); ++i)
	{
		GenParticle const& gen = (*_gencol)[i];
		_genPdgId   = gen.pdgId();
		_genCrCode  = gen.generatorId().id();
		_genmomX    = gen.momentum().vect().x();
		_genmomY    = gen.momentum().vect().y();
		_genmomZ    = gen.momentum().vect().z();
		_genStartX  = gen.position().x()+ 3904;
		_genStartY  = gen.position().y();
		_genStartZ  = gen.position().z();
		_genStartT  = gen.time();
	} 


//=======================Get Cluster Info ==============//
	//std:cout<<"[In Analyze()] Getting Cluster Info..."<<std::endl;
        for (unsigned int tclu=0; tclu<_calclustercol->size();++tclu){
		CaloCluster const& cluster = (*_calclustercol)[tclu];
		const CaloCluster::CaloCrystalHitPtrVector caloClusterHits = 		
						 cluster.caloCrystalHitsPtrVector();

		CLHEP::Hep3Vector crystalPos   	= cal.geomUtil().mu2eToDiskFF(cluster.diskId(), 	 						cluster.cog3Vector()); 
		_clusterNHits 			=  caloClusterHits.size();
		_clusterdiskId  		= cluster.diskId();
		_clustertime      		= cluster.time();
		_clusterEdep      		= cluster.energyDep(); 
		_clusterEDepErr     		= cluster.energyDepErr();
		_clusterangle      		= cluster.angle();
		_clustercog3VectorX      	= cluster.cog3Vector().x();
		_clustercog3VectorY      	= cluster.cog3Vector().y();
		_clustercog3VectorZ      	= cluster.cog3Vector().z();
		_clusterR  			= sqrt(cluster.cog3Vector().x()* cluster.cog3Vector().x()
								+ cluster.cog3Vector().y()*cluster.cog3Vector().y());
		outputfile<<_evt<<","<<_run<<","<<_calcryhitcol->size()<<","<<_clusterEdep <<std::endl;
        	_nClusters++;
	}
		
//=====================Crystal Hits Info =======================//
	//std:cout<<"[In Analyze()] Getting Crystal Info..."<<std::endl;
	_nHits = _calcryhitcol->size();
  
	for (unsigned int ic=0; ic<_calcryhitcol->size();++ic) 
	{	   
		   CaloCrystalHit const& hit      = (*_calcryhitcol)[ic];
		   int diskId                     = cal.crystal(hit.id()).diskId();
		   CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit.id()).position());  
		   _cryId[ic] 	 	= hit.id();
		   _cryTime[ic]      	= hit.time();
		   _cryEdep[ic]      	= hit.energyDep();
		   _cryTotE[ic] 	= hit.energyDepTot();
		   _cryTotEErr[ic] 	= hit.energyDepTotErr();
		   _cryPosX[ic]     	= crystalPos.x();
		   _cryPosY[ic]      	= crystalPos.y();
		   _cryPosZ[ic]      	= crystalPos.z();
		   _cryRadius[ic] 	= sqrt(crystalPos.x()*crystalPos.x() +  
						    crystalPos.y()*crystalPos.y());
	}
	_Ntup->Fill();
	_nEvents++;
}


bool CaloClusterAna::findData(const art::Event& evt){

	_calcryhitcol =0;
	_calclustercol=0;
	_gencol=0;
	
	auto genpart = evt.getValidHandle<GenParticleCollection>(_genTag);
	_gencol = genpart.product();
	auto cryhit = evt.getValidHandle<CaloCrystalHitCollection>(_calocrysTag);
	_calcryhitcol =cryhit.product();
	auto cluster= evt.getValidHandle<CaloClusterCollection>(_caloclusterTag);
	_calclustercol =cluster.product();
	
        
	return  _gencol!=0 and _calcryhitcol!=0 && _calclustercol !=0;
       }

 void CaloClusterAna::endJob(){} 

}

DEFINE_ART_MODULE(mu2e::CaloClusterAna);
