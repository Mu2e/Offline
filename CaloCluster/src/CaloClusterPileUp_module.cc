//
// Visualization of pile up on the calorimeter clusters
//
// $Id: CaloClusterPileUp_module.cc,v 1.6 2013/03/15 15:52:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:03 $
//
// Original author G. Pezzullo & G. Tassielli
//


// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"

#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"


// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/VisibleGenElTrack.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"

//calorimeter packages
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
//#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "CaloCluster/inc/CaloClusterUtilities.hh"
//#include "CaloCluster/inc/CaloClusterer.hh"

//#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
//#include "CaloCluster/inc/CaloClusterTools.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib/exception.h"
#include "BaBar/BaBar/include/Constants.hh"

//Root includes
#include "TFile.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TTree.h"

//c++ includes
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "TMath.h"



using namespace std;

namespace mu2e {

//struct electronData{
//        double _cluEnergy;
//        double _cluTime;
//        double _impTime;
//        double _impEnergy;
//        CLHEP::Hep3Vector _cluCog;
//        CLHEP::Hep3Vector _impPos;
//        CLHEP::Hep3Vector _impPosCryFrame;
//        unsigned int _row;
//        unsigned int _column;
//        unsigned int _Nseed;
//        CLHEP::Hep3Vector _cryOrigin;
//        CLHEP::Hep3Vector _impMom3Vec;
//        double _cryTime;
//
//        bool operator<( const electronData other) const{
//                return ( _impTime< other._impTime);
//        }
//        electronData & operator=(const electronData& other) {
//                _cluEnergy = other._cluEnergy;
//                _cluTime   = other._cluTime;
//                _impTime = other._impTime;
//                _impEnergy = other._impEnergy;
//                _cluCog    = other._cluCog;
//                _impPos = other._impPos;
//                _impPosCryFrame = other._impPosCryFrame;
//                _row       = other._row;
//                _column    = other._column;
//                _Nseed     = other._Nseed;
//                _cryOrigin = other._cryOrigin;
//                _impMom3Vec =other._impMom3Vec;
//                _cryTime    = other._cryTime;
//                return *this;
//        }
//        electronData():
//                _impTime(1e10),
//                _impEnergy(0.0),
//                _Nseed(0){
//        }
//};

struct cryInfo {
        CLHEP::Hep3Vector position;
        CLHEP::Hep3Vector cryPosition;
        CLHEP::Hep3Vector momentum;
        double energy;
        double  time;
        double  totEnergyDep;
        int pdgId;
        int isGen;
        int isSignal;
        int isDIO;
        int isPionCapture;
        int isPionCaptureComb;
        int isMuonCapture;
        int isMuonDecayInFlight;
        int isEPlusfromStoppedPi;
        int vane;
        int qualityCuts;
};
//the key is the trackId
typedef std::map<int , cryInfo> CryMap;


//the key is the the vane
//typedef std::map<unsigned int,std::map<unsigned int, electronData > > ElectronMap;

static int ncalls(0);

class CaloClusterPileUp : public art::EDAnalyzer {
public:
        explicit CaloClusterPileUp(fhicl::ParameterSet const& pset):
        _diagLevel(pset.get<int>("diagLevel",1)),
        _qualityCuts(pset.get<int>("qualityCuts",2)),
        _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
        _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
        _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
        _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
        _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel", "makeCaloCluster")),
        _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
        _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "time")),
        _producerName("Algo"+TOUpper(_caloClusterAlgorithm)+"SeededBy"+TOUpper(_caloClusterSeeding)),
        _elextractModuleLabel(pset.get<std::string>("elextractModuleLabel", "extractElData")),
        _extractElectronsData(pset.get<string>("elextractModuleLabel")),
        _rowToCanc(pset.get<int>("rowToCancel",-1)),
        _columnToCanc(pset.get<int>("columnToCancel",-1)),
        _nAnalyzed(0),
        _nBadG4Status(0),
        _nOverflow(0),
        _nKilled(0),
        _totalcputime(0),
        _totalrealtime(0),
        _Ntup(0),
        _EnergyClusterCut(pset.get<double>("energyClusterCut",60.)),//MeV
        _application(nullptr),
        _directory(0)
        {
        }
        virtual ~CaloClusterPileUp() {
                // if (_fakeCanvas)        delete _fakeCanvas;
        }
        virtual void beginJob();
        virtual void endJob();

        void analyze(art::Event const& e );
        // void produce(art::Event /*const*/& e );

private:

        //void doTracker(art::Event const& evt, bool skip);

        void doCalorimeter(art::Event const& evt, bool skip);

        // Diagnostic level
        int _diagLevel;

        int _qualityCuts;

        // Label of the module that made the hits.
        //std::string _makerModuleLabel;

        // Label of the generator.
        std::string _generatorModuleLabel;

        // Label of the G4 module
        std::string _g4ModuleLabel;

        // Label of the calo readout hits maker
        std::string _caloReadoutModuleLabel;

        // Label of the calo crystal hists maker
        std::string _caloCrystalModuleLabel;

        // Label of the calo clusters  maker
        std::string _caloClusterModuleLabel;

        string _caloClusterAlgorithm;
        string _caloClusterSeeding;
        string _producerName;
        string _elextractModuleLabel;
        string _extractElectronsData;

        int _rowToCanc, _columnToCanc;

        //number of analyzed events
        int _nAnalyzed;

        int _nBadG4Status, _nOverflow, _nKilled;
        float _totalcputime, _totalrealtime;

        TTree* _Ntup;

        double _EnergyClusterCut;

        bool _skipEvent;

        // The job needs exactly one instance of TApplication.  See note 1.
        unique_ptr<TApplication> _application;

        // Save directory from beginJob so that we can go there in endJob. See note 3.
        TDirectory* _directory;

        Int_t _clNo,
        _nCryCl,
        _clSeeds,
        _clVane;

        Float_t _evt,
        _clE,
        _clT,
        _clCOGu,
        _clCOGv,
        _clCOGw,
        _clCryEnergyMaxRow,
        _clCryEnergyMaxColumn;

        Int_t _clSeedTrackId[10000],
        _clSeedPdgId[10000],
        _clSeedQC[10000],
        _clSeedIsGen[10000],
        _clSeedIsSignal[10000],
        _clSeedIsDIO[10000],
        _clSeedIsPionCapture[10000],
        _clSeedIsPionCaptureComb[10000],
        _clSeedIsMuonCapture[10000],
        _clSeedIsMuonDecayInFlight[10000],
        _clSeedIsEPlusfromStoppedPi[10000];

        Float_t _clSeedTime[10000],
        _clSeedTotEnergyDep[10000],
        _clSeedPx[10000],
        _clSeedPy[10000],
        _clSeedPz[10000],
        _clSeedPu[10000],
        _clSeedPv[10000],
        _clSeedPw[10000],
        _clSeedCryFramePu[10000],
        _clSeedCryFramePv[10000],
        _clSeedCryFramePw[10000],
        _clSeedVaneFramePu[10000],
        _clSeedVaneFramePv[10000],
        _clSeedVaneFramePw[10000],
        _clSeedPpx[10000],
        _clSeedPpy[10000],
        _clSeedPpz[10000],
        _clSeedE[10000],
        _clSeedPp[10000],
        _clSeedPpu[10000],
        _clSeedPpv[10000],
        _clSeedPpw[10000],
        _clSeedThetaW[10000],
        _clSeedThetaV[10000];




};

bool findTrkId(std::vector<unsigned int> vec, unsigned int t){
        bool res = false;

        unsigned int size = vec.size();
        if(size!=0){
                unsigned int cont = 0;
                while(cont!=size){
                        if(vec[cont] == t) {
                                res = true;
                        }
                        ++cont;

                }
        }
        return res;
}


void CaloClusterPileUp::beginJob( ) {

        cout << "start CaloClusterEnergyResolMap..."<<endl;
}

void CaloClusterPileUp::analyze(art::Event const & evt ) {

        ++_nAnalyzed;
        ++ncalls;

        //CaloClusterer c;

        art::ServiceHandle<GeometryService> geom;
        GeomHandle<VaneCalorimeter> cg;
        if (ncalls == 1) {

                // cout << "This should be done only in the first event" << endl;
                art::ServiceHandle<art::TFileService> tfs;

                _Ntup        = tfs->make<TTree>("ClusterSeedMap", "Cluster seeds info");

                _Ntup->Branch("evt", &_evt , "evt/F");
                _Ntup->Branch("clNo",&_clNo , "clNo/I");
                _Ntup->Branch("clE",&_clE , "clE/F");
                _Ntup->Branch("clT",&_clT , "clT/F");
                _Ntup->Branch("clCOGu",&_clCOGu , "clCOGu/F");
                _Ntup->Branch("clCOGv",&_clCOGv , "clCOGv/F");
                _Ntup->Branch("clCOGw",&_clCOGw , "clCOGw/F");
                _Ntup->Branch("clCryEnergyMaxRow", &_clCryEnergyMaxRow , "clCryEnergyMaxRow/F");
                _Ntup->Branch("clCryEnergyMaxColumn", &_clCryEnergyMaxColumn , "clCryEnergyMaxColumn/F");
                _Ntup->Branch("nCryCl",&_nCryCl , "nCryCl/I");
                _Ntup->Branch("clVane",&_clVane , "clVane/I");
                _Ntup->Branch("clSeeds",&_clSeeds , "clSeeds/I");

                _Ntup->Branch("clSeedTrakId[clSeeds]",_clSeedTrackId , "clSeedTrakId[clSeeds]/I");
                _Ntup->Branch("clSeedPdgId[clSeeds]", _clSeedPdgId , "clSeedPdgId[clSeeds]/I");
                _Ntup->Branch("clSeedQC[clSeeds]", _clSeedQC , "clSeedQC[clSeeds]/I");
                _Ntup->Branch("clSeedIsGen[clSeeds]", _clSeedIsGen , "clSeedIsGen[clSeeds]/I");
                _Ntup->Branch("clSeedIsPionCapture[clSeeds]", _clSeedIsPionCapture , "clSeedIsPionCapture[clSeeds]/I");
                _Ntup->Branch("clSeedIsPionCaptureComb[clSeeds]", _clSeedIsPionCaptureComb , "clSeedIsPionCaptureComb[clSeeds]/I");
                _Ntup->Branch("clSeedIsDIO[clSeeds]", _clSeedIsDIO , "clSeedIsDIO[clSeeds]/I");
                _Ntup->Branch("clSeedIsMuonCapture[clSeeds]", _clSeedIsMuonCapture , "clSeedIsMuonCapture[clSeeds]/I");
                _Ntup->Branch("clSeedIsMuonDecayInFlight[clSeeds]", _clSeedIsMuonDecayInFlight , "clSeedIsMuonDecayInFlight[clSeeds]/I");
                _Ntup->Branch("clSeedIsEPlusfromStoppedPi[clSeeds]", _clSeedIsEPlusfromStoppedPi, "clSeedIsEPlusfromStoppedPi[clSeeds]/I");
                _Ntup->Branch("clSeedIsSignal[clSeeds]", _clSeedIsSignal , "clSeedIsSignal[clSeeds]/I");
                _Ntup->Branch("clSeedTime[clSeeds]",_clSeedTime , "clSeedTime[clSeeds]/F");
                _Ntup->Branch("clSeedTotEnergyDep[clSeeds]",_clSeedTotEnergyDep , "clSeedTotEnergyDep[clSeeds]/F");
                _Ntup->Branch("clSeedPx[clSeeds]",_clSeedPx , "clSeedPx[clSeeds]/F");
                _Ntup->Branch("clSeedPy[clSeeds]",_clSeedPy , "clSeedPy[clSeeds]/F");
                _Ntup->Branch("clSeedPz[clSeeds]",_clSeedPz , "clSeedPz[clSeeds]/F");
                _Ntup->Branch("clSeedPu[clSeeds]",_clSeedPu , "clSeedPu[clSeeds]/F");
                _Ntup->Branch("clSeedPv[clSeeds]",_clSeedPv , "clSeedPv[clSeeds]/F");
                _Ntup->Branch("clSeedPw[clSeeds]",_clSeedPw , "clSeedPw[clSeeds]/F");
                _Ntup->Branch("clSeedCryFramePu[clSeeds]",_clSeedCryFramePu , "clSeedCryFramePu[clSeeds]/F");
                _Ntup->Branch("clSeedCryFramePv[clSeeds]",_clSeedCryFramePv , "clSeedCryFramePv[clSeeds]/F");
                _Ntup->Branch("clSeedCryFramePw[clSeeds]",_clSeedCryFramePw , "clSeedCryFramePw[clSeeds]/F");
                _Ntup->Branch("clSeedVaneFramePu[clSeeds]",_clSeedVaneFramePu , "clSeedVaneFramePu[clSeeds]/F");
                _Ntup->Branch("clSeedVaneFramePv[clSeeds]",_clSeedVaneFramePv , "clSeedVaneFramePv[clSeeds]/F");
                _Ntup->Branch("clSeedVaneFramePw[clSeeds]",_clSeedVaneFramePw , "clSeedVaneFramePw[clSeeds]/F");
                _Ntup->Branch("clSeedPpx[clSeeds]",_clSeedPpx , "clSeedPpx[clSeeds]/F");
                _Ntup->Branch("clSeedPpy[clSeeds]",_clSeedPpy , "clSeedPpy[clSeeds]/F");
                _Ntup->Branch("clSeedPpz[clSeeds]",_clSeedPpz , "clSeedPpz[clSeeds]/F");
                _Ntup->Branch("clSeedE[clSeeds]",_clSeedE , "clSeedE[clSeeds]/F");
                _Ntup->Branch("clSeedPp[clSeeds]",_clSeedPp , "clSeedPp[clSeeds]/F");
                _Ntup->Branch("clSeedPpu[clSeeds]",_clSeedPpu , "clSeedPpu[clSeeds]/F");
                _Ntup->Branch("clSeedPpv[clSeeds]",_clSeedPpv , "clSeedPpv[clSeeds]/F");
                _Ntup->Branch("clSeedPpw[clSeeds]",_clSeedPpw , "clSeedPpw[clSeeds]/F");
                _Ntup->Branch("clSeedThetaW[clSeeds]",_clSeedThetaW , "clSeedThetaW[clSeeds]/F");
                _Ntup->Branch("clSeedThetaV[clSeeds]",_clSeedThetaV , "clSeedThetaV[clSeeds]/F");

        }

        doCalorimeter(evt, _skipEvent);

} // end of analyze


void CaloClusterPileUp::endJob() { }


void CaloClusterPileUp::doCalorimeter(art::Event const& evt, bool skip){

        cout << "Event number : " << evt.id().event()<<
                        ", CaloClusterPileUp: begin" << endl;

        GlobalConstantsHandle<ParticleDataTable> pdt;

        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;
        if(! geom->hasElement<VaneCalorimeter>() ) return;
        GeomHandle<VaneCalorimeter> cg;

        // Get handles to calorimeter collections
        art::Handle<CaloHitCollection> caloHits;
        evt.getByLabel(_caloReadoutModuleLabel, caloHits);

        art::Handle<CaloCrystalHitCollection>  caloCrystalHits;
        evt.getByLabel(_caloCrystalModuleLabel, caloCrystalHits);

        // Get the persistent data about pointers to StepPointMCs
        art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
        evt.getByLabel(_caloReadoutModuleLabel,"CaloHitMCCrystalPtr",mcptrHandle);

        PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();
        if (!( caloHits.isValid())) {
                return;
        }

        if (!caloCrystalHits.isValid()) {
                cout << "NO CaloCrystalHits" << endl;
                return;
        }

        // Get handles to the generated and simulated particles.
        art::Handle<GenParticleCollection> genParticles;
        evt.getByLabel(_generatorModuleLabel, genParticles);

        art::Handle<SimParticleCollection> simParticles;
        evt.getByLabel(_g4ModuleLabel, simParticles);

        art::Handle<CaloClusterCollection> caloClusters;
        evt.getByLabel(_caloClusterModuleLabel,_producerName,caloClusters );

        art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
        evt.getByLabel(_extractElectronsData,genEltrksHandle);

        VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();
        std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it;

        double trkMomCut = 100.0;//MeV
        std::vector<unsigned int> tmpV, tmpVTot, tmpVsignal;
        trkIdVector qualityCollection;


        int NtrkCut =0;
        int NtrkTot = 0;

        if ( _diagLevel < 0 ){
                cout<<"-------------------- 1 -----------------------"<<endl;
                cout<<"start counting & mapping generated electrons..."<<endl;
        }
        //mapping &counting the electrons with quality cuts in the TRK
        for ( genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
                ++NtrkTot;
                VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);
                GenElHitData& hdil = iEltrk.getithLoopHit(0);
                //GenElHitData& ldil = iEltrk.getHit((int)(iEltrk.getNumOfHit() - 1) );

                if(!findTrkId(tmpVTot, iEltrk.getTrkID().asUint() ) ){

                        tmpVTot.push_back( iEltrk.getTrkID().asUint() );

                }

                double cosTheta = iEltrk.getTrkLrntzVec().cosTheta() ;

                double cosPitch = hdil._hitMomentum.cosTheta();

                //double lcosPitch = ldil._hitMomentum.cosTheta();

                bool condition = true;
                condition &= ( iEltrk.getNumOfHit() >= 20 );
                condition &= ( hdil._hitMomentum.mag() >= trkMomCut );
                condition &= ( cosTheta >= -0.5 );
                condition &= ( cosTheta <=  0.5 );
                condition &= ( cosPitch > 0.5 );
                condition &= ( cosPitch < 0.70710678118655 );// 1 / sqrt(2)
                //condition &= ( lcosPitch / cosPitch ) >= (1.0 - iEltrk.getNumOfHit()*8.e-4);
                //condition &= (ldil._hitMomentum.mag() /hdil._hitMomentum.mag() ) >= (1.0 - iEltrk.getNumOfHit()*2.e-4);


                if( condition){
                        NtrkCut++;// # of electrons at the entrance of the TRK which have trkMomCut MeV of kinetic energy and which will light at least 20 straws
                        if(!findTrkId(tmpV, iEltrk.getTrkID().asUint() ) ){
                                qualityCollection.push_back(iEltrk.getTrkID().asUint());
                                tmpV.push_back( iEltrk.getTrkID().asUint() );
                                if ( iEltrk.isConversionEl() ) {
                                        tmpVsignal.push_back( iEltrk.getTrkID().asUint() );
                                }

                        }

                }



        }//end TRK mapping
        if ( _diagLevel < 0 ){
                cout<<"-------------------- 2 -----------------------"<<endl;

                cout<<"end counting & mapping generated electrons..."<<endl;
                cout <<"NtrkCut = " <<NtrkCut<<endl;

        }

        if (NtrkTot==0) return;

        std::vector<unsigned int> trkVec, trkVecTot;
        for(unsigned int f = 0; f != tmpV.size(); ++f){
                trkVec.push_back(tmpV[f]);
        }
        for(unsigned int f = 0; f != tmpVTot.size(); ++f){
                trkVecTot.push_back(tmpVTot[f]);
        }

        if(caloClusters->size()>0 ){
                if ( _diagLevel < 0 ){
                        cout<<"-------------------- 3 -----------------------"<<endl;
                }
                int iVane =-1;
                for(size_t icl=0; icl<caloClusters->size(); ++icl){

                        double eDepClu = 0.;

                        _evt = evt.id().event();

                        CaloCluster const& clu = (*caloClusters).at(icl);
			CaloClusterTools cluTool(clu);

                        _clNo = icl;
                        _clE  = clu.energyDep();

                        _clT  = cluTool.timeFasterCrystal();//time();

                        _clCOGu = clu.cog3Vector().x();
                        _clCOGv = clu.cog3Vector().y();
                        _clCOGw = clu.cog3Vector().z();

                        _nCryCl = clu.size();//clusterSize;

                        eDepClu = clu.energyDep();
                        iVane = clu.vaneId();
                        _clVane = iVane;

                       
                        _clCryEnergyMaxRow = cluTool.cryEnergydepMaxRow(),
                          _clCryEnergyMaxColumn = cluTool.cryEnergydepMaxColumn();
                        CaloCrystalHitPtrVector caloClusterHits = clu.caloCrystalHitsPtrVector();

                        if(eDepClu >= _EnergyClusterCut){

                                if(_diagLevel < 0){
                                        cout <<"eDepClu >energyDepClusterCut"<< endl;
                                }

                                std::map<unsigned int, unsigned int> seedMap, signalMap;
                                CryMap cryMap;
				int totGen = 0;
				
                                for(size_t i=0; i<caloClusterHits.size(); ++i){

                                        if ( _diagLevel < 0 ){
                                                cout<<"-------------------- 3."<<i <<" -----------------------"<<endl;
                                        }
                                        CaloCrystalHit const& hit = *(caloClusterHits.at(i));

                                        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                                        if(ROIds.size()<1 ) continue;

                                        CaloHit const& thehit = *ROIds.at(0);

                                        size_t collectionPosition = ROIds.at(0).key();

                                        PtrStepPointMCVector const & mcptr(hits_mcptr->at(collectionPosition));

                                        size_t nHitsPerCrystal = mcptr.size();
                                        CLHEP::Hep3Vector crystalFrame(0.0, 0.0, 0.0);

                                        for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {

                                                if(_diagLevel < 0){
                                                        cout<<"-------------------- 3."<<i <<"."<< j2<<"-----------------------"<<
                                                                        "mcptr.size() = "<< mcptr.size()<<endl;

                                                }

                                                StepPointMC const& mchit = *mcptr[j2];

                                                if(_diagLevel < 0){
                                                        cout<<"stored mchit..."<<endl;
                                                        mchit.print(cout);
                                                        cout<<endl;
                                                }
                                                art::Ptr<SimParticle> const& simptr = mchit.simParticle();
                                                SimParticle const& sim = *simptr;


                                                if(_diagLevel < 0){
                                                        cout<<"stored sim..."<<endl;
                                                }

                                                SimParticleCollection::key_type trackId(mchit.trackId());

                                                if(_diagLevel < 0){
                                                        cout<<"trackId = "<< trackId.asUint()<<endl;
                                                }

                                                int isConv = 0
                                                                , isDIO = 0
                                                                , isPionCapture = 0
                                                                , isPionCaptureComb = 0
                                                                , isMuonCapture = 0
                                                                , isMuonDecayInFlight = 0
                                                                , isEPlusfromStoppedPi = 0;

                                                if(_diagLevel < 0){
                                                        cout<<"sim info : "<<
                                                                        "sim.nSteps() = "<< sim.nSteps()<<
                                                                        "sim.startMomentum() = "<<sim.startMomentum()<<
                                                                        "sim.startPosition() = " <<sim.startPosition()<<endl;
                                                }

                                                if(_diagLevel < 0){
                                                        cout<<"Conditon yes"<<endl;
                                                }
						

						double mass = 0.0;
						double trackKine = 0.0;
						if (!( sim.madeInG4()) ){
						
                                                const HepPDT::ParticleData& data = pdt->particle(sim.pdgId() ).ref();
                                                mass = data.mass().value();

                                                //add non-linearity effect
                                                trackKine = std::sqrt(mchit.momentum().mag2() + std::pow(mass, 2) ) - mass;
						}else{
						  if(_diagLevel < 0){
						    cout<<"ATTENTION :: !( sim.madeInG4())..."<<endl;
						  }
						}
                                                if(cryMap.find(mchit.trackId().asInt() )==cryMap.end() ){
                                                        if(_diagLevel < 0){
                                                                cout<<"Add new particle to cryMap..."<<endl;
                                                        }
                                                        cryMap[mchit.trackId().asInt()].time         = mchit.time();
                                                        cryMap[mchit.trackId().asInt()].totEnergyDep = mchit.totalEDep();
                                                        cryMap[mchit.trackId().asInt()].pdgId        = sim.pdgId();
                                                        cryMap[mchit.trackId().asInt()].isGen        = sim.fromGenerator();
							if(sim.fromGenerator()==1){
							  ++totGen;
							}
                                                        cryMap[mchit.trackId().asInt()].position     = mchit.position();
                                                        cryMap[mchit.trackId().asInt()].momentum     = mchit.momentum();
                                                        cryMap[mchit.trackId().asInt()].energy       = trackKine;// std::sqrt(std::pow(mchit.momentum().mag(), 2) + std::pow(gen.momentum().m(), 2) ) - gen.momentum().m();
                                                        crystalFrame = cg->toCrystalFrame(thehit.id(), mchit.position());
                                                        cryMap[mchit.trackId().asInt()].cryPosition  = crystalFrame;
                                                        cryMap[mchit.trackId().asInt()].vane         = iVane;
                                                        cryMap[mchit.trackId().asInt()].qualityCuts  = qualityCollection.find(mchit.trackId().asInt());



                                                        if(sim.genParticle().isNonnull()){

                                                                GenParticle const& gen = *sim.genParticle();
                                                                if(_diagLevel < 0){
                                                                        cout<<"gen stored..."<<endl;
                                                                }
                                                                GenId genId = gen.generatorId();
                                                                if(_diagLevel < 0){
                                                                        cout<<"genId stored..."<<endl;
                                                                }

                                                                if(_diagLevel < 0){
                                                                        string name = genId.name();
                                                                        cout<< name.c_str()<<endl;
                                                                }
                                                                if(genId==GenId::conversionGun){
                                                                        isConv = 1;
                                                                }


                                                                if(genId.isDio()){
                                                                        isDIO = 1;
                                                                }

                                                                if(genId==GenId::pionCapture){
                                                                        isPionCapture = 1;
                                                                }

                                                                if(genId==GenId::PiCaptureCombined){
                                                                        isPionCaptureComb = 1;
                                                                }

                                                                if(genId==GenId::muonCapture){
                                                                        isMuonCapture = 1;
                                                                }

                                                                if(genId==GenId::muonCapture){
                                                                        isMuonDecayInFlight = 1;
                                                                }

                                                                if(genId==GenId::muonCapture){
                                                                        isEPlusfromStoppedPi = 1;
                                                                }
                                                        }
                                                        cryMap[mchit.trackId().asInt()].isSignal     = isConv;
                                                        cryMap[mchit.trackId().asInt()].isDIO     = isDIO;

                                                        cryMap[mchit.trackId().asInt()].isPionCapture    = isPionCapture;

                                                        cryMap[mchit.trackId().asInt()].isPionCaptureComb     = isPionCaptureComb;

                                                        cryMap[mchit.trackId().asInt()].isMuonCapture     = isMuonCapture;

                                                        cryMap[mchit.trackId().asInt()].isMuonDecayInFlight    = isMuonDecayInFlight;

                                                        cryMap[mchit.trackId().asInt()].isEPlusfromStoppedPi    = isEPlusfromStoppedPi;
                                                }else if(mchit.time() < cryMap[mchit.trackId().asInt()].time ){
                                                        if(_diagLevel < 0){
                                                                cout<<"Update cryMap..."<<endl;
                                                        }
                                                        cryMap[mchit.trackId().asInt()].time         = mchit.time();
                                                        cryMap[mchit.trackId().asInt()].totEnergyDep = mchit.totalEDep();
                                                        cryMap[mchit.trackId().asInt()].pdgId        = sim.pdgId();
                                                        cryMap[mchit.trackId().asInt()].isGen        = sim.fromGenerator();
                                                        cryMap[mchit.trackId().asInt()].position     = mchit.position();
                                                        cryMap[mchit.trackId().asInt()].momentum     = mchit.momentum();
                                                        cryMap[mchit.trackId().asInt()].energy       = trackKine;//std::sqrt(std::pow(mchit.momentum().mag(), 2) + std::pow(gen.momentum().m(), 2) ) - gen.momentum().m();

                                                        crystalFrame = cg->toCrystalFrame(thehit.id(), mchit.position());

                                                        cryMap[mchit.trackId().asInt()].cryPosition  = crystalFrame;
                                                        cryMap[mchit.trackId().asInt()].vane         = iVane;//cg->getVaneByRO(thehit.id());

                                                        if(sim.genParticle().isNonnull())  {

                                                                GenParticle const& gen = *sim.genParticle();
                                                                if(_diagLevel < 0){
                                                                        cout<<"gen stored..."<<endl;
                                                                }
                                                                GenId genId = gen.generatorId();
                                                                if(_diagLevel < 0){
                                                                        cout<<"genId stored..."<<endl;
                                                                }

                                                                if(_diagLevel < 0){
                                                                        string name = genId.name();
                                                                        cout<< "GenId = "<< name.c_str()<<endl;
                                                                }

                                                                if(genId==GenId::conversionGun){
                                                                        isConv = 1;
                                                                }


                                                                if(genId.isDio()){
                                                                        isDIO = 1;
                                                                }


                                                                if(genId==GenId::pionCapture){
                                                                        isPionCapture = 1;
                                                                }


                                                                if(genId==GenId::PiCaptureCombined){
                                                                        isPionCaptureComb = 1;
                                                                }


                                                                if(genId==GenId::muonCapture){
                                                                        isMuonCapture = 1;
                                                                }

                                                                if(genId==GenId::muonCapture){
                                                                        isMuonDecayInFlight = 1;
                                                                }

                                                                if(genId==GenId::muonCapture){
                                                                        isEPlusfromStoppedPi = 1;
                                                                }
                                                        }
                                                        cryMap[mchit.trackId().asInt()].isSignal     = isConv;
                                                        cryMap[mchit.trackId().asInt()].isDIO     = isDIO;
                                                        cryMap[mchit.trackId().asInt()].isPionCapture    = isPionCapture;
                                                        cryMap[mchit.trackId().asInt()].isPionCaptureComb   = isPionCaptureComb;
                                                        cryMap[mchit.trackId().asInt()].isMuonCapture     = isMuonCapture;
                                                        cryMap[mchit.trackId().asInt()].isMuonDecayInFlight    = isMuonDecayInFlight;
                                                        cryMap[mchit.trackId().asInt()].isEPlusfromStoppedPi    = isEPlusfromStoppedPi;
                                                }


                                                CLHEP::Hep3Vector cryFrame = cg->toCrystalFrame(thehit.id(), mchit.position());

                                        }//end loop for(nHitsPerCrystal)


                                }//end loop for(caloClusterHits.size())

                                _clSeeds = totGen;//cryMap.size();

                                if(_clSeeds > 0){
                                        int j=0, tmpVane=-1;
                                        for(CryMap::iterator it = cryMap.begin(); it != cryMap.end(); ++it){
                                                _clSeedTrackId[j]      = it->first;
                                                tmpVane                = it->second.vane;
                                                _clSeedPdgId[j]        = it->second.pdgId;
                                                _clSeedQC[j]           = it->second.qualityCuts;
                                                _clSeedIsGen[j]        = it->second.isGen;

                                                if(_diagLevel < 0){
                                                        cout<<"it->second.isSignal = "<< it->second.isSignal<<endl;
                                                }

                                                _clSeedIsSignal[j]     = it->second.isSignal;
                                                _clSeedIsDIO[j]        = it->second.isDIO;
                                                _clSeedIsMuonCapture[j]= it->second.isMuonCapture;
                                                _clSeedIsMuonDecayInFlight[j]     = it->second.isMuonDecayInFlight;
                                                _clSeedIsEPlusfromStoppedPi[j]     = it->second.isEPlusfromStoppedPi;
                                                _clSeedIsPionCapture[j]     = it->second.isPionCapture;
                                                _clSeedIsPionCaptureComb[j]     = it->second.isPionCaptureComb;
                                                _clSeedTime[j]         = it->second.time ;
                                                _clSeedTotEnergyDep[j] = it->second.totEnergyDep;
                                                _clSeedPx[j]           = it->second.position.x();
                                                _clSeedPy[j]           = it->second.position.y();
                                                _clSeedPz[j]           = it->second.position.z();

                                                _clSeedCryFramePu[j]   = it->second.cryPosition.x();
                                                _clSeedCryFramePv[j]   = it->second.cryPosition.y();
                                                _clSeedCryFramePw[j]   = it->second.cryPosition.z();

                                                CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(tmpVane, it->second.position);
                                                _clSeedVaneFramePu[j]  = vaneFrame.x();
                                                _clSeedVaneFramePv[j]  = vaneFrame.y();
                                                _clSeedVaneFramePw[j]  = vaneFrame.z();
                                                _clSeedPpx[j]          = it->second.momentum.x();
                                                _clSeedPpy[j]          = it->second.momentum.y();
                                                _clSeedPpz[j]          = it->second.momentum.z();
                                                _clSeedE[j]            = it->second.momentum.mag();
                                                _clSeedPp[j]           = it->second.energy;

                                                Vane const &vane = cg->vane(tmpVane);
                                                CLHEP::Hep3Vector Mom_rotated = (vane.rotation())*(it->second.momentum);
                                                _clSeedPpu[j]          = Mom_rotated.x();
                                                _clSeedPpv[j]          = Mom_rotated.y();
                                                _clSeedPpw[j]          = Mom_rotated.z();
                                                double thetaWimpact = std::atan(-1.0*Mom_rotated.getZ() / Mom_rotated.getX() ) ;
                                                _clSeedThetaW[j]       = thetaWimpact*Constants::radToDegrees;
                                                double thetaVimpact = std::atan(Mom_rotated.getY() /  Mom_rotated.getX() ) ;
                                                _clSeedThetaV[j]       = thetaVimpact*Constants::radToDegrees;
                                                ++j;
                                        }

                                }else {
                                        _clSeeds                        = 1;
                                        _clSeedTrackId[0]               = -1;
                                        _clSeedPdgId[0]                 = 0;
                                        _clSeedQC[0]                    = -1;
                                        _clSeedIsGen[0]                 = 0.0;
                                        _clSeedIsSignal[0]              = 0;
                                        _clSeedIsMuonCapture[0]         = 0;
                                        _clSeedIsMuonDecayInFlight[0]   = 0;
                                        _clSeedIsPionCapture[0]         = 0;
                                        _clSeedIsPionCaptureComb[0]     = 0;
                                        _clSeedIsEPlusfromStoppedPi[0]  = 0;
                                        _clSeedIsDIO[0]                 = 0;
                                        _clSeedTime[0]                  = 0.0;
                                        _clSeedTotEnergyDep[0]          = 0.0;
                                        _clSeedPx[0]                    = 0.0;
                                        _clSeedPy[0]                    = 0.0;
                                        _clSeedPz[0]                    = 0.0;
                                        _clSeedCryFramePu[0]            = 0.0;
                                        _clSeedCryFramePv[0]            = 0.0;
                                        _clSeedCryFramePw[0]            = 0.0;
                                        _clSeedVaneFramePu[0]           = 0.0;
                                        _clSeedVaneFramePv[0]           = 0.0;
                                        _clSeedVaneFramePw[0]           = 0.0;
                                        _clSeedPpx[0]                   = 0.0;
                                        _clSeedPpy[0]                   = 0.0;
                                        _clSeedPpz[0]                   = 0.0;
                                        _clSeedE[0]                     = 0.0;
                                        _clSeedPp[0]                    = 0.0;
                                        _clSeedPpu[0]                   = 0.0;
                                        _clSeedPpv[0]                   = 0.0;
                                        _clSeedPpw[0]                   = 0.0;
                                        _clSeedThetaW[0]                = 0.0;
                                        _clSeedThetaV[0]                = 0.0;

                                }



                                if(_diagLevel < 0){
                                        cout<<"-------------------- Filling _Ntup... -----------------------"<<endl;
                                }

                                if(_diagLevel < 0){
                                        cout<<"-------------------- Filled _Ntup -----------------------"<<endl;
                                }

                        }else {
                                _clSeeds                                  = 1;
                                _clSeedTrackId[0]                         = -1;
                                _clSeedPdgId[0]                           = 0;
                                _clSeedIsGen[0]                           = 0.0;
                                _clSeedIsSignal[0]                        = 0;
                                _clSeedTime[0]                            = 0.0;
                                _clSeedTotEnergyDep[0]                    = 0.0;
                                _clSeedPx[0]                              = 0.0;
                                _clSeedPy[0]                              = 0.0;
                                _clSeedPz[0]                              = 0.0;
                                _clSeedCryFramePu[0]                      = 0.0;
                                _clSeedCryFramePv[0]                      = 0.0;
                                _clSeedCryFramePw[0]                      = 0.0;
                                _clSeedVaneFramePu[0]                     = 0.0;
                                _clSeedVaneFramePv[0]                     = 0.0;
                                _clSeedVaneFramePw[0]                     = 0.0;
                                _clSeedPpx[0]                             = 0.0;
                                _clSeedPpy[0]                             = 0.0;
                                _clSeedPpz[0]                             = 0.0;
                                _clSeedE[0]                               = 0.0;
                                _clSeedPpu[0]                             = 0.0;
                                _clSeedPpv[0]                             = 0.0;
                                _clSeedPpw[0]                             = 0.0;
                                _clSeedThetaW[0]                          = 0.0;
                                _clSeedThetaV[0]                          = 0.0;

                        }
                        _Ntup->Fill();

                }//end for(caloClsuters.size)
        }//end if(caloClsuters.size() > 0)

        if(_diagLevel < 0){
                cout <<"-------------------> number of generated electrons ( trkVecTot.size() ) = "<< trkVecTot.size() <<endl;
        }

        cout << "Event "<<evt.id().event()<<" CaloClusterPileUp done..."<<endl;


}

}


using mu2e::CaloClusterPileUp;
DEFINE_ART_MODULE(CaloClusterPileUp);


