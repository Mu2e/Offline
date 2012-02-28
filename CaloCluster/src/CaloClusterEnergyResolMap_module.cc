//
// Visualization of the energy resolution  on the rows and on the columns
//
// $Id: CaloClusterEnergyResolMap_module.cc,v 1.1 2012/02/28 22:24:48 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/02/28 22:24:48 $
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
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "Mu2eUtilities/inc/sort_functors.hh"
//#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "CaloCluster/inc/CaloClusterer.hh"
//#include "CaloCluster/inc/CaloClusterFinder.hh"
//#include "CaloCluster/inc/ClosestCaloClusterFinder.hh"
//#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

// Other includes.
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib/exception.h"
//#include "Analyses/inc/MCCaloUtilities.hh"

//Root includes
#include "TFile.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"

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


struct elecData{
        double _cluEnergy;
        double _cluTime;
        double _impTime;
        double _impEnergy;
        CLHEP::Hep3Vector _cluCog;
        CLHEP::Hep3Vector _impPos;
        CLHEP::Hep3Vector _impPosCryFrame;
        unsigned int _row;
        unsigned int _column;
        unsigned int _Nseed;
        CLHEP::Hep3Vector _cryOrigin;
        CLHEP::Hep3Vector _impMom3Vec;

        bool operator<( const elecData other) const{
                return ( _impTime< other._impTime);
        }
        elecData & operator=(const elecData& other) {
                _cluEnergy = other._cluEnergy;
                _cluTime   = other._cluTime;
                _impTime = other._impTime;
                _impEnergy = other._impEnergy;
                _cluCog    = other._cluCog;
                _impPos = other._impPos;
                _impPosCryFrame = other._impPosCryFrame;
                _row       = other._row;
                _column    = other._column;
                _Nseed     = other._Nseed;
                _cryOrigin = other._cryOrigin;
                _impMom3Vec =other._impMom3Vec;
                return *this;
        }
        elecData():
                _impTime(1e10),
                _impEnergy(0.0),
                _Nseed(0){
        }
};

struct cluData{
        double _cluEnergy;
        double _cluTime;
        double _impTime;
        double _impMomentum;
        //CLHEP::Hep3Vector _cluCog;
        CLHEP::Hep3Vector _impPos;
       // CLHEP::Hep3Vector _impPosCryFrame;
        //CLHEP::Hep3Vector _cryOrigin;
        CLHEP::Hep3Vector _impMom3Vec;

        bool operator<( const cluData other) const{
                return ( _impTime< other._impTime);
        }
        cluData & operator=(const cluData& other) {
                _cluEnergy = other._cluEnergy;
                _cluTime   = other._cluTime;
                _impTime = other._impTime;
                _impMomentum = other._impMomentum;
                //_cluCog    = other._cluCog;
                _impPos = other._impPos;
                //_impPosCryFrame = other._impPosCryFrame;
                //_cryOrigin = other._cryOrigin;
                _impMom3Vec =other._impMom3Vec;
                return *this;
        }
        cluData():
                _impTime(1e10),
                _impMomentum(0.0){
        }
};




//the key is the the vane
typedef std::map<unsigned int,std::map<unsigned int, elecData > > ElecMap;

//the first key is the row index, the second is the column index
typedef std::map<unsigned int, std::map<unsigned int, std::map<unsigned int, cluData > > > inCluMap;

typedef std::map<unsigned int, inCluMap > CluMap;

static int ncalls(0);

class CaloClusterEnergyResolMap : public art::EDAnalyzer {
public:
        explicit CaloClusterEnergyResolMap(fhicl::ParameterSet const& pset):
        _diagLevel(pset.get<int>("diagLevel",0)),
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
        _hTHistDeltaEnergy_row(0),
        _hTHistDeltaEnergyRec_row(0),
        _hTHistDeltaEnergy_column(0),
        _hTHistDeltaEnergyRec_column(0),
        _hTHistDeltaPmag(0),
        _hTHistDeltaPitch(0),
        _hTHistDeltaPfirst(0),
        _hTHistDeltaPlast (0),
        _hTHistDeltaEnergyClu(0),
        _hTHistDeltaEnergyCluRec(0),
        _hTHistThetaZ(0),
        _hTHistThetaZRec(0),
        _hTHistPhi(0),
        _hTHistPhiRec(0),
        _hTHistPhi_v0(0),
        _hTHistPhiRec_v0(0),
        _hTHistPhi_v1(0),
        _hTHistPhiRec_v1(0),
        _hTHistPhi_v2(0),
        _hTHistPhiRec_v2(0),
        _hTHistPhi_v3(0),
        _hTHistPhiRec_v3(0),
        _hTHistMulti(0),
        _hTHistEnergyCluMultiDIOdio(0),
        _hTHistEnergyMultiDIOdio(0),
        _hTHistTimeCluMultiDIOdio(0),
        _hTHistTimeMultiDIOdio(0),
        _hTHistEnergyCluMultiDIOsignal(0),
        _hTHistEnergyMultiDIOsignal(0),
        _hTHistTimeCluMultiDIOsignal(0),
        _hTHistTimeMultiDIOsignal(0),
        _hTHistEnergyMultiDIO(0),
        _hTHistTimeMultiDIO(0),
        _hTHistEnergyMultiSignal(0),
        _hTHistTimeMultiSignal(0),
//        _hTHistEnergyCluMulti(0),
//        _hTHistEnergyMulti(0),
//        _hTHistTimeCluMulti(0),
//        _hTHistTimeMulti(0),
        //_hTHistTest(0),
        _EnergyClusterCut(pset.get<double>("energyClusterCut",60.)),//MeV
        _application(0),
        _directory(0)
        {
        }
        virtual ~CaloClusterEnergyResolMap() {
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

        TH2D* _hTHistDeltaEnergy_row;
        TH2D* _hTHistDeltaEnergyRec_row;
        TH2D* _hTHistDeltaEnergy_column;
        TH2D* _hTHistDeltaEnergyRec_column;

        TH1D* _hTHistDeltaPmag;
        TH1D* _hTHistDeltaPitch;
        TH1D* _hTHistDeltaPfirst;
        TH1D* _hTHistDeltaPlast;

        TH1D* _hTHistDeltaEnergyClu;
        TH1D* _hTHistDeltaEnergyCluRec;

        TH1D* _hTHistThetaZ;
        TH1D* _hTHistThetaZRec;
        TH1D* _hTHistPhi;
        TH1D* _hTHistPhiRec;
        TH1D* _hTHistPhi_v0;
        TH1D* _hTHistPhiRec_v0;
        TH1D* _hTHistPhi_v1;
        TH1D* _hTHistPhiRec_v1;
        TH1D* _hTHistPhi_v2;
        TH1D* _hTHistPhiRec_v2;
        TH1D* _hTHistPhi_v3;
        TH1D* _hTHistPhiRec_v3;
        TH2D* _hTHistMulti;
        TH1D* _hTHistEnergyCluMultiDIOdio;
        TH1D* _hTHistEnergyMultiDIOdio;
        TH1D* _hTHistTimeCluMultiDIOdio;
        TH1D* _hTHistTimeMultiDIOdio;
        TH1D* _hTHistEnergyCluMultiDIOsignal;
        TH1D* _hTHistEnergyMultiDIOsignal;
        TH1D* _hTHistTimeCluMultiDIOsignal;
        TH1D* _hTHistTimeMultiDIOsignal;
        TH1D* _hTHistEnergyMultiDIO;
        TH1D* _hTHistTimeMultiDIO;
        TH1D* _hTHistEnergyMultiSignal;
        TH1D* _hTHistTimeMultiSignal;


        //TH1D* _hTHistTest;
        double _EnergyClusterCut;

        std::auto_ptr<MCCaloUtilities> CaloManager;

        bool _skipEvent;

        //TCanvas*      _fakeCanvas;

        // The job needs exactly one instance of TApplication.  See note 1.
        auto_ptr<TApplication> _application;

        // Save directory from beginJob so that we can go there in endJob. See note 3.
        TDirectory* _directory;


};

bool findTrkId(std::vector<unsigned int> vec, unsigned int t){
        bool res = false;

        unsigned int size = vec.size();
        if(size!=0){
                unsigned int cont = 0;
                while(/*!res ||*/ cont!=size){
                        if(vec[cont] == t) {
                                res = true;
                        }
                        ++cont;

                }
        }
        return res;
}


void CaloClusterEnergyResolMap::beginJob( ) {
        //        ++_nAnalyzed;
        //        ++mu2e::ncalls;

        cout << "start CaloClusterEnergyResolMap..."<<endl;

        CaloManager = auto_ptr<MCCaloUtilities>(new MCCaloUtilities());

        //art::ServiceHandle<art::TFileService> tfs;

        // If needed, create the ROOT interactive environment. See note 1.
        if ( !gApplication ){
                int    tmp_argc(0);
                char** tmp_argv(0);
                _application = auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
        }

        gStyle->SetPalette(1);
        gROOT->SetStyle("Plain");

        // _fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);

        // See note 3.
        _directory = gDirectory;


}

void CaloClusterEnergyResolMap::analyze(art::Event const & evt ) {

        ++_nAnalyzed;
        ++ncalls;

        CaloClusterer c;

        art::ServiceHandle<GeometryService> geom;
        GeomHandle<Calorimeter> cg;
        if (ncalls == 1) {

                // cout << "This should be done only in the first event" << endl;
                art::ServiceHandle<art::TFileService> tfs;

                _hTHistDeltaEnergy_row   = tfs->make<TH2D>( "DeltaEnergy_row", "DeltaEnergy;row index;Eseed-Eclu [MeV]", cg->nCrystalR(), 0, cg->nCrystalR() , 320, -50., 110.);
                _hTHistDeltaEnergyRec_row= tfs->make<TH2D>( "DeltaEnergyRec_row", "DeltaEnergyRec;row index;Eseed-Eclu [MeV]", cg->nCrystalR(), 0, cg->nCrystalR(), 320, -50., 110.);
                _hTHistDeltaEnergy_column   = tfs->make<TH2D>( "DeltaEnergy_column", "DeltaEnergy;column index;Eseed-Eclu [MeV]", cg->nCrystalZ(), 0, cg->nCrystalZ() , 320, -50., 110.);
                _hTHistDeltaEnergyRec_column= tfs->make<TH2D>( "DeltaEnergyRec_column", "DeltaEnergyRec;column index;Eseed-Eclu [MeV]", cg->nCrystalZ(), 0, cg->nCrystalZ() , 320, -50., 110.);

                _hTHistDeltaPmag  = tfs->make<TH1D>( "TrkDeltaPmag", "TrkDeltaPmag;deltaP/Pfirst [%];entries", 2000 , 0., 0.01);
                _hTHistDeltaPitch = tfs->make<TH1D>( "TrkDeltaPitch","TrkDeltaPitch;deltaPitch/pitchFirst [%];entries", 2000 , 0., 0.01);

                _hTHistDeltaPfirst = tfs->make<TH1D>( "TrkPfirst", "momFirst; momentum [MeV];entries", 2000 , 50., 110.0);
                _hTHistDeltaPlast  = tfs->make<TH1D>( "TrkPlast","momLast ; momentum [MeV];entries", 2000 , 50., 110.0);

                _hTHistDeltaEnergyClu   = tfs->make<TH1D>( "DeltaEnergyClu", "DeltaEnergyClu; E_{seed} - Energy_{cluster} [MeV];entries", 320,  -50.0, 110.);
                _hTHistDeltaEnergyCluRec= tfs->make<TH1D>( "DeltaEnergyCluRec","DeltaEnergyCluRec ; E_{seed} - Energy_{cluster} [MeV];entries", 320, -50.0, 110.);

                _hTHistThetaZ      = tfs->make<TH1D>( "CosThetaZ", "CosThetaZ; cos#theta [MeV];entries", 720 , -180., 180.0);
                _hTHistThetaZRec      = tfs->make<TH1D>( "CosThetaZRec", "CosThetaZRec; cos#theta [MeV];entries",720 , -180., 180.0);

                _hTHistPhi         = tfs->make<TH1D>( "CosPhi","CosPhi ; cos#phi ;entries", 720 , -180., 180.0);
                _hTHistPhiRec         = tfs->make<TH1D>( "CosPhiRec","CosPhiRec ; cos#phi ;entries", 720 , -180., 180.0);
                _hTHistPhi_v0         = tfs->make<TH1D>( "CosPhi_v0","CosPhi_v0 ; cos#phi ;entries", 720 , -180., 180.0);
                _hTHistPhiRec_v0         = tfs->make<TH1D>( "CosPhiRec_v0","CosPhiRec_v0 ; cos#phi ;entries", 720 , -180., 180.0);
                _hTHistPhi_v1         = tfs->make<TH1D>( "CosPhi_v1","CosPhi_v1 ; cos#phi ;entries", 720 , -180., 180.0);
                _hTHistPhiRec_v1      = tfs->make<TH1D>( "CosPhiRec_v1","CosPhiRec_v1 ; cos#phi ;entries", 720 , -180., 180.0);
                _hTHistPhi_v2         = tfs->make<TH1D>( "CosPhi_v2","CosPhi_v2 ; cos#phi ;entries", 720 , -180., 180.0);
                _hTHistPhiRec_v2      = tfs->make<TH1D>( "CosPhiRec_v2","CosPhiRec_v2 ; cos#phi ;entries", 720 , -180., 180.0);
                _hTHistPhi_v3         = tfs->make<TH1D>( "CosPhi_v3","CosPhi_v3 ; cos#phi ;entries", 720 , -180., 180.0);
                _hTHistPhiRec_v3      = tfs->make<TH1D>( "CosPhiRec_v3","CosPhiRec_v3 ; cos#phi ;entries", 720 , -180., 180.0);

                _hTHistMulti          = tfs->make<TH2D>( "MultiSeedDistr","MultiSeedDistr ; row index ; column index",  cg->nCrystalR(), 0, cg->nCrystalR(), cg->nCrystalZ(), 0, cg->nCrystalZ() );

                _hTHistEnergyCluMultiDIOdio = tfs->make<TH1D>( "EnergyCluMultiDIOdio","EnergyCluMultiDIOdio ; E_{cluster} [MeV];entries", 1000 , 0., 250.0);
                _hTHistEnergyMultiDIOdio    = tfs->make<TH1D>( "EnergyMultiDIOdio","EnergyMultiDIOdio ; E_{seed} ;entries", 440 , 0., 110.0);
                _hTHistTimeCluMultiDIOdio   = tfs->make<TH1D>( "TimeCluMultiDIOdio","TimeCluMultiDIOdio ; Time_{cluster} [ns] ;entries", 300 , -0., 3000.0);
                _hTHistTimeMultiDIOdio      = tfs->make<TH1D>( "TimeMultiDIOdio","TimeMultiDIOdio ; Time_{seed} [ns] ;entries", 300 , 0., 3000.0);

                _hTHistEnergyCluMultiDIOsignal = tfs->make<TH1D>( "EnergyCluMultiDIOsignal","EnergyCluMultiDIOsignal ; E_{cluster} [MeV];entries", 1000 , 0., 250.0);
                _hTHistEnergyMultiDIOsignal    = tfs->make<TH1D>( "EnergyMultiDIOsignal","EnergyMultiDIOsignal ; E_{seed} ;entries", 440 , 0., 110.0);
                _hTHistTimeCluMultiDIOsignal   = tfs->make<TH1D>( "TimeCluMultiDIOsignal","TimeCluMultiDIOsignal ; Time_{cluster} [ns] ;entries", 300 , -0., 3000.0);
                _hTHistTimeMultiDIOsignal      = tfs->make<TH1D>( "TimeMultiDIOsignal","TimeMultiDIOsignal ; Time_{seed} [ns] ;entries", 300 , 0., 3000.0);

                _hTHistEnergyMultiDIO    = tfs->make<TH1D>( "EnergyMultiDIO","EnergyMultiDIO ; E_{seed} ;entries", 440 , 0., 110.0);
                _hTHistTimeMultiDIO      = tfs->make<TH1D>( "TimeMultiDIO","TimeMultiDIO ; Time_{seed} [ns] ;entries", 300 , 0., 3000.0);

                _hTHistEnergyMultiSignal    = tfs->make<TH1D>( "EnergyMultiSignal","EnergyMultiSignal ; E_{seed} ;entries", 440 , 0., 110.0);
                _hTHistTimeMultiSignal      = tfs->make<TH1D>( "TimeMultiSignal","TimeMultiSignal ; Time_{seed} [ns] ;entries", 300 , 0., 3000.0);
               // _hTHistTest        = tfs->make<TH1D>( "CosPhiRecProva","CosPhiRecProva ; cos#phi ;entries", 720 , -180., 180.0);
        }

        doCalorimeter(evt, _skipEvent);

} // end of analyze

void CaloClusterEnergyResolMap::endJob() {

}



void CaloClusterEnergyResolMap::doCalorimeter(art::Event const& evt, bool skip){


        if ( _diagLevel > 0 ) cout << "MakeCaloCluster: produce() begin" << endl;

        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;
        if(! geom->hasElement<Calorimeter>() ) return;
        GeomHandle<Calorimeter> cg;

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
        int NtrkCut =0;
        int NtrkTot = 0;

        //mapping &counting the electrons with quality cuts in the TRK
        for ( genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
                ++NtrkTot;
                VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);
                GenElHitData& hdil = iEltrk.getithLoopHit(0);
                GenElHitData& ldil = iEltrk.getHit((int)(iEltrk.getNumOfHit() - 1) );

                if(!findTrkId(tmpVTot, iEltrk.getTrkID().asUint() ) ){

                        //cout<<"faccio il puschback in tmpV..."<<endl;
                        tmpVTot.push_back( iEltrk.getTrkID().asUint() );
                        //cout<<"--------> puschback in tmpV eseguito..."<<endl;

                }

                double cosTheta = iEltrk.getTrkLrntzVec().cosTheta() ;

                double cosPitch = hdil._hitMomentum.cosTheta();

                double lcosPitch = ldil._hitMomentum.cosTheta();

//
//                cout << "cosThetafirst = "<<cosPitch <<", "<<"costhetaLast = "<<lcosPitch<<endl;
//                cout<< "ldil._hitMomentum.mag() = "<<ldil._hitMomentum.mag()<<"hdil._hitMomentum.mag() = "<< hdil._hitMomentum.mag()<<endl;


                bool condition = true;
                condition &= ( iEltrk.getNumOfHit() >= 20 );
                condition &= ( hdil._hitMomentum.mag() >= trkMomCut );
                condition &= ( cosTheta >= -0.5 );
                condition &= ( cosTheta <=  0.5 );
                condition &= ( cosPitch > 0.5 );
                condition &= ( cosPitch < 0.70710678118655 );// 1 / sqrt(2)
                //condition &= ( lcosPitch / cosPitch ) >= (1.0 - iEltrk.getNumOfHit()*8.e-4);
                //condition &= (ldil._hitMomentum.mag() /hdil._hitMomentum.mag() ) >= (1.0 - iEltrk.getNumOfHit()*2.e-4);

                _hTHistDeltaPfirst->Fill(hdil._hitMomentum.mag());
                _hTHistDeltaPlast->Fill( ldil._hitMomentum.mag());
                _hTHistDeltaPmag->Fill((hdil._hitMomentum.mag() - ldil._hitMomentum.mag()) /hdil._hitMomentum.mag()/ iEltrk.getNumOfHit());
                _hTHistDeltaPitch->Fill((cosPitch - lcosPitch)/ cosPitch /iEltrk.getNumOfHit());

                if( condition){
                        NtrkCut++;// # of electrons at the entrance of the TRK which have trkMomCut MeV of kinetic energy and which will light at least 20 straws
                        if(!findTrkId(tmpV, iEltrk.getTrkID().asUint() ) ){

                                //cout<<"faccio il puschback in tmpV..."<<endl;
                                tmpV.push_back( iEltrk.getTrkID().asUint() );
                                if ( iEltrk.isConversionEl() ) {
                                        tmpVsignal.push_back( iEltrk.getTrkID().asUint() );
                                }
                                //cout<<"--------> puschback in tmpV eseguito..."<<endl;

                        }

                }



        }//end TRK mapping
        //cout <<"NtrkCut = " <<NtrkCut<<endl;

        if (NtrkTot==0) return;

        std::vector<unsigned int> trkVec, trkVecTot;
        for(unsigned int f = 0; f != tmpV.size(); ++f){
                trkVec.push_back(tmpV[f]);
        }
        for(unsigned int f = 0; f != tmpVTot.size(); ++f){
                trkVecTot.push_back(tmpVTot[f]);
        }

        ElecMap elecMap;
        CluMap elecMapMulti;

        //adesso leggo i cluster che il mio evento ha generato
        if(caloClusters->size()>0 ){
                int iVane;
                for(size_t icl=0; icl<caloClusters->size(); ++icl){

                        double eDepClu = 0.;

                        CaloCluster const& clu = (*caloClusters).at(icl);
                        eDepClu = clu._energy;
                        iVane = clu._iVane;

                        CaloCrystalHitPtrVector caloClusterHits = clu._caloCrystalHits;
                        if(eDepClu >= _EnergyClusterCut){
                                //cout <<"eDepClu >_energyClusterCut"<< endl;

                                // std::vector<elecData> vecElecData;
                                std::map<unsigned int, unsigned int> seedMap, signalMap;

                                for(size_t i=0; i<caloClusterHits.size(); ++i){
                                        CaloCrystalHit const& hit = *(caloClusterHits.at(i));

                                        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                                        if(ROIds.size()<1 ) continue;

                                        CaloHit const& thehit = *ROIds.at(0);
                                        size_t collectionPosition = ROIds.at(0).key();

                                        PtrStepPointMCVector const & mcptr(hits_mcptr->at(collectionPosition));
                                        size_t nHitsPerCrystal = mcptr.size();
                                        //size_t tmpJ2 = -1;

                                        for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {

                                                //cout<< "Start loop..."<< "j2 = "<< j2<<endl;

                                                StepPointMC const& mchit = *mcptr[j2];

                                                // The simulated particle that made this hit.
                                                SimParticleCollection::key_type trackId(mchit.trackId());

                                                CLHEP::Hep3Vector cryFrame = cg->toCrystalFrame(thehit.id(), mchit.position());

                                                if( (cg->getCrystalRByRO(thehit.id()) !=(_rowToCanc - 1) && cg->getCrystalZByRO(thehit.id()) != (_columnToCanc-1) ) ){
                                                        unsigned int row = cg->getCrystalRByRO(thehit.id() );
                                                        unsigned int column = cg->getCrystalZByRO(thehit.id() );

                                                        if(elecMap[iVane][trackId.asUint()]._impTime > mchit.time() ){
                                                                bool search = false, searchSignal = false;
                                                                for(unsigned int i=0; i!= trkVecTot.size(); ++i){
                                                                        if(trkVecTot.at(i) == trackId.asUint() ){
                                                                                search = true;
                                                                        }
                                                                }
                                                                for(unsigned int i=0; i!= tmpVsignal.size(); ++i){
                                                                        if(tmpVsignal.at(i) == trackId.asUint() ){
                                                                                searchSignal = true;
                                                                        }
                                                                }
                                                                if(search){
                                                                        elecMap[iVane][trackId.asUint()]._Nseed = elecMap[iVane][trackId.asUint()]._Nseed + 1;
                                                                        seedMap[trackId.asUint()] = 1;
                                                                }
                                                                if(searchSignal){
                                                                        signalMap[trackId.asUint()] = 1;
                                                                }
                                                                elecMap[iVane][trackId.asUint()]._cluEnergy = clu._energy;
                                                                elecMap[iVane][trackId.asUint()]._cluTime = clu._time;
                                                                elecMap[iVane][trackId.asUint()]._impTime = mchit.time();
                                                                elecMap[iVane][trackId.asUint()]._impEnergy = mchit.momentum().mag();
                                                                elecMap[iVane][trackId.asUint()]._cluCog = clu._impactPoint;
                                                                elecMap[iVane][trackId.asUint()]._impPos = mchit.position();
                                                                elecMap[iVane][trackId.asUint()]._impPosCryFrame = cryFrame;
                                                                elecMap[iVane][trackId.asUint()]._row    = row;
                                                                elecMap[iVane][trackId.asUint()]._column = column;
                                                                elecMap[iVane][trackId.asUint()]._cryOrigin  = cg->getCrystalOriginByRO( thehit.id() );
                                                                elecMap[iVane][trackId.asUint()]._impMom3Vec = mchit.momentum();

                                                                //                                                                cout<< "###################"<<endl;
                                                                //                                                                cout<< "idVande = "<< iVane<<endl;
                                                                //                                                                cout << "cluX = "<<elecMap[iVane][trackId.asUint()]._impPos.getX()<<endl;
                                                                //                                                                cout << "cluY = "<<elecMap[iVane][trackId.asUint()]._impPos.getY()<<endl;
                                                                //                                                                cout << "cluZ = "<<elecMap[iVane][trackId.asUint()]._impPos.getZ()<<endl;
                                                                //                                                                cout<< "###################"<<endl;
                                                        }
//                                                        cout<<"----------- step 1 --------------"<<endl;
//
//                                                        cout<<"----------- step 2:--------------"<<endl;
//
//                                                        cout<<"----------- step 2_bis: seed found --------------"<<endl;

//                                                        bool search = false;
//                                                        for(unsigned int i=0; i!= trkVecTot.size(); ++i){
//                                                                if(trkVecTot.at(i) == trackId.asUint() ){
//                                                                        search = true;
//                                                                }
//                                                        }
//
//                                                        if( search ){
//                                                                elecMapMulti[iVane][row][column][trackId.asUint()]._cluEnergy  = clu._energy;
//                                                                elecMapMulti[iVane][row][column][trackId.asUint()]._cluTime    = clu._time;
//                                                                elecMapMulti[iVane][row][column][trackId.asUint()]._impTime    = mchit.time();
//                                                                elecMapMulti[iVane][row][column][trackId.asUint()]._impMomentum = mchit.momentum().mag();
//                                                                elecMapMulti[iVane][row][column][trackId.asUint()]._impPos   = mchit.position();
//                                                                elecMapMulti[iVane][row][column][trackId.asUint()]._impMom3Vec  = mchit.momentum();
//                                                                cout<<"cluE = "<<elecMapMulti[iVane][row][column][trackId.asUint()]._cluEnergy<<endl;
//                                                                cout<<"cluT = "<<elecMapMulti[iVane][row][column][trackId.asUint()]._cluTime<<endl;
//                                                                cout<<"elecMapMulti Fillata"<<endl;
//                                                        }


                                                }
                                        }

                                }// esco da questo ciclo con: time e Id della particella che per prima ha dato origine al cluster
                                if(seedMap.size() > 1){
                                        if(signalMap.size() > 0){
                                                _hTHistEnergyCluMultiDIOsignal->Fill( eDepClu);
                                                _hTHistTimeCluMultiDIOsignal->Fill( clu._time);

                                                for(std::map<unsigned int, unsigned int>::iterator its = seedMap.begin(); its != seedMap.end(); ++its){
                                                        //std::map<unsigned int, unsigned int>::iterator its = seedMap.begin();
                                                        //its = its + i;
                                                        std::map<unsigned int, elecData >::iterator it = elecMap[iVane].find(its->first);
                                                        _hTHistEnergyMultiDIOsignal->Fill( (it)->second._impEnergy );
                                                        _hTHistTimeMultiDIOsignal->Fill( (it)->second._impTime);

                                                        if(signalMap.find(its->first) != signalMap.end() ){
                                                                _hTHistEnergyMultiSignal->Fill( (it)->second._impEnergy );
                                                                _hTHistTimeMultiSignal->Fill( (it)->second._impTime);
                                                        }else{
                                                                _hTHistEnergyMultiDIO->Fill( (it)->second._impEnergy );
                                                                _hTHistTimeMultiDIO->Fill( (it)->second._impTime);
                                                        }
                                                }
                                        }else{
                                                _hTHistEnergyCluMultiDIOdio->Fill( eDepClu);
                                                _hTHistTimeCluMultiDIOdio->Fill( clu._time);


                                                for(std::map<unsigned int, unsigned int>::iterator its = seedMap.begin(); its != seedMap.end(); ++its){
                                                        //std::map<unsigned int, unsigned int>::iterator its = seedMap.begin();
                                                        //its = its + i;
                                                        std::map<unsigned int, elecData >::iterator it = elecMap[iVane].find(its->first);
                                                        _hTHistEnergyMultiDIOdio->Fill( (it)->second._impEnergy );
                                                        _hTHistTimeMultiDIOdio->Fill( (it)->second._impTime);
                                                }
                                        }
                                }



                        }// end if(eDepClu >= (_hTHistEff->GetXaxis()->GetBinCenter(bin) -  0.5*_hTHistEff->GetXaxis()->GetBinWidth(bin) ) ){

                }//end for(caloClsuters.size)
        }//end if(caloClsuters.size() > 0)

        //counting how many reconstructed clusters match with the generated particles
        unsigned int canc2 = 0;

        //cout <<"-------------------> numero di elettroni generati = "<< trkVecTot.size() <<endl;
        unsigned int size2 = trkVecTot.size();
        unsigned int it2=0;
        while( it2 < size2){
                ElecMap::iterator ite = elecMap.begin();
                bool trovato = false;
                while(!trovato && ite!=elecMap.end() ){
                        if(ite->second.find(trkVecTot[it2]) != ite->second.end()){
                                //                                        cout<< "$$$$$$$$$$$$----> vane = "<< ite->first<<endl;
                                //                                        cout<<"------------------------------------------------------"<<endl;
                                //                                        cout << "impX = "<<ite->second[trkVec[it2]]._impPos.getX()<<endl;
                                //                                        cout << "impY = "<<ite->second[trkVec[it2]]._impPos.getY()<<endl;
                                //                                        cout << "impZ = "<<ite->second[trkVec[it2]]._impPos.getZ()<<endl;
                                //                                        cout << "impXcryFRame = "<<ite->second[trkVec[it2]]._impPosCryFrame.getX()<<endl;
                                //                                        cout << "impYcryFRame = "<<ite->second[trkVec[it2]]._impPosCryFrame.getY()<<endl;
                                //                                        cout << "impZcryFRame = "<<ite->second[trkVec[it2]]._impPosCryFrame.getZ()<<endl;
                                //                                        cout<< "@@@@@@@@@@@@@@@@@@@"<<endl;
                                //cout << "cluX = "<<ite->second[trkVec[it2]]._cluCog.getX()<<endl;
                                //                                        cout << "cluY = "<<ite->second[trkVec[it2]]._cluCog.getY()<<endl;
                                //                                        cout << "cluZ = "<<ite->second[trkVec[it2]]._cluCog.getZ()<<endl;
                                //                                        cout<< "@@@@@@@@@@@@@@@@@@@"<<endl;
                                //                                        cout << "delta X = " << ite->second[trkVec[it2]]._impPos.getX() - ite->second[trkVec[it2]]._cluCog.getX()<<endl;
                                //                                        cout << "delta Y = " << ite->second[trkVec[it2]]._impPos.getY() - ite->second[trkVec[it2]]._cluCog.getY()<<endl;
                                //                                        cout << "delta Z = " << ite->second[trkVec[it2]]._impPos.getZ() - ite->second[trkVec[it2]]._cluCog.getZ()<<endl;
                                _hTHistDeltaEnergyCluRec->Fill( ite->second[trkVecTot[it2]]._impEnergy - ite->second[trkVecTot[it2]]._cluEnergy);

                                _hTHistDeltaEnergyRec_row->Fill(_hTHistDeltaEnergyRec_row->GetXaxis()->GetBinCenter(ite->second[trkVecTot[it2]]._row + 1), ite->second[trkVecTot[it2]]._impEnergy - ite->second[trkVecTot[it2]]._cluEnergy);
                                _hTHistDeltaEnergyRec_column->Fill(_hTHistDeltaEnergyRec_column->GetXaxis()->GetBinCenter(ite->second[trkVecTot[it2]]._column+ 1), ite->second[trkVecTot[it2]]._impEnergy - ite->second[trkVecTot[it2]]._cluEnergy);

                                _hTHistThetaZRec->Fill(ite->second[trkVecTot[it2]]._impMom3Vec.getTheta()*180./TMath::Pi());
                                _hTHistPhiRec->Fill(ite->second[trkVecTot[it2]]._impMom3Vec.getPhi()*180./TMath::Pi());
//                                _hTHistTest->Fill( std::acos(ite->second[trkVecTot[it2]]._impMom3Vec.getZ()/ite->second[trkVecTot[it2]]._impMom3Vec.mag() )*180./TMath::Pi() );
//                                cout<<"ite->second[trkVecTot[it2]]._Nseed = "<<ite->second[trkVecTot[it2]]._Nseed<<endl;
//                                if(ite->second[trkVecTot[it2]]._Nseed > 1){
//                                        _hTHistMulti->Fill(ite->second[trkVecTot[it2]]._row + 1, ite->second[trkVecTot[it2]]._column + 1);
//                                        _hTHistEnergyCluMulti->Fill( ite->second[trkVecTot[it2]]._cluEnergy);
//                                        _hTHistEnergyMulti->Fill( ite->second[trkVecTot[it2]]._impEnergy);
//                                        _hTHistTimeCluMulti->Fill( ite->second[trkVecTot[it2]]._cluTime);
//                                        _hTHistTimeMulti->Fill( ite->second[trkVecTot[it2]]._impTime);
//                                }

                                if(ite->first == 0 ){
                                        _hTHistPhiRec_v0->Fill(ite->second[trkVecTot[it2]]._impMom3Vec.getPhi()*180./TMath::Pi());
                                }
                                if(ite->first == 1 ){
                                        _hTHistPhiRec_v1->Fill(ite->second[trkVecTot[it2]]._impMom3Vec.getPhi()*180./TMath::Pi());
                                }
                                if(ite->first == 2 ){
                                        _hTHistPhiRec_v2->Fill(ite->second[trkVecTot[it2]]._impMom3Vec.getPhi()*180./TMath::Pi());
                                }
                                if(ite->first == 3 ){
                                        _hTHistPhiRec_v3->Fill(ite->second[trkVecTot[it2]]._impMom3Vec.getPhi()*180./TMath::Pi());
                                }
                                trovato = true;
                                canc2 = it2;
                                //cout<<"$$$ yes! the calo works! $$$"<<endl;
                        }
                        ++ite;
                }

                if(trovato){
                        std::vector<unsigned int>::iterator er = trkVecTot.begin();
                        er +=canc2;
                        trkVecTot.erase(er);
                        size2 = trkVecTot.size();
                }else{
                        ++it2;
                }

        }
        //                cout <<"after.... trkVecTot.size = "<< trkVecTot.size() <<endl;
        //
        //                cout<<"...\n...\n "<<endl;

        //-----------------------------------
        //counting how many reconstructed clusters match with the generated particles which pass the quality cuts of the tracker

        unsigned int canc = 0;

        //cout <<"-------------------> numero di elettroni di qualita' generati = "<< trkVec.size() <<endl;
        //for(unsigned int it = 0; it != trkVec.size(); ++it){
        unsigned int size = trkVec.size();
        unsigned int it=0;
        while( it < size){
                ElecMap::iterator ite = elecMap.begin();
                bool trovato = false;
                while(!trovato && ite!=elecMap.end() ){
                        if(ite->second.find(trkVec[it]) != ite->second.end()){
                                //                                        cout << "impX = "<<ite->second[trkVec[it]]._impPosCryFrame.getX()<<endl;
                                //                                        cout << "impY = "<<ite->second[trkVec[it]]._impPosCryFrame.getY()<<endl;
                                //                                        cout << "impZ = "<<ite->second[trkVec[it]]._impPosCryFrame.getZ()<<endl;

                                _hTHistDeltaEnergyClu->Fill( ite->second[trkVec[it]]._impEnergy - ite->second[trkVec[it]]._cluEnergy);

                                _hTHistDeltaEnergy_row->Fill(_hTHistDeltaEnergy_row->GetXaxis()->GetBinCenter(ite->second[trkVec[it]]._row + 1), ite->second[trkVec[it]]._impEnergy - ite->second[trkVec[it]]._cluEnergy);
                                _hTHistDeltaEnergy_column->Fill(_hTHistDeltaEnergy_column->GetXaxis()->GetBinCenter(ite->second[trkVec[it]]._column + 1), ite->second[trkVec[it]]._impEnergy - ite->second[trkVec[it]]._cluEnergy);

                                _hTHistThetaZ->Fill(ite->second[trkVecTot[it]]._impMom3Vec.getTheta()*180./TMath::Pi());
                                _hTHistPhi->Fill(ite->second[trkVecTot[it]]._impMom3Vec.getPhi()*180./TMath::Pi());

                                if(ite->first == 0 ){
                                        _hTHistPhi_v0->Fill(ite->second[trkVecTot[it]]._impMom3Vec.getPhi()*180./TMath::Pi());
                                }
                                if(ite->first == 1 ){
                                        _hTHistPhi_v1->Fill(ite->second[trkVecTot[it]]._impMom3Vec.getPhi()*180./TMath::Pi());
                                }
                                if(ite->first == 2 ){
                                        _hTHistPhi_v2->Fill(ite->second[trkVecTot[it]]._impMom3Vec.getPhi()*180./TMath::Pi());
                                }
                                if(ite->first == 3 ){
                                        _hTHistPhi_v3->Fill(ite->second[trkVecTot[it]]._impMom3Vec.getPhi()*180./TMath::Pi());
                                }

                                trovato = true;
                                canc = it;

                        }
                        ++ite;
                }

                if(trovato){
                        std::vector<unsigned int>::iterator er = trkVec.begin();
                        er +=canc;
                        trkVec.erase(er);
                        size = trkVec.size();
                }else{
                        ++it;
                }

        }


        //                cout <<"after.... trkVec.size = "<< trkVec.size() <<endl;
        cout << "Event "<<evt.id().event()<<" CaloClusterEnergyResolMap done..."<<endl;


}

}


using mu2e::CaloClusterEnergyResolMap;
DEFINE_ART_MODULE(CaloClusterEnergyResolMap);

