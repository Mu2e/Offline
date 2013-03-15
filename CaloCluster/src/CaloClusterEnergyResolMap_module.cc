//
// Visualization of the energy resolution  on the rows and on the columns
//
// $Id: CaloClusterEnergyResolMap_module.cc,v 1.14 2013/03/15 15:52:03 kutschke Exp $
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
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "CaloCluster/inc/CaloClusterer.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

// Other includes.
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib/exception.h"

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


struct electronData{
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
        double _cryTime;

        bool operator<( const electronData other) const{
                return ( _impTime< other._impTime);
        }
        electronData & operator=(const electronData& other) {
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
                _cryTime    = other._cryTime;
                return *this;
        }
        electronData():
                _impTime(1e10),
                _impEnergy(0.0),
                _Nseed(0){
        }
};

struct cryInfo {
        CLHEP::Hep3Vector position;
        CLHEP::Hep3Vector cryPosition;
        CLHEP::Hep3Vector momentum;
        double                time;
        double        totEnergyDep;
        int pdgId;
        int isGen;
        int isSignal;
        int vane;
};
//the key is the trackId
typedef std::map<int , cryInfo> CryMap;


//the key is the the vane
typedef std::map<unsigned int,std::map<unsigned int, electronData > > ElectronMap;

static int ncalls(0);

class CaloClusterEnergyResolMap : public art::EDAnalyzer {
public:
        explicit CaloClusterEnergyResolMap(fhicl::ParameterSet const& pset):
        _diagLevel(pset.get<int>("diagLevel",1)),
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
        _hTHistDeltaEnergyRowProj(0),
        _hTHistDeltaEnergyColumnProj(0),
        _hTHistDeltaEnergyRecRowProj(0),
        _hTHistDeltaEnergyRecColumnProj(0),
        _sigmaVsRow(0),
        _sigmaVsColumn(0),
        _EpeakVsRow(0),
        _EpeakVsColumn(0),
        _hTHistThetaW(0),
        _hTHistThetaW_v0(0),
        _hTHistThetaW_v1(0),
        _hTHistThetaW_v2(0),
        _hTHistThetaW_v3(0),
        _hTHistThetaV(0),
        _hTHistThetaV_v0(0),
        _hTHistThetaV_v1(0),
        _hTHistThetaV_v2(0),
        _hTHistThetaV_v3(0),
        _hTHistCryVsStepPointTime(0),
        _cTCanvSigmaVsRow(0),
        _TGr(0),
        _Ntup(0),

        _EnergyClusterCut(pset.get<double>("energyClusterCut",60.)),//MeV
        _binWidth(0.5),
        _application(nullptr),
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
        TH1D* _hTHistDeltaEnergyRowProj;
        TH1D* _hTHistDeltaEnergyColumnProj;
        TH1D* _hTHistDeltaEnergyRecRowProj;
        TH1D* _hTHistDeltaEnergyRecColumnProj;
        TH1D* _sigmaVsRow;
        TH1D* _sigmaVsColumn;
        TH1D* _EpeakVsRow;
        TH1D* _EpeakVsColumn;
        TH1D* _hTHistMomDotVaneNorm;
        TH1D* _hTHistMomRecDotVaneNorm;
        TH1D* _hTHistThetaW;
        TH1D* _hTHistThetaW_v0;
        TH1D* _hTHistThetaW_v1;
        TH1D* _hTHistThetaW_v2;
        TH1D* _hTHistThetaW_v3;
        TH1D* _hTHistThetaV;
        TH1D* _hTHistThetaV_v0;
        TH1D* _hTHistThetaV_v1;
        TH1D* _hTHistThetaV_v2;
        TH1D* _hTHistThetaV_v3;
        TH1D*  _hTHistCryVsStepPointTime;
        TCanvas *_cTCanvSigmaVsRow;
        TGraphErrors *_TGr;
        TTree* _Ntup;

        //TH1D* _hTHistTest;
        double _EnergyClusterCut;
        double _binWidth;//bin / MeV

        std::unique_ptr<MCCaloUtilities> CaloManager;

        bool _skipEvent;

        //TCanvas*      _fakeCanvas;

        // The job needs exactly one instance of TApplication.  See note 1.
        unique_ptr<TApplication> _application;

        // Save directory from beginJob so that we can go there in endJob. See note 3.
        TDirectory* _directory;

        Int_t _clNo,
        _nCryCl,
        _clSeeds,
        _clVane;
        //          _cryId,
        //          _vane,
        //          _cryPdgId,
        //          _cryIsGen,
        //          _crytrkId,
        //          _crynSteps;

        Float_t _evt,
        _clE,
        _clT,
        _clCOGu,
        _clCOGv,
        _clCOGw;
        //            _cryEdep,
        //            _cryEdepTot,
        //            _cryT;

        Int_t _clSeedTrackId[10000],
        _clSeedPdgId[10000],
        _clSeedIsGen[10000],
        _clSeedIsSignal[10000];

        //        _SPpdgId[10000],
        //          _SPIsGen[10000],
        //          _SPTrkId[10000],
        //          _SPIsSignal[10000];

        Float_t         _clSeedTime[10000],
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

        cout << "start CaloClusterEnergyResolMap..."<<endl;

        CaloManager = unique_ptr<MCCaloUtilities>(new MCCaloUtilities());


        // If needed, create the ROOT interactive environment. See note 1.
        if ( !gApplication ){
                int    tmp_argc(0);
                char** tmp_argv(0);
                _application = unique_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
        }

        gStyle->SetPalette(1);
        gROOT->SetStyle("Plain");

       _directory = gDirectory;


}

void CaloClusterEnergyResolMap::analyze(art::Event const & evt ) {

        ++_nAnalyzed;
        ++ncalls;

        CaloClusterer c;

        art::ServiceHandle<GeometryService> geom;
        GeomHandle<VaneCalorimeter> cg;
        if (ncalls == 1) {

                art::ServiceHandle<art::TFileService> tfs;

                art::TFileDirectory impAngles = tfs->mkdir("ImpactAngles");
                art::TFileDirectory multClu = tfs->mkdir("MultClu");
                art::TFileDirectory trkInfo = tfs->mkdir("TrkInfo");


                _hTHistDeltaEnergy_row   = tfs->make<TH2D>( "DeltaEnergy_row", "DeltaEnergy;row index;Eseed-Eclu [MeV]", cg->nCrystalR(), 0, cg->nCrystalR() , 160 / _binWidth, -50., 110.);
                _hTHistDeltaEnergyRec_row= tfs->make<TH2D>( "DeltaEnergyRec_row", "DeltaEnergyRec;row index;Eseed-Eclu [MeV]", cg->nCrystalR(), 0, cg->nCrystalR(),  160 / _binWidth, -50., 110.);
                _hTHistDeltaEnergy_column   = tfs->make<TH2D>( "DeltaEnergy_column", "DeltaEnergy;column index;Eseed-Eclu [MeV]", cg->nCrystalZ(), 0, cg->nCrystalZ() ,  160 / _binWidth, -50., 110.);
                _hTHistDeltaEnergyRec_column= tfs->make<TH2D>( "DeltaEnergyRec_column", "DeltaEnergyRec;column index;Eseed-Eclu [MeV]", cg->nCrystalZ(), 0, cg->nCrystalZ() ,  160 / _binWidth, -50., 110.);

                _hTHistDeltaPmag  = trkInfo.make<TH1D>( "TrkDeltaPmag", "TrkDeltaPmag;deltaP/Pfirst [%];Entries [#]", 2000 , 0., 0.01);
                _hTHistDeltaPitch = trkInfo.make<TH1D>( "TrkDeltaPitch","TrkDeltaPitch;deltaPitch/pitchFirst [%];Entries [#]", 2000 , 0., 0.01);

                _hTHistDeltaPfirst = trkInfo.make<TH1D>( "TrkPfirst", "momFirst; momentum [MeV];Entries [#]", 2000 , 50., 110.0);
                _hTHistDeltaPlast  = trkInfo.make<TH1D>( "TrkPlast","momLast ; momentum [MeV];Entries [#]", 2000 , 50., 110.0);

                _hTHistDeltaEnergyClu   = tfs->make<TH1D>( "DeltaEnergyClu", "DeltaEnergyClu; E_{seed} - Energy_{cluster} [MeV];Entries [#]", 320,  -50.0, 110.);
                _hTHistDeltaEnergyCluRec= tfs->make<TH1D>( "DeltaEnergyCluRec","DeltaEnergyCluRec ; E_{seed} - Energy_{cluster} [MeV];Entries [#]", 320, -50.0, 110.);

                _hTHistThetaZ         = impAngles.make<TH1D>( "CosThetaZ", "CosThetaZ; cos#theta [MeV];Entries [#]", 720 , -180., 180.0);
                _hTHistThetaZRec      = impAngles.make<TH1D>( "CosThetaZRec", "CosThetaZRec; cos#theta [MeV];Entries [#]",720 , -180., 180.0);

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

                _hTHistEnergyCluMultiDIOdio = multClu.make<TH1D>( "EnergyCluMultiDIOdio","EnergyCluMultiDIOdio ; E_{cluster} [MeV];Entries [#]", 1000 , 0., 250.0);
                _hTHistEnergyMultiDIOdio    = multClu.make<TH1D>( "EnergyMultiDIOdio","EnergyMultiDIOdio ; E_{seed} ;Entries [#]", 440 , 0., 110.0);
                _hTHistTimeCluMultiDIOdio   = multClu.make<TH1D>( "TimeCluMultiDIOdio","TimeCluMultiDIOdio ; Time_{cluster} [ns] ;Entries [#]", 300 , -0., 3000.0);
                _hTHistTimeMultiDIOdio      = multClu.make<TH1D>( "TimeMultiDIOdio","TimeMultiDIOdio ; Time_{seed} [ns] ;Entries [#]", 300 , 0., 3000.0);

                _hTHistEnergyCluMultiDIOsignal = multClu.make<TH1D>( "EnergyCluMultiDIOsignal","EnergyCluMultiDIOsignal ; E_{cluster} [MeV];Entries [#]", 1000 , 0., 250.0);
                _hTHistEnergyMultiDIOsignal    = multClu.make<TH1D>( "EnergyMultiDIOsignal","EnergyMultiDIOsignal ; E_{seed} ;Entries [#]", 440 , 0., 110.0);
                _hTHistTimeCluMultiDIOsignal   = multClu.make<TH1D>( "TimeCluMultiDIOsignal","TimeCluMultiDIOsignal ; Time_{cluster} [ns] ;Entries [#]", 300 , -0., 3000.0);
                _hTHistTimeMultiDIOsignal      = multClu.make<TH1D>( "TimeMultiDIOsignal","TimeMultiDIOsignal ; Time_{seed} [ns] ;Entries [#]", 300 , 0., 3000.0);

                _hTHistEnergyMultiDIO          = multClu.make<TH1D>( "EnergyMultiDIO","EnergyMultiDIO ; E_{seed} ;Entries [#]", 440 , 0., 110.0);
                _hTHistTimeMultiDIO            = multClu.make<TH1D>( "TimeMultiDIO","TimeMultiDIO ; Time_{seed} [ns] ;Entries [#]", 300 , 0., 3000.0);

                _hTHistEnergyMultiSignal       = multClu.make<TH1D>( "EnergyMultiSignal","EnergyMultiSignal ; E_{seed} ;Entries [#]", 440 , 0., 110.0);
                _hTHistTimeMultiSignal         = multClu.make<TH1D>( "TimeMultiSignal","TimeMultiSignal ; Time_{seed} [ns] ;Entries [#]", 300 , 0., 3000.0);

                _hTHistMomDotVaneNorm          = impAngles.make<TH1D>( "MomDotVaneNorm","MomDotVaneNorm ; (#vec{p}, #hat{n}) / |#vec{p}| ;Entries [#]", 720 , -180., 180.0);
                _hTHistMomRecDotVaneNorm       = impAngles.make<TH1D>( "MomRecDotVaneNorm","MomRecDotVaneNorm ; (#vec{p}, #hat{n}) / |#vec{p}| ;Entries [#]", 720 , -180., 180.0);

                _hTHistThetaW                  = impAngles.make<TH1D>( "ThetaW","ThetaW ; arcTan(#vec{p_{u}}/(#vec{p}, #hat{n}) ) [deg];Entries [#]", 720 , -180., 180.0);
                _hTHistThetaW_v0               = impAngles.make<TH1D>( "ThetaW_v0","ThetaW_v0 ;  arcTan(#vec{p_{u}}/(#vec{p}, #hat{n}) ) [deg];Entries [#]", 720 , -180., 180.0);
                _hTHistThetaW_v1               = impAngles.make<TH1D>( "ThetaW_v1","ThetaW_v1 ; arcTan(#vec{p_{u}}/(#vec{p}, #hat{n}) ) [deg];Entries [#]", 720 , -180., 180.0);
                _hTHistThetaW_v2               = impAngles.make<TH1D>( "ThetaW_v2","ThetaW_v2 ;  arcTan(#vec{p_{u}}/(#vec{p}, #hat{n}) ) [deg];Entries [#]", 720 , -180., 180.0);
                _hTHistThetaW_v3               = impAngles.make<TH1D>( "ThetaW_v3","ThetaW_v3 ;  arcTan(#vec{p_{u}}/(#vec{p}, #hat{n}) ) [deg];Entries [#]", 720 , -180., 180.0);
                _hTHistThetaV                  = impAngles.make<TH1D>( "ThetaV","ThetaV ; arcTan(#vec{p_{v}}/(#vec{p}, #hat{n}) ) [deg];Entries [#]", 720 , -180., 180.0);
                _hTHistThetaV_v0               = impAngles.make<TH1D>( "ThetaV_v0","ThetaV_v0 ; arcTan(#vec{p_{v}}/(#vec{p}, #hat{n}) ) [deg];Entries [#]", 720 , -180., 180.0);
                _hTHistThetaV_v1               = impAngles.make<TH1D>( "ThetaV_v1","ThetaV_v1 ;arcTan(#vec{p_{v}}/(#vec{p}, #hat{n}) ) [deg];Entries [#]", 720 , -180., 180.0);
                _hTHistThetaV_v2               = impAngles.make<TH1D>( "ThetaV_v2","ThetaV_v2 ;arcTan(#vec{p_{v}}/(#vec{p}, #hat{n}) ) [deg];Entries [#]", 720 , -180., 180.0);
                _hTHistThetaV_v3               = impAngles.make<TH1D>( "ThetaV_v3","ThetaV_v3 ; arcTan(#vec{p_{v}}/(#vec{p}, #hat{n}) ) [deg];Entries [#]", 720 , -180., 180.0);
                _hTHistCryVsStepPointTime      = tfs->make<TH1D>( "CryVsStepPointTime","CryVsStepPointTime ; t_{Cry} - t_{Step} [ns];Entries", 200000 , -1000.0, 1000.0);

                _Ntup        = tfs->make<TTree>("ClusterSeedMap", "Cluster seeds info");

                _Ntup->Branch("evt", &_evt , "evt/F");//1
                _Ntup->Branch("clNo",&_clNo , "clNo/I");//2
                _Ntup->Branch("clE",&_clE , "clE/F");//3
                _Ntup->Branch("clT",&_clT , "clT/F");//4
                _Ntup->Branch("clCOGu",&_clCOGu , "clCOGu/F");//5
                _Ntup->Branch("clCOGv",&_clCOGv , "clCOGv/F");//6
                _Ntup->Branch("clCOGw",&_clCOGw , "clCOGw/F");//7
                _Ntup->Branch("nCryCl",&_nCryCl , "nCryCl/I");//8
                _Ntup->Branch("clVane",&_clVane , "clVane/I");//9
                _Ntup->Branch("clSeeds",&_clSeeds , "clSeeds/I");//10

                _Ntup->Branch("clSeedTrakId[clSeeds]",_clSeedTrackId , "clSeedTrakId[clSeeds]/I");//11
                _Ntup->Branch("clSeedPdgId[clSeeds]", _clSeedPdgId , "clSeedPdgId[clSeeds]/I");//12
                _Ntup->Branch("clSeedIsGen[clSeeds]", _clSeedIsGen , "clSeedIsGen[clSeeds]/I");
                _Ntup->Branch("clSeedIsSignal[clSeeds]", _clSeedIsSignal , "clSeedIsSignal[clSeeds]/I");
                _Ntup->Branch("clSeedTime[clSeeds]",_clSeedTime , "clSeedTime[clSeeds]/F");//13
                _Ntup->Branch("clSeedTotEnergyDep[clSeeds]",_clSeedTotEnergyDep , "clSeedTotEnergyDep[clSeeds]/F");//14
                _Ntup->Branch("clSeedPx[clSeeds]",_clSeedPx , "clSeedPx[clSeeds]/F");//15
                _Ntup->Branch("clSeedPy[clSeeds]",_clSeedPy , "clSeedPy[clSeeds]/F");//16
                _Ntup->Branch("clSeedPz[clSeeds]",_clSeedPz , "clSeedPz[clSeeds]/F");//17
                _Ntup->Branch("clSeedPu[clSeeds]",_clSeedPu , "clSeedPu[clSeeds]/F");//18
                _Ntup->Branch("clSeedPv[clSeeds]",_clSeedPv , "clSeedPv[clSeeds]/F");//29
                _Ntup->Branch("clSeedPw[clSeeds]",_clSeedPw , "clSeedPw[clSeeds]/F");//20
                _Ntup->Branch("clSeedCryFramePu[clSeeds]",_clSeedCryFramePu , "clSeedCryFramePu[clSeeds]/F");//21
                _Ntup->Branch("clSeedCryFramePv[clSeeds]",_clSeedCryFramePv , "clSeedCryFramePv[clSeeds]/F");//22
                _Ntup->Branch("clSeedCryFramePw[clSeeds]",_clSeedCryFramePw , "clSeedCryFramePw[clSeeds]/F");//23
                _Ntup->Branch("clSeedVaneFramePu[clSeeds]",_clSeedVaneFramePu , "clSeedVaneFramePu[clSeeds]/F");//24
                _Ntup->Branch("clSeedVaneFramePv[clSeeds]",_clSeedVaneFramePv , "clSeedVaneFramePv[clSeeds]/F");//25
                _Ntup->Branch("clSeedVaneFramePw[clSeeds]",_clSeedVaneFramePw , "clSeedVaneFramePw[clSeeds]/F");//26
                _Ntup->Branch("clSeedPpx[clSeeds]",_clSeedPpx , "clSeedPpx[clSeeds]/F");//27
                _Ntup->Branch("clSeedPpy[clSeeds]",_clSeedPpy , "clSeedPpy[clSeeds]/F");//28
                _Ntup->Branch("clSeedPpz[clSeeds]",_clSeedPpz , "clSeedPpz[clSeeds]/F");//29
                _Ntup->Branch("clSeedPpu[clSeeds]",_clSeedPpu , "clSeedPpu[clSeeds]/F");//30
                _Ntup->Branch("clSeedPpv[clSeeds]",_clSeedPpv , "clSeedPpv[clSeeds]/F");//31
                _Ntup->Branch("clSeedPpw[clSeeds]",_clSeedPpw , "clSeedPpw[clSeeds]/F");//32
                _Ntup->Branch("clSeedThetaW[clSeeds]",_clSeedThetaW , "clSeedThetaW[clSeeds]/F");//33
                _Ntup->Branch("clSeedThetaV[clSeeds]",_clSeedThetaV , "clSeedThetaV[clSeeds]/F");//34

        }

        doCalorimeter(evt, _skipEvent);

} // end of analyze


double logNormale(double *x, double *par){
        double y = 0.0;

        //if par[4]=-1 the logNormal function is defined with the higher tail on the left side
        double tempx = (-par[1]+x[0]);

        if(tempx > 0){
                y = TMath::Gaus(TMath::Log(tempx), par[2], par[3], kTRUE) / tempx ;
                y += par[4];
        }
        y*=par[0];
        y *= par[5];
        return y;
}

void CaloClusterEnergyResolMap::endJob() {
        art::ServiceHandle<GeometryService> geom;
        GeomHandle<VaneCalorimeter> cg;
        const int n = cg->nCrystalR();

        float x[n] ;
        float y[n] ;
        float ex[n];
        float ey[n];

        art::ServiceHandle<art::TFileService> tfs;
        art::TFileDirectory dirRow = tfs->mkdir("DirRow");
        art::TFileDirectory dirColumn = tfs->mkdir("DirColumn");

        _sigmaVsRow            = dirRow.make<TH1D>( "sigmaVsRow", "sigmaVsRow; row index; #sigma [MeV]", cg->nCrystalR(), 0.,cg->nCrystalR() );
        _sigmaVsColumn         = dirColumn.make<TH1D>( "sigmaVsColumn","sigmaVsColumn ; column index;  #sigma [MeV]", cg->nCrystalZ(), 0. ,cg->nCrystalZ() );
        _EpeakVsRow            = dirRow.make<TH1D>( "EpeakVsRow", "EpeakVsRow; row index; E_{peak} [MeV]", cg->nCrystalR(), 0.,cg->nCrystalR() );
        _EpeakVsColumn         = dirColumn.make<TH1D>( "EpeakVsColum","EpeakVsColum ; column index;  E_{peak} [MeV]", cg->nCrystalZ(), 0. ,cg->nCrystalZ() );
        _cTCanvSigmaVsRow      = dirRow.make<TCanvas>( "SigmaVsRow", "SigmaVsRow; row index; #sigma [MeV]" );


        TF1 *fitLog = new TF1("fitLog", logNormale, -50., 110., 6);
        fitLog->SetParameters(500., 0.1, 1., 2.);



        for( int i = 1; i<= _hTHistDeltaEnergy_row->GetNbinsX(); ++i){

                _hTHistDeltaEnergyRowProj= dirRow.make<TH1D>(Form("row_%d",i),Form("row_%d",i),_hTHistDeltaEnergy_row->GetNbinsY(),-50,110);

                _hTHistDeltaEnergyRowProj = _hTHistDeltaEnergy_row->ProjectionY(Form("row_%d",i), i, i);
                fitLog->SetParameters(_hTHistDeltaEnergyRowProj->GetEntries(), 0.1, 1., 2., 0., _binWidth);
                fitLog->FixParameter(5, _binWidth);
                fitLog->SetParLimits(4, 0., _hTHistDeltaEnergyRowProj->GetMaximum() / _binWidth );
                fitLog->SetParLimits(3, 0., 6.);

                _hTHistDeltaEnergyRowProj->Fit("fitLog", "L");
                _hTHistDeltaEnergyRowProj->Fit("fitLog", "L");
                _hTHistDeltaEnergyRowProj->Fit("fitLog", "L");
                _sigmaVsRow->SetBinContent(i, fitLog->GetParameter(3) );
                _sigmaVsRow->SetBinError(i, fitLog->GetParError(3) );
                double sigma = fitLog->GetParameter(3);
                double mean = fitLog->GetParameter(2);
                double shift = fitLog->GetParameter(1);

                double Dshift = fitLog->GetParError(1);
                double Dsigma = fitLog->GetParError(3);
                double Dmean = fitLog->GetParError(2);

                double Ep = exp(-pow(sigma, 2) )*(exp(mean) +  shift*exp(pow(sigma, 2) ) );
                double deltaEp = Dshift + exp(mean)*exp(-pow(sigma, 2) )*(Dmean + 2.0*sigma*Dsigma);

                _EpeakVsRow->SetBinContent(i, Ep);
                _EpeakVsRow->SetBinError(i, deltaEp);

                x[i] = i;
                y[i] = sigma;
                ex[i] = 0.;
                ey[i] = Dsigma;

        }
        _cTCanvSigmaVsRow->cd();
        _cTCanvSigmaVsRow->SetGrid();
        TGraphErrors *_TGr = new TGraphErrors(n,x,y,ex,ey);
        _TGr->SetTitle("Sigma Vs Row");
        _TGr->SetMarkerColor(4);
        _TGr->SetMarkerStyle(21);
        _TGr->Draw("AP");

        for( int i = 1; i<= _hTHistDeltaEnergy_column->GetNbinsX(); ++i){

                _hTHistDeltaEnergyColumnProj= dirColumn.make<TH1D>(Form("column_%d",i),Form("column_%d",i),_hTHistDeltaEnergy_column->GetNbinsY(),-50,110);

                _hTHistDeltaEnergyColumnProj = _hTHistDeltaEnergy_column->ProjectionY(Form("column_%d",i), i, i);
                fitLog->SetParameters(_hTHistDeltaEnergyColumnProj->GetEntries(), 0.1, 1., 2.);

                fitLog->SetParameters(_hTHistDeltaEnergyColumnProj->GetEntries(), 0.1, 1., 2., 0., _binWidth);
                fitLog->FixParameter(5, _binWidth);
                fitLog->SetParLimits(4, 0., _hTHistDeltaEnergyColumnProj->GetMaximum() / _binWidth );
                fitLog->SetParLimits(3, 0., 6.);


                _hTHistDeltaEnergyColumnProj->Fit("fitLog", "L");
                _hTHistDeltaEnergyColumnProj->Fit("fitLog", "L");
                _hTHistDeltaEnergyColumnProj->Fit("fitLog", "L");
                _sigmaVsColumn->SetBinContent(i, fitLog->GetParameter(3) );
                _sigmaVsColumn->SetBinError(i, fitLog->GetParError(3) );

                double sigma = fitLog->GetParameter(3);
                double mean = fitLog->GetParameter(2);
                double shift = fitLog->GetParameter(1);

                double Dshift = fitLog->GetParError(1);
                double Dsigma = fitLog->GetParError(3);
                double Dmean = fitLog->GetParError(2);

                double Ep = exp(-pow(sigma, 2) )*(exp(mean) +  shift*exp(pow(sigma, 2) ) );
                double deltaEp = Dshift + exp(mean)*exp(-pow(sigma, 2) )*(Dmean + 2.0*sigma*Dsigma);

                _EpeakVsColumn->SetBinContent(i, Ep);
                _EpeakVsColumn->SetBinError(i, deltaEp);

        }
}



void CaloClusterEnergyResolMap::doCalorimeter(art::Event const& evt, bool skip){


        if ( _diagLevel < 0 ) cout << "CaloClusterEnergyResolMap: begin" << endl;

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
                GenElHitData& ldil = iEltrk.getHit((int)(iEltrk.getNumOfHit() - 1) );

                if(!findTrkId(tmpVTot, iEltrk.getTrkID().asUint() ) ){

                        tmpVTot.push_back( iEltrk.getTrkID().asUint() );
                }

                double cosTheta = iEltrk.getTrkLrntzVec().cosTheta() ;

                double cosPitch = hdil._hitMomentum.cosTheta();

                double lcosPitch = ldil._hitMomentum.cosTheta();


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

        ElectronMap elecMap;

        if(caloClusters->size()>0 ){
                if ( _diagLevel < 0 ){
                        cout<<"-------------------- 3 -----------------------"<<endl;
                }
                int iVane =-1;
                for(size_t icl=0; icl<caloClusters->size(); ++icl){

                        double eDepClu = 0.;

                        _evt = evt.id().event();

                        CaloCluster const& clu = (*caloClusters).at(icl);

                        _clNo = icl;
                        _clE = clu.energyDep();
                        _clT = clu.time();

                        _clCOGu = clu.cog3Vector().x();
                        _clCOGv = clu.cog3Vector().y();
                        _clCOGw = clu.cog3Vector().z();
                        _nCryCl = clu.size();

                        eDepClu = clu.energyDep();
                        iVane = clu.vaneId();
                        _clVane = iVane;

                        CaloCrystalHitPtrVector caloClusterHits = clu.caloCrystalHitsPtrVector();
                        if(eDepClu >= _EnergyClusterCut){

                                if(_diagLevel < 0){
                                        cout <<"eDepClu >energyDepClusterCut"<< endl;
                                }

                                std::map<unsigned int, unsigned int> seedMap, signalMap;
                                CryMap cryMap;

                                for(size_t i=0; i<caloClusterHits.size(); ++i){



                                        if ( _diagLevel < 0 ){
                                                cout<<"-------------------- 3."<<i <<" -----------------------"<<endl;
                                        }
                                        CaloCrystalHit const& hit = *(caloClusterHits.at(i));


                                        double CryTime = hit.time();

                                        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                                        if(ROIds.size()<1 ) continue;

                                        CaloHit const& thehit = *ROIds.at(0);

                                        size_t collectionPosition = ROIds.at(0).key();

                                        PtrStepPointMCVector const & mcptr(hits_mcptr->at(collectionPosition));

                                        size_t nHitsPerCrystal = mcptr.size();



//                                        float earliest = 100000;

                                        for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {

                                                if(_diagLevel < 0){
                                                        cout<<"-------------------- 3."<<i <<"."<< j2<<"-----------------------"<<endl;
                                                }

                                                StepPointMC const& mchit = *mcptr[j2];

                                                art::Ptr<SimParticle> const& simptr = mchit.simParticle();
                                                SimParticle const& sim = *simptr;

                                                // The simulated particle that made this hit.
                                                SimParticleCollection::key_type trackId(mchit.trackId());

                                                if(sim.fromGenerator() ){
                                                        GenParticle const& gen = *sim.genParticle();
                                                        GenId genId = gen.generatorId();

                                                        if(cryMap.find(mchit.trackId().asInt() )==cryMap.end() ){
                                                                cryMap[mchit.trackId().asInt()].time         = mchit.time();
                                                                cryMap[mchit.trackId().asInt()].totEnergyDep = mchit.totalEDep();
                                                                cryMap[mchit.trackId().asInt()].pdgId        = sim.pdgId();
                                                                cryMap[mchit.trackId().asInt()].isGen        = sim.fromGenerator();
                                                                cryMap[mchit.trackId().asInt()].position     = mchit.position();
                                                                cryMap[mchit.trackId().asInt()].momentum     = mchit.momentum();
                                                                CLHEP::Hep3Vector crystalFrame = cg->toCrystalFrame(thehit.id(), mchit.position());
                                                                cryMap[mchit.trackId().asInt()].cryPosition  = crystalFrame;
                                                                cryMap[mchit.trackId().asInt()].vane         = iVane;
                                                                cryMap[mchit.trackId().asInt()].isSignal = genId==GenId::conversionGun;

                                                        }else if(mchit.time() < cryMap[mchit.trackId().asInt()].time ){
                                                                cryMap[mchit.trackId().asInt()].time         = mchit.time();
                                                                cryMap[mchit.trackId().asInt()].totEnergyDep = mchit.totalEDep();
                                                                cryMap[mchit.trackId().asInt()].pdgId        = sim.pdgId();
                                                                cryMap[mchit.trackId().asInt()].isGen        = sim.fromGenerator();
                                                                cryMap[mchit.trackId().asInt()].position     = mchit.position();
                                                                cryMap[mchit.trackId().asInt()].momentum     = mchit.momentum();
                                                                CLHEP::Hep3Vector crystalFrame = cg->toCrystalFrame(thehit.id(), mchit.position());
                                                                cryMap[mchit.trackId().asInt()].cryPosition  = crystalFrame;
                                                                cryMap[mchit.trackId().asInt()].vane         = iVane;
                                                                cryMap[mchit.trackId().asInt()].isSignal = genId==GenId::conversionGun;


                                                        }
                                                }


                                                _hTHistCryVsStepPointTime->Fill(CryTime- mchit.time());

                                                CLHEP::Hep3Vector cryFrame = cg->toCrystalFrame(thehit.id(), mchit.position());

                                                if( (cg->crystalRByRO(thehit.id()) !=(_rowToCanc - 1) && cg->crystalZByRO(thehit.id()) != (_columnToCanc-1) ) ){
                                                        unsigned int row = cg->crystalRByRO(thehit.id() );
                                                        unsigned int column = cg->crystalZByRO(thehit.id() );

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
                                                                elecMap[iVane][trackId.asUint()]._cluEnergy = clu.energyDep();
                                                                elecMap[iVane][trackId.asUint()]._cluTime = clu.time();
                                                                elecMap[iVane][trackId.asUint()]._impTime = mchit.time();
                                                                elecMap[iVane][trackId.asUint()]._impEnergy = mchit.momentum().mag();
                                                                elecMap[iVane][trackId.asUint()]._cluCog = clu.cog3Vector();
                                                                elecMap[iVane][trackId.asUint()]._impPos = mchit.position();
                                                                elecMap[iVane][trackId.asUint()]._impPosCryFrame = cryFrame;
                                                                elecMap[iVane][trackId.asUint()]._row    = row;
                                                                elecMap[iVane][trackId.asUint()]._column = column;
                                                                elecMap[iVane][trackId.asUint()]._cryOrigin  = cg->crystalOriginByRO( thehit.id() );
                                                                elecMap[iVane][trackId.asUint()]._impMom3Vec = mchit.momentum();
                                                                elecMap[iVane][trackId.asUint()]._cryTime = CryTime;
                                                                if(_diagLevel < 0){
                                                                        cout<< "###################"<<endl;
                                                                        cout<< "idVane = "<< iVane<<endl;
                                                                        cout << "cluX = "<<elecMap[iVane][trackId.asUint()]._impPos.getX()<<endl;
                                                                        cout << "cluY = "<<elecMap[iVane][trackId.asUint()]._impPos.getY()<<endl;
                                                                        cout << "cluZ = "<<elecMap[iVane][trackId.asUint()]._impPos.getZ()<<endl;
                                                                        cout<< "###################"<<endl;
                                                                }
                                                        }

                                                }
                                        }//end loop for(nHitsPerCrystal)


                                }//end loop for(caloClusterHits.size())
                                //--------------------------------------------
                                _clSeeds = cryMap.size();
                                //_clVane = -1;
                                if(_clSeeds > 0){
                                        int j=0, tmpVane=-1;
                                        for(CryMap::iterator it = cryMap.begin(); it != cryMap.end(); ++it){
                                                _clSeedTrackId[j]      = it->first;
                                                tmpVane                = it->second.vane;
                                                _clSeedPdgId[j]        = it->second.pdgId;
                                                _clSeedIsGen[j]        = it->second.isGen;
                                                _clSeedIsSignal[j]     = it->second.isSignal;
                                                _clSeedTime[j]         = it->second.time ;
                                                _clSeedTotEnergyDep[j] = it->second.totEnergyDep;
                                                _clSeedPx[j]           = it->second.position.x();
                                                _clSeedPy[j]           = it->second.position.y();
                                                _clSeedPz[j]           = it->second.position.z();
//                                                CLHEP::Hep3Vector crystalFrame = cg->toCrystalFrame(thehit.id(), it->second.position);
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
                                                Vane const &vane = cg->vane(tmpVane);
                                                CLHEP::Hep3Vector Mom_rotated = (vane.rotation())*(it->second.momentum);
                                                _clSeedPpu[j]          = Mom_rotated.x();
                                                _clSeedPpv[j]          = Mom_rotated.y();
                                                _clSeedPpw[j]          = Mom_rotated.z();
                                                double thetaWimpact = std::atan(-1.0*Mom_rotated.getZ() / Mom_rotated.getX() ) ;
                                                _clSeedThetaW[j]       = thetaWimpact*180./TMath::Pi();
                                                double thetaVimpact = std::atan(Mom_rotated.getY() /  Mom_rotated.getX() ) ;
                                                _clSeedThetaV[j]       = thetaVimpact*180./TMath::Pi();
                                                ++j;
                                        }

                                }else {
                                        _clSeeds = 1;
                                        _clSeedTrackId[0]      = -1;
                                        _clSeedPdgId[0]        = 0;
                                        _clSeedIsGen[0]        = 0.0;
                                        _clSeedIsSignal[0]     = 0;
                                        _clSeedTime[0]         = 0.0;
                                        _clSeedTotEnergyDep[0] = 0.0;
                                        _clSeedPx[0]           = 0.0;
                                        _clSeedPy[0]           = 0.0;
                                        _clSeedPz[0]           = 0.0;
                                        _clSeedCryFramePu[0]   = 0.0;
                                        _clSeedCryFramePv[0]   = 0.0;
                                        _clSeedCryFramePw[0]   = 0.0;
                                        _clSeedVaneFramePu[0]  = 0.0;
                                        _clSeedVaneFramePv[0]  = 0.0;
                                        _clSeedVaneFramePw[0]  = 0.0;
                                        _clSeedPpx[0]          = 0.0;
                                        _clSeedPpy[0]          = 0.0;
                                        _clSeedPpz[0]          = 0.0;
                                        _clSeedPpu[0]          = 0.0;
                                        _clSeedPpv[0]          = 0.0;
                                        _clSeedPpw[0]          = 0.0;
                                        _clSeedThetaW[0]       = 0.0;
                                        _clSeedThetaV[0]       = 0.0;

                                }



                                if(_diagLevel < 0){
                                        cout<<"-------------------- Filling _Ntup... -----------------------"<<endl;
                                }
                                _Ntup->Fill();
                                if(_diagLevel < 0){
                                        cout<<"-------------------- Filled _Ntup -----------------------"<<endl;
                                }

                                //---------------------------------------------
                                if(seedMap.size() > 1){
                                        if(signalMap.size() > 0){
                                                _hTHistEnergyCluMultiDIOsignal->Fill( eDepClu);
                                                _hTHistTimeCluMultiDIOsignal->Fill( clu.time());

                                                for(std::map<unsigned int, unsigned int>::iterator its = seedMap.begin(); its != seedMap.end(); ++its){
                                                        //std::map<unsigned int, unsigned int>::iterator its = seedMap.begin();
                                                        //its = its + i;
                                                        std::map<unsigned int, electronData >::iterator it = elecMap[iVane].find(its->first);
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
                                                _hTHistTimeCluMultiDIOdio->Fill( clu.time());


                                                for(std::map<unsigned int, unsigned int>::iterator its = seedMap.begin(); its != seedMap.end(); ++its){

                                                        std::map<unsigned int, electronData >::iterator it = elecMap[iVane].find(its->first);
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
        if(_diagLevel < 0){
                cout <<"-------------------> number of generated electrons ( trkVecTot.size() ) = "<< trkVecTot.size() <<endl;
        }
        unsigned int size2 = trkVecTot.size();
        unsigned int it2=0;
        while( it2 < size2){
                ElectronMap::iterator ite = elecMap.begin();
                bool trovato = false;
                while(!trovato && ite!=elecMap.end() ){
                        if(ite->second.find(trkVecTot[it2]) != ite->second.end()){

                                if(_diagLevel < 0){
                                        cout<< "$$$$$$$$$$$$----> vane = "<< ite->first<<endl;
                                        cout<<"------------------------------------------------------"<<endl;
                                        cout << "impX = "<<ite->second[trkVecTot[it2]]._impPos.getX()<<endl;
                                        cout << "impY = "<<ite->second[trkVecTot[it2]]._impPos.getY()<<endl;
                                        cout << "impZ = "<<ite->second[trkVecTot[it2]]._impPos.getZ()<<endl;
                                        cout << "impXcryFRame = "<<ite->second[trkVecTot[it2]]._impPosCryFrame.getX()<<endl;
                                        cout << "impYcryFRame = "<<ite->second[trkVecTot[it2]]._impPosCryFrame.getY()<<endl;
                                        cout << "impZcryFRame = "<<ite->second[trkVecTot[it2]]._impPosCryFrame.getZ()<<endl;
                                        cout<< "@@@@@@@@@@@@@@@@@@@"<<endl;
                                        cout << "cluX = "<<ite->second[trkVecTot[it2]]._cluCog.getX()<<endl;
                                        cout << "cluY = "<<ite->second[trkVecTot[it2]]._cluCog.getY()<<endl;
                                        cout << "cluZ = "<<ite->second[trkVecTot[it2]]._cluCog.getZ()<<endl;
                                        cout<< "@@@@@@@@@@@@@@@@@@@"<<endl;
                                        cout << "delta X = " << ite->second[trkVecTot[it2]]._impPos.getX() - ite->second[trkVecTot[it2]]._cluCog.getX()<<endl;
                                        cout << "delta Y = " << ite->second[trkVecTot[it2]]._impPos.getY() - ite->second[trkVecTot[it2]]._cluCog.getY()<<endl;
                                        cout << "delta Z = " << ite->second[trkVecTot[it2]]._impPos.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ()<<endl;
                                }

                                _hTHistDeltaEnergyCluRec->Fill( ite->second[trkVecTot[it2]]._impEnergy - ite->second[trkVecTot[it2]]._cluEnergy);

                                _hTHistDeltaEnergyRec_row->Fill(_hTHistDeltaEnergyRec_row->GetXaxis()->GetBinCenter(ite->second[trkVecTot[it2]]._row + 1), ite->second[trkVecTot[it2]]._impEnergy - ite->second[trkVecTot[it2]]._cluEnergy);
                                _hTHistDeltaEnergyRec_column->Fill(_hTHistDeltaEnergyRec_column->GetXaxis()->GetBinCenter(ite->second[trkVecTot[it2]]._column+ 1), ite->second[trkVecTot[it2]]._impEnergy - ite->second[trkVecTot[it2]]._cluEnergy);

                                _hTHistThetaZRec->Fill(ite->second[trkVecTot[it2]]._impMom3Vec.getTheta()*180./TMath::Pi());
                                _hTHistPhiRec->Fill(ite->second[trkVecTot[it2]]._impMom3Vec.getPhi()*180./TMath::Pi());

                                Vane const &vane = cg->vane(ite->first);
                                CLHEP::Hep3Vector dirMom = ite->second[trkVecTot[it2]]._impMom3Vec.unit();
                                CLHEP::Hep3Vector dirMom_rotated = (vane.rotation())*dirMom;
                                _hTHistMomRecDotVaneNorm->Fill(acos(-1.0*dirMom_rotated.getX())*180./TMath::Pi() );

                                if(ite->first == 0 ){
                                        _hTHistPhiRec_v0->Fill(ite->second[trkVecTot[it2]]._impMom3Vec.getPhi()*180./TMath::Pi());
                                        _hTHistMomRecDotVaneNorm->Fill(acos(-1.0*ite->second[trkVecTot[it2]]._impMom3Vec.getY() /ite->second[trkVecTot[it2]]._impMom3Vec.mag() )*180./TMath::Pi() );
                                }
                                if(ite->first == 1 ){
                                        _hTHistPhiRec_v1->Fill(ite->second[trkVecTot[it2]]._impMom3Vec.getPhi()*180./TMath::Pi());
                                        _hTHistMomRecDotVaneNorm->Fill(acos(ite->second[trkVecTot[it2]]._impMom3Vec.getX() /ite->second[trkVecTot[it2]]._impMom3Vec.mag() )*180./TMath::Pi() );
                                }
                                if(ite->first == 2 ){
                                        _hTHistPhiRec_v2->Fill(ite->second[trkVecTot[it2]]._impMom3Vec.getPhi()*180./TMath::Pi());
                                        _hTHistMomRecDotVaneNorm->Fill( acos(ite->second[trkVecTot[it2]]._impMom3Vec.getY() /ite->second[trkVecTot[it2]]._impMom3Vec.mag() )*180./TMath::Pi() );
                                }
                                if(ite->first == 3 ){
                                        _hTHistPhiRec_v3->Fill(ite->second[trkVecTot[it2]]._impMom3Vec.getPhi()*180./TMath::Pi());
                                        _hTHistMomRecDotVaneNorm->Fill( acos(-1.0*ite->second[trkVecTot[it2]]._impMom3Vec.getX() /ite->second[trkVecTot[it2]]._impMom3Vec.mag() )*180./TMath::Pi() );
                                }
                                trovato = true;
                                canc2 = it2;
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

        if(_diagLevel < 0){
                cout <<"after.... trkVecTot.size = "<< trkVecTot.size() <<endl;
                cout<<"...\n...\n "<<endl;
                cout <<"-------------------> number of quality electrons generated ( trkVec.size() ) = "<< trkVec.size() <<endl;
        }

        //counting how many reconstructed clusters match with the generated particles which pass the quality cuts of the tracker
        unsigned int canc = 0;

        //for(unsigned int it = 0; it != trkVec.size(); ++it){
        unsigned int size = trkVec.size();
        unsigned int it=0;
        while( it < size){
                ElectronMap::iterator ite = elecMap.begin();
                bool trovato = false;
                while(!trovato && ite!=elecMap.end() ){
                        if(ite->second.find(trkVec[it]) != ite->second.end()){
                                if(_diagLevel < 0){
                                        cout << "impX = "<<ite->second[trkVec[it]]._impPosCryFrame.getX()<<endl;
                                        cout << "impY = "<<ite->second[trkVec[it]]._impPosCryFrame.getY()<<endl;
                                        cout << "impZ = "<<ite->second[trkVec[it]]._impPosCryFrame.getZ()<<endl;
                                }

                                _hTHistDeltaEnergyClu->Fill( ite->second[trkVec[it]]._impEnergy - ite->second[trkVec[it]]._cluEnergy);

                                _hTHistDeltaEnergy_row->Fill(_hTHistDeltaEnergy_row->GetXaxis()->GetBinCenter(ite->second[trkVec[it]]._row + 1), ite->second[trkVec[it]]._impEnergy - ite->second[trkVec[it]]._cluEnergy);
                                _hTHistDeltaEnergy_column->Fill(_hTHistDeltaEnergy_column->GetXaxis()->GetBinCenter(ite->second[trkVec[it]]._column + 1), ite->second[trkVec[it]]._impEnergy - ite->second[trkVec[it]]._cluEnergy);

                                _hTHistThetaZ->Fill(ite->second[trkVecTot[it]]._impMom3Vec.getTheta()*180./TMath::Pi());
                                _hTHistPhi->Fill(ite->second[trkVecTot[it]]._impMom3Vec.getPhi()*180./TMath::Pi());


                                CLHEP::Hep3Vector dirMom = ite->second[trkVecTot[it2]]._impMom3Vec;//.unit();

                                Vane const &vane = cg->vane(ite->first);
                                dirMom = (vane.rotation())*dirMom;

                                _hTHistMomDotVaneNorm->Fill(acos(-1.0*dirMom.getX() / dirMom.mag() )*180./TMath::Pi() );
                                _hTHistThetaW->Fill(std::atan(-1.0*dirMom.getZ() / dirMom.getX() )*180./TMath::Pi() );

                                _hTHistThetaV->Fill(std::atan(dirMom.getY() /  dirMom.getX() )*180./TMath::Pi() );

                                if(ite->first == 0 ){
                                        _hTHistPhi_v0->Fill(ite->second[trkVec[it]]._impMom3Vec.getPhi()*180./TMath::Pi());

                                        _hTHistThetaW_v0->Fill(std::atan(-1.0*ite->second[trkVec[it]]._impMom3Vec.getZ() /  ite->second[trkVec[it]]._impMom3Vec.getY()  )*180./TMath::Pi() );
                                        if(_diagLevel < 0){
                                                cout<<"vaneIndex = "<<ite->first <<endl;
                                                cout<<"Vx = "<<ite->second[trkVec[it]]._impMom3Vec.getX()<<", Vy = "<<ite->second[trkVec[it]]._impMom3Vec.getY() <<", Vz = "<< ite->second[trkVec[it]]._impMom3Vec.getZ()<<endl;
                                        }
                                        _hTHistThetaV_v0->Fill(std::atan(ite->second[trkVec[it]]._impMom3Vec.getX() /  ite->second[trkVec[it]]._impMom3Vec.getY() )*180./TMath::Pi() );

                                }
                                if(ite->first == 1 ){
                                        _hTHistPhi_v1->Fill(ite->second[trkVec[it]]._impMom3Vec.getPhi()*180./TMath::Pi());

                                        _hTHistThetaW_v1->Fill(std::atan(ite->second[trkVec[it]]._impMom3Vec.getZ() /  ite->second[trkVec[it]]._impMom3Vec.getX()  )*180./TMath::Pi());
                                        if(_diagLevel < 0){
                                                cout<<"vaneIndex = "<<ite->first <<endl;
                                                cout<<"Vx = "<<ite->second[trkVec[it]]._impMom3Vec.getX()<<", Vy = "<<ite->second[trkVec[it]]._impMom3Vec.getY() <<", Vz = "<< ite->second[trkVec[it]]._impMom3Vec.getZ()<<endl;
                                        }
                                        _hTHistThetaV_v1->Fill(std::atan(ite->second[trkVec[it]]._impMom3Vec.getY() / ite->second[trkVec[it]]._impMom3Vec.getX()  )*180./TMath::Pi());

                                }
                                if(ite->first == 2 ){
                                        _hTHistPhi_v2->Fill(ite->second[trkVec[it]]._impMom3Vec.getPhi()*180./TMath::Pi());

                                        _hTHistThetaW_v2->Fill(std::atan(ite->second[trkVec[it]]._impMom3Vec.getZ() / ite->second[trkVec[it]]._impMom3Vec.getY() )*180./TMath::Pi() );
                                        if(_diagLevel < 0){
                                                cout<<"vaneIndex = "<<ite->first <<endl;
                                                cout<<"Vx = "<<ite->second[trkVec[it]]._impMom3Vec.getX()<<", Vy = "<<ite->second[trkVec[it]]._impMom3Vec.getY() <<", Vz = "<< ite->second[trkVec[it]]._impMom3Vec.getZ()<<endl;
                                        }
                                        _hTHistThetaV_v2->Fill(std::atan(-1.0*ite->second[trkVec[it]]._impMom3Vec.getX() /  ite->second[trkVec[it]]._impMom3Vec.getY()  )*180./TMath::Pi() );


                                }
                                if(ite->first == 3 ){
                                        _hTHistPhi_v3->Fill(ite->second[trkVec[it]]._impMom3Vec.getPhi()*180./TMath::Pi());

                                        _hTHistThetaW_v3->Fill(std::atan(-1.0*ite->second[trkVec[it]]._impMom3Vec.getZ() / ite->second[trkVec[it]]._impMom3Vec.getX()  )*180./TMath::Pi());
                                        if(_diagLevel < 0){
                                                cout<<"vaneIndex = "<<ite->first <<endl;
                                                cout<<"Vx = "<<ite->second[trkVec[it]]._impMom3Vec.getX()<<", Vy = "<<ite->second[trkVec[it]]._impMom3Vec.getY() <<", Vz = "<< ite->second[trkVec[it]]._impMom3Vec.getZ()<<endl;
                                        }
                                        _hTHistThetaV_v3->Fill(std::atan(ite->second[trkVec[it]]._impMom3Vec.getY() /  ite->second[trkVec[it]]._impMom3Vec.getX()  )*180./TMath::Pi());
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

        cout << "Event "<<evt.id().event()<<" CaloClusterEnergyResolMap done..."<<endl;


}

}


using mu2e::CaloClusterEnergyResolMap;
DEFINE_ART_MODULE(CaloClusterEnergyResolMap);


