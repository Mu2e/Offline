//
//
//
// $Id: CaloClusterLogCog_module.cc,v 1.6 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
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
#include "Mu2eUtilities/inc/LinePointPCA.hh"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//calorimeter packages
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
//#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "CaloCluster/inc/CaloClusterer.hh"
//#include "CaloCluster/inc/CaloClusterFinder.hh"
//#include "CaloCluster/inc/ClosestCaloClusterFinder.hh"

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"


// Other includes.
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib/exception.h"

//#include "Analyses/inc/MCCaloUtilities.hh"

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

//root includes
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

// From the art tool-chain
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "TMath.h"

//#include "CaloCluster/inc/CaloClusterUtilities.hh"

using namespace std;

namespace mu2e {


struct electronData{
        double _cluEnergy;
        double _cluTime;
        int    _cluSize;
        double _cryEnergyDep;
        double _cryEnergyDepTotal;
        double _impTime;
        double _impEnergy;
        CLHEP::Hep3Vector _cluCog;
        CLHEP::Hep3Vector _impPos;
        CLHEP::Hep3Vector _impPosCryFrame;
        unsigned int _vane;
        unsigned int _row;
        unsigned int _colum;
        CLHEP::Hep3Vector _impMom3Vec;
        int _impPdgId;
        int _impIsGen;
        ClusterMap _clusterMap;


        bool operator<( const electronData other) const{
                return ( _impTime< other._impTime);
        }
        electronData & operator=(const electronData& other) {
                _cluEnergy = other._cluEnergy;
                _cluTime   = other._cluTime;
                _cluSize   = other._cluSize;
                _impTime = other._impTime;
                _impEnergy = other._impEnergy;
                _cluCog    = other._cluCog;
                _impPos = other._impPos;
                _impPosCryFrame = other._impPosCryFrame;
                _vane     = other._vane;
                _row      = other._row;
                _colum    = other._colum;
                _impMom3Vec =other._impMom3Vec;
                _impPdgId  = other._impPdgId;
                _clusterMap = other._clusterMap;


                return *this;
        }
        electronData():
                _impTime(1e10),
                _impEnergy(0.0){
        }
};

//the key is the the vane
typedef std::map<unsigned int,std::map<unsigned int, electronData > > ElectronMap;


static int ncalls(0);

class CaloClusterLogCog : public art::EDAnalyzer {
public:
        explicit CaloClusterLogCog(fhicl::ParameterSet const& pset):
        _diagLevel(pset.get<int>("diagLevel",0)),
        _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
        _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
        _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
        _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
        _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel", "makeCaloCluster")),
        _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
        _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "time")),
        _producerName("Algo"+mu2e::TOUpper(_caloClusterAlgorithm)+"SeededBy"+mu2e::TOUpper(_caloClusterSeeding)),
        _elextractModuleLabel(pset.get<std::string>("elextractModuleLabel", "extractElData")),
        _extractElectronsData(pset.get<string>("elextractModuleLabel")),
        _rowToCanc(pset.get<int>("rowToCancel",-1)),//row index running from 1, 2, ... nCryR
        _columnToCanc(pset.get<int>("columnToCancel",-1)),//column index running from 1, 2, ... nCryZ
        _depthPitch(pset.get<double>("depthPitch",0.5)),//[mm]
        _thetaVimpact(pset.get<double>("thetaVimpact",0.785)),//45 [grad]
        _thetaWimpact(pset.get<double>("thetaWimpact",0.785)),//45 [grad]
        _nAnalyzed(0),
        //_hTHistGlobalEffNorm(0),
        _hTHistDeltaEnergy(0),
        _hTHistDeltaEnergVRec(0),
        //_hTHistGlobalEffRec(0),
        _hTHistDeltaUquality(0),
        _hTHistDeltaURec(0),
        _hTHistDeltaVquality(0),
        _hTHistDeltaVRec(0),
        _hTHistDeltaWquality(0),
        _hTHistDeltaWRec(0),
        _hTHistImpactParam(0),
        _hTHistDistanceToCog(0),
        _hTHistMomDotVaneNorm(0),
        _hTHistRecImpactParam(0),
        _hTHistRecDistanceToCog(0),
        _hTHistMomRecDotVaneNorm(0),
        _hTHistOffSet(0),
        _Ntup(0),
        _MaxCogOffSet(pset.get<double>("MaxCogOffset",5.0)),
        _MinCogOffSet(pset.get<double>("MinCogOffset",0.8)),
        _PitchCogOffSet(pset.get<double>("PitchCogOffset",0.2)),
        _MaxDepth(pset.get<double>("maxDepth",11.)),
        _application(0),
        _directory(0)
        {
        }
        virtual ~CaloClusterLogCog() {
        }
        virtual void beginJob();
        virtual void endJob();

        void analyze(art::Event const& e );

private:

        void doCalorimeter(art::Event const& evt, bool skip);

        // Diagnostic level
        int _diagLevel;

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

        double _minimumEnergy; //minimum energy deposition of hits
        int _rowToCanc, _columnToCanc;
        double _depthPitch;
        double _thetaVimpact;
        double _thetaWimpact;
        //number of analyzed events
        int _nAnalyzed;

        //TH1D* _hTHistGlobalEffNorm;
        TH1D* _hTHistDeltaEnergy;
        TH1D* _hTHistDeltaEnergVRec;
        //TH1D* _hTHistGlobalEffRec;

        TH2D* _hTHistDeltaUquality;
        TH2D* _hTHistDeltaURec;
        TH2D* _hTHistDeltaVquality;
        TH2D* _hTHistDeltaVRec;
        TH2D* _hTHistDeltaWquality;
        TH2D* _hTHistDeltaWRec;
        TH2D* _hTHistImpactParam;
        TH2D* _hTHistDistanceToCog;
        TH2D* _hTHistMomDotVaneNorm;
        TH2D*_hTHistRecImpactParam;
        TH2D*_hTHistRecDistanceToCog;
        TH2D*_hTHistMomRecDotVaneNorm;
        TH2D* _hTHistOffSet;
        TTree* _Ntup;



        double _MaxCogOffSet, _MinCogOffSet,_PitchCogOffSet;
        double _MaxDepth;

        double globalNtrkCut;
        double globalNtrkTot;
        double *globalCaloCut;
        double *globalRecCaloCut;

        std::auto_ptr<MCCaloUtilities> CaloManager;

        bool _skipEvent;

        Int_t _nDepthSteps,
        _seedId[10000],
        _seedPdgId[10000],
        _seedIsGen[10000],
        _clVane[10000],
        _clSize[10000],
        _clCogCrySize[10000];

        Float_t _evt,
        _cogOffSet,
        _clE[10000],
        _clT[10000],
        _clCOGx[10000],
        _clCOGy[10000],
        _clCOGz[10000],
        _clCOGu[10000],
        _clCOGv[10000],
        _clCOGw[10000],
        _clCOGrow[10000],
        _clCOGcolumn[10000],
        _clRows[10000][10000],
        _clColumns[10000][10000],
        _clCogRows[10000][10000],
        _clCogColumns[10000][10000],
        _cryEdep[10000],
        _cryEdepTot[10000],
        _DCAu[10000],
        _DCAv[10000],
        _DCAw[10000],
        _depth[10000],
        _seedPx[10000],
        _seedPy[10000],
        _seedPz[10000] ,
        _seedE[10000] ,
        _seedTime[10000],
        _seedPpx[10000] ,
        _seedPpy[10000] ,
        _seedPpz[10000] ,
        _seedPpu[10000] ,
        _seedPpv[10000] ,
        _seedPpw[10000] ,
        _seedRow[10000],
        _seedColumn[10000],
        _seedPu[10000] ,
        _seedPv[10000] ,
        _seedPw[10000] ,
        _seedPpCosTh[10000],
        _seedPpPhi[10000] ,
        _seedPpRotCosTh[10000],
        _seedPpRotPhi[10000] ,
        _seedThetaW[10000],
        _seedThetaV[10000],
        _seedDeltaW[10000],
        _seedDeltaV[10000],
        _seedDeltaCorrW[10000],
        _seedDeltaCorrV[10000],
        _distToCog[10000],
        _cogDca[10000],
        _prova1[10000],
        _prova2[10000];


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

double triangoloVar(double x, string d){
        double y =0.;
        double par[6] = {0.};
        if(d.compare("V")==1){
                double p1[6] = { 0.1689, 29.92, -1.344, 4.136, 0.6882, 0.004642};
                for(int y=0; y<6; ++y){
                        par[y] = p1[y];
                }
        } else if (d.compare("W")==1){
                double p2[6] = { 0.343, 30., -4.25, -1.441, 0.5152, 0.0001403};
                for(int y=0; y<6; ++y){
                        par[y] = p2[y];
                }
        }

        if(par[1] == 0.) return y;

        double m = par[0], xp = par[1], q = par[2], s = par[3], p1 = par[4], p0 = par[5];

        double amp = p0*x;
        amp = TMath::Exp(amp);

        double tempx = x - s;
        double temp;
        double p2 = 1.- p1;
        //if(tempx < 0.){
        //  temp = tempx/xp +1.;
        //}else {
        temp = tempx/xp;
        //}
        double cut = (temp - ( (int) temp)+( tempx <0 ? +1. : 0)  );

        if(cut < p1){
                y = m*xp*cut + q;
        } else {
                double mp = -p1*m/p2, qp = q - mp*xp;
                y= mp*xp*cut + qp ;
        }

        y *= amp;
        return y;

}



void CaloClusterLogCog::beginJob( ) {

        cout << "start CaloClusterLogCog..."<<endl;

        CaloManager = auto_ptr<MCCaloUtilities>(new MCCaloUtilities());

        // If needed, create the ROOT interactive environment. See note 1.
        if ( !gApplication ){
                int    tmp_argc(0);
                char** tmp_argv(0);
                _application = auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
        }

        gStyle->SetPalette(1);
        gROOT->SetStyle("Plain");

        _directory = gDirectory;


}

void CaloClusterLogCog::analyze(art::Event const & evt ) {

        ++_nAnalyzed;
        ++ncalls;

        CaloClusterer c;

        art::ServiceHandle<GeometryService> geom;
        if (ncalls == 1) {

                // cout << "This should be done only in the first event" << endl;

                art::ServiceHandle<art::TFileService> tfs;
                art::TFileDirectory cogResol = tfs->mkdir("CogResol");
                art::TFileDirectory energyResol = tfs->mkdir("EnergyResol");
                art::TFileDirectory dcaCog = tfs->mkdir("DcaCog");

                _hTHistDeltaEnergy   = energyResol.make<TH1D>( "DeltaEnergy", "DeltaEnergy;Eseed-Eclu [MeV]", 100, -20., 110.);
                _hTHistDeltaEnergVRec= energyResol.make<TH1D>( "DeltaEnergVRec", "DeltaEnergVRec;Eseed-Eclu [MeV]", 100, -20., 110.);

                _hTHistDeltaUquality = cogResol.make<TH2D>( "DeltaUquality", "DeltaUquality; depth_{MaxShower} [mm];Xseed-Xclu [mm]", (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth, 200, -100., 100. );
                _hTHistDeltaURec     = cogResol.make<TH2D>( "DeltaURec", "DeltaURec; depth_{MaxShower} [mm];Xseed-Xclu [mm]", (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth, 200, -100., 100.);
                _hTHistDeltaVquality = cogResol.make<TH2D>( "DeltaVquality", "DeltaVquality; depth_{MaxShower} [mm];Yseed-Yclu [mm]", (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth, 200, -100., 100.);
                _hTHistDeltaVRec     = cogResol.make<TH2D>( "DeltaVRec", "DeltaVRec; depth_{MaxShower} [mm]; Yseed-Yclu [mm]", (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth, 200, -100., 100.);
                _hTHistDeltaWquality = cogResol.make<TH2D>( "DeltaWquality", "DeltaWquality; depth_{MaxShower} [mm];vZseed-Zclu [mm]", (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth, 200, -100., 100.);
                _hTHistDeltaWRec     = cogResol.make<TH2D>( "DeltaWRec", "DeltaWRec; depth_{MaxShower} [mm];Zseed-Zclu [mm]", (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth,200, -100., 100.);

                _hTHistImpactParam   = dcaCog.make<TH2D>( "ImpactParam", "ImpactParam;  depth_{MaxShower} [mm];ImpactParam [mm]", (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth, 600,  0., 150.);
                _hTHistDistanceToCog = dcaCog.make<TH2D>( "DistanceToCog", "DistanceToCog; depth_{MaxShower} [mm]; |#vec{x}_{seed} - #vec{x}_{cog}| [mm]", (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth,  600,  0., 150.);
                _hTHistMomDotVaneNorm= dcaCog.make<TH2D>( "MomDotVaneNorm","MomDotVaneNorm ; depth_{MaxShower} [mm] ;ArcCos(#vec{p}, #hat{n}) / |#vec{p}| ", (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth, 720 , -180., 180.0);

                _hTHistRecImpactParam   = dcaCog.make<TH2D>( "RecImpactParam", "RecImpactParam;  depth_{MaxShower} [mm];ImpactParam [mm]", (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth, 400,  0., 100.);
                _hTHistRecDistanceToCog = dcaCog.make<TH2D>( "RecDistanceToCog", "RecDistanceToCog; depth_{MaxShower} [mm]; |#vec{x}_{seed} - #vec{x}_{cog}| [mm]", (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth, 400,  0., 100.);
                _hTHistMomRecDotVaneNorm= dcaCog.make<TH2D>( "MomRecDotVaneNorm","MomRecDotVaneNorm ; depth_{MaxShower} [mm]; ArcCos(#vec{p}, #hat{n}) / |#vec{p}| ", (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth, 720 , -180., 180.0);

                _hTHistOffSet           = tfs->make<TH2D>("OffSetVsDepth","OffSetVsDepth; OffSet [#]; Depth [mm]", (_MaxCogOffSet  - _MinCogOffSet) / _PitchCogOffSet, _MinCogOffSet, _MaxCogOffSet, (_MaxDepth - 10.) / _depthPitch, 10., _MaxDepth) ;
                _Ntup        = tfs->make<TTree>("ClusterTrj", "Cluster trajectory info");

                _Ntup->Branch("evt", &_evt , "evt/F");
                _Ntup->Branch("nDepthSteps",    &_nDepthSteps , "nDepthSteps/I");
                _Ntup->Branch("clVane[nDepthSteps]",    _clVane , "clVane[nDepthSteps]/I");
                _Ntup->Branch("clE[nDepthSteps]",        _clE , "clE[nDepthSteps]/F");
                _Ntup->Branch("clT[nDepthSteps]",        _clT , "clT[nDepthSteps]/F");
                _Ntup->Branch("clSize[nDepthSteps]",     _clSize , "clSize[nDepthSteps]/I");
                _Ntup->Branch("clCOGx[nDepthSteps]",     _clCOGx , "clCOGx[nDepthSteps]/F");
                _Ntup->Branch("clCOGy[nDepthSteps]",     _clCOGy , "clCOGy[nDepthSteps]/F");
                _Ntup->Branch("clCOGz[nDepthSteps]",     _clCOGz , "clCOGz[nDepthSteps]/F");
                _Ntup->Branch("clCOGu[nDepthSteps]",     _clCOGu , "clCOGu[nDepthSteps]/F");
                _Ntup->Branch("clCOGv[nDepthSteps]",     _clCOGv , "clCOGv[nDepthSteps]/F");
                _Ntup->Branch("clCOGw[nDepthSteps]",     _clCOGw , "clCOGw[nDepthSteps]/F");

                _Ntup->Branch("clCOGrow[nDepthSteps]",                        _clCOGrow , "clCOGrow[nDepthSteps]/F");
                _Ntup->Branch("clCOGcolumn[nDepthSteps]",                  _clCOGcolumn , "clCOGcolumn[nDepthSteps]/F");
                _Ntup->Branch("clCogCrySize[nDepthSteps]",                _clCogCrySize , "clCogCrySize[nDepthSteps]/I");
                _Ntup->Branch("clRows[nDepthSteps][*clSize]",                   _clRows , "clRows[nDepthSteps][*clSize]/F");
                _Ntup->Branch("clColumns[nDepthSteps][*clSize]",             _clColumns , "clColumns[nDepthSteps][*clSize]/F");
                _Ntup->Branch("clCogRows[nDepthSteps][*clCogCrySize]",       _clCogRows , "clCogRows[nDepthSteps][*clCogCrySize]/F");
                _Ntup->Branch("clCogColumns[nDepthSteps][*clCogCrySize]", _clCogColumns , "clCogColumns[nDepthSteps][*clCogCrySize]/F");

                _Ntup->Branch("cryEdep[nDepthSteps]",    _cryEdep , "cryEdep[nDepthSteps]/F");
                _Ntup->Branch("cryEdepTot[nDepthSteps]", _cryEdepTot , "cryEdepTot[nDepthSteps]/F");
                _Ntup->Branch("DCAu[nDepthSteps]",       _DCAu , "DCAu[nDepthSteps]/F");
                _Ntup->Branch("DCAv[nDepthSteps]",       _DCAv , "DCAv[nDepthSteps]/F");
                _Ntup->Branch("DCAw[nDepthSteps]",       _DCAw , "DCAw[nDepthSteps]/F");
                _Ntup->Branch("depth[nDepthSteps]",      _depth , "depth[nDepthSteps]/F");
                _Ntup->Branch("seedPx[nDepthSteps]",     _seedPx , "seedPx[nDepthSteps]/F");
                _Ntup->Branch("seedPy[nDepthSteps]",     _seedPy , "seedPy[nDepthSteps]/F");
                _Ntup->Branch("seedPz[nDepthSteps]",     _seedPz , "seedPz[nDepthSteps]/F");
                _Ntup->Branch("seedE[nDepthSteps]",      _seedE , "seedE[nDepthSteps]/F");
                _Ntup->Branch("seedTime[nDepthSteps]",   _seedTime , "seedTime[nDepthSteps]/F");
                _Ntup->Branch("seedPpx[nDepthSteps]",    _seedPpx , "seedPpx[nDepthSteps]/F");
                _Ntup->Branch("seedPpy[nDepthSteps]",    _seedPpy , "seedPpy[nDepthSteps]/F");
                _Ntup->Branch("seedPpz[nDepthSteps]",    _seedPpz , "seedPpz[nDepthSteps]/F");
                _Ntup->Branch("seedPpu[nDepthSteps]",    _seedPpu , "seedPpu[nDepthSteps]/F");
                _Ntup->Branch("seedPpv[nDepthSteps]",    _seedPpv , "seedPpv[nDepthSteps]/F");
                _Ntup->Branch("seedPpw[nDepthSteps]",    _seedPpw , "seedPpw[nDepthSteps]/F");
                _Ntup->Branch("seedRow[nDepthSteps]",    _seedRow, "seedRow[nDepthSteps]/F");
                _Ntup->Branch("seedColumn[nDepthSteps]", _seedColumn , "seedColumn[nDepthSteps]/F");
                _Ntup->Branch("seedPu[nDepthSteps]",     _seedPu , "seedPu[nDepthSteps]/F");
                _Ntup->Branch("seedPv[nDepthSteps]",     _seedPv , "seedPv[nDepthSteps]/F");
                _Ntup->Branch("seedPw[nDepthSteps]",     _seedPw , "seedPw[nDepthSteps]/F");
                _Ntup->Branch("seedPpCosTh[nDepthSteps]",_seedPpCosTh , "seedPpCosTh[nDepthSteps]/F");
                _Ntup->Branch("seedPpPhi[nDepthSteps]",  _seedPpPhi , "seedPpPhi[nDepthSteps]/F");
                _Ntup->Branch("seedPpRotCosTh[nDepthSteps]",_seedPpRotCosTh , "seedPpRotCosTh[nDepthSteps]/F");
                _Ntup->Branch("seedPpRotPhi[nDepthSteps]",  _seedPpRotPhi , "seedPpRotPhi[nDepthSteps]/F");
                _Ntup->Branch("seedPdgId[nDepthSteps]",  _seedPdgId , "seedPdgId[nDepthSteps]/I");
                _Ntup->Branch("seedIsGen[nDepthSteps]",  _seedIsGen , "seedIsGen[nDepthSteps]/I");
                _Ntup->Branch("seedThetaW[nDepthSteps]", _seedThetaW , "seedThetaW[nDepthSteps]/F");
                _Ntup->Branch("seedThetaV[nDepthSteps]", _seedThetaV , "seedThetaV[nDepthSteps]/F");
                _Ntup->Branch("seedDeltaW[nDepthSteps]", _seedDeltaW , "seedDeltaW[nDepthSteps]/F");
                _Ntup->Branch("seedDeltaV[nDepthSteps]", _seedDeltaV , "seedDeltaV[nDepthSteps]/F");
                _Ntup->Branch("seedDeltaCorrW[nDepthSteps]", _seedDeltaCorrW , "seedDeltaCorrW[nDepthSteps]/F");
                _Ntup->Branch("seedDeltaCorrV[nDepthSteps]", _seedDeltaCorrV , "seedDeltaCorrV[nDepthSteps]/F");
                _Ntup->Branch("cogDca[nDepthSteps]", _cogDca , "cogDca[nDepthSteps]/F");
                _Ntup->Branch("distToCog[nDepthSteps]", _distToCog , "distToCog[nDepthSteps]/F");
                _Ntup->Branch("prova1[nDepthSteps]", _prova1 , "prova1[nDepthSteps]/F");
                _Ntup->Branch("prova2[nDepthSteps]", _prova2 , "prova2[nDepthSteps]/F");
        }

        doCalorimeter(evt, _skipEvent);

} // end of analyze

void CaloClusterLogCog::endJob() {

}



void CaloClusterLogCog::doCalorimeter(art::Event const& evt, bool skip){


        if ( _diagLevel > 0 ) cout << "MakeCaloCluster: produce() begin" << endl;

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

        auto_ptr<CaloClusterCollection>             caloClustersPointer(new CaloClusterCollection);

        art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
        evt.getByLabel(_extractElectronsData,genEltrksHandle);

        VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();
        std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it;

        double trkMomCut = 100.0;//MeV
        std::vector<unsigned int> tmpV, tmpVTot;
        int NtrkCut =0;
        int NtrkTot = 0;

        //mapping &counting the electrons with quality cuts in the TRK
        for ( genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
                ++NtrkTot;
                VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);
                GenElHitData& hdil = iEltrk.getithLoopHit(0);

                if(!findTrkId(tmpVTot, iEltrk.getTrkID().asUint() ) ){

                        //cout<<"faccio il puschback in tmpV..."<<endl;
                        tmpVTot.push_back( iEltrk.getTrkID().asUint() );
                        //cout<<"--------> puschback in tmpV eseguito..."<<endl;

                }

                if (iEltrk.getNumOfHit()>=20){



                        if(hdil._hitMomentum.mag() >= trkMomCut){
                                NtrkCut++;// # of electrons at the entrance of the TRK which have trkMomCut MeV of kinetic energy and which will light at least 20 straws
                                if(!findTrkId(tmpV, iEltrk.getTrkID().asUint() ) ){

                                        //cout<<"faccio il puschback in tmpV..."<<endl;
                                        tmpV.push_back( iEltrk.getTrkID().asUint() );
                                        //cout<<"--------> puschback in tmpV eseguito..."<<endl;

                                }

                                //                                for(unsigned int b=1; b <=(unsigned int) _hTHistEff->GetNbinsX(); ++b){
                                //
                                //                                        if (hdil._hitMomentum.mag() >= (_hTHistEff->GetXaxis()->GetBinCenter(b) -  0.5*_hTHistEff->GetXaxis()->GetBinWidth(b) )){
                                //                                               // NGoodElec[b-1]++;
                                //                                                if ( iEltrk.isConversionEl() ) {
                                //                                                     //   NGoodConvEl[b-1]++;
                                //                                                }
                                //                                        }
                                //
                                //                                }
                        }

                }

        }//end TRK mapping
        //cout <<"NtrkCut = " <<NtrkCut<<endl;

        if (NtrkCut==0) return;

        globalNtrkCut += NtrkCut;
        globalNtrkTot += NtrkTot;

        for(int b =1; b <= _hTHistOffSet->GetNbinsX(); ++b){

                _evt = evt.id().event();
                _cogOffSet = _hTHistOffSet->GetXaxis()->GetBinCenter(b);


                for(int bin =1; bin <= _hTHistDeltaURec->GetNbinsX(); ++bin){
                        std::vector<unsigned int> trkVec, trkVecTot;
                        for(unsigned int f = 0; f != tmpV.size(); ++f){
                                trkVec.push_back(tmpV[f]);
                        }
                        for(unsigned int f = 0; f != tmpVTot.size(); ++f){
                                trkVecTot.push_back(tmpVTot[f]);
                        }
                        ElectronMap elecMap;

                        if(caloClusters->size()>0 ){
                                int iVane;
                                for(size_t icl=0; icl<caloClusters->size(); ++icl){
                                        CaloCluster const& clu = (*caloClusters).at(icl);

                                        caloClustersPointer->push_back(clu);
                                        CaloClusterCollection::iterator tmpCluster = caloClustersPointer->end();
                                        tmpCluster--;


                                        //                                        CLHEP::Hep3Vector cogDepth = cog_depth(*tmpCluster,_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) );//- 0.5*_hTHistDeltaURec->GetXaxis()->GetBinWidth(bin) );
                                        ClusterMap clusterMap;
                                        LOGcogMap( *tmpCluster, _hTHistOffSet->GetXaxis()->GetBinCenter(b), _hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) , clusterMap );
                                        CLHEP::Hep3Vector cogDepth = clusterMap._cluCOG;

                                        iVane = clu.vaneId();

                                        CaloCrystalHitPtrVector caloClusterHits = clu.caloCrystalHitsPtrVector();

                                        for(size_t i=0; i<caloClusterHits.size(); ++i){
                                                CaloCrystalHit const& hit = *(caloClusterHits.at(i));

                                                std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                                                if(ROIds.size()<1 ) continue;

                                                CaloHit const& thehit = *ROIds.at(0);
                                                size_t collectionPosition = ROIds.at(0).key();

                                                PtrStepPointMCVector const & mcptr(hits_mcptr->at(collectionPosition));
                                                if(mcptr.size() <=0) continue;
                                                size_t nHitsPerCrystal = mcptr.size();

                                                for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {

                                                        //cout<< "Start loop..."<< "j2 = "<< j2<<endl;

                                                        StepPointMC const& mchit = *mcptr[j2];

                                                        // The simulated particle that made this hit.
                                                        SimParticleCollection::key_type trackId(mchit.trackId());
                                                        SimParticle const& sim = *(simParticles->getOrNull(mchit.trackId()));

                                                        CLHEP::Hep3Vector cryFrame = cg->toCrystalFrame(thehit.id(), mchit.position());

                                                        if( (cg->crystalRByRO(thehit.id()) != ( _rowToCanc - 1 ) && cg->crystalZByRO(thehit.id()) != ( _columnToCanc - 1 ) ) ){

                                                                if(elecMap[iVane][trackId.asUint()]._impTime > mchit.time() ){
                                                                        elecMap[iVane][trackId.asUint()]._cluEnergy = clu.energyDep();
                                                                        elecMap[iVane][trackId.asUint()]._cluTime = clu.time();
                                                                        elecMap[iVane][trackId.asUint()]._cluSize = clu.size();
                                                                        elecMap[iVane][trackId.asUint()]._cryEnergyDep = hit.energyDep();
                                                                        elecMap[iVane][trackId.asUint()]._cryEnergyDepTotal = hit.energyDepTotal();
                                                                        elecMap[iVane][trackId.asUint()]._impTime = mchit.time();
                                                                        elecMap[iVane][trackId.asUint()]._impEnergy = mchit.momentum().mag();
                                                                        elecMap[iVane][trackId.asUint()]._cluCog = cogDepth;// clu._impactPoint;//logCog;
                                                                        elecMap[iVane][trackId.asUint()]._impPos = mchit.position();
                                                                        elecMap[iVane][trackId.asUint()]._impPosCryFrame = cryFrame;
                                                                        elecMap[iVane][trackId.asUint()]._vane = iVane;
                                                                        elecMap[iVane][trackId.asUint()]._row = cg->crystalRByRO(thehit.id());
                                                                        elecMap[iVane][trackId.asUint()]._colum = cg->crystalZByRO(thehit.id());
                                                                        elecMap[iVane][trackId.asUint()]._impMom3Vec = mchit.momentum();
                                                                        elecMap[iVane][trackId.asUint()]._impPdgId = sim.pdgId();
                                                                        elecMap[iVane][trackId.asUint()]._impIsGen = sim.fromGenerator();
                                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluSize        = clusterMap._cluSize      ;
                                                                        elecMap[iVane][trackId.asUint()]._clusterMap._COGcrySize     = clusterMap._COGcrySize   ;
                                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluSize        = clusterMap._cluSize      ;
                                                                        elecMap[iVane][trackId.asUint()]._clusterMap._vaneId         = clusterMap._vaneId       ;
                                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluCogRow      = clusterMap._cluCogRow    ;
                                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluCogColumn   = clusterMap._cluCogColumn ;
                                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluCOG         = clusterMap._cluCOG       ;
                                                                        elecMap[iVane][trackId.asUint()]._clusterMap._rowVec         = clusterMap._rowVec       ;
                                                                        elecMap[iVane][trackId.asUint()]._clusterMap._columnVec      = clusterMap._columnVec    ;
                                                                        elecMap[iVane][trackId.asUint()]._clusterMap._COGrowVec      = clusterMap._COGrowVec    ;
                                                                        elecMap[iVane][trackId.asUint()]._clusterMap._COGcolumnVec   = clusterMap._COGcolumnVec ;



                                                                        if(_diagLevel < 0){
                                                                                cout<< "###################"<<endl;
                                                                                cout<< "idVande = "<< iVane<<endl;
                                                                                cout << "cluU = "<<elecMap[iVane][trackId.asUint()]._impPos.getX()<<endl;
                                                                                cout << "cluV = "<<elecMap[iVane][trackId.asUint()]._impPos.getY()<<endl;
                                                                                cout << "cluW = "<<elecMap[iVane][trackId.asUint()]._impPos.getZ()<<endl;
                                                                                cout<< "###################"<<endl;
                                                                        }

                                                                }
                                                        }
                                                }

                                        }


                                }//end for(caloClsuters.size)
                        }//end if(caloClsuters.size() > 0)

                        //counting how many reconstructed clusters match with the generated particles
                        unsigned int canc2 = 0;

                        //cout <<"-------------------> numero di elettroni generati = "<< trkVecTot.size() <<endl;
                        unsigned int size2 = trkVecTot.size();
                        unsigned int it2=0;

                        _nDepthSteps = size2;

                        while( it2 < size2){
                                ElectronMap::iterator ite = elecMap.begin();
                                bool trovato = false;
                                while(!trovato && ite!=elecMap.end() ){
                                        if(ite->second.find(trkVecTot[it2]) != ite->second.end()){
                                                _clVane[it2] = ite->second[trkVec[it2]]._vane;
                                                _clE[it2]    = ite->second[trkVec[it2]]._cluEnergy;
                                                _clT[it2]    = ite->second[trkVec[it2]]._cluTime;
                                                _clSize[it2] = ite->second[trkVec[it2]]._cluSize;
                                                _clCOGu[it2] = ite->second[trkVec[it2]]._cluCog.x();
                                                _clCOGv[it2] = ite->second[trkVec[it2]]._cluCog.y();
                                                _clCOGw[it2] = ite->second[trkVec[it2]]._cluCog.z();
                                                CLHEP::Hep3Vector Mu2eFrame = cg->fromVaneFrame(ite->first, ite->second[trkVec[it2]]._cluCog);
                                                _clCOGx[it2] = Mu2eFrame.x();
                                                _clCOGy[it2] = Mu2eFrame.y();
                                                _clCOGz[it2] = Mu2eFrame.z();
                                                _cryEdep[it2] = ite->second[trkVec[it2]]._cryEnergyDep;
                                                _cryEdepTot[it2] = ite->second[trkVec[it2]]._cryEnergyDepTotal;
                                                _clCogCrySize[it2]     = ite->second[trkVec[it2]]._clusterMap._COGcrySize;
                                                _clCOGrow[it2]         = ite->second[trkVec[it2]]._clusterMap._cluCogRow;
                                                _clCOGcolumn[it2]      = ite->second[trkVec[it2]]._clusterMap._cluCogColumn;
                                                for(int i = 0; i<_clSize[it2]; ++i ){
                                                        _clRows[it2][i]         = ite->second[trkVec[it2]]._clusterMap._rowVec[i];
                                                        _clColumns[it2][i]      = ite->second[trkVec[it2]]._clusterMap._columnVec[i];
                                                }
                                                for(int j=0; j<_clCogCrySize[it2]; ++j ){
                                                        _clCogRows[it2][j]      = ite->second[trkVec[it2]]._clusterMap._COGrowVec[j];
                                                        _clCogColumns[it2][j]   = ite->second[trkVec[it2]]._clusterMap._COGcolumnVec[j];
                                                }
                                                _seedPx[it2]    = ite->second[trkVec[it2]]._impPos.x();
                                                _seedPy[it2]    = ite->second[trkVec[it2]]._impPos.y();
                                                _seedPz[it2]    =ite->second[trkVec[it2]]._impPos.z();
                                                _seedTime[it2]  = ite->second[trkVec[it2]]._impTime;
                                                _seedE[it2]     =ite->second[trkVec[it2]]._impEnergy;
                                                _seedPpx[it2]   = ite->second[trkVec[it2]]._impMom3Vec.x();
                                                _seedPpy[it2]   = ite->second[trkVec[it2]]._impMom3Vec.y();
                                                _seedPpz[it2]   = ite->second[trkVec[it2]]._impMom3Vec.z();
                                                _seedRow[it2]   = ite->second[trkVec[it2]]._row;
                                                _seedColumn[it2]= ite->second[trkVec[it2]]._colum;

                                                _seedPpCosTh[it2] = ite->second[trkVec[it2]]._impMom3Vec.cosTheta()*180./TMath::Pi();
                                                _seedPpPhi[it2] = ite->second[trkVec[it2]]._impMom3Vec.phi()*180./TMath::Pi();

                                                _seedPdgId[it2] = ite->second[trkVec[it2]]._impPdgId;
                                                _seedIsGen[it2] = ite->second[trkVec[it2]]._impIsGen;

                                                //                                        cout << "delta X = " << ite->second[trkVec[it2]]._impPos.getX() - ite->second[trkVec[it2]]._cluCog.getX()<<endl;
                                                //                                        cout << "delta Y = " << ite->second[trkVec[it2]]._impPos.getY() - ite->second[trkVec[it2]]._cluCog.getY()<<endl;
                                                //                                        cout << "delta Z = " << ite->second[trkVec[it2]]._impPos.getZ() - ite->second[trkVec[it2]]._cluCog.getZ()<<endl;



                                                CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, ite->second[trkVecTot[it2]]._impPos);
                                                _seedPu[it2] = vaneFrame.x();
                                                _seedPv[it2] = vaneFrame.y();
                                                _seedPw[it2] = vaneFrame.z();

                                                CLHEP::Hep3Vector dirMom = ite->second[trkVecTot[it2]]._impMom3Vec.unit();

                                                Vane const &vane = cg->vane(ite->first);
                                                CLHEP::Hep3Vector dirMom_rotated = (vane.rotation())*dirMom;
                                                _seedPpu[it2] = dirMom_rotated.x()*ite->second[trkVecTot[it2]]._impMom3Vec.mag();
                                                _seedPpv[it2] = dirMom_rotated.y()*ite->second[trkVecTot[it2]]._impMom3Vec.mag();
                                                _seedPpw[it2] = dirMom_rotated.z()*ite->second[trkVecTot[it2]]._impMom3Vec.mag();

                                                _seedPpRotCosTh[it2] = dirMom_rotated.cosTheta()*180./TMath::Pi();
                                                _seedPpRotPhi[it2]   = dirMom_rotated.phi()*180./TMath::Pi();

                                                LinePointPCA lppca(vaneFrame,dirMom_rotated,  ite->second[trkVecTot[it2]]._cluCog);
                                                //                                        LinePointPCA lppca(ite->second[trkVecTot[it2]]._impPos,dirMom,  ite->second[trkVecTot[it2]]._cluCog);

                                                CLHEP::Hep3Vector dcaV = lppca.unit();
                                                double impactParam = lppca.dca();
                                                _cogDca[it2] = impactParam;
                                                _DCAu[it2]   = dcaV.x()*impactParam;
                                                _DCAv[it2]   = dcaV.y()*impactParam;
                                                _DCAw[it2]   = dcaV.z()*impactParam;

                                                _prova1[it2] = sqrt( pow(dcaV.z(),2) + pow(dcaV.x(),2))*impactParam;

                                                _prova2[it2] = sqrt( pow(dcaV.y(),2) + pow(dcaV.x(),2))*impactParam;

                                                double distanceToCog = sqrt( pow(vaneFrame.getX() -ite->second[trkVecTot[it2]]._cluCog.getX() ,2) + pow(vaneFrame.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() ,2) + pow(vaneFrame.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() ,2) );
                                                _distToCog[it2] = distanceToCog;

                                                //                                        double distanceToCog = sqrt( pow(ite->second[trkVecTot[it2]]._impPos.getX() -ite->second[trkVecTot[it2]]._cluCog.getX() ,2) + pow(ite->second[trkVecTot[it2]]._impPos.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() ,2) + pow(ite->second[trkVecTot[it2]]._impPos.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() ,2) );
                                                double tmpDepth = _hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaURec->GetXaxis()->GetBinWidth(bin)*0.5;
                                                _depth[it2] = tmpDepth;

                                                //                                                _hTHistRecDistanceToCog->Fill(tmpDepth, distanceToCog);
                                                //                                                _hTHistRecImpactParam->Fill(tmpDepth, impactParam);
                                                //                                                _hTHistMomRecDotVaneNorm->Fill(tmpDepth, asin(impactParam / distanceToCog)*180./TMath::Pi());

                                                //                                                if(bin == 1){
                                                //                                                        _hTHistDeltaEnergVRec->Fill( ite->second[trkVecTot[it2]]._impEnergy - ite->second[trkVecTot[it2]]._cluEnergy);
                                                //
                                                //
                                                //                                                }
                                                //CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, vane.getOrigin() );
                                                double thetaWimpact = std::atan(-1.0*dirMom_rotated.getZ() / dirMom_rotated.getX() ) ;
                                                _seedThetaW[it2] = thetaWimpact*180./TMath::Pi();

                                                double deltaZ = (tmpDepth /*- impactParam/cos(thetaWimpact) */)*tan(thetaWimpact);
                                                deltaZ =  vaneFrame.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() + deltaZ ;
                                                _seedDeltaW[it2] = deltaZ;
                                                std::string W="W";
                                                double corrW = triangoloVar(vaneFrame.z(), W);
                                                deltaZ -= corrW;
                                                _seedDeltaCorrW[it2] = deltaZ ;

                                                double thetaVimpact = std::atan(dirMom_rotated.getY() /  dirMom_rotated.getX() ) ;
                                                _seedThetaV[it2] = thetaVimpact*180./TMath::Pi();

                                                double deltaY = (tmpDepth )*tan(thetaVimpact);
                                                deltaY = vaneFrame.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() - deltaY ;
                                                _seedDeltaV[it2] = deltaY;
                                                std::string V="V";
                                                double corrV = triangoloVar(vaneFrame.y(), V);
                                                deltaY -= corrV;
                                                _seedDeltaCorrV[it2] = deltaY ;

                                                if(_diagLevel < 0){
                                                        cout<<"vaneId = "<<ite->first<<endl;
                                                        cout<<"posX = "<<ite->second[trkVecTot[it2]]._impPos.getX()<<", posY = "<<ite->second[trkVecTot[it2]]._impPos.getY()<<", posZ = "<< ite->second[trkVecTot[it2]]._impPos.getZ()<<endl;
                                                        cout<<"on vane ref..."<<endl;
                                                        cout<<"posU = "<<vaneFrame.getX()<<", posV = "<<vaneFrame.getY()<<", posW = "<< vaneFrame.getZ()<<endl;
                                                        cout << "dirMomX = "<< dirMom.getX()<<", dirMomY = "<< dirMom.getY()<<", dirMomZ = "<<dirMom.getZ() <<endl;
                                                        cout<<"after rotation..."<<endl;
                                                        cout << "dirMomU = "<< dirMom_rotated.getX()<<", dirMomV = "<< dirMom_rotated.getY()<<", dirMomW = "<<dirMom_rotated.getZ() <<endl;
                                                        cout<<"cog coordinates... on vane ref!"<<endl;
                                                        cout<<"cluU = "<< ite->second[trkVecTot[it2]]._cluCog.getX()<< ", cluV = "<<ite->second[trkVecTot[it2]]._cluCog.getY() <<", cluW = "<<ite->second[trkVecTot[it2]]._cluCog.getZ()<<endl;
                                                        cout<<"------------------------------------------------------"<<endl;
                                                        cout<< "dcaV_U = "<<dcaV.getX()<< ", dcaV_V = "<<dcaV.getY()<<", dcaV_W = "<<dcaV.getZ()<<endl;
                                                        cout << "impactParam = "<< impactParam<<endl;
                                                        cout<< "distance ImpPointToCog = "<< distanceToCog << endl;
                                                        cout<< "-----------------------------> vane = "<< ite->first<<endl;
                                                        cout << "impX = "<< vaneFrame.getX()<<endl;
                                                        cout << "impY = "<< vaneFrame.getY()<<endl;
                                                        cout << "impZ = "<< vaneFrame.getZ()<<endl;
                                                        cout<<"------------------------------------------------------"<<endl;
                                                        cout << "cluX = "<<ite->second[trkVec[it2]]._cluCog.getX()<<endl;
                                                        cout << "cluY = "<<ite->second[trkVec[it2]]._cluCog.getY()<<endl;
                                                        cout << "cluZ = "<<ite->second[trkVec[it2]]._cluCog.getZ()<<endl;
                                                        cout<<"------------------------------------------------------"<<endl;
                                                        cout<<"seedDeltaCorrW = "<< deltaZ<<endl;
                                                        cout<<"seedDeltaCorrV = "<< deltaY<<endl;
                                                }

                                                //                                                _hTHistDeltaURec->Fill(tmpDepth, vaneFrame.getX() - ite->second[trkVecTot[it2]]._cluCog.getX() );
                                                //                                                _hTHistDeltaWRec->Fill(tmpDepth,  deltaZ );
                                                //                                                _hTHistDeltaVRec->Fill(tmpDepth,  deltaY );


                                                //                                        if(ite->first == 0 || ite->first == 2){
                                                //
                                                //                                                _hTHistDeltaURec->Fill(_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._impPos.getX() - ite->second[trkVecTot[it2]]._cluCog.getX() );
                                                //                                        }else{
                                                //                                                _hTHistDeltaVRec->Fill(_hTHistDeltaVRec->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._impPos.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() );
                                                //                                        }
                                                //                                        _hTHistDeltaWRec->Fill(_hTHistDeltaWRec->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._impPos.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() );

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

                        //counting how many reconstructed clusters match with the generated particles which pass the quality cuts of the tracker
                        //                        unsigned int canc = 0;
                        if(_diagLevel < 0){
                                cout <<"after.... trkVecTot.size = "<< trkVecTot.size() <<endl;
                                cout <<"-------------------> number of quality generated electrons = "<< trkVec.size() <<endl;
                        }
                        //                        unsigned int size = trkVec.size();
                        //                        unsigned int it=0;
                        //                        while( it < size){
                        //                                ElectronMap::iterator ite = elecMap.begin();
                        //                                bool trovato = false;
                        //                                while(!trovato && ite!=elecMap.end() ){
                        //                                        if(ite->second.find(trkVec[it]) != ite->second.end()){
                        //
                        //                                                if(_diagLevel < 0){
                        //                                                        cout << "impX = "<<ite->second[trkVec[it]]._impPosCryFrame.getX()<<endl;
                        //                                                        cout << "impY = "<<ite->second[trkVec[it]]._impPosCryFrame.getY()<<endl;
                        //                                                        cout << "impZ = "<<ite->second[trkVec[it]]._impPosCryFrame.getZ()<<endl;
                        //                                                }
                        //
                        //                                                //globalCaloCut[bin - 1] += 1.0;
                        //
                        //                                                trovato = true;
                        //                                                canc = it;
                        //
                        //                                                CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, ite->second[trkVec[it]]._impPos);
                        //
                        //                                                CLHEP::Hep3Vector dirMom = ite->second[trkVec[it]]._impMom3Vec.unit();
                        //
                        //                                                Vane const &vane = cg->getVane(ite->first);
                        //                                                CLHEP::Hep3Vector dirMom_rotated = *(vane.getRotation())*dirMom;
                        //
                        //                                                LinePointPCA lppca(vaneFrame,dirMom_rotated,  ite->second[trkVec[it]]._cluCog);
                        //
                        //                                                CLHEP::Hep3Vector dcaV = lppca.unit();
                        //                                                double impactParam = lppca.dca();
                        //
                        //                                                double distanceToCog = sqrt( pow(vaneFrame.getX() -ite->second[trkVec[it]]._cluCog.getX() ,2) + pow(vaneFrame.getY() - ite->second[trkVec[it]]._cluCog.getY() ,2) + pow(vaneFrame.getZ() - ite->second[trkVec[it]]._cluCog.getZ() ,2) );
                        //
                        ////                                                _hTHistDistanceToCog->Fill(_hTHistDistanceToCog->GetXaxis()->GetBinCenter(bin) -_hTHistDistanceToCog->GetXaxis()->GetBinWidth(bin)*0.5, distanceToCog);
                        ////                                                _hTHistImpactParam->Fill(_hTHistImpactParam->GetXaxis()->GetBinCenter(bin) -_hTHistImpactParam->GetXaxis()->GetBinWidth(bin)*0.5, impactParam);
                        ////                                                _hTHistMomDotVaneNorm->Fill(_hTHistMomDotVaneNorm->GetXaxis()->GetBinCenter(bin) -_hTHistMomDotVaneNorm->GetXaxis()->GetBinWidth(bin)*0.5, asin(impactParam / distanceToCog)*180./TMath::Pi());
                        //
                        //                                                if(bin ==1){
                        //                                           //             _hTHistDeltaEnergy->Fill(ite->second[trkVec[it]]._impEnergy - ite->second[trkVec[it]]._cluEnergy);
                        //                                                        //_hTHistDeltaWquality->Fill(vaneFrame.getY() - ite->second[trkVec[it]]._cluCog.getY() );
                        //                                                }
                        //                                                //CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, vane.getOrigin() );
                        //                                                double tmpDepth = _hTHistDeltaUquality->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaUquality->GetXaxis()->GetBinWidth(bin)*0.5;
                        //
                        //                                                double impactParamPrjW = impactParam*(pow(dcaV.getX(), 2) + pow(dcaV.getZ(), 2));
                        //                                                double impactParamPrjV = impactParam*(pow(dcaV.getX(), 2) + pow(dcaV.getY(), 2));
                        //
                        //                                                double deltaW = (tmpDepth - impactParamPrjW/cos(_thetaWimpact) )*tan(_thetaWimpact);
                        //                                                int sign = -1;
                        //                                                if(dcaV.getY() < 0){
                        //                                                        sign = 1;
                        //                                                }
                        //                                                double deltaV = sign*(tmpDepth - impactParamPrjV/cos(_thetaVimpact) )*tan(_thetaVimpact);
                        //
                        ////                                                _hTHistDeltaUquality->Fill(tmpDepth, vaneFrame.getX() - ite->second[trkVec[it]]._cluCog.getX() );
                        ////                                                _hTHistDeltaWquality->Fill(tmpDepth, vaneFrame.getZ() - ite->second[trkVec[it]]._cluCog.getZ() + deltaW );
                        ////                                                _hTHistDeltaVquality->Fill(tmpDepth, vaneFrame.getY() - ite->second[trkVec[it]]._cluCog.getY() + deltaV );
                        //
                        //                                                //  _hTHistDeltaUquality->Fill(_hTHistDeltaUquality->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaUquality->GetXaxis()->GetBinWidth(bin)*0.5, vaneFrame.getX() - ite->second[trkVec[it]]._cluCog.getX() );
                        //                                                // _hTHistDeltaWquality->Fill(_hTHistDeltaWquality->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaWquality->GetXaxis()->GetBinWidth(bin)*0.5, vaneFrame.getZ() - ite->second[trkVec[it]]._cluCog.getZ() );
                        //
                        //                                        }
                        //                                        ++ite;
                        //                                }
                        //
                        //                                if(trovato){
                        //                                        std::vector<unsigned int>::iterator er = trkVec.begin();
                        //                                        er +=canc;
                        //                                        trkVec.erase(er);
                        //                                        size = trkVec.size();
                        //                                }else{
                        //                                        ++it;
                        //                                }
                        //
                        //                        }

                        //                        if(_diagLevel < 0){
                        //                                cout <<"after.... trkVec.size = "<< trkVec.size() <<endl;
                        //                                cout<<"...\n...\n "<<endl;
                        //                        }

                }//end for(_hTHistEff->GetNbinsX())


                _Ntup->Fill();
        }

        cout << "Event "<<evt.id().event()<<" CaloClusterLogCog done..."<<endl;

}

}


using mu2e::CaloClusterLogCog;
DEFINE_ART_MODULE(CaloClusterLogCog);

