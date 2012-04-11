//
// implementation of different algorithm to reconstruct the impact position
//
// $Id: CaloClusterCogCorrFunc_module.cc,v 1.6 2012/04/11 10:35:13 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/04/11 10:35:13 $
//
// Original author G. Pezzullo
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
#include "CalorimeterGeom/inc/Calorimeter.hh"
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

struct elecData{
        double _cluEnergy;
        double _cluTime;
        int _cluSize;
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

        bool operator<( const elecData other) const{
                return ( _impTime< other._impTime);
        }
        elecData & operator=(const elecData& other) {
                _cluEnergy = other._cluEnergy;
                _cluTime   = other._cluTime;
                _cluSize   = other._cluSize;
                _cryEnergyDep = other._cryEnergyDep;
                _cryEnergyDepTotal = other._cryEnergyDepTotal;
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
                _impIsGen  = other._impIsGen;
                _clusterMap = other._clusterMap;


                return *this;
        }
        elecData():
                _impTime(1e10),
                _impEnergy(0.0){
        }
};

//the key is the the vane
typedef std::map<unsigned int,std::map<unsigned int, elecData > > ElecMap;


static int ncalls(0);

class CaloClusterCogCorrFunc : public art::EDAnalyzer {
public:
        explicit CaloClusterCogCorrFunc(fhicl::ParameterSet const& pset):
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
        _Ntup(0),
        _CogOffSet(pset.get<double>("cogOffset",3.)),
        _Depth(pset.get<double>("depth",11.)),
        _application(0),
        _directory(0)
        {
        }
        virtual ~CaloClusterCogCorrFunc() {
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

        TH1D* _hTHistDeltaUquality;
        TH1D* _hTHistDeltaURec;
        TH1D* _hTHistDeltaVquality;
        TH1D* _hTHistDeltaVRec;
        TH1D* _hTHistDeltaWquality;
        TH1D* _hTHistDeltaWRec;
        TH1D* _hTHistImpactParam;
        TH1D* _hTHistDistanceToCog;
        TH1D* _hTHistMomDotVaneNorm;
        TH1D*_hTHistRecImpactParam;
        TH1D*_hTHistRecDistanceToCog;
        TH1D*_hTHistMomRecDotVaneNorm;
        TTree* _Ntup;



        double _CogOffSet;
        double _Depth;

        double globalNtrkCut;
        double globalNtrkTot;
        double *globalCaloCut;
        double *globalRecCaloCut;

        std::auto_ptr<MCCaloUtilities> CaloManager;

        bool _skipEvent;




        Int_t _seedId,
        _clSize,
        _seedPdgId ,
        _seedIsGen,
        _clVane,
        _clCogCrySize;

        Float_t _evt,
        _clE ,
        _clT ,
        _clCOGx ,
        _clCOGy ,
        _clCOGz ,
        _clCOGu ,
        _clCOGv ,
        _clCOGw ,
        _clCOGrow,
        _clCOGcolumn,
        _clRows[10000],
        _clColumns[10000],
        _clCogRows[10000],
        _clCogColumns[10000],
        _clShowerDir,
        _clErrShowerDir,
        _cryEdep[10000],
        _cryEdepTot,
        _DCAu,
        _DCAv,
        _DCAw,
        _seedPx ,
        _seedPy ,
        _seedPz ,
        _seedPpu ,
        _seedPpv ,
        _seedPpw ,
        _seedE ,
        _seedTime,
        _seedPpx ,
        _seedPpy ,
        _seedPpz ,
        _seedRow,
        _seedColumn,
        _seedPu ,
        _seedPv ,
        _seedPw ,
        _seedPpCosTh,
        _seedPpPhi ,
        _seedPpRotCosTh,
        _seedPpRotPhi ,
        _seedThetaW,
        _seedThetaV,
        _seedDeltaW,
        _seedDeltaV,
        _seedDeltaCorrW,
        _seedDeltaCorrV,
        _distToCog,
        _cogDca,
        _depth,
        _prova1,
        _prova2;


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

double triangoloVar(double x, char d){
        double y =0.;
        double par[6] = {0.};
        if(d == 'V' ){
                //                double p1[6] = { 0.1689, 29.92, -1.344, 4.136, 0.6882, 0.004642};//linear corr
                double p1[6] = {0.3822, 29.75, -7.595, 5.384, 0.6845, 4.219e-15};//pol4 corr
                for(int y=0; y<6; ++y){
                        par[y] = p1[y];
                }
        } else if (d == 'W'){
                //                double p2[6] = { 0.343, 30., -4.25, -1.441, 0.5152, 0.0001403};//linear corr
                double p2[6] = {0.3201, 30., -24.76, -1.257, 0.5072, 9.911e-6};//pol4 corr
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

double corrFunc(double x, double *par){
        double res=0.0;

        double y = par[0],
                        c = par[1],
                        d = par[2],
                        e = par[3],
                        f = par[4],
                        g = par[5],
                        h = par[6],
                        i = par[7],
                        l = par[8],
                        y1 = par[9],
                        o = par[10],
                        p = par[11],
                        q = par[12];


        double a = -1.*(-3.*e + 3.*l + -2.*d*y - c*y*y + h*y*y -f*pow(y,4) )/pow(y,4),
                        b = -1.*(4.*e - 4.*l +3.*d*y -3.*i*y +2.*c*y*y - 2.*h*y*y - g*y*y*y) / pow(y,3),
                        m = -1.*(3.*l -3.*q + 2*i*y1 - 2.*p*y1 + h*y1*y1 - o*y1*y1 - f*pow(y1, 4) )/pow(y1, 4),
                        n = -1.*(-4.*l +4.*q -3*i*y1 + 3.*p*y1 -2.*h*y1*y1 +2.*o*y1*y1 - g*pow(y1, 3))/pow(y1, 3);

        if(x >= y){
                res += a*pow(x, 4);
                res += b*pow(x, 3);
                res += c*pow(x, 2);
                res += d*x;
                res += e;
        } else if(x < y && x > y1){
                res += f*pow(x, 4);
                res += g*pow(x, 3);
                res += h*pow(x, 2);
                res += i*x;
                res += l;
        } else if(x<=y1){
                res += m*pow(x, 4);
                res += n*pow(x, 3);
                res += o*pow(x, 2);
                res += p*x;
                res += q;
        }

        return res;
}

double cogVcorrFunc(double x){
        //        double par[13] = {-0.5943, -3.393, -26.8, 0.9915, -7.886, -30.87, -35.87, -0.3304, -6.319, -1.594, 171.8, 250.5, 114.3};//function of tan(seedThetaV)
        //        double par[13] = {4139., 0.8839, -33.3, 441.9, -2.449e-7, -4.651e-5, -0.0006987, 0.4832, 0.9593, -2.156e4, 0., 0., 0.};
        double par[13] = {1.06e4, 0.8839, -33.3, 441.9, -0.0002132, 0.03959, -2.759, 85.47, -975.8, 40., -0.0007464, 0.4805, 0.9861};//func of seedThetaV
        return corrFunc(x, par);
}

double cogWcorrFunc(double x){
        //        double par[13] = {4.263e4, 0., 0., 1., 7.466, -56.24, 159.3, -200.3, 71.76, 1.487, -32.87, -21.72, -0.8187};
        //        double par[13] = {44.03, 0.8839, -33.3, 441.9, -6.497e-6, 0.001035, -0.04497, -0.00266, 0.6607, -219.7, 0., 0., 0.};//func of seedThetaW
        double par[13] = {51.86, -859.6, 9.515e-7, 0.0002434, -0.02172, -0.002659, -4.938, 11.84, 5.795, -29.08, 48.75};//func of seedThetaW
        return corrFunc(x, par);
}


void CaloClusterCogCorrFunc::beginJob( ) {

        cout << "start CaloClusterCogCorrFunc..."<<endl;

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

void CaloClusterCogCorrFunc::analyze(art::Event const & evt ) {

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

                _hTHistDeltaUquality = cogResol.make<TH1D>( "DeltaUquality", "DeltaUquality ;Xseed-Xclu [mm]", 200, -100., 100. );
                _hTHistDeltaURec     = cogResol.make<TH1D>( "DeltaURec", "DeltaURec ;Xseed-Xclu [mm]",  200, -100., 100.);
                _hTHistDeltaVquality = cogResol.make<TH1D>( "DeltaVquality", "DeltaVquality ;Yseed-Yclu [mm]",  200, -100., 100.);
                _hTHistDeltaVRec     = cogResol.make<TH1D>( "DeltaVRec", "DeltaVRec ; Yseed-Yclu [mm]",  200, -100., 100.);
                _hTHistDeltaWquality = cogResol.make<TH1D>( "DeltaWquality", "DeltaWquality ;vZseed-Zclu [mm]",  200, -100., 100.);
                _hTHistDeltaWRec     = cogResol.make<TH1D>( "DeltaWRec", "DeltaWRec ;Zseed-Zclu [mm]", 200, -100., 100.);

                _hTHistImpactParam   = dcaCog.make<TH1D>( "ImpactParam", "ImpactParam;  depth_{MaxShower} [mm];ImpactParam [mm]",  600,  0., 150.);
                _hTHistDistanceToCog = dcaCog.make<TH1D>( "DistanceToCog", "DistanceToCog ; |#vec{x}_{seed} - #vec{x}_{cog}| [mm]",   600,  0., 150.);
                _hTHistMomDotVaneNorm= dcaCog.make<TH1D>( "MomDotVaneNorm","MomDotVaneNorm   ;ArcCos(#vec{p}, #hat{n}) / |#vec{p}| ",  720 , -180., 180.0);

                _hTHistRecImpactParam   = dcaCog.make<TH1D>( "RecImpactParam", "RecImpactParam;  depth_{MaxShower} [mm];ImpactParam [mm]",  400,  0., 100.);
                _hTHistRecDistanceToCog = dcaCog.make<TH1D>( "RecDistanceToCog", "RecDistanceToCog ; |#vec{x}_{seed} - #vec{x}_{cog}| [mm]",  400,  0., 100.);
                _hTHistMomRecDotVaneNorm= dcaCog.make<TH1D>( "MomRecDotVaneNorm","MomRecDotVaneNorm  ; ArcCos(#vec{p}, #hat{n}) / |#vec{p}| ",  720 , -180., 180.0);
                _Ntup        = tfs->make<TTree>("ClusterTrj", "Cluster trajectory info");

                _Ntup->Branch("evt", &_evt , "evt/F");
                _Ntup->Branch("clVane",    &_clVane , "clVane/I");
                _Ntup->Branch("clE",        &_clE , "clE/F");
                _Ntup->Branch("clT",        &_clT , "clT/F");
                _Ntup->Branch("clSize",    &_clSize , "clSize/I");
                _Ntup->Branch("clCOGx",     &_clCOGx , "clCOGx/F");
                _Ntup->Branch("clCOGy",     &_clCOGy , "clCOGy/F");
                _Ntup->Branch("clCOGz",     &_clCOGz , "clCOGz/F");
                _Ntup->Branch("clCOGu",     &_clCOGu , "clCOGu/F");
                _Ntup->Branch("clCOGv",     &_clCOGv , "clCOGv/F");
                _Ntup->Branch("clCOGw",     &_clCOGw , "clCOGw/F");

                _Ntup->Branch("clCOGrow",    &_clCOGrow , "clCOGrow/F");
                _Ntup->Branch("clCOGcolumn",    &_clCOGcolumn , "clCOGcolumn/F");
                _Ntup->Branch("clCogCrySize",    &_clCogCrySize , "clCogCrySize/I");
                _Ntup->Branch("clRows[clSize]",_clRows , "clRows[clSize]/F");
                _Ntup->Branch("clColumns[clSize]",_clColumns , "clColumns[clSize]/F");
                _Ntup->Branch("clCogRows[clCogCrySize]",_clCogRows , "clCogRows[clCogCrySize]/F");
                _Ntup->Branch("clCogColumns[clCogCrySize]",_clCogColumns , "clCogColumns[clCogCrySize]/F");


                _Ntup->Branch("clShowerDir", &_clShowerDir , "clShowerDir/F");
                _Ntup->Branch("clErrShowerDir", &_clErrShowerDir , "clErrShowerDir/F");
                _Ntup->Branch("cryEdep[clSize]",     _cryEdep , "cryEdep[clSize]/F");
                _Ntup->Branch("cryEdepTot",     &_cryEdepTot , "cryEdepTot/F");
                _Ntup->Branch("DCAu",     &_DCAu , "DCAu/F");
                _Ntup->Branch("DCAv",     &_DCAv , "DCAv/F");
                _Ntup->Branch("DCAw",     &_DCAw , "DCAw/F");
                _Ntup->Branch("seedPx",     &_seedPx , "seedPx/F");
                _Ntup->Branch("seedPy",     &_seedPy , "seedPy/F");
                _Ntup->Branch("seedPz",     &_seedPz , "seedPz/F");
                _Ntup->Branch("seedE",      &_seedE , "seedE/F");
                _Ntup->Branch("seedTime",   &_seedTime , "seedTime/F");
                _Ntup->Branch("seedPpx",    &_seedPpx , "seedPpx/F");
                _Ntup->Branch("seedPpy",    &_seedPpy , "seedPpy/F");
                _Ntup->Branch("seedPpz",    &_seedPpz , "seedPpz/F");
                _Ntup->Branch("seedPpu",    &_seedPpu , "seedPpu/F");
                _Ntup->Branch("seedPpv",    &_seedPpv , "seedPpv/F");
                _Ntup->Branch("seedPpw",    &_seedPpw , "seedPpw/F");
                _Ntup->Branch("seedRow",    &_seedRow, "seedRow/F");
                _Ntup->Branch("seedColumn", &_seedColumn , "seedColumn/F");
                _Ntup->Branch("seedPu",     &_seedPu , "seedPu/F");
                _Ntup->Branch("seedPv",     &_seedPv , "seedPv/F");
                _Ntup->Branch("seedPw",     &_seedPw , "seedPw/F");
                _Ntup->Branch("seedPpCosTh",&_seedPpCosTh , "seedPpCosTh/F");
                _Ntup->Branch("seedPpPhi",  &_seedPpPhi , "seedPpPhi/F");
                _Ntup->Branch("seedPpRotCosTh",&_seedPpRotCosTh , "seedPpRotCosTh/F");
                _Ntup->Branch("seedPpRotPhi",  &_seedPpRotPhi , "seedPpRotPhi/F");
                _Ntup->Branch("seedPdgId",  &_seedPdgId , "seedPdgId/I");
                _Ntup->Branch("seedIsGen",  &_seedIsGen , "seedIsGen/I");
                _Ntup->Branch("seedThetaW", &_seedThetaW , "seedThetaW/F");
                _Ntup->Branch("seedThetaV", &_seedThetaV , "seedThetaV/F");
                _Ntup->Branch("seedDeltaW", &_seedDeltaW , "seedDeltaW/F");
                _Ntup->Branch("seedDeltaV", &_seedDeltaV , "seedDeltaV/F");
                _Ntup->Branch("seedDeltaCorrW", &_seedDeltaCorrW , "seedDeltaCorrW/F");
                _Ntup->Branch("seedDeltaCorrV", &_seedDeltaCorrV , "seedDeltaCorrV/F");

                _Ntup->Branch("cogDca", &_cogDca , "cogDca/F");
                _Ntup->Branch("depth", &_depth , "depth/F");
                _Ntup->Branch("distToCog", &_distToCog , "distToCog/F");
                _Ntup->Branch("prova1", &_prova1 , "prova1/F");
                _Ntup->Branch("prova2", &_prova2 , "prova2/F");
        }

        doCalorimeter(evt, _skipEvent);

} // end of analyze

void CaloClusterCogCorrFunc::endJob() {

}



void CaloClusterCogCorrFunc::doCalorimeter(art::Event const& evt, bool skip){


        /*if ( _diagLevel > 0 )*/ cout << "CaloClusterCogCorrFunc: analyze() begin" << endl;

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


        //        for(int bin =1; bin <= _hTHistDeltaURec->GetNbinsX(); ++bin){
        std::vector<unsigned int> trkVec, trkVecTot;
        for(unsigned int f = 0; f != tmpV.size(); ++f){
                trkVec.push_back(tmpV[f]);
        }
        for(unsigned int f = 0; f != tmpVTot.size(); ++f){
                trkVecTot.push_back(tmpVTot[f]);
        }
        ElecMap elecMap;

        if(caloClusters->size()>0 ){
                int iVane;
                for(size_t icl=0; icl<caloClusters->size(); ++icl){
                        CaloCluster const& clu = (*caloClusters).at(icl);

                        caloClustersPointer->push_back(clu);
                        CaloClusterCollection::iterator tmpCluster = caloClustersPointer->end();
                        tmpCluster--;

                        if(_diagLevel <0) cout<<"calculating cog()..."<<endl;
                        //CLHEP::Hep3Vector cogDepth = cog_depth(*tmpCluster,_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) );//- 0.5*_hTHistDeltaURec->GetXaxis()->GetBinWidth(bin) );
                        //CLHEP::Hep3Vector cogDepth = LOGcog(*tmpCluster, _CogOffSet, _hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) );//- 0.5*_hTHistDeltaURec->GetXaxis()->GetBinWidth(bin) );
                        ClusterMap clusterMap;
                        //LOGcogMap( *tmpCluster, _CogOffSet, _hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) , clusterMap );
                        cog_depth( *tmpCluster, _Depth , clusterMap );
                        CLHEP::Hep3Vector cogDepth = clusterMap._cluCOG;
                        if(_diagLevel <0) cout<<"cog() calculated!!"<<endl;
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

                                        if( (cg->getCrystalRByRO(thehit.id()) != ( _rowToCanc - 1 ) && cg->getCrystalZByRO(thehit.id()) != ( _columnToCanc - 1 ) ) ){

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
                                                        elecMap[iVane][trackId.asUint()]._row = cg->getCrystalRByRO(thehit.id());
                                                        elecMap[iVane][trackId.asUint()]._colum = cg->getCrystalZByRO(thehit.id());
                                                        elecMap[iVane][trackId.asUint()]._impMom3Vec = mchit.momentum();
                                                        elecMap[iVane][trackId.asUint()]._impPdgId = sim.pdgId();
                                                        elecMap[iVane][trackId.asUint()]._impIsGen = sim.fromGenerator();
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluSize        = clusterMap._cluSize      ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._COGcrySize     = clusterMap._COGcrySize   ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluSize        = clusterMap._cluSize      ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._vane           = clusterMap._vane         ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluCogRow      = clusterMap._cluCogRow    ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluCogColumn   = clusterMap._cluCogColumn ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluCOG         = clusterMap._cluCOG       ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._rowVec         = clusterMap._rowVec       ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._columnVec      = clusterMap._columnVec    ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cryEdepVec     = clusterMap._cryEdepVec  ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._COGrowVec      = clusterMap._COGrowVec    ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._COGcolumnVec   = clusterMap._COGcolumnVec ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._showerDir      = clusterMap._showerDir ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._errShowerDir   = clusterMap._errShowerDir;


                                                        if(_diagLevel < 0){
                                                                cout<< "###################"<<endl;
                                                                cout<< "idVande = "<< iVane<<endl;
                                                                cout << "cluU = "<<elecMap[iVane][trackId.asUint()]._impPos.getX()<<endl;
                                                                cout << "cluV = "<<elecMap[iVane][trackId.asUint()]._impPos.getY()<<endl;
                                                                cout << "cluW = "<<elecMap[iVane][trackId.asUint()]._impPos.getZ()<<endl;
                                                                cout<< "elecMapclusterMar.CogRow = "<<elecMap[iVane][trackId.asUint()]._clusterMap._cluCogRow<<endl;
                                                                cout<< "elecMapclusterMar.CogColumn = "<<elecMap[iVane][trackId.asUint()]._clusterMap._cluCogColumn<<endl;
                                                                cout<< "elecMapclusterMar.COGcrySize ="<<elecMap[iVane][trackId.asUint()]._clusterMap._COGcrySize<<endl;
                                                                cout<< "clusterMar.CogRow = "<<clusterMap._cluCogRow<<endl;
                                                                cout<< "clusterMar.CogColumn = "<<clusterMap._cluCogColumn<<endl;
                                                                cout<< "clusterMar.COGcrySize ="<<clusterMap._COGcrySize<<endl;
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
        while( it2 < size2){
                ElecMap::iterator ite = elecMap.begin();
                bool trovato = false;
                while(!trovato && ite!=elecMap.end() ){
                        if(ite->second.find(trkVecTot[it2]) != ite->second.end()){
                                _evt = evt.id().event();
                                _seedPdgId = 0;
                                _seedIsGen= 0;

                                _clVane = ite->second[trkVec[it2]]._vane;
                                _clE = ite->second[trkVec[it2]]._cluEnergy;
                                _clT = ite->second[trkVec[it2]]._cluTime;
                                _clSize = ite->second[trkVec[it2]]._cluSize;
                                _clCOGu = ite->second[trkVec[it2]]._cluCog.x();
                                _clCOGv = ite->second[trkVec[it2]]._cluCog.y();
                                _clCOGw = ite->second[trkVec[it2]]._cluCog.z();
                                CLHEP::Hep3Vector Mu2eFrame = cg->fromVaneFrame(ite->first, ite->second[trkVec[it2]]._cluCog);
                                _clCOGx = Mu2eFrame.x();
                                _clCOGy = Mu2eFrame.y();
                                _clCOGz = Mu2eFrame.z();

                                _cryEdepTot = ite->second[trkVec[it2]]._cryEnergyDepTotal;
                                _clCogCrySize     = ite->second[trkVec[it2]]._clusterMap._COGcrySize;
                                _clCOGrow         = ite->second[trkVec[it2]]._clusterMap._cluCogRow;
                                _clCOGcolumn      = ite->second[trkVec[it2]]._clusterMap._cluCogColumn;
                                for(int i = 0; i<_clSize; ++i ){
                                        _cryEdep[i] = ite->second[trkVec[it2]]._cryEnergyDep;
                                        _clRows[i]         = ite->second[trkVec[it2]]._clusterMap._rowVec[i];
                                        _clColumns[i]      = ite->second[trkVec[it2]]._clusterMap._columnVec[i];
                                }
                                for(int j=0; j<_clCogCrySize; ++j ){
                                        _clCogRows[j]      = ite->second[trkVec[it2]]._clusterMap._COGrowVec[j];
                                        _clCogColumns[j]   = ite->second[trkVec[it2]]._clusterMap._COGcolumnVec[j];
                                }
                                _clShowerDir =  ite->second[trkVec[it2]]._clusterMap._showerDir;
                                _clErrShowerDir = ite->second[trkVec[it2]]._clusterMap._errShowerDir;

                                _seedPx = ite->second[trkVec[it2]]._impPos.x();
                                _seedPy = ite->second[trkVec[it2]]._impPos.y();
                                _seedPz =ite->second[trkVec[it2]]._impPos.z();
                                _seedTime  = ite->second[trkVec[it2]]._impTime;
                                _seedE =ite->second[trkVec[it2]]._impEnergy;
                                _seedPpx = ite->second[trkVec[it2]]._impMom3Vec.x();
                                _seedPpy = ite->second[trkVec[it2]]._impMom3Vec.y();
                                _seedPpz = ite->second[trkVec[it2]]._impMom3Vec.z();
                                _seedRow = ite->second[trkVec[it2]]._row;
                                _seedColumn = ite->second[trkVec[it2]]._colum;

                                _seedPpCosTh = ite->second[trkVec[it2]]._impMom3Vec.cosTheta()*180./TMath::Pi();
                                _seedPpPhi = ite->second[trkVec[it2]]._impMom3Vec.phi()*180./TMath::Pi();

                                _seedPdgId = ite->second[trkVec[it2]]._impPdgId;
                                _seedIsGen = ite->second[trkVec[it2]]._impIsGen;
                                //                                        cout << "delta X = " << ite->second[trkVec[it2]]._impPos.getX() - ite->second[trkVec[it2]]._cluCog.getX()<<endl;
                                //                                        cout << "delta Y = " << ite->second[trkVec[it2]]._impPos.getY() - ite->second[trkVec[it2]]._cluCog.getY()<<endl;
                                //                                        cout << "delta Z = " << ite->second[trkVec[it2]]._impPos.getZ() - ite->second[trkVec[it2]]._cluCog.getZ()<<endl;



                                CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, ite->second[trkVecTot[it2]]._impPos);
                                _seedPu = vaneFrame.x();
                                _seedPv = vaneFrame.y();
                                _seedPw = vaneFrame.z();

                                CLHEP::Hep3Vector dirMom = ite->second[trkVecTot[it2]]._impMom3Vec.unit();
                                if(_diagLevel < 1){
                                        if(cry(ite->second[trkVec[it2]]._cluCog.y()) ==15.0 || cry( ite->second[trkVec[it2]]._cluCog.z() ) ==15.0 ){
                                                cout<<"ecco il mukkio...."<<endl;
                                                cout<<"cogDepth.v() = " <<cry(ite->second[trkVec[it2]]._cluCog.y())<<", cogDepth.w() = "<< cry( ite->second[trkVec[it2]]._cluCog.z() )<<", cogDepth.u() = "<<cry( ite->second[trkVec[it2]]._cluCog.x() )<<endl;
                                                cout <<"iVane = "<<ite->second[trkVec[it2]]._vane<<", cluSize = "<<ite->second[trkVec[it2]]._cluSize<<", cluEdep = "<< ite->second[trkVec[it2]]._cluEnergy<<", depth = "<<_Depth<<endl;
                                                //                                                        CLHEP::Hep3Vector cogDepth2 = cog_depth(*tmpCluster,_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) );
                                                //                                                        cout<<"instead, the linear COG give..."<<endl;
                                                //                                                        cout<<"cogDepth2.v() = " <<cry( cogDepth2.y() )<<", cogDepth2.w() = "<< cry( cogDepth2.z() )<<", cogDepth2.u() = "<< cry( cogDepth2.x() )<<endl;
                                                cout<<"fine del mukkio...."<<endl;


                                        }
                                        cout<<"vaneId = "<<ite->first<<endl;
                                        cout<<"posX = "<<ite->second[trkVecTot[it2]]._impPos.getX()<<", posY = "<<ite->second[trkVecTot[it2]]._impPos.getY()<<", posZ = "<< ite->second[trkVecTot[it2]]._impPos.getZ()<<endl;
                                        cout<<"on vane ref..."<<endl;
                                        cout<<"posU = "<<vaneFrame.getX()<<", posV = "<<vaneFrame.getY()<<", posW = "<< vaneFrame.getZ()<<endl;
                                        cout << "dirMomX = "<< dirMom.getX()<<", dirMomY = "<< dirMom.getY()<<", dirMomZ = "<<dirMom.getZ() <<endl;
                                }
                                Vane const &vane = cg->getVane(ite->first);
                                CLHEP::Hep3Vector dirMom_rotated = *(vane.getRotation())*dirMom;
                                _seedPpu = dirMom_rotated.x()*ite->second[trkVecTot[it2]]._impMom3Vec.mag();
                                _seedPpv = dirMom_rotated.y()*ite->second[trkVecTot[it2]]._impMom3Vec.mag();
                                _seedPpw = dirMom_rotated.z()*ite->second[trkVecTot[it2]]._impMom3Vec.mag();
                                _seedPpRotCosTh = dirMom_rotated.cosTheta()*180./TMath::Pi();
                                _seedPpRotPhi = dirMom_rotated.phi()*180./TMath::Pi();

                                if(_diagLevel < 0){
                                        cout<<"after rotation..."<<endl;
                                        cout << "dirMomU = "<< dirMom_rotated.getX()<<", dirMomV = "<< dirMom_rotated.getY()<<", dirMomW = "<<dirMom_rotated.getZ() <<endl;
                                        cout<<"cog coordinates... on vane ref!"<<endl;
                                        cout<<"cluU = "<< ite->second[trkVecTot[it2]]._cluCog.getX()<< ", cluV = "<<ite->second[trkVecTot[it2]]._cluCog.getY() <<", cluW = "<<ite->second[trkVecTot[it2]]._cluCog.getZ()<<endl;
                                }

                                LinePointPCA lppca(vaneFrame,dirMom_rotated,  ite->second[trkVecTot[it2]]._cluCog);
                                //                                        LinePointPCA lppca(ite->second[trkVecTot[it2]]._impPos,dirMom,  ite->second[trkVecTot[it2]]._cluCog);

                                CLHEP::Hep3Vector dcaV = lppca.unit();
                                double impactParam = lppca.dca();
                                _cogDca = impactParam;
                                _DCAu = dcaV.x()*impactParam;
                                _DCAv = dcaV.y()*impactParam;
                                _DCAw = dcaV.z()*impactParam;

                                //                                        _prova1 = sqrt( pow(dcaV.z(),2) + pow(dcaV.x(),2))*impactParam;

                                _prova2 = sqrt( pow(dcaV.y(),2) + pow(dcaV.x(),2))*impactParam;

                                double distanceToCog = sqrt( pow(vaneFrame.getX() -ite->second[trkVecTot[it2]]._cluCog.getX() ,2) + pow(vaneFrame.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() ,2) + pow(vaneFrame.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() ,2) );
                                _distToCog = distanceToCog;

                                //                                        double distanceToCog = sqrt( pow(ite->second[trkVecTot[it2]]._impPos.getX() -ite->second[trkVecTot[it2]]._cluCog.getX() ,2) + pow(ite->second[trkVecTot[it2]]._impPos.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() ,2) + pow(ite->second[trkVecTot[it2]]._impPos.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() ,2) );
                                double tmpDepth = _Depth;//_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaURec->GetXaxis()->GetBinWidth(bin)*0.5;

                                _depth = tmpDepth;

                                _hTHistRecDistanceToCog->Fill(  distanceToCog);
                                _hTHistRecImpactParam->Fill(  impactParam);
                                _hTHistMomRecDotVaneNorm->Fill(  asin(impactParam / distanceToCog)*180./TMath::Pi());
                                if(_diagLevel < 0){
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
                                }
                                // if(bin == 1){
                                _hTHistDeltaEnergVRec->Fill( ite->second[trkVecTot[it2]]._impEnergy - ite->second[trkVecTot[it2]]._cluEnergy);


                                //}
                                //CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, vane.getOrigin() );
                                double thetaWimpact = std::atan(-1.0*dirMom_rotated.getZ() / dirMom_rotated.getX() ) ;
                                _seedThetaW = thetaWimpact*180./TMath::Pi();

                                double deltaZ = cogWcorrFunc(thetaWimpact );
                                deltaZ =  vaneFrame.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() - deltaZ ;
                                _seedDeltaW = deltaZ;


                                double corrW = triangoloVar(vaneFrame.z(), 'W');
                                deltaZ -= corrW;
                                _seedDeltaCorrW = deltaZ ;


                                double thetaVimpact = std::atan(dirMom_rotated.getY() /  dirMom_rotated.getX() ) ;
                                _seedThetaV = thetaVimpact*180./TMath::Pi();

                                double deltaY = cogVcorrFunc( thetaVimpact );
                                _prova1     = vaneFrame.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() - deltaY;

                                deltaY = vaneFrame.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() - deltaY ;

                                _seedDeltaV = deltaY;

                                double corrV = triangoloVar(vaneFrame.y(), 'V');
                                deltaY -= corrV;
                                _seedDeltaCorrV = deltaY ;

                                _hTHistDeltaURec->Fill(  vaneFrame.getX() - ite->second[trkVecTot[it2]]._cluCog.getX() );
                                _hTHistDeltaWRec->Fill(   deltaZ );
                                _hTHistDeltaVRec->Fill(   deltaY );


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

        _Ntup->Fill();

        //counting how many reconstructed clusters match with the generated particles which pass the quality cuts of the tracker
        unsigned int canc = 0;
        if(_diagLevel < 0){
                cout <<"after.... trkVecTot.size = "<< trkVecTot.size() <<endl;
                cout <<"-------------------> number of quality generated electrons = "<< trkVec.size() <<endl;
        }
        unsigned int size = trkVec.size();
        unsigned int it=0;
        while( it < size){
                ElecMap::iterator ite = elecMap.begin();
                bool trovato = false;
                while(!trovato && ite!=elecMap.end() ){
                        if(ite->second.find(trkVec[it]) != ite->second.end()){

                                if(_diagLevel < 0){
                                        cout << "impX = "<<ite->second[trkVec[it]]._impPosCryFrame.getX()<<endl;
                                        cout << "impY = "<<ite->second[trkVec[it]]._impPosCryFrame.getY()<<endl;
                                        cout << "impZ = "<<ite->second[trkVec[it]]._impPosCryFrame.getZ()<<endl;
                                }

                                //globalCaloCut[bin - 1] += 1.0;

                                trovato = true;
                                canc = it;

                                CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, ite->second[trkVec[it]]._impPos);

                                CLHEP::Hep3Vector dirMom = ite->second[trkVec[it]]._impMom3Vec.unit();

                                Vane const &vane = cg->getVane(ite->first);
                                CLHEP::Hep3Vector dirMom_rotated = *(vane.getRotation())*dirMom;

                                LinePointPCA lppca(vaneFrame,dirMom_rotated,  ite->second[trkVec[it]]._cluCog);

                                CLHEP::Hep3Vector dcaV = lppca.unit();
                                double impactParam = lppca.dca();

                                double distanceToCog = sqrt( pow(vaneFrame.getX() -ite->second[trkVec[it]]._cluCog.getX() ,2) + pow(vaneFrame.getY() - ite->second[trkVec[it]]._cluCog.getY() ,2) + pow(vaneFrame.getZ() - ite->second[trkVec[it]]._cluCog.getZ() ,2) );

                                _hTHistDistanceToCog->Fill( distanceToCog);
                                _hTHistImpactParam->Fill(  impactParam);
                                _hTHistMomDotVaneNorm->Fill(  asin(impactParam / distanceToCog)*180./TMath::Pi() );

                                //if(bin ==1){
                                _hTHistDeltaEnergy->Fill(ite->second[trkVec[it]]._impEnergy - ite->second[trkVec[it]]._cluEnergy);
                                //_hTHistDeltaWquality->Fill(vaneFrame.getY() - ite->second[trkVec[it]]._cluCog.getY() );
                                //}
                                //CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, vane.getOrigin() );
                                double tmpDepth = _Depth;//_hTHistDeltaUquality->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaUquality->GetXaxis()->GetBinWidth(bin)*0.5;

                                double impactParamPrjW = impactParam*(pow(dcaV.getX(), 2) + pow(dcaV.getZ(), 2));
                                double impactParamPrjV = impactParam*(pow(dcaV.getX(), 2) + pow(dcaV.getY(), 2));

                                double deltaW = (tmpDepth - impactParamPrjW/cos(_thetaWimpact) )*tan(_thetaWimpact);
                                int sign = -1;
                                if(dcaV.getY() < 0){
                                        sign = 1;
                                }
                                double deltaV = sign*(tmpDepth - impactParamPrjV/cos(_thetaVimpact) )*tan(_thetaVimpact);

                                _hTHistDeltaUquality->Fill( vaneFrame.getX() - ite->second[trkVec[it]]._cluCog.getX() );
                                _hTHistDeltaWquality->Fill(  vaneFrame.getZ() - ite->second[trkVec[it]]._cluCog.getZ() + deltaW );
                                _hTHistDeltaVquality->Fill(  vaneFrame.getY() - ite->second[trkVec[it]]._cluCog.getY() + deltaV );

                                //  _hTHistDeltaUquality->Fill(_hTHistDeltaUquality->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaUquality->GetXaxis()->GetBinWidth(bin)*0.5, vaneFrame.getX() - ite->second[trkVec[it]]._cluCog.getX() );
                                // _hTHistDeltaWquality->Fill(_hTHistDeltaWquality->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaWquality->GetXaxis()->GetBinWidth(bin)*0.5, vaneFrame.getZ() - ite->second[trkVec[it]]._cluCog.getZ() );

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

        if(_diagLevel < 0){
                cout <<"after.... trkVec.size = "<< trkVec.size() <<endl;
                cout<<"...\n...\n "<<endl;
        }

        // }//end for(_hTHistEff->GetNbinsX())
        cout << "Event "<<evt.id().event()<<" CaloClusterCogCorrFunc done..."<<endl;

}

}


using mu2e::CaloClusterCogCorrFunc;
DEFINE_ART_MODULE(CaloClusterCogCorrFunc);

