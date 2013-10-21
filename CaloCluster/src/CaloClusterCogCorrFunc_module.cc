//
// implementation of different algorithm to reconstruct the impact position
//
// $Id: CaloClusterCogCorrFunc_module.cc,v 1.17 2013/10/21 20:34:14 gandr Exp $
// $Author: gandr $
// $Date: 2013/10/21 20:34:14 $
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
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "CaloCluster/inc/CaloClusterTools.hh"
#include "CaloCluster/inc/CaloClusterer.hh"
//#include "CaloCluster/inc/CaloClusterFinder.hh"
//#include "CaloCluster/inc/ClosestCaloClusterFinder.hh"

//#include "RecoDataProducts/inc/CaloCluster.hh"
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

struct simData{
        double entranceEnergy, exitEnergy;
        double cosTheta, cosPitch, lcosPitch, entranceTime, exitTime;
        size_t nHits;
        int isQC;
        simData(){};
};

typedef std::map<unsigned int, simData> SimMap;



struct electronData{
        double _gentime;
        double _genx;
        double _geny;
        double _genz;
        double _gencosth;
        double _genphi;
        double _genp;
        double _gene;
        double _cluEnergy;
        double _cluTime;
        double _fastCryTime;
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
  int _impIsConv;
        ClusterMap _clusterMap;

        bool operator<( const electronData other) const{
                return ( _impTime< other._impTime);
        }
        electronData & operator=(const electronData& other) {
                _cluEnergy = other._cluEnergy;
                _cluTime   = other._cluTime;
                _fastCryTime = other._fastCryTime;
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
		_impIsConv = other._impIsConv;
                _clusterMap = other._clusterMap;


                return *this;
        }
        electronData():
                _impTime(1e10),
		_impEnergy(0.0),
		_impIsConv(0){
        }
};

//the key is the the vane
typedef std::map<unsigned int,std::map<unsigned int, electronData > > ElectronMap;


static int ncalls(0);

class CaloClusterCogCorrFunc : public art::EDAnalyzer {
public:
  explicit CaloClusterCogCorrFunc(fhicl::ParameterSet const& pset):
        art::EDAnalyzer(pset),
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
        _application(nullptr),
        _directory(0)
        {
        }
        virtual ~CaloClusterCogCorrFunc() {
        }
        void beginJob();
        void endJob();

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

        std::unique_ptr<MCCaloUtilities> CaloManager;

        bool _skipEvent;




        Int_t _seedId,
	  _mcNHits,
        _clSize,
        _seedPdgId ,
        _seedIsGen,
	  _seedIsConv,
        _clVane,
        _clCogCrySize;

        Float_t _evt,
	  _mcTrkEntrE,
	  _mcTrkExitE,
	  _mcTrkEntrT0,
	  _mcTrkExitT0,
	  _mcCosTh,
	  _mcCosPitch,
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
        _clCryTimes[10000],
        _clShowerDir,
        _clErrShowerDir,
        _clCryEnergyMaxRow,
        _clCryEnergyMaxColumn,
        _clCryMaxEdep,
        _clWsize,
        _clVsize,
        _cryEdep[10000],
        _cryEdepTot,
        _DCAu,
        _DCAv,
        _DCAw,
        _seedIsGoodTrk,
        _seedNGoodTrkPerEvent,
        _seedGenp,
        _seedGene,
        _seedGenx,
        _seedGeny,
        _seedGenz,
        _seedGencosth,
        _seedGenphi,
        _seedGentime,
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
        unique_ptr<TApplication> _application;

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
                //                double p1[6] = {0.3822, 29.75, -7.595, 5.384, 0.6845, 4.219e-15};//pol4 corr
                double p1[6] = {0.4825, 30.06, -9.645, 3.767, 0.6561, -0.001441};
                for(int y=0; y<6; ++y){
                        par[y] = p1[y];
                }
        } else if (d == 'W'){
                //                double p2[6] = { 0.343, 30., -4.25, -1.441, 0.5152, 0.0001403};//linear corr
                //                double p2[6] = {0.3201, 30., -24.76, -1.257, 0.5072, 9.911e-6};//pol4 corr
                double p2[6] = {0.3203, 30., -35.87, -1.29, 0.51, 6.046e-6};
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

double triangoloVarSpread(double x, char d){
        double y =0.;
        double par[6] = {0.};
        if(d == 'V' ){
                //                double p1[6] = { 0.1689, 29.92, -1.344, 4.136, 0.6882, 0.004642};//linear corr
                //                double p1[6] = {0.3822, 29.75, -7.595, 5.384, 0.6845, 4.219e-15};//pol4 corr
                //                double p1[6] = {0.5092, 29.22, -5.829, 11.64, 0.54, 0.00167};
                double p1[6] = {0.2165, 30.11, -1.985, 93.42, 0.6901, 0.002987};//1may
                for(int y=0; y<6; ++y){
                        par[y] = p1[y];
                }
        } else if (d == 'W'){
                //                double p2[6] = { 0.343, 30., -4.25, -1.441, 0.5152, 0.0001403};//linear corr
                //                double p2[6] = {0.3201, 30., -24.76, -1.257, 0.5072, 9.911e-6};//pol4 corr
                //                double p2[6] = {0.3148, 30., -35.87, -1.33, 0.5083, 2.726e-6};
                double p2[6] = {0.269, 30.01, -2.213, -2.27, 0.5647, -1.135e-5};//1may
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

        if(d == 'V'&& x<98. ){
                double p3[9] = {-5.824, -1.428, 0.4556, -0.03528, 0.001294, -2.604e-5, 2.946e-7, -1.758e-9, 4.307e-12};
                double res = p3[0];
                res += p3[1]*x; res += p3[2]*x*x; res += p3[3]*x*x*x; res += p3[4]*pow(x,4);res += p3[5]*pow(x,5); res += p3[6]*pow(x,6);res += p3[7]*pow(x,7);res += p3[8]*pow(x,8);
                y = res;
        }else if (d == 'W' && x<30.){
                double p4[9] = {-13.59, -3.906, 1.826, -0.3313, 0.03293, -0.001898, 6.333e-5, -1.137e-6, 8.527e-9};
                y += p4[0] + p4[1]*x; y += p4[2]*x*x; y += p4[3]*x*x*x; y += p4[4]*pow(x,4);y += p4[5]*pow(x,5); y += p4[6]*pow(x,6);y += p4[7]*pow(x,7);y += p4[8]*pow(x,8);
        }
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
        //        double par[13] = {51.86, -859.6, 9.515e-7, 0.0002434, -0.02172, -0.002659, -4.938, 11.84, 5.795, -29.08, 48.75};//func of seedThetaW
        double par[13] = {7.562e05, -859.6, 9.515e-07, 0.0002434, -8.12e-06, -0.0137, 1.58e08, 4.244e16, 1340., 2.168e10, -0.04434, -0.167, 3.187 };
        return corrFunc(x, par);
}


void CaloClusterCogCorrFunc::beginJob( ) {
  //if( evt.id().event() %1000 ==0){
  //cout << "start CaloClusterCogCorrFunc..."<<endl;
	//}
        //CaloManager = unique_ptr<MCCaloUtilities>(new MCCaloUtilities());

        // If needed, create the ROOT interactive environment. See note 1.
        //if ( !gApplication ){
  //      int    tmp_argc(0);
  //            char** tmp_argv(0);
  //            _application = unique_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
  //    }

  //        gStyle->SetPalette(1);
  //    gROOT->SetStyle("Plain");

  //    _directory = gDirectory;


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
		_Ntup->Branch("mcNHits",     &_mcNHits , "mcNHits/I");
                _Ntup->Branch("mcTrkEntrE",     &_mcTrkEntrE , "mcTrkEntrE/F");
                _Ntup->Branch("mcTrkExitE",     &_mcTrkExitE , "mcTrkExitE/F");
                _Ntup->Branch("mcTrkEntrT0",     &_mcTrkEntrT0 , "mcTrkEntrT0/F");
                _Ntup->Branch("mcTrkExitT0",     &_mcTrkExitT0 , "mcTrkExitT0/F");
                _Ntup->Branch("mcCosTh",     &_mcCosTh , "mcCosTh/F");
                _Ntup->Branch("mcCosPitch",     &_mcCosPitch , "mcCosPitch/F");
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
                _Ntup->Branch("clCryTimes[clSize]",_clCryTimes , "clCryTimes[clSize]/F");


                _Ntup->Branch("clShowerDir", &_clShowerDir , "clShowerDir/F");
                _Ntup->Branch("clErrShowerDir", &_clErrShowerDir , "clErrShowerDir/F");
                _Ntup->Branch("clCryEnergyMaxRow", &_clCryEnergyMaxRow , "clCryEnergyMaxRow/F");
                _Ntup->Branch("clCryEnergyMaxColumn", &_clCryEnergyMaxColumn , "clCryEnergyMaxColumn/F");
                _Ntup->Branch("clCryMaxEdep", &_clCryMaxEdep , "clCryMaxEdep/F");
                _Ntup->Branch("clWsize", &_clWsize , "clWsize/F");
                _Ntup->Branch("clVsize", &_clVsize , "clVsize/F");
                _Ntup->Branch("cryEdep[clSize]",     _cryEdep , "cryEdep[clSize]/F");
                _Ntup->Branch("cryEdepTot",     &_cryEdepTot , "cryEdepTot/F");
                _Ntup->Branch("DCAu",     &_DCAu , "DCAu/F");
                _Ntup->Branch("DCAv",     &_DCAv , "DCAv/F");
                _Ntup->Branch("DCAw",     &_DCAw , "DCAw/F");
                _Ntup->Branch("seedIsGoodTrk",     &_seedIsGoodTrk , "seedIsGoodTrk/F");
                _Ntup->Branch("seedNGoodTrkPerEvent",     &_seedNGoodTrkPerEvent , "seedNGoodTrkPerEvent/F");
                _Ntup->Branch("seedGenp",     &_seedGenp , "seedGenp/F");
                _Ntup->Branch("seedGene",     &_seedGene , "seedGene/F");
                _Ntup->Branch("seedGenx",     &_seedGenx , "seedGenx/F");
                _Ntup->Branch("seedGeny",     &_seedGeny , "seedGeny/F");
                _Ntup->Branch("seedGenz",     &_seedGenz , "seedGenz/F");
                _Ntup->Branch("seedGentime",     &_seedGentime , "seedGentime/F");
                _Ntup->Branch("seedGencosth",     &_seedGencosth , "seedGencosth/F");
                _Ntup->Branch("seedGenphi",     &_seedGenphi , "seedGenphi/F");
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
                _Ntup->Branch("seedIsConv",  &_seedIsConv , "seedIsConv/I");
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

  	if( evt.id().event() %1000 ==0){
        /*if ( _diagLevel > 0 )*/ cout << "CaloClusterCogCorrFunc: analyze() begin" << endl;
	}
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
        evt.getByLabel(_caloClusterModuleLabel,_producerName, caloClusters );

        unique_ptr<CaloClusterCollection> caloClustersPointer(new CaloClusterCollection);

        art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
        evt.getByLabel(_extractElectronsData,genEltrksHandle);

        VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();
        std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it;

        double trkMomCut = 100.0;//MeV
        std::vector<unsigned int> tmpV, tmpVTot;
        int NtrkCut =0;
        int NtrkTot = 0;


	SimMap simMap;
        //mapping &counting the electrons with quality cuts in the TRK
        for ( genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
                ++NtrkTot;
                VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);
                GenElHitData& hdil = iEltrk.getithLoopHit(0);
		GenElHitData& ldil = iEltrk.getHit((int)(iEltrk.getNumOfHit() - 1) );

                double cosTheta = iEltrk.getTrkLrntzVec().cosTheta() ;
                double cosPitch = hdil._hitMomentum.cosTheta();
                double lcosPitch = ldil._hitMomentum.cosTheta();

                if(_diagLevel>2){
                        cout << "cosThetafirst = "<<cosPitch <<", "<<"costhetaLast = "<<lcosPitch<<endl;
                        cout<< "ldil._hitMomentum.mag() = "<<ldil._hitMomentum.mag()<<"hdil._hitMomentum.mag() = "<< hdil._hitMomentum.mag()<<endl;
                }

		bool condition = true;
		simData sData;
                sData.nHits = iEltrk.getNumOfHit() ;
                sData.cosTheta = cosTheta ;
                sData.cosPitch = cosPitch;
                sData.lcosPitch = lcosPitch;
                sData.entranceEnergy = hdil._hitMomentum.mag();
                sData.exitEnergy = ldil._hitMomentum.mag();
                sData.isQC = 0;
                sData.entranceTime = hdil._mcHitTime;
                sData.exitTime = ldil._mcHitTime;
	
                //the following condition are the same used by Dave Brown for TTracker studies
                condition &= ( iEltrk.getNumOfHit() >= 20 );
                condition &= ( hdil._hitMomentum.mag() >= trkMomCut );
                condition &= ( cosTheta >= -0.5 );
                condition &= ( cosTheta <=  0.5 );
                condition &= ( cosPitch > 0.5 );
                condition &= ( cosPitch < 0.70710678118655 );// 1 / sqrt(2)

                if(!findTrkId(tmpVTot, iEltrk.getTrkID().asUint() ) ){

                        //cout<<"faccio il puschback in tmpV..."<<endl;
                        tmpVTot.push_back( iEltrk.getTrkID().asUint() );
                        //cout<<"--------> puschback in tmpV eseguito..."<<endl;

                }

		//                if (iEltrk.getNumOfHit()>=20){
		  


		  //                        if(hdil._hitMomentum.mag() >= trkMomCut){
		NtrkCut++;// # of electrons at the entrance of the TRK which have trkMomCut MeV of kinetic energy and which will light at least 20 straws
		if( condition ){
		  sData.isQC = 1;
		  if(!findTrkId(tmpV, iEltrk.getTrkID().asUint() ) ){
		    //cout<<"faccio il puschback in tmpV..."<<endl;
		    tmpV.push_back( iEltrk.getTrkID().asUint() );
		    //cout<<"--------> puschback in tmpV eseguito..."<<endl;
		    
		  }
		}		
				//                        }
		simMap[iEltrk.getTrkID().asUint()] =  sData ;
			//  }

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
        ElectronMap elecMap;

        if(caloClusters->size()>0 ){
                int iVane;
                for(size_t icl=0; icl<caloClusters->size(); ++icl){
                        CaloCluster const& clu = (*caloClusters).at(icl);

                        caloClustersPointer->push_back(clu);
                        CaloClusterCollection::iterator tmpCluster = caloClustersPointer->end();
                        tmpCluster--;

                        if(_diagLevel >0) cout<<"calculating cog()..."<<endl;
                        //CLHEP::Hep3Vector cogDepth = cog_depth(*tmpCluster,_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) );//- 0.5*_hTHistDeltaURec->GetXaxis()->GetBinWidth(bin) );
                        //CLHEP::Hep3Vector cogDepth = LOGcog(*tmpCluster, _CogOffSet, _hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) );//- 0.5*_hTHistDeltaURec->GetXaxis()->GetBinWidth(bin) );
                        ClusterMap clusterMap;
                        //LOGcogMap( *tmpCluster, _CogOffSet, _hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) , clusterMap );
                        cog_depth( *tmpCluster, _Depth , clusterMap );
                        CLHEP::Hep3Vector cogDepth = clusterMap._cluCOG;
                        if(_diagLevel >0) cout<<"cog() calculated!!"<<endl;
                        iVane = clu.vaneId();

                        CaloCrystalHitPtrVector caloClusterHits = clu.caloCrystalHitsPtrVector();

                        for(size_t i=0; i<caloClusterHits.size(); ++i){
                                CaloCrystalHit const& hit = *(caloClusterHits.at(i));

                                std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                                if(ROIds.size()<1 ) continue;

                                CaloHit const& thehit = *ROIds.at(0);
                                size_t collectionPosition = ROIds.at(0).key();
				int crystalId = cg->crystalByRO(thehit.id());
				
                                PtrStepPointMCVector const & mcptr(hits_mcptr->at(collectionPosition));
                                if(mcptr.size() <=0) continue;

                                size_t nHitsPerCrystal = mcptr.size();

                                for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {

                                        //cout<< "Start loop..."<< "j2 = "<< j2<<endl;

                                        StepPointMC const& mchit = *mcptr[j2];

                                        // The simulated particle that made this hit.
                                        SimParticleCollection::key_type trackId(mchit.trackId());
                                        SimParticle const& sim = *(simParticles->getOrNull(mchit.trackId()));


                                        CLHEP::Hep3Vector cryFrame = cg->toCrystalFrame(crystalId, mchit.position());

                                        if( (cg->crystalRByRO(thehit.id()) != ( _rowToCanc - 1 ) && cg->crystalZByRO(thehit.id()) != ( _columnToCanc - 1 ) ) ){

                                                if(elecMap[iVane][trackId.asUint()]._impTime > mchit.time() ){
                                                        if(sim.fromGenerator() ){
                                                                //GenParticle const& gen = *sim.genParticle();
                                                                GenParticle const& gen = genParticles->at(sim.generatorIndex());
                                                                //                                                                cout <<"gen = "<< gen << endl;
                                                                //                                                                cout <<" gen.momentum().mag() = "<<gen.momentum().vect().mag()<<endl;
                                                                //                                                                cout <<" gen.momentum().e() = "<<gen.momentum().e()<<endl;
                                                                //                                                                cout <<" gen.pdgId() = "<<gen.pdgId()<<endl;
                                                                //                                                                cout <<" gen.position().x() = "<<gen.position().x()<<endl;
                                                                //                                                                cout <<" gen.position().y() = "<<gen.position().y()<<endl;
                                                                //                                                                cout <<" gen.position().z() = "<<gen.position().z()<<endl;

                                                                elecMap[iVane][trackId.asUint()]._genp      = gen.momentum().vect().mag();//sim.startMomentum().vect().mag();//
                                                                elecMap[iVane][trackId.asUint()]._gene      = gen.momentum().e();//sim.startMomentum().e();//
                                                                elecMap[iVane][trackId.asUint()]._gentime   = gen.time();//sim.startGlobalTime();//
                                                                elecMap[iVane][trackId.asUint()]._genx      = gen.position().x();//sim.startPosition().x();//
                                                                elecMap[iVane][trackId.asUint()]._geny      = gen.position().y();//sim.startPosition().y();//
                                                                elecMap[iVane][trackId.asUint()]._genz      = gen.position().z();//sim.startPosition().z();//
                                                                elecMap[iVane][trackId.asUint()]._gencosth  = gen.momentum().cosTheta();//sim.startMomentum().cosTheta();//
                                                                elecMap[iVane][trackId.asUint()]._genphi    = gen.momentum().phi();//sim.startMomentum().phi();//
                                                        }else {
                                                                elecMap[iVane][trackId.asUint()]._genp      = -1e6;
                                                                elecMap[iVane][trackId.asUint()]._gene      = -1e6;
                                                                elecMap[iVane][trackId.asUint()]._gentime   = -1e6;
                                                                elecMap[iVane][trackId.asUint()]._genx      = -1e6;
                                                                elecMap[iVane][trackId.asUint()]._geny      = -1e6;
                                                                elecMap[iVane][trackId.asUint()]._genz      = -1e6;
                                                                elecMap[iVane][trackId.asUint()]._gencosth  = -1e6;
                                                                elecMap[iVane][trackId.asUint()]._genphi    = -1e6;
                                                        }
                                                        CaloClusterTools cluTool(clu);
                                                        elecMap[iVane][trackId.asUint()]._cluEnergy                  = clu.energyDep();
                                                        elecMap[iVane][trackId.asUint()]._cluTime                    = clu.time();
                                                        elecMap[iVane][trackId.asUint()]._fastCryTime                = cluTool.timeFasterCrystal();
                                                        elecMap[iVane][trackId.asUint()]._cluSize                    = clu.size();
                                                        elecMap[iVane][trackId.asUint()]._cryEnergyDep               = hit.energyDep();
                                                        elecMap[iVane][trackId.asUint()]._cryEnergyDepTotal          = hit.energyDepTotal();
                                                        elecMap[iVane][trackId.asUint()]._impTime                    = mchit.time();
                                                        elecMap[iVane][trackId.asUint()]._impEnergy                  = mchit.momentum().mag();
                                                        elecMap[iVane][trackId.asUint()]._cluCog                     = cogDepth;// clu._impactPoint;//logCog;
                                                        elecMap[iVane][trackId.asUint()]._impPos                     = mchit.position();
                                                        elecMap[iVane][trackId.asUint()]._impPosCryFrame             = cryFrame;
                                                        elecMap[iVane][trackId.asUint()]._vane                       = iVane;
                                                        elecMap[iVane][trackId.asUint()]._row                        = cg->crystalRByRO(thehit.id());
                                                        elecMap[iVane][trackId.asUint()]._colum                      = cg->crystalZByRO(thehit.id());
                                                        elecMap[iVane][trackId.asUint()]._impMom3Vec                 = mchit.momentum();
                                                        elecMap[iVane][trackId.asUint()]._impPdgId                   = sim.pdgId();
                                                        elecMap[iVane][trackId.asUint()]._impIsGen                   = sim.fromGenerator();
							if(sim.fromGenerator() ){
							  GenParticle const& gen = *sim.genParticle();
							  GenId genId = gen.generatorId();
							  if(genId==GenId::conversionGun){
							    elecMap[iVane][trackId.asUint()]._impIsConv        =  1;
							  }
							}
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluSize        = clusterMap._cluSize      ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._COGcrySize     = clusterMap._COGcrySize   ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluSize        = clusterMap._cluSize      ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._vaneId           = clusterMap._vaneId         ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluCogRow      = clusterMap._cluCogRow    ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluCogColumn   = clusterMap._cluCogColumn ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cluCOG         = clusterMap._cluCOG       ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._rowVec         = clusterMap._rowVec       ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._columnVec      = clusterMap._columnVec    ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._cryEdepVec     = clusterMap._cryEdepVec  ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._COGrowVec      = clusterMap._COGrowVec    ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._COGcolumnVec   = clusterMap._COGcolumnVec ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._timeVec         = clusterMap._timeVec       ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._showerDir      = clusterMap._showerDir ;
                                                        elecMap[iVane][trackId.asUint()]._clusterMap._errShowerDir   = clusterMap._errShowerDir;


                                                        if(_diagLevel > 0){
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

        //counting how many reconstructed clusters match with the generated particles which pass the quality cuts of the tracker
        unsigned int canc2 = 0, canc1;

        //cout <<"-------------------> numero di elettroni generati = "<< trkVecTot.size() <<endl;
        unsigned int size2 = trkVec.size();
        unsigned int it2=0;
        while( it2 < size2){
                ElectronMap::iterator ite = elecMap.begin();
                bool trovato = false;
                while(!trovato && ite!=elecMap.end() ){
                        if(ite->second.find(trkVec[it2]) != ite->second.end()){
                                _evt                  = evt.id().event();
                                _seedPdgId            = 0;
                                _seedIsGen            = 0;
				_seedIsConv           = 0;
                                _clCryEnergyMaxRow    = 0.0;
                                _clCryEnergyMaxColumn = 0.0;
				_mcCosPitch     = 0.0;
				_mcCosTh        = 0.0;
				_mcNHits        = 0;
				_mcTrkEntrE     = 0.0;
				_mcTrkExitE     = 0.0;
				_mcTrkEntrT0    = 0.0;
				_mcTrkExitT0    = 0.0;

                                _clVane = ite->second[trkVec[it2]]._vane;
                                _clE    = ite->second[trkVec[it2]]._cluEnergy;
                                _clT    = ite->second[trkVec[it2]]._fastCryTime;//_cluTime;
                                _clSize = ite->second[trkVec[it2]]._cluSize;
                                _clCOGu = ite->second[trkVec[it2]]._cluCog.x();
                                _clCOGv = ite->second[trkVec[it2]]._cluCog.y();
                                _clCOGw = ite->second[trkVec[it2]]._cluCog.z();
                                CLHEP::Hep3Vector Mu2eFrame = cg->fromVaneFrame(ite->first, ite->second[trkVec[it2]]._cluCog);
                                _clCOGx = Mu2eFrame.x();
                                _clCOGy = Mu2eFrame.y();
                                _clCOGz = Mu2eFrame.z();

                                _cryEdepTot       = ite->second[trkVec[it2]]._cryEnergyDepTotal;
                                _clCogCrySize     = ite->second[trkVec[it2]]._clusterMap._COGcrySize;
                                _clCOGrow         = ite->second[trkVec[it2]]._clusterMap._cluCogRow;
                                _clCOGcolumn      = ite->second[trkVec[it2]]._clusterMap._cluCogColumn;

                                _clCryEnergyMaxRow    = ite->second[trkVec[it2]]._clusterMap._rowVec[0];
                                _clCryEnergyMaxColumn = ite->second[trkVec[it2]]._clusterMap._columnVec[0];

                                double tempEmaxCry    = ite->second[trkVec[it2]]._clusterMap._cryEdepVec[0];

                                _clCryMaxEdep = 0.0;

                                double tmpClVmin = cg->nCrystalR(), tmpClVmax=0., tmpClWmin = cg->nCrystalZ(), tmpClWmax = 0.;

                                for(int i = 0; i<_clSize; ++i ){
                                        _cryEdep[i] = ite->second[trkVec[it2]]._clusterMap._cryEdepVec[i];//_cryEnergyDep;
                                        if(_cryEdep[i] > _clCryMaxEdep){
                                                _clCryMaxEdep = _cryEdep[i];
                                        }

                                        _clRows[i]         = ite->second[trkVec[it2]]._clusterMap._rowVec[i];

                                        if(_clRows[i] < tmpClVmin) {
                                                tmpClVmin = _clRows[i];
                                        }
                                        if(_clRows[i] > tmpClVmax) {
                                                tmpClVmax = _clRows[i];
                                        }

                                        _clColumns[i]      = ite->second[trkVec[it2]]._clusterMap._columnVec[i];

					_clCryTimes[i]     = ite->second[trkVec[it2]]._clusterMap._timeVec[i];

                                        if(_clColumns[i] < tmpClWmin) {
                                                tmpClWmin = _clColumns[i];
                                        }
                                        if(_clColumns[i] > tmpClWmax) {
                                                tmpClWmax = _clColumns[i];
                                        }



                                        if(_cryEdep[i] > tempEmaxCry){
                                                _clCryEnergyMaxRow =      _clRows[i]   ;
                                                _clCryEnergyMaxColumn =   _clColumns[i];
                                                tempEmaxCry = _cryEdep[i];
                                        }
                                }

                                _clWsize = tmpClWmax - tmpClWmin + 1.;

                                _clVsize = tmpClVmax - tmpClVmin + 1.;

                                for(int j=0; j<_clCogCrySize; ++j ){
                                        _clCogRows[j]      = ite->second[trkVec[it2]]._clusterMap._COGrowVec[j];
                                        _clCogColumns[j]   = ite->second[trkVec[it2]]._clusterMap._COGcolumnVec[j];
                                }

                                _clShowerDir    =  ite->second[trkVec[it2]]._clusterMap._showerDir;

                                _clErrShowerDir = ite->second[trkVec[it2]]._clusterMap._errShowerDir;

                                _seedIsGoodTrk  = 0.0;
                                _seedNGoodTrkPerEvent = tmpV.size();

				  if(simMap.find((unsigned int)trkVec[it2]) != simMap.end() ){
				    SimMap::iterator it = simMap.find((unsigned int)trkVec[it2]);
				    _mcCosPitch     = it->second.cosPitch;
				    _mcCosTh        = it->second.cosTheta;
				    _mcNHits        = it->second.nHits;
				    _mcTrkEntrE     = it->second.entranceEnergy;
				    _mcTrkExitE     = it->second.exitEnergy;
				    _mcTrkEntrT0    = it->second.entranceTime;
				    _mcTrkExitT0    = it->second.exitTime;
				  }
                                _seedGenp = ite->second[trkVec[it2]]._genp;
                                _seedGene = ite->second[trkVec[it2]]._gene;
                                _seedGenx = ite->second[trkVec[it2]]._genx;
                                _seedGeny = ite->second[trkVec[it2]]._geny;
                                _seedGenz = ite->second[trkVec[it2]]._genz;
                                _seedGencosth = ite->second[trkVec[it2]]._gencosth;
                                _seedGenphi = ite->second[trkVec[it2]]._genphi;
                                _seedGentime = ite->second[trkVec[it2]]._gentime;

                                _seedPx = ite->second[trkVec[it2]]._impPos.x();
                                _seedPy = ite->second[trkVec[it2]]._impPos.y();
                                _seedPz = ite->second[trkVec[it2]]._impPos.z();
                                _seedTime  = ite->second[trkVec[it2]]._impTime;
                                _seedE = ite->second[trkVec[it2]]._impEnergy;
                                _seedPpx = ite->second[trkVec[it2]]._impMom3Vec.x();
                                _seedPpy = ite->second[trkVec[it2]]._impMom3Vec.y();
                                _seedPpz = ite->second[trkVec[it2]]._impMom3Vec.z();
                                _seedRow = ite->second[trkVec[it2]]._row;
                                _seedColumn = ite->second[trkVec[it2]]._colum;

                                _seedPpCosTh = ite->second[trkVec[it2]]._impMom3Vec.cosTheta()*180./TMath::Pi();
                                _seedPpPhi   = ite->second[trkVec[it2]]._impMom3Vec.phi()*180./TMath::Pi();

                                _seedPdgId = ite->second[trkVec[it2]]._impPdgId;
                                _seedIsGen = ite->second[trkVec[it2]]._impIsGen;
                                _seedIsConv = ite->second[trkVec[it2]]._impIsConv;
                                //                                        cout << "delta X = " << ite->second[trkVec[it2]]._impPos.getX() - ite->second[trkVec[it2]]._cluCog.getX()<<endl;
                                //                                        cout << "delta Y = " << ite->second[trkVec[it2]]._impPos.getY() - ite->second[trkVec[it2]]._cluCog.getY()<<endl;
                                //                                        cout << "delta Z = " << ite->second[trkVec[it2]]._impPos.getZ() - ite->second[trkVec[it2]]._cluCog.getZ()<<endl;



                                CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, ite->second[trkVec[it2]]._impPos);
                                _seedPu = vaneFrame.x();
                                _seedPv = vaneFrame.y();
                                _seedPw = vaneFrame.z();

                                CLHEP::Hep3Vector dirMom = ite->second[trkVec[it2]]._impMom3Vec.unit();
                                if(_diagLevel > 1){
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
                                        cout<<"posX = "<<ite->second[trkVec[it2]]._impPos.getX()<<", posY = "<<ite->second[trkVec[it2]]._impPos.getY()<<", posZ = "<< ite->second[trkVec[it2]]._impPos.getZ()<<endl;
                                        cout<<"on vane ref..."<<endl;
                                        cout<<"posU = "<<vaneFrame.getX()<<", posV = "<<vaneFrame.getY()<<", posW = "<< vaneFrame.getZ()<<endl;
                                        cout << "dirMomX = "<< dirMom.getX()<<", dirMomY = "<< dirMom.getY()<<", dirMomZ = "<<dirMom.getZ() <<endl;
                                }
                                Vane const &van = cg->vane(ite->first);
                                CLHEP::Hep3Vector const dirMom_rotated = (van.rotation())*dirMom;
                                _seedPpu = dirMom_rotated.x()*ite->second[trkVec[it2]]._impMom3Vec.mag();
                                _seedPpv = dirMom_rotated.y()*ite->second[trkVec[it2]]._impMom3Vec.mag();
                                _seedPpw = dirMom_rotated.z()*ite->second[trkVec[it2]]._impMom3Vec.mag();
                                _seedPpRotCosTh = dirMom_rotated.cosTheta()*180./TMath::Pi();
                                _seedPpRotPhi = dirMom_rotated.phi()*180./TMath::Pi();

                                if(_diagLevel > 0){
                                        cout<<"after rotation..."<<endl;
                                        cout << "dirMomU = "<< dirMom_rotated.getX()<<", dirMomV = "<< dirMom_rotated.getY()<<", dirMomW = "<<dirMom_rotated.getZ() <<endl;
                                        cout<<"cog coordinates... on vane ref!"<<endl;
                                        cout<<"cluU = "<< ite->second[trkVec[it2]]._cluCog.getX()<< ", cluV = "<<ite->second[trkVec[it2]]._cluCog.getY() <<", cluW = "<<ite->second[trkVec[it2]]._cluCog.getZ()<<endl;
                                }

                                LinePointPCA lppca(vaneFrame,dirMom_rotated,  ite->second[trkVec[it2]]._cluCog);
                                //                                        LinePointPCA lppca(ite->second[trkVec[it2]]._impPos,dirMom,  ite->second[trkVec[it2]]._cluCog);

                                CLHEP::Hep3Vector dcaV = lppca.unit();
                                double impactParam = lppca.dca();
                                _cogDca = impactParam;
                                _DCAu = dcaV.x()*impactParam;
                                _DCAv = dcaV.y()*impactParam;
                                _DCAw = dcaV.z()*impactParam;

                                _prova1 = ite->second[trkVec[it2]]._fastCryTime;

                                //                                _prova2 = sqrt( pow(dcaV.y(),2) + pow(dcaV.x(),2))*impactParam;

                                double distanceToCog = sqrt( pow(vaneFrame.getX() -ite->second[trkVec[it2]]._cluCog.getX() ,2) + pow(vaneFrame.getY() - ite->second[trkVec[it2]]._cluCog.getY() ,2) + pow(vaneFrame.getZ() - ite->second[trkVec[it2]]._cluCog.getZ() ,2) );
                                _distToCog = distanceToCog;

                                //                                        double distanceToCog = sqrt( pow(ite->second[trkVec[it2]]._impPos.getX() -ite->second[trkVec[it2]]._cluCog.getX() ,2) + pow(ite->second[trkVec[it2]]._impPos.getY() - ite->second[trkVec[it2]]._cluCog.getY() ,2) + pow(ite->second[trkVec[it2]]._impPos.getZ() - ite->second[trkVec[it2]]._cluCog.getZ() ,2) );
                                double tmpDepth = _Depth;//_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaURec->GetXaxis()->GetBinWidth(bin)*0.5;

                                _depth = tmpDepth;

                                _hTHistRecDistanceToCog->Fill(  distanceToCog);
                                _hTHistRecImpactParam->Fill(  impactParam);
                                _hTHistMomRecDotVaneNorm->Fill(  asin(impactParam / distanceToCog)*180./TMath::Pi());
                                if(_diagLevel > 0){
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
                                _hTHistDeltaEnergVRec->Fill( ite->second[trkVec[it2]]._impEnergy - ite->second[trkVec[it2]]._cluEnergy);


                                //}
                                //CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, van.getOrigin() );
                                double thetaWimpact = std::atan(-1.0*dirMom_rotated.getZ() / dirMom_rotated.getX() ) ;
                                _seedThetaW = thetaWimpact*180./TMath::Pi();

                                double deltaZ = cogWcorrFunc( _seedThetaW );
                                deltaZ =  vaneFrame.getZ() - ite->second[trkVec[it2]]._cluCog.getZ() - deltaZ ;
                                _seedDeltaW = deltaZ;


                                double corrW = triangoloVarSpread(vaneFrame.z(), 'W');
                                deltaZ -= corrW;
                                _seedDeltaCorrW = deltaZ ;

                                //double corrSpreadW = triangoloVarSpread(vaneFrame.z(), 'W');
                                deltaZ += corrW;
                                deltaZ += corrW;
                                //deltaZ -= corrSpreadW;
                                //_prova1 = deltaZ;

                                double thetaVimpact = std::atan(dirMom_rotated.getY() /  dirMom_rotated.getX() ) ;
                                _seedThetaV = thetaVimpact*180./TMath::Pi();

                                double deltaY = cogVcorrFunc( _seedThetaV );


                                deltaY = vaneFrame.getY() - ite->second[trkVec[it2]]._cluCog.getY() - deltaY ;

                                _seedDeltaV = deltaY;

                                double corrV = triangoloVarSpread(vaneFrame.y(), 'V');
                                deltaY -= corrV;
                                _seedDeltaCorrV = deltaY ;

                                //double corrSpreadV = triangoloVarSpread(vaneFrame.y(), 'V');
                                deltaY += corrV;
                                deltaY -= corrV;//corrSpreadV;
                                _prova2 = deltaY;


                                _hTHistDeltaURec->Fill(  vaneFrame.getX() - ite->second[trkVec[it2]]._cluCog.getX() );
                                _hTHistDeltaWRec->Fill(   deltaZ );
                                _hTHistDeltaVRec->Fill(   deltaY );


                                //                                        if(ite->first == 0 || ite->first == 2){
                                //
                                //                                                _hTHistDeltaURec->Fill(_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin), ite->second[trkVec[it2]]._impPos.getX() - ite->second[trkVec[it2]]._cluCog.getX() );
                                //                                        }else{
                                //                                                _hTHistDeltaVRec->Fill(_hTHistDeltaVRec->GetXaxis()->GetBinCenter(bin), ite->second[trkVec[it2]]._impPos.getY() - ite->second[trkVec[it2]]._cluCog.getY() );
                                //                                        }
                                //                                        _hTHistDeltaWRec->Fill(_hTHistDeltaWRec->GetXaxis()->GetBinCenter(bin), ite->second[trkVec[it2]]._impPos.getZ() - ite->second[trkVec[it2]]._cluCog.getZ() );

                                _Ntup->Fill();

                                trovato = true;
                                canc2 = it2;
                                canc1 = ite->first;
                        }
                        ++ite;
                }

                if(trovato){
                        std::vector<unsigned int>::iterator er = trkVec.begin();
                        er +=canc2;
                        trkVec.erase(er);
                        size2 = trkVec.size();
                        //  cout<< "before erase trkVec[] size =  "<< size2 <<endl;
                        elecMap[canc1].erase( trkVec[it2]);
                        //cout<< "after erase trkVec[] size =  "<< trkVec.size() <<endl;

                        if(elecMap[canc1].size() == 0){
                                //      cout<< "elecMap[ "<< canc1 <<" ] == 0"<<endl;
                                elecMap.erase(elecMap.find(canc1) );
                                //    cout<< "elecMap[ "<< canc1 <<" ] == erased..."<<endl;
                        }
                        ite = elecMap.begin();

                }else{
                        ++it2;
                }


        }

        //        for(ElectronMap::iterator ite = elecMap.begin(); ite != elecMap.end(); ++ite ){
        //                if(ite->second.size() == 0) continue;
        //                cout<< "start filling non good tracks"<<endl;
        //                cout<< "elecMap[.] size =  "<< ite->second.size() <<endl;
        //                for(std::map<unsigned int, electronData >::iterator it = ite->second.begin();it != ite->second.end() ;++it){
        //                        if(it == ite->second.end() ) {
        //                                cout<<"krekko!!!"<<endl;
        //                                break;
        //                        }
        //                        cout<<"start filling Ntupla..."<<endl;
        //
        //                        _evt = evt.id().event();
        //                        _seedPdgId = 0;
        //                        _seedIsGen= 0;
        //                        _clCryEnergyMaxRow = 0.0;
        //                        _clCryEnergyMaxColumn = 0.0;
        //                        cout<<"-------0.1--------"<<endl;
        //
        //                        _clVane =it->second._vane;
        //                        _clE =it->second._cluEnergy;
        //                        _clT =it->second._cluTime;
        //                        cout<<"-------0.2--------"<<endl;
        //                        _clSize =it->second._cluSize;
        //                        _clCOGu =it->second._cluCog.x();
        //                        _clCOGv =it->second._cluCog.y();
        //                        _clCOGw =it->second._cluCog.z();
        //                        cout<<"-------0.3--------"<<endl;
        //                        CLHEP::Hep3Vector Mu2eFrame = cg->fromVaneFrame(ite->first,it->second._cluCog);
        //                        _clCOGx = Mu2eFrame.x();
        //                        _clCOGy = Mu2eFrame.y();
        //                        _clCOGz = Mu2eFrame.z();
        //                        cout<<"-------0.4--------"<<endl;
        //                        _cryEdepTot =it->second._cryEnergyDepTotal;
        //                        _clCogCrySize     =it->second._clusterMap._COGcrySize;
        //                        _clCOGrow         =it->second._clusterMap._cluCogRow;
        //                        _clCOGcolumn      =it->second._clusterMap._cluCogColumn;
        //                        cout<<"-------0.5--------"<<endl;
        //                        _clCryEnergyMaxRow = -1.;//it->second._clusterMap._rowVec[0];
        //                        _clCryEnergyMaxColumn = -1.;//it->second._clusterMap._columnVec[0];
        //                        double tempEmaxCry = -1.;//it->second._clusterMap._cryEdepVec[0];
        //                        cout<<"-------1--------"<<endl;
        //                        cout<<"_clSize = "<< _clSize <<endl;
        //                        for(int i = 0; i<_clSize; ++i ){
        //                                _cryEdep[i] =it->second._clusterMap._cryEdepVec[i];//_cryEnergyDep;
        //
        //                                _clRows[i]         =it->second._clusterMap._rowVec[i];
        //                                _clColumns[i]      =it->second._clusterMap._columnVec[i];
        //
        //                                if(_cryEdep[i] > tempEmaxCry){
        //                                        _clCryEnergyMaxRow =      _clRows[i]   ;
        //                                        _clCryEnergyMaxColumn =   _clColumns[i];
        //                                        tempEmaxCry = _cryEdep[i];
        //                                }
        //                        }
        //                        cout<<"-------2--------"<<endl;
        //                        for(int j=0; j<_clCogCrySize; ++j ){
        //                                _clCogRows[j]      =it->second._clusterMap._COGrowVec[j];
        //                                _clCogColumns[j]   =it->second._clusterMap._COGcolumnVec[j];
        //                        }
        //                        _clShowerDir = it->second._clusterMap._showerDir;
        //                        _clErrShowerDir =it->second._clusterMap._errShowerDir;
        //
        //                        _seedIsGoodTrk = 1.0;
        //                        _seedNGoodTrkPerEvent = tmpV.size();
        //
        //                        _seedGenp =it->second._genp;
        //                        _seedGene =it->second._gene;
        //                        _seedGenx =it->second._genx;
        //                        _seedGeny =it->second._geny;
        //                        _seedGenz =it->second._genz;
        //                        _seedGencosth =it->second._gencosth;
        //                        _seedGenphi =it->second._genphi;
        //                        _seedGentime =it->second._gentime;
        //                        _seedPx =it->second._impPos.x();
        //                        _seedPy =it->second._impPos.y();
        //                        _seedPz =it->second._impPos.z();
        //                        _seedTime  =it->second._impTime;
        //                        _seedE =it->second._impEnergy;
        //                        _seedPpx =it->second._impMom3Vec.x();
        //                        _seedPpy =it->second._impMom3Vec.y();
        //                        _seedPpz =it->second._impMom3Vec.z();
        //                        _seedRow =it->second._row;
        //                        _seedColumn =it->second._colum;
        //                        cout<<"-------3--------"<<endl;
        //                        _seedPpCosTh =it->second._impMom3Vec.cosTheta()*180./TMath::Pi();
        //                        _seedPpPhi =it->second._impMom3Vec.phi()*180./TMath::Pi();
        //
        //                        _seedPdgId =it->second._impPdgId;
        //                        _seedIsGen =it->second._impIsGen;
        //                        //                                        cout << "delta X = " <<it->second._impPos.getX() -it->second._cluCog.getX()<<endl;
        //                        //                                        cout << "delta Y = " <<it->second._impPos.getY() -it->second._cluCog.getY()<<endl;
        //                        //                                        cout << "delta Z = " <<it->second._impPos.getZ() -it->second._cluCog.getZ()<<endl;
        //
        //
        //
        //                        CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first,it->second._impPos);
        //                        _seedPu = vaneFrame.x();
        //                        _seedPv = vaneFrame.y();
        //                        _seedPw = vaneFrame.z();
        //
        //                        CLHEP::Hep3Vector dirMom =it->second._impMom3Vec.unit();
        //                        if(_diagLevel > 1){
        //                                if(cry(ite->second[trkVec[it2]]._cluCog.y()) ==15.0 || cry(it->second._cluCog.z() ) ==15.0 ){
        //                                        cout<<"ecco il mukkio...."<<endl;
        //                                        cout<<"cogDepth.v() = " <<cry(ite->second[trkVec[it2]]._cluCog.y())<<", cogDepth.w() = "<< cry(it->second._cluCog.z() )<<", cogDepth.u() = "<<cry(it->second._cluCog.x() )<<endl;
        //                                        cout <<"iVane = "<<ite->second[trkVec[it2]]._vane<<", cluSize = "<<ite->second[trkVec[it2]]._cluSize<<", cluEdep = "<<it->second._cluEnergy<<", depth = "<<_Depth<<endl;
        //                                        //                                                        CLHEP::Hep3Vector cogDepth2 = cog_depth(*tmpCluster,_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) );
        //                                        //                                                        cout<<"instead, the linear COG give..."<<endl;
        //                                        //                                                        cout<<"cogDepth2.v() = " <<cry( cogDepth2.y() )<<", cogDepth2.w() = "<< cry( cogDepth2.z() )<<", cogDepth2.u() = "<< cry( cogDepth2.x() )<<endl;
        //                                        cout<<"fine del mukkio...."<<endl;
        //
        //
        //                                }
        //                                cout<<"vaneId = "<<ite->first<<endl;
        //                                cout<<"posX = "<<ite->second[trkVec[it2]]._impPos.getX()<<", posY = "<<ite->second[trkVec[it2]]._impPos.getY()<<", posZ = "<<it->second._impPos.getZ()<<endl;
        //                                cout<<"on vane ref..."<<endl;
        //                                cout<<"posU = "<<vaneFrame.getX()<<", posV = "<<vaneFrame.getY()<<", posW = "<< vaneFrame.getZ()<<endl;
        //                                cout << "dirMomX = "<< dirMom.getX()<<", dirMomY = "<< dirMom.getY()<<", dirMomZ = "<<dirMom.getZ() <<endl;
        //                        }
        //                        Vane const &vane = cg->getVane(ite->first);
        //                        CLHEP::Hep3Vector dirMom_rotated = *(van.getRotation())*dirMom;
        //                        _seedPpu = dirMom_rotated.x()*ite->second[trkVec[it2]]._impMom3Vec.mag();
        //                        _seedPpv = dirMom_rotated.y()*ite->second[trkVec[it2]]._impMom3Vec.mag();
        //                        _seedPpw = dirMom_rotated.z()*ite->second[trkVec[it2]]._impMom3Vec.mag();
        //                        _seedPpRotCosTh = dirMom_rotated.cosTheta()*180./TMath::Pi();
        //                        _seedPpRotPhi = dirMom_rotated.phi()*180./TMath::Pi();
        //                        cout<<"-------4--------"<<endl;
        //                        if(_diagLevel > 0){
        //                                cout<<"after rotation..."<<endl;
        //                                cout << "dirMomU = "<< dirMom_rotated.getX()<<", dirMomV = "<< dirMom_rotated.getY()<<", dirMomW = "<<dirMom_rotated.getZ() <<endl;
        //                                cout<<"cog coordinates... on vane ref!"<<endl;
        //                                cout<<"cluU = "<<it->second._cluCog.getX()<< ", cluV = "<<ite->second[trkVec[it2]]._cluCog.getY() <<", cluW = "<<ite->second[trkVec[it2]]._cluCog.getZ()<<endl;
        //                        }
        //
        //                        LinePointPCA lppca(vaneFrame,dirMom_rotated, it->second._cluCog);
        //                        //                                        LinePointPCA lppca(ite->second[trkVec[it2]]._impPos,dirMom, it->second._cluCog);
        //
        //                        CLHEP::Hep3Vector dcaV = lppca.unit();
        //                        double impactParam = lppca.dca();
        //                        _cogDca = impactParam;
        //                        _DCAu = dcaV.x()*impactParam;
        //                        _DCAv = dcaV.y()*impactParam;
        //                        _DCAw = dcaV.z()*impactParam;
        //
        //                        //                                        _prova1 = sqrt( pow(dcaV.z(),2) + pow(dcaV.x(),2))*impactParam;
        //
        //                        //                                _prova2 = sqrt( pow(dcaV.y(),2) + pow(dcaV.x(),2))*impactParam;
        //
        //                        double distanceToCog = sqrt( pow(vaneFrame.getX() -ite->second[trkVec[it2]]._cluCog.getX() ,2) + pow(vaneFrame.getY() -it->second._cluCog.getY() ,2) + pow(vaneFrame.getZ() -it->second._cluCog.getZ() ,2) );
        //                        _distToCog = distanceToCog;
        //
        //                        //                                        double distanceToCog = sqrt( pow(ite->second[trkVec[it2]]._impPos.getX() -ite->second[trkVec[it2]]._cluCog.getX() ,2) + pow(ite->second[trkVec[it2]]._impPos.getY() -it->second._cluCog.getY() ,2) + pow(ite->second[trkVec[it2]]._impPos.getZ() -it->second._cluCog.getZ() ,2) );
        //                        double tmpDepth = _Depth;//_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin) -_hTHistDeltaURec->GetXaxis()->GetBinWidth(bin)*0.5;
        //
        //                        _depth = tmpDepth;
        //
        //                        _hTHistRecDistanceToCog->Fill(  distanceToCog);
        //                        _hTHistRecImpactParam->Fill(  impactParam);
        //                        _hTHistMomRecDotVaneNorm->Fill(  asin(impactParam / distanceToCog)*180./TMath::Pi());
        //                        if(_diagLevel > 0){
        //                                cout<<"------------------------------------------------------"<<endl;
        //                                cout<< "dcaV_U = "<<dcaV.getX()<< ", dcaV_V = "<<dcaV.getY()<<", dcaV_W = "<<dcaV.getZ()<<endl;
        //                                cout << "impactParam = "<< impactParam<<endl;
        //                                cout<< "distance ImpPointToCog = "<< distanceToCog << endl;
        //                                cout<< "-----------------------------> vane = "<< ite->first<<endl;
        //                                cout << "impX = "<< vaneFrame.getX()<<endl;
        //                                cout << "impY = "<< vaneFrame.getY()<<endl;
        //                                cout << "impZ = "<< vaneFrame.getZ()<<endl;
        //                                cout<<"------------------------------------------------------"<<endl;
        //                                cout << "cluX = "<<ite->second[trkVec[it2]]._cluCog.getX()<<endl;
        //                                cout << "cluY = "<<ite->second[trkVec[it2]]._cluCog.getY()<<endl;
        //                                cout << "cluZ = "<<ite->second[trkVec[it2]]._cluCog.getZ()<<endl;
        //                                cout<<"------------------------------------------------------"<<endl;
        //                        }
        //                        // if(bin == 1){
        //                        _hTHistDeltaEnergVRec->Fill(it->second._impEnergy -it->second._cluEnergy);
        //
        //
        //                        //}
        //                        //CLHEP::Hep3Vector vaneFrame = cg->toVaneFrame(ite->first, van.getOrigin() );
        //                        double thetaWimpact = std::atan(-1.0*dirMom_rotated.getZ() / dirMom_rotated.getX() ) ;
        //                        _seedThetaW = thetaWimpact*180./TMath::Pi();
        //
        //                        double deltaZ = cogWcorrFunc(thetaWimpact );
        //                        deltaZ =  vaneFrame.getZ() -it->second._cluCog.getZ() - deltaZ ;
        //                        _seedDeltaW = deltaZ;
        //
        //
        //                        double corrW = triangoloVar(vaneFrame.z(), 'W');
        //                        deltaZ -= corrW;
        //                        _seedDeltaCorrW = deltaZ ;
        //
        //                        double corrSpreadW = triangoloVarSpread(vaneFrame.z(), 'W');
        //                        deltaZ += corrW;
        //                        deltaZ -= corrSpreadW;
        //                        _prova1 = deltaZ;
        //
        //                        double thetaVimpact = std::atan(dirMom_rotated.getY() /  dirMom_rotated.getX() ) ;
        //                        _seedThetaV = thetaVimpact*180./TMath::Pi();
        //
        //                        double deltaY = cogVcorrFunc( thetaVimpact );
        //
        //
        //                        deltaY = vaneFrame.getY() -it->second._cluCog.getY() - deltaY ;
        //
        //                        _seedDeltaV = deltaY;
        //
        //                        double corrV = triangoloVar(vaneFrame.y(), 'V');
        //                        deltaY -= corrV;
        //                        _seedDeltaCorrV = deltaY ;
        //
        //                        double corrSpreadV = triangoloVarSpread(vaneFrame.y(), 'V');
        //                        deltaY += corrV;
        //                        deltaY -= corrSpreadV;
        //                        _prova2 = deltaY;
        //                        cout<<"-------5--------"<<endl;
        //
        //                        _hTHistDeltaURec->Fill(  vaneFrame.getX() -it->second._cluCog.getX() );
        //                        _hTHistDeltaWRec->Fill(   deltaZ );
        //                        _hTHistDeltaVRec->Fill(   deltaY );
        //
        //
        //                        //                                        if(ite->first == 0 || ite->first == 2){
        //                                //
        //                                //                                                _hTHistDeltaURec->Fill(_hTHistDeltaURec->GetXaxis()->GetBinCenter(bin),it->second._impPos.getX() -it->second._cluCog.getX() );
        //                        //                                        }else{
        //                        //                                                _hTHistDeltaVRec->Fill(_hTHistDeltaVRec->GetXaxis()->GetBinCenter(bin),it->second._impPos.getY() -it->second._cluCog.getY() );
        //                        //                                        }
        //                        //                                        _hTHistDeltaWRec->Fill(_hTHistDeltaWRec->GetXaxis()->GetBinCenter(bin),it->second._impPos.getZ() -it->second._cluCog.getZ() );
        //
        //                        _Ntup->Fill();
        //                        cout<< "Ntupla filled..."<<endl;
        //
        //                }
        //        }


        //        _Ntup->Fill();

        //counting how many reconstructed clusters match with the generated particles which pass the quality cuts of the tracker


        // }//end for(_hTHistEff->GetNbinsX())
	if( evt.id().event() %1000 ==0){
	  cout << "Event "<<evt.id().event()<<" CaloClusterCogCorrFunc done..."<<endl;
	}
}

}


using mu2e::CaloClusterCogCorrFunc;
DEFINE_ART_MODULE(CaloClusterCogCorrFunc);

