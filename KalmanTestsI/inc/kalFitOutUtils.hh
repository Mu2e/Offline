//
// output utilities for reco modules
//
// $Id: kalFitOutUtils.hh,v 1.2 2012/12/04 00:51:26 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:26 $
//
#ifndef kalFitOutUtils_HH
#define kalFitOutUtils_HH

#include <iostream>
#include <string>
#include <memory>

//#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
//#include "art/Framework/Principal/Run.h"
//#include "art/Framework/Principal/SubRun.h"
//#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "ITrackerGeom/inc/ITracker.hh"
//#include "GeometryService/inc/WorldG4.hh"
//#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"

#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalRep.hh"
//#include "TrkBase/TrkRecoTrk.hh"
//#include "BaBar/PdtPid.hh"
#include "TrkBase/HelixTraj.hh"
//#include "TrkBase/TrkHotListFull.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/HelixParams.hh"
#include "TrajGeom/TrkLineTraj.hh"
#include "BField/BField.hh"

#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "KalmanTestsI/inc/kalFitDataOuts.hh"
#include "KalmanTestsI/inc/TrkCellHit.hh"
#include "KalmanTests/inc/KalFitResult.hh"

// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "MCDataProducts/inc/GenId.hh"

using namespace std; 

namespace mu2e 
{

class kalFitOutUtils{
public:
        kalFitOutUtils(std::string g4ModuleLabel="", std::string generatorModuleLabel="",
                       std::string targetStepPoints="", string hitMakerModuleLabel="",
                       int debugLvl=0, int oneturn=0,std::string vdStepPoints="",GenId genidcompare=GenId::conversionGun) :
                                _debugLvl(debugLvl),
                                _oneturn(oneturn),
                                nbadfit(0),
                                startpar(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1)),
                                recopar(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1))
        {
                _g4ModuleLabel=g4ModuleLabel;
                _generatorModuleLabel=generatorModuleLabel;
                _targetStepPoints=targetStepPoints;
                _hitMakerModuleLabel=hitMakerModuleLabel;
                _vdStepPoints=vdStepPoints;
                _genidcompare=genidcompare;
        }
        ~kalFitOutUtils(){}

        void bookHitos();
        void finalizeHistos();

        //bool FillMCInfo(art::Event& event, std::vector<hitIndex> &strawhits, HelixTraj &seed);
        bool FillMCInfo(art::Event const& event, std::vector<hitIndex> &strawhits, HelixTraj &seed);
        bool findMCData(const art::Event& evt);
        void hitsDiag(std::vector<const TrkCellHit*> const& hits);
        //void hitDiag(const TrkCellHit* strawhit);
        //void findArcs(std::vector<const TrkCellHit*> const& straws, std::vector<TrkArc>&  arcs) const;
        //static int findArc(size_t itsh,std::vector<TrkArc>& arcs );

        // DIO spectrum
        static double DIOspectrum(double ee);

        void FillHistos(KalFitResult& myfit, HelixTraj &seed, int iseed=0);

        // Label of the module that created the data products.
        std::string _g4ModuleLabel;
        // Label of the generator.
        std::string _generatorModuleLabel;
        // Instance names of data products
        std::string _targetStepPoints;
        string _hitMakerModuleLabel;
        std::string _vdStepPoints;
        GenId _genidcompare;

        int _debugLvl;
        int _oneturn;

        //std::vector<hitIndex> strawhits;
        std::vector<double> hitflt;
        double time0, time0_max;
        double s0;
        CLHEP::Hep3Vector gen_pos;
        CLHEP::Hep3Vector gen_mom;
        const BField* _bfield;

        CLHEP::Hep3Vector firstcell_pos;
        // Pointers to histograms, ntuples, TGraphs.

        int nbadfit;
        TH2D *hpull;
        TH2D *hdev;
        TH1D *hchi2;
        TH1D *hprob;
        TH1D *hchi2_0;
        TH1D *hmom;
        TH1D *hpars[5];
        TH1D *hnhits;
        TH1D *hnhitstotal;
        TH1D *hnturns;
        TH1D *heloss,*helossIWall,*helossIWallNorm;
        TH1D *hdtheta,*hdthetaturn,*hdthetaIWall,*hdthetaIWallNorm;
        TH1D *hlenturn;
        TH1D *hdedx,*hdedxturn;
        TH1D *hdistres;
        TH1D *htargeteloss;
        TTree *treefit;
        HelixParams startpar;
        HelixParams recopar;

        // cache of event data
        MCEvtData _mcdata;

        std::vector<TrkCellHitInfo> _tcellinfo;
        std::vector<HitInfo> _mchitinfo;

        RecoInfo recoinfo;
        EventInfo eventinfo;

};

}  // end namespace mu2e
#endif
