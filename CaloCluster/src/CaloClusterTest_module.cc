//
// Visualization of the clusters on each vane with the relative distributions of time, energy and size
//
// $Id: CaloClusterTest_module.cc,v 1.2 2012/03/07 18:00:38 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/03/07 18:00:38 $
//
// Original author G. Pezzullo & G. Tassielli
//


//Mu2e includes
#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "CaloCluster/inc/CaloClusterer.hh"

//Root includes
#include "TFile.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"

//art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"

//other includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"

//c++ includes
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include "TMath.h"



using namespace std;

namespace mu2e {

static int ncalls(0);

class CaloClusterTest : public art::EDAnalyzer {
public:
        explicit CaloClusterTest(fhicl::ParameterSet const& pset):
        _diagLevel(pset.get<int>("diagLevel",0)),
        _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
        _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
        _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
        _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
        _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel", "makeCaloCluster")),
        _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
        _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "time")),
        _producerName("Algo"+mu2e::TOUpper(_caloClusterAlgorithm)+"SeededBy"+mu2e::TOUpper(_caloClusterSeeding)),
        _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)), // MeV
        _EnergyClusterCut(pset.get<double>("energyClusterCut",60.)),//MeV
        _nAnalyzed(0),
        _hTHStackMap(0),
        _hTHStackE(0),
        _hTHStackT(0),
        _hTCanvasMap(0),
        _hTCanvasE(0),
        _hTCanvasT(0),
        _fakeCanvas(0),
        _application(0),
        _directory(0)
        {
        }
        virtual ~CaloClusterTest() {
                if (_fakeCanvas)        delete _fakeCanvas;
        }
        virtual void beginJob();
        virtual void endJob();

        void analyze(art::Event const& e );

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
        const string _producerName;

        double _minimumEnergy;
        double _EnergyClusterCut;


        //number of analyzed events
        int _nAnalyzed;

        THStack* _hTHStackMap;
        THStack* _hTHStackE;
        THStack* _hTHStackT;
        TCanvas* _hTCanvasMap;
        TCanvas* _hTCanvasE;
        TCanvas* _hTCanvasT;

        std::auto_ptr<MCCaloUtilities> CaloManager;

        bool _skipEvent;

        TCanvas*      _fakeCanvas;

        // The job needs exactly one instance of TApplication.  See note 1.
        auto_ptr<TApplication> _application;

        // Save directory from beginJob so that we can go there in endJob. See note 3.
        TDirectory* _directory;


        std::vector<TCanvas*> vecTCAN;
        std::vector<TCanvas*> vecTCANclu;

};



void CaloClusterTest::beginJob( ) {

        CaloManager = auto_ptr<MCCaloUtilities>(new MCCaloUtilities());

        // If needed, create the ROOT interactive environment. See note 1.
        if ( !gApplication ){
          int    tmp_argc(0);
          char** tmp_argv(0);
          _application = auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
        }

        gStyle->SetPalette(1);
        gROOT->SetStyle("Plain");

        _fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);

        _directory = gDirectory;

}

void CaloClusterTest::analyze(art::Event const& evt ) {

        ++_nAnalyzed;
        ++ncalls;

        art::ServiceHandle<GeometryService> geom;

        doCalorimeter(evt, _skipEvent);

} // end of analyze

void CaloClusterTest::endJob() {


}

void CaloClusterTest::doCalorimeter(art::Event const& evt, bool skip){

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

        art::Handle<CaloClusterCollection> caloClusters;
        evt.getByLabel(_caloClusterModuleLabel,_producerName,caloClusters );

        int nVane = cg->nVane();
        int nCryR = cg->nCrystalR();   //starts from 1
        int nCryZ = cg->nCrystalZ();   //starts from 1

        std::vector<TH2D*> vecMap;
        std::vector<TH1D*> vecE;
        std::vector<TH1D*> vecT;

        std::vector<TH1D*> vecSizeC;//gio
        std::vector<TH1D*> vecEneC;//gio

        std::vector<THStack*> vecTHSMap;
        std::vector<THStack*> vecTHSE;
        std::vector<THStack*> vecTHST;

        std::vector<THStack*> vecTHSSizeC;//gio
        std::vector<THStack*> vecTHSEneC;//gio

        std::vector<TH2D*> vecVaneMap;
        std::vector<TH1D*> vecVaneE;
        std::vector<TH1D*> vecVaneT;

        for(int p=0; p<nVane; ++p){

                vecTHSMap.push_back(new THStack(Form("clustStackMap_vane%i",p), "ClusterHitMap"));
                vecTHSE.push_back(new THStack(Form("clustStackE_vane%i",p), "histE"));
                vecTHST.push_back(new THStack(Form("clustStackT_vane%i",p), "histT"));

                vecTHSSizeC.push_back(new THStack(Form("clustStackSizeC_vane%i",p), "histSizeC"));//gio
                vecTHSEneC.push_back(new THStack(Form("clustStackEneC_vane%i",p), "histEneC"));//gio


                if (ncalls == 1) {
                        vecTCAN.push_back(new TCanvas(Form("canv_vane%i",p),Form("canv_vane%i",p)));
                        vecTCAN.at(p)->Divide(3,2);
                }
                vecVaneMap.push_back(new TH2D(Form("histMap_vane%i",p), "histMap", nCryZ+1, 0., nCryZ, nCryR, 0., nCryR));
                vecVaneE.push_back(new TH1D(Form("histE_vane%i",p), "histE", 500, 0., 110.));
                vecVaneT.push_back(new TH1D(Form("histT_vane%i",p), "histT", 500, 0., 3000.));

        }


        if(caloClusters->size()>0 ){
                int colorCode = 2;
                int iVane;

                for(size_t icl=0; icl<caloClusters->size(); ++icl){
                        vecMap.push_back( new TH2D(Form("histMap_cluster%i",(int)icl), "ClusterHitMap", nCryZ+1, 0., nCryZ, nCryR, 0., nCryR));
                        vecE.push_back( new TH1D(Form("histE_cluster%i",(int)icl), "histE", 500, 0., 110.));
                        vecT.push_back( new TH1D(Form("histT_cluster%i",(int)icl), "histT", 500, 0., 3000.));

                        vecEneC.push_back( new TH1D(Form("histEneC_cluster%i",(int)icl), "histEneC",500,  0., 110.));//gio
                        vecSizeC.push_back( new TH1D(Form("histSizeC_cluster%i",(int)icl), "histSizeC", 500, 0., 100.));//gio

                        std::vector<TH2D*>::iterator itMap = vecMap.end();
                        itMap--;

                        std::vector<TH1D*>::iterator itE = vecE.end();
                        itE--;

                        std::vector<TH1D*>::iterator itT = vecT.end();
                        itT--;

                        std::vector<TH1D*>::iterator itEneC = vecEneC.end();
                        itEneC--;

                        std::vector<TH1D*>::iterator itSizeC = vecSizeC.end();
                        itSizeC--;

                        int cryZ, cryR;

                        double eDepClu = 0.;
                        double eDepCry = 0.;
                        double eDepTotCry = 0.;
                        double timeCry = 0.;

                        double eDepClu_check = 0.;
                        double eDepTotClu_check = 0.;
                        double timeClu_check = 0.;

                        CaloCluster const& clu = (*caloClusters).at(icl);
                        eDepClu = clu.energyDep;
                        iVane = clu.vaneId;

                        if(eDepClu < _EnergyClusterCut) continue;

                        for(size_t icry = 0; icry <  clu.caloCrystalHitsPtrVector.size(); ++icry){

                                eDepCry = clu.caloCrystalHitsPtrVector.at(icry)->energyDep();
                                eDepTotCry = clu.caloCrystalHitsPtrVector.at(icry)->energyDepTotal();
                                timeCry = clu.caloCrystalHitsPtrVector.at(icry)->time();

                                eDepClu_check += eDepCry;
                                eDepTotClu_check += eDepTotCry;
                                timeClu_check += timeCry;

                                CaloCrystalHit const& hitClu = *(clu.caloCrystalHitsPtrVector.at(icry));
                                std::vector<art::Ptr<CaloHit> > const& ROIdsClu = hitClu.readouts();
                                CaloHit const& thehitClu = *ROIdsClu.at(0);

                                cryZ = cg->getCrystalZByRO(thehitClu.id());//get z-coordinate (from 0 to nCryZ-1)
                                cryR = cg->getCrystalRByRO(thehitClu.id());//get r-coordinate (from 0 to nCryR-1)

                                (*itMap)->Fill(cryZ, cryR);
                                (*itE)->Fill(eDepCry);
                                (*itT)->Fill(timeCry);
                        }
                        timeClu_check /= clu.caloCrystalHitsPtrVector.size();


                        (*itEneC)->Fill(eDepClu);
                        (*itSizeC)->Fill(clu.caloCrystalHitsPtrVector.size());


                        if(colorCode == 10){
                                ++colorCode;
                        }
                        (*itMap)->SetFillColor(colorCode);
                        (*itE)->SetFillColor(colorCode);
                        (*itT)->SetFillColor(colorCode);

                        (*itEneC)->SetFillColor(colorCode);
                        (*itSizeC)->SetFillColor(colorCode);

                        ++colorCode;

                        vecTHSMap.at(iVane)->Add(*itMap,"lego1");
                        vecTHSE.at(iVane)->Add(*itE);
                        vecTHST.at(iVane)->Add(*itT);

                        vecTHSEneC.at(iVane)->Add(*itEneC);
                        vecTHSSizeC.at(iVane)->Add(*itSizeC);



                }//end for(caloClsuters.size())
        }//end if( caloClsuters.size() > 0 )

        if (!( caloHits.isValid())) {
                return;
        }

        if (!caloCrystalHits.isValid()) {
                // cout << "NO CaloCrystalHits" << endl;
                return;
        }

        if(caloCrystalHits->size()>0 ){
                for(size_t i=0; i<caloCrystalHits->size(); ++i){
                        CaloCrystalHit const& hit = (*caloCrystalHits).at(i);

                        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                        if(ROIds.size()<1 ) continue;

                        if(hit.energyDep() < _minimumEnergy) continue;

                        CaloHit const& thehit = *ROIds.at(0);
                        int idVane = cg->getVaneByRO(thehit.id());
                        int cryZ = cg->getCrystalZByRO(thehit.id());//get z-coordinate (from 0 to nCryZ-1)
                        int cryR = cg->getCrystalRByRO(thehit.id());//get r-coordinate (from 0 to nCryR-1)
                        double edep    = hit.energyDepTotal();
                        double Time = hit.time();
                        vecVaneMap.at(idVane)->Fill(cryZ, cryR);
                        vecVaneE.at(idVane)->Fill(edep);
                        vecVaneT.at(idVane)->Fill(Time);
                }//end caloCrystalHits loop
        }

        for(int p=0; p<nVane; ++p){
                vecTCAN.at(p)->cd(1);
                vecTHSMap.at(p)->Draw("lego1");
                vecTCAN.at(p)->cd(2);
                vecVaneMap.at(p)->SetFillColor(0);
                vecVaneMap.at(p)->GetZaxis()->SetRangeUser(0,vecVaneMap.at(p)->GetMaximumStored());
                vecVaneMap.at(p)->Draw("lego1");
                vecTCAN.at(p)->cd(3);
                vecVaneE.at(p)->Draw();
                vecTHSE.at(p)->Draw("same");
                vecTCAN.at(p)->cd(4);
                vecVaneT.at(p)->Draw();
                vecTHST.at(p)->Draw("same");

                vecTCAN.at(p)->cd(5);
                vecTHSEneC.at(p)->Draw();
                vecTCAN.at(p)->cd(6);
                vecTHSSizeC.at(p)->Draw();

                vecTCAN.at(p)->Modified();
                vecTCAN.at(p)->Update();
        }

        cerr << "Double click in the canvas_Fake to continue:" ;
        _fakeCanvas->cd();
        TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d",evt.id().event()));
        printEvN->SetTextFont(62);
        printEvN->SetTextSizePixels(180);
        printEvN->Draw();
        _fakeCanvas->Update();
        _fakeCanvas->WaitPrimitive();
        cerr << endl;
        delete printEvN;

        for(int p=0; p<nVane; ++p){

                if (vecTHSMap.at(p)!=0x0) delete vecTHSMap.at(p);
                if (vecTHSE.at(p)!=0x0) delete vecTHSE.at(p);
                if (vecTHST.at(p)!=0x0) delete vecTHST.at(p);
                if (vecVaneMap.at(p)!=0x0) delete vecVaneMap.at(p);
                if (vecVaneE.at(p)!=0x0) delete vecVaneE.at(p);
                if (vecVaneT.at(p)!=0x0) delete vecVaneT.at(p);

                if (vecTHSEneC.at(p)!=0x0) delete vecTHSEneC.at(p);
                if (vecTHSSizeC.at(p)!=0x0) delete vecTHSSizeC.at(p);

        }

}

}


using mu2e::CaloClusterTest;
DEFINE_ART_MODULE(CaloClusterTest);

