//
// module for the calculation of the efficiency Vs energy cluster cut and other distributions related to the efficiency
//
// $Id: CaloClusterEff_module.cc,v 1.5 2012/07/10 00:02:19 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:19 $
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
#include "art/Framework/Principal/DataViewImpl.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/VisibleGenElTrack.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

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
//#include "CaloCluster/inc/CaloClusterUtilities.hh"

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
        CLHEP::Hep3Vector _cryOrigin;

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
                _cryOrigin = other._cryOrigin;
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

class CaloClusterEff : public art::EDAnalyzer {
public:
        explicit CaloClusterEff(fhicl::ParameterSet const& pset):
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
        _rowToCanc(pset.get<int>("rowToCancel",-1)),//row index running from 1, 2, ... nCryR
        _columnToCanc(pset.get<int>("columnToCancel",-1)),//column index running from 1, 2, ... nCryZ
        _nAnalyzed(0),
        _totalcputime(0),
        _totalrealtime(0),
        _hTHistEff(0),
        _hTHistEffTrk(0),
        _hTHistGlobalEff(0),
        _hTHistGlobalEffNorm(0),
        _hTHistGlobalEffTot(0),
        _hTHistDeltaEnergy(0),
        _hTHistDeltaEnergyRec(0),
        _hTHistGlobalEffRec(0),
        _hTHistDeltaXquality(0),
        _hTHistDeltaXRec(0),
        _hTHistDeltaYquality(0),
        _hTHistDeltaYRec(0),
        _hTHistDeltaZquality(0),
        _hTHistDeltaZRec(0),
        _hTHistDeltaPmag(0),
        _hTHistDeltaPitch(0),
        _hTHistDeltaPfirst(0),
        _hTHistDeltaPlast (0),
        _hTHistEnergyClu(0),
        _hTHistEnergyCluRec(0),
        _hTHistDistrRow(0),
        _hTHistDistrColumn(0),
        _hTHistDistrRecRow(0),
        _hTHistDistrRecColumn(0),
        _hTHistDistrRowSlice(0),
        _hTHistDistrColumnSlice(0),
        _hTHistDistrRecRowSlice(0),
        _hTHistDistrRecColumnSlice(0),
        _EnergyClusterCut(pset.get<double>("energyClusterCut",60.)),//MeV
        _application(0),
        _directory(0)
        {
        }
        virtual ~CaloClusterEff() {
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

        int _rowToCanc, _columnToCanc;

        //number of analyzed events
        int _nAnalyzed;

        float _totalcputime, _totalrealtime;

        TH2D* _hTHistEff;
        TH2D* _hTHistEffTrk;
        TH1D* _hTHistGlobalEff;
        TH1D* _hTHistGlobalEffNorm;
        TH1D* _hTHistGlobalEffTot;
        TH2D* _hTHistDeltaEnergy;
        TH2D* _hTHistDeltaEnergyRec;
        TH1D* _hTHistGlobalEffRec;

        TH2D* _hTHistDeltaXquality;
        TH2D* _hTHistDeltaXRec;
        TH2D* _hTHistDeltaYquality;
        TH2D* _hTHistDeltaYRec;
        TH2D* _hTHistDeltaZquality;
        TH2D* _hTHistDeltaZRec;
        TH1D* _hTHistDeltaPmag;
        TH1D* _hTHistDeltaPitch;
        TH1D* _hTHistDeltaPfirst;
        TH1D* _hTHistDeltaPlast;
        TH2D* _hTHistEnergyClu;
        TH2D* _hTHistEnergyCluRec;

        TH2D*  _hTHistDistrRow;
        TH2D*  _hTHistDistrColumn;
        TH2D*  _hTHistDistrRecRow;
        TH2D*  _hTHistDistrRecColumn;

        TH1D*  _hTHistDistrRowSlice;
        TH1D*  _hTHistDistrColumnSlice;
        TH1D*  _hTHistDistrRecRowSlice;
        TH1D*  _hTHistDistrRecColumnSlice;

        double _EnergyClusterCut;

        double globalNtrkCut;
        double globalNtrkTot;
        double *globalVecRec;
        double *globalCaloCut;
        double *globalRecCaloCut;

        std::auto_ptr<MCCaloUtilities> CaloManager;

        bool _skipEvent;

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


void CaloClusterEff::beginJob( ) {

        cout<<"start CaloClusterEff..."<<endl;

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

void CaloClusterEff::analyze(art::Event const & evt ) {

        ++_nAnalyzed;
        ++ncalls;

        CaloClusterer c;

        art::ServiceHandle<GeometryService> geom;
        GeomHandle<Calorimeter> cg;

        if (ncalls == 1) {

                // cout << "This should be done only in the first event" << endl;
                art::ServiceHandle<art::TFileService> tfs;
                art::TFileDirectory seedInfo = tfs->mkdir("SeedInfo");
                art::TFileDirectory cog = tfs->mkdir("Cog");
                art::TFileDirectory trkInfo = tfs->mkdir("TrkInfo");
                art::TFileDirectory clusterInfo = tfs->mkdir("ClusterInfo");


                _hTHistEff           = tfs->make<TH2D>( "CaloEff", "CaloEff;EnergyLowCut [MeV];N_{rec} / N_{imp}", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., 100, 0., 6.);
                _hTHistEffTrk        = tfs->make<TH2D>( "CaloEffTrk", "CaloEffTrk;EnergyLowCut [MeV];N_{rec} / N_{imp}", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., 100, 0., 2.);
                _hTHistGlobalEff     = tfs->make<TH1D>( "CaloGlobalEffTrk", "CaloGlobalEffTrk;EnergyLowCut [MeV];#epsilon", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110.);
                _hTHistGlobalEffNorm = tfs->make<TH1D>( "CaloGlobalEffTrkNorm", "CaloGlobalEffTrkNorm;EnergyLowCut [MeV];#epsilon", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110.);
                _hTHistGlobalEffTot  = tfs->make<TH1D>( "CaloGlobalEffTrkTot", "CaloGlobalEffTrkTot;EnergyLowCut [MeV];#epsilon", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110.);

                _hTHistDeltaEnergy     = clusterInfo.make<TH2D>( "DeltaEnergy", "DeltaEnergy;EnergyLowCut [MeV];Eseed-Eclu [MeV]", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., 100, -20., 110.);
                _hTHistDeltaEnergyRec  = clusterInfo.make<TH2D>( "DeltaEnergyRec", "DeltaEnergyRec;EnergyLowCut [MeV];Eseed-Eclu [MeV]", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., 100, -20., 110.);

                _hTHistGlobalEffRec    = tfs->make<TH1D>( "CaloGlobalEffRec", "CaloGlobalEffRec;EnergyLowCut [MeV];#epsilon", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110.);


                _hTHistDeltaXquality   = cog.make<TH2D>( "DeltaXquality", "DeltaXquality;EnergyLowCut [MeV];Xseed-Xclu [mm]", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., 200, -100., 100.);
                _hTHistDeltaXRec       = cog.make<TH2D>( "DeltaXRec", "DeltaXRec;EnergyLowCut [MeV];Xseed-Xclu [mm]", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., 200, -100., 100.);
                _hTHistDeltaYquality   = cog.make<TH2D>( "DeltaYquality", "DeltaYquality;EnergyLowCut [MeV];Yseed-Yclu [mm]", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., 200, -100., 100.);
                _hTHistDeltaYRec       = cog.make<TH2D>( "DeltaYRec", "DeltaYRec;EnergyLowCut [MeV];Yseed-Yclu [mm]", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., 200, -100., 100.);
                _hTHistDeltaZquality   = cog.make<TH2D>( "DeltaZquality", "DeltaZquality;EnergyLowCut [MeV];Zseed-Zclu [mm]", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., 200, -100., 100.);
                _hTHistDeltaZRec       = cog.make<TH2D>( "DeltaZRec", "DeltaZRec;EnergyLowCut [MeV];Zseed-Zclu [mm]", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., 200, -100., 100.);

                _hTHistDeltaPmag       = trkInfo.make<TH1D>( "TrkDeltaPmag", "TrkDeltaPmag;deltaP/Pfirst [%];entries", 2000 , 0., 0.01);
                _hTHistDeltaPitch      = trkInfo.make<TH1D>( "TrkDeltaPitch","TrkDeltaPitch;deltaPitch/pitchFirst [%];entries", 2000 , 0., 0.01);

                _hTHistDeltaPfirst     = trkInfo.make<TH1D>( "TrkPfirst", "momFirst; momentum [MeV];entries", 2000 , 50., 110.0);
                _hTHistDeltaPlast      = trkInfo.make<TH1D>( "TrkPlast","momLast ; momentum [MeV];entries", 2000 , 50., 110.0);

                _hTHistEnergyClu       = clusterInfo.make<TH2D>( "EnergyClu", "EnergyClu; E_{cluster}^{cut};Energy_{cluster} [MeV]", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., 450 , 0., 150.0);
                _hTHistEnergyCluRec    = clusterInfo.make<TH2D>( "EnergyCluRec","EnergyCluRec ; E_{cluster}^{cut}; Energy_{cluster} [MeV]", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., 450 ,0., 150.0);


                _hTHistDistrRow        = seedInfo.make<TH2D>( "DistrRow", "DistrRow; E_{cluster}^{cut};row index", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., cg->nCrystalR() , 0., cg->nCrystalR() );
                _hTHistDistrColumn     = seedInfo.make<TH2D>( "DistrColumn","DistrColumn ; E_{cluster}^{cut}; column index", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110.,cg->nCrystalZ() ,0., cg->nCrystalZ() );

                _hTHistDistrRecRow     = seedInfo.make<TH2D>( "DistrRecRow", "DistrRecRow; E_{cluster}^{cut};row index", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110., cg->nCrystalR() , 0., cg->nCrystalR() );
                _hTHistDistrRecColumn  = seedInfo.make<TH2D>( "DistrRecColumn","DistrRecColumn ; E_{cluster}^{cut}; column index", TMath::Nint((110.-_EnergyClusterCut)/2.0), _EnergyClusterCut, 110.,cg->nCrystalZ() ,0., cg->nCrystalZ() );


                globalNtrkCut = 0.0;
                globalNtrkTot = 0.0;
                globalVecRec =new double [_hTHistEffTrk->GetNbinsX()];
                globalCaloCut = new double [_hTHistGlobalEffNorm->GetNbinsX()];
                globalRecCaloCut= new double [_hTHistGlobalEffRec->GetNbinsX()];
                for(int j= 0; j< _hTHistEffTrk->GetNbinsX(); ++j ){
                        globalVecRec[j] = 0.0;
                        globalCaloCut[j] = 0.0;
                        globalRecCaloCut[j] = 0.0;
                }

        }

        doCalorimeter(evt, _skipEvent);

} // end of analyze

void CaloClusterEff::endJob() {
        art::ServiceHandle<art::TFileService> tfs;
        art::TFileDirectory seedInfoSlice = tfs->mkdir("SeedInfoSlice");

        //cout<< "globalNtrkCut = "<< globalNtrkCut << ", globalNtrkTot = "<<globalNtrkTot<<" ==> globalNtrkCut / globalNtrkTot = "<<globalNtrkCut / globalNtrkTot<<endl;
        for(unsigned int b=1; b <=(unsigned int) _hTHistGlobalEff->GetNbinsX(); ++b){
                double tmpEff = globalVecRec[b-1] / globalNtrkCut;
                _hTHistGlobalEff->SetBinContent(b,  tmpEff);
                _hTHistGlobalEff->SetBinError(b, TMath::Sqrt(tmpEff*fabs(1.0-tmpEff) / globalNtrkCut ));


                double tmpEff2 = globalCaloCut[b-1] / globalNtrkCut;
                _hTHistGlobalEffNorm->SetBinContent(b,  tmpEff2);
                _hTHistGlobalEffNorm->SetBinError(b, TMath::Sqrt(tmpEff2*fabs(1.0-tmpEff2) / globalNtrkCut ));

                double tmpEff3 = globalVecRec[b-1] / globalNtrkTot;
                _hTHistGlobalEffTot->SetBinContent(b,  tmpEff3);
                _hTHistGlobalEffTot->SetBinError(b, TMath::Sqrt(tmpEff3*fabs(1.0-tmpEff3) / globalNtrkTot ));

                double tmpEff4 = globalRecCaloCut[b-1] / globalNtrkTot;
                _hTHistGlobalEffRec->SetBinContent(b,  tmpEff4);
                _hTHistGlobalEffRec->SetBinError(b, TMath::Sqrt(tmpEff4*fabs(1.0-tmpEff4) / globalNtrkTot ));

        }
        int i = 1;
        _hTHistDistrRowSlice        = seedInfoSlice.make<TH1D>(Form("DistrRowSlice_%d",i), Form("DistrRowSlice_%d",i), _hTHistDistrRow->GetNbinsY(), 0.,  _hTHistDistrRow->GetNbinsY());
        _hTHistDistrColumnSlice     = seedInfoSlice.make<TH1D>(Form("DistrColumnSlice_%d",i), Form("DistrColumnSlice_%d",i), _hTHistDistrColumn->GetNbinsY(), 0.,  _hTHistDistrColumn->GetNbinsY());

       _hTHistDistrRecRowSlice      = seedInfoSlice.make<TH1D>(Form("DistrRecRowSlice_%d",i), Form("DistrRecRowSlice_%d",i), _hTHistDistrRecRow->GetNbinsY(), 0.,  _hTHistDistrRecRow->GetNbinsY());
        _hTHistDistrRecColumnSlice  = seedInfoSlice.make<TH1D>(Form("DistrRecColumnSlice_%d",i), Form("DistrRecColumnSlice_%d",i), _hTHistDistrRecColumn->GetNbinsY(), 0.,  _hTHistDistrRecColumn->GetNbinsY());

        _hTHistDistrRowSlice        =          _hTHistDistrRow->ProjectionY(Form("DistrRowSlice_%d",i), i, i);
        _hTHistDistrColumnSlice     =         _hTHistDistrColumn->ProjectionY(Form("DistrColumnSlice_%d",i), i, i);

        _hTHistDistrRecRowSlice     =         _hTHistDistrRecRow->ProjectionY(Form("DistrRecRowSlice_%d",i), i, i);
        _hTHistDistrRecColumnSlice  =         _hTHistDistrRecColumn->ProjectionY(Form("DistrRecColumnSlice_%d",i), i, i);

}

void CaloClusterEff::doCalorimeter(art::Event const& evt, bool skip){
        //cout << "start CaloClusterEff..."<<endl;

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

        int NGoodElec [ _hTHistEff->GetNbinsX()];
        int NGoodConvEl[ _hTHistEff->GetNbinsX()];

        for(unsigned int b=1; b <=(unsigned int) _hTHistEff->GetNbinsX(); ++b){
                NGoodElec[b-1] = 0;
                NGoodConvEl[b-1] = 0;
        }

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
                GenElHitData& ldil = iEltrk.getHit((int)(iEltrk.getNumOfHit() - 1) );

                if(!findTrkId(tmpVTot, iEltrk.getTrkID().asUint() ) ){

                        //cout<<"faccio il puschback in tmpV..."<<endl;
                        tmpVTot.push_back( iEltrk.getTrkID().asUint() );
                        //cout<<"--------> puschback in tmpV eseguito..."<<endl;

                }



                double cosTheta = iEltrk.getTrkLrntzVec().cosTheta() ;

                double cosPitch = hdil._hitMomentum.cosTheta();

                double lcosPitch = ldil._hitMomentum.cosTheta();


//                cout << "cosThetafirst = "<<cosPitch <<", "<<"costhetaLast = "<<lcosPitch<<endl;
//                cout<< "ldil._hitMomentum.mag() = "<<ldil._hitMomentum.mag()<<"hdil._hitMomentum.mag() = "<< hdil._hitMomentum.mag()<<endl;


                bool condition = true;
                //the following condition are the same used by Dave Brown for TTracker studies
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

                if( condition /*&& (iEltrk.getNumOfHit() >= 20) &&  (cosTheta >= -0.5) && (cosTheta <=  0.5) && ( cosPitch > 0.5) && (cosPitch < 0.70710678118655) && (hdil._hitMomentum.mag() >= trkMomCut)*/ ){
                        NtrkCut++;// # of electrons at the entrance of the TRK which have trkMomCut MeV of kinetic energy and which will light at least 20 straws
                        if(!findTrkId(tmpV, iEltrk.getTrkID().asUint() ) ){

                                //cout<<"faccio il puschback in tmpV..."<<endl;
                                tmpV.push_back( iEltrk.getTrkID().asUint() );
                                //cout<<"--------> puschback in tmpV eseguito..."<<endl;

                        }

                        for(unsigned int b=1; b <=(unsigned int) _hTHistEff->GetNbinsX(); ++b){

                                if (hdil._hitMomentum.mag() >= (_hTHistEff->GetXaxis()->GetBinCenter(b) -  0.5*_hTHistEff->GetXaxis()->GetBinWidth(b) )){
                                        NGoodElec[b-1]++;
                                        if ( iEltrk.isConversionEl() ) {
                                                NGoodConvEl[b-1]++;
                                        }
                                }

                        }
                }



        }//end TRK mapping
        //cout <<"NtrkCut = " <<NtrkCut<<endl;

        if (NtrkTot==0) return;

        globalNtrkCut += NtrkCut;
        globalNtrkTot += NtrkTot;

        int nVane = cg->nVane();
        //cout<<"\n ----------------------> cg->roHalfThickness() = "<< cg->roHalfThickness() <<endl;
        double vecRec[nVane][_hTHistEff->GetNbinsX()];
        double vecRecNorm[nVane][_hTHistGlobalEffNorm->GetNbinsX()];

        for(int l=0; l<nVane; ++l){

                for (int ibin=1; ibin<=_hTHistEff->GetNbinsX();ibin++){
                        vecRec[l][ibin-1]=0.0;
                        vecRecNorm[l][ibin-1] = 0.0;
                }
        }

        for(int bin =1; bin <= _hTHistEff->GetNbinsX(); ++bin ){
                std::vector<unsigned int> trkVec, trkVecTot;
                for(unsigned int f = 0; f != tmpV.size(); ++f){
                        trkVec.push_back(tmpV[f]);
                }
                for(unsigned int f = 0; f != tmpVTot.size(); ++f){
                        trkVecTot.push_back(tmpVTot[f]);
                }
                ElectronMap elecMap;

                //start reading the clusters that I reconstructed from the generated particles
                if(caloClusters->size()>0 ){
                        int iVane;
                        for(size_t icl=0; icl<caloClusters->size(); ++icl){

                                double eDepClu = 0.;

                                CaloCluster const& clu = (*caloClusters).at(icl);
                                eDepClu = clu.energyDep();
                                iVane = clu.vaneId();

                                CaloCrystalHitPtrVector caloClusterHits = clu.caloCrystalHitsPtrVector();

                                if(eDepClu >= (_hTHistEff->GetXaxis()->GetBinCenter(bin) -  0.5*_hTHistEff->GetXaxis()->GetBinWidth(bin) ) ){
                                        vecRec[iVane][bin - 1]+=1.0;

                                        for(size_t i=0; i<caloClusterHits.size(); ++i){
                                                CaloCrystalHit const& hit = *(caloClusterHits.at(i));

                                                std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                                                if(ROIds.size()<1 ) continue;

                                                CaloHit const& thehit = *ROIds.at(0);
                                                size_t collectionPosition = ROIds.at(0).key();

                                                PtrStepPointMCVector const & mcptr(hits_mcptr->at(collectionPosition));
                                                size_t nHitsPerCrystal = mcptr.size();


                                                for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {

                                                        //cout<< "Start loop..."<< "j2 = "<< j2<<endl;

                                                        StepPointMC const& mchit = *mcptr[j2];

                                                        // The simulated particle that made this hit.
                                                        SimParticleCollection::key_type trackId(mchit.trackId());

                                                        CLHEP::Hep3Vector cryFrame = cg->toCrystalFrame(thehit.id(), mchit.position());

                                                        if( (cg->getCrystalRByRO(thehit.id()) != ( _rowToCanc - 1.0 ) && cg->getCrystalZByRO(thehit.id()) != ( _columnToCanc - 1.0 ) ) ){
                                                                if(elecMap[iVane][trackId.asUint()]._impTime > mchit.time() ){
                                                                        elecMap[iVane][trackId.asUint()]._cluEnergy = clu.energyDep();
                                                                        elecMap[iVane][trackId.asUint()]._cluTime = clu.time();
                                                                        elecMap[iVane][trackId.asUint()]._impTime = mchit.time();
                                                                        elecMap[iVane][trackId.asUint()]._impEnergy = mchit.momentum().mag();
                                                                        elecMap[iVane][trackId.asUint()]._cluCog = clu.cog3Vector();
                                                                        elecMap[iVane][trackId.asUint()]._impPos = mchit.position();
                                                                        elecMap[iVane][trackId.asUint()]._impPosCryFrame = cryFrame;
                                                                        elecMap[iVane][trackId.asUint()]._row    = cg->getCrystalRByRO(thehit.id() );
                                                                        elecMap[iVane][trackId.asUint()]._column    = cg->getCrystalZByRO(thehit.id() );
                                                                        elecMap[iVane][trackId.asUint()]._cryOrigin  = cg->getCrystalOriginByRO( thehit.id() );

                                                                        //                                                                cout<< "###################"<<endl;
                                                                        //                                                                cout<< "idVande = "<< iVane<<endl;
                                                                        //                                                                cout << "cluX = "<<elecMap[iVane][trackId.asUint()]._impPos.getX()<<endl;
                                                                        //                                                                cout << "cluY = "<<elecMap[iVane][trackId.asUint()]._impPos.getY()<<endl;
                                                                        //                                                                cout << "cluZ = "<<elecMap[iVane][trackId.asUint()]._impPos.getZ()<<endl;
                                                                        //                                                                cout<< "###################"<<endl;


                                                                }
                                                        }
                                                }

                                        }// esco da questo ciclo con: time e Id della particella che per prima ha dato origine al cluster

                                }// end if(eDepClu >= (_hTHistEff->GetXaxis()->GetBinCenter(bin) -  0.5*_hTHistEff->GetXaxis()->GetBinWidth(bin) ) ){

                        }//end for(caloClsuters.size)
                }//end if(caloClsuters.size() > 0)

                //counting how many reconstructed clusters match with the generated particles
                unsigned int canc2 = 0;
                //cout <<"-------------------> numero di elettroni generati = "<< trkVecTot.size() <<endl;
                unsigned int size2 = trkVecTot.size();
                unsigned int it2=0;
                while( it2 < size2){
                        ElectronMap::iterator ite = elecMap.begin();
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
                                        _hTHistEnergyCluRec->Fill(_hTHistEnergyCluRec->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._cluEnergy);

                                        globalRecCaloCut[bin - 1] += 1.0;

                                        _hTHistDistrRecRow->Fill(_hTHistDistrRecRow->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._row );
                                        _hTHistDistrRecColumn->Fill(_hTHistDistrRecColumn->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._column );


                                        _hTHistDeltaEnergyRec->Fill(_hTHistDeltaEnergyRec->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._impEnergy - ite->second[trkVecTot[it2]]._cluEnergy);

                                        //the following lines were used ti verify that I used the correct translations in the file "CaloClustertilities.hh"
                                        //double cryHalfLength = cg->crystalHalfLength();
                                        //double cryHalfSize = cg->crystalHalfSize();
                                        //double deltaYc = 0.0, deltaXc = 0.0;
                                        //                                        if( ite->first == 3){
                                        //                                                // = (deltaR )*cryHalfSize*2.;
                                        //                                                //deltaZc = (deltaZ )*cryHalfSize*2.;
                                        //                                                deltaXc = cryHalfLength + 2.0*cg->roHalfThickness();
                                        //
                                        //                                                //                                                       resError.setX(defaultError);
                                        //                                                //                                                       resError.setY( RMScryHalfSize/TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);
                                        //                                        }
                                        //                                        else if(ite->first == 2){
                                        //                                                //deltaXc = (deltaR )*cryHalfSize*2.;
                                        //                                                //deltaZc = (deltaZ)*cryHalfSize*2.;
                                        //                                                deltaYc = -1.*cryHalfLength - 2.0*cg->roHalfThickness();
                                        //
                                        //                                                //                                                       resError.setY(defaultError);
                                        //                                                //                                                       resError.setX(RMScryHalfSize/TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);
                                        //
                                        //                                        }else if(ite->first == 1){
                                        //                                                //deltaYc = -1.*(deltaR)*cryHalfSize*2. ;
                                        //                                                //deltaZc = (deltaZ)*cryHalfSize*2.;
                                        //                                                deltaXc = -1.*cryHalfLength - 2.0*cg->roHalfThickness();
                                        //
                                        //                                                //                                                       resError.setX(defaultError);
                                        //                                                //                                                       resError.setY(RMScryHalfSize/TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);
                                        //                                        }else if(ite->first == 0){
                                        //                                                //deltaXc = -1.0*(deltaR)*cryHalfSize*2.;
                                        //                                                //deltaZc = (deltaZ)*cryHalfSize*2.;
                                        //                                                deltaYc = cryHalfLength + 2.0*cg->roHalfThickness();
                                        //
                                        //                                                //                                                       resError.setY(defaultError);
                                        //                                                //                                                       resError.setX(RMScryHalfSize/TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);
                                        //                                        }

                                        if(ite->first == 0 || ite->first == 2){

                                                _hTHistDeltaXRec->Fill(_hTHistDeltaXRec->GetXaxis()->GetBinCenter(bin),/*ite->second[trkVecTot[it2]]._impPos.getX() - ite->second[trkVecTot[it2]]._cryOrigin.getX() - deltaXc);*/ ite->second[trkVecTot[it2]]._impPos.getX() - ite->second[trkVecTot[it2]]._cluCog.getX() );
                                        }else{
                                                _hTHistDeltaYRec->Fill(_hTHistDeltaYRec->GetXaxis()->GetBinCenter(bin),/*ite->second[trkVecTot[it2]]._impPos.getY() - ite->second[trkVecTot[it2]]._cryOrigin.getY() -deltaYc);*/ ite->second[trkVecTot[it2]]._impPos.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() );
                                        }
                                        _hTHistDeltaZRec->Fill(_hTHistDeltaZRec->GetXaxis()->GetBinCenter(bin),/*ite->second[trkVecTot[it2]]._impPos.getZ() - ite->second[trkVecTot[it2]]._cryOrigin.getZ());*/ite->second[trkVecTot[it2]]._impPos.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() );
                                        //}

                                        //std::map<unsigned int, electronData >::iterator jo = ite->second  ;
                                        //cout<< "cervavo "<< trkVec[it]<<", trovato "<< jo->first <<endl;
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
                                if(size2 == 0) break;
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
                        ElectronMap::iterator ite = elecMap.begin();
                        bool trovato = false;
                        while(!trovato && ite!=elecMap.end() ){
                                if(ite->second.find(trkVec[it]) != ite->second.end()){
                                        //                                        cout << "impX = "<<ite->second[trkVec[it]]._impPosCryFrame.getX()<<endl;
                                        //                                        cout << "impY = "<<ite->second[trkVec[it]]._impPosCryFrame.getY()<<endl;
                                        //                                        cout << "impZ = "<<ite->second[trkVec[it]]._impPosCryFrame.getZ()<<endl;
                                        vecRecNorm[(ite->first)][bin - 1]+=1.0;
                                        globalCaloCut[bin - 1] += 1.0;

                                        _hTHistDistrRow->Fill(_hTHistDistrRow->GetXaxis()->GetBinCenter(bin), ite->second[trkVec[it]]._row );
                                        _hTHistDistrColumn->Fill(_hTHistDistrColumn->GetXaxis()->GetBinCenter(bin), ite->second[trkVec[it]]._column );


                                        _hTHistEnergyClu->Fill(_hTHistEnergyClu->GetXaxis()->GetBinCenter(bin), ite->second[trkVec[it]]._cluEnergy);

                                        _hTHistDeltaEnergy->Fill(_hTHistDeltaEnergy->GetXaxis()->GetBinCenter(bin), ite->second[trkVec[it]]._impEnergy - ite->second[trkVec[it]]._cluEnergy);

                                        trovato = true;
                                        canc = it;
                                        //the following lines were used ti verify that I used the correct translations in the file "CaloClustertilities.hh"

                                        //cout<<"$$$ yes! the calo works! $$$"<<endl;
                                        //double cryHalfLength = cg->crystalHalfLength();
                                        //double cryHalfSize = cg->crystalHalfSize();
                                        //double deltaYc = 0.0, deltaXc = 0.0;
                                        //                                        if( ite->first == 3){
                                        //                                                // = (deltaR )*cryHalfSize*2.;
                                        //                                                //deltaZc = (deltaZ )*cryHalfSize*2.;
                                        //                                                deltaXc = cryHalfLength + 2.0*cg->roHalfThickness();
                                        //
                                        //                                                //                                                       resError.setX(defaultError);
                                        //                                                //                                                       resError.setY( RMScryHalfSize/TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);
                                        //                                        }
                                        //                                        else if(ite->first == 2){
                                        //                                                //deltaXc = (deltaR )*cryHalfSize*2.;
                                        //                                                //deltaZc = (deltaZ)*cryHalfSize*2.;
                                        //                                                deltaYc = -1.*cryHalfLength - 2.0*cg->roHalfThickness();
                                        //
                                        //                                                //                                                       resError.setY(defaultError);
                                        //                                                //                                                       resError.setX(RMScryHalfSize/TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);
                                        //
                                        //                                        }else if(ite->first == 1){
                                        //                                                //deltaYc = -1.*(deltaR)*cryHalfSize*2. ;
                                        //                                                //deltaZc = (deltaZ)*cryHalfSize*2.;
                                        //                                                deltaXc = -1.*cryHalfLength - 2.0*cg->roHalfThickness();
                                        //
                                        //                                                //                                                       resError.setX(defaultError);
                                        //                                                //                                                       resError.setY(RMScryHalfSize/TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);
                                        //                                        }else if(ite->first == 0){
                                        //                                                //deltaXc = -1.0*(deltaR)*cryHalfSize*2.;
                                        //                                                //deltaZc = (deltaZ)*cryHalfSize*2.;
                                        //                                                deltaYc = cryHalfLength + 2.0*cg->roHalfThickness();
                                        //
                                        //                                                //                                                       resError.setY(defaultError);
                                        //                                                //                                                       resError.setX(RMScryHalfSize/TMath::Sqrt(cluster._nCrystal) * TMath::Sqrt(tmpEq)/cluster.energyDep);
                                        //                                        }
                                        if(ite->first == 0 || ite->first == 2){
                                                _hTHistDeltaXquality->Fill(_hTHistDeltaXquality->GetXaxis()->GetBinCenter(bin),/* ite->second[trkVec[it]]._impPos.getX() - ite->second[trkVec[it]]._cryOrigin.getX() - deltaXc);*/ite->second[trkVec[it]]._impPos.getX() - ite->second[trkVec[it]]._cluCog.getX() );
                                        }else{
                                                _hTHistDeltaYquality->Fill(_hTHistDeltaYquality->GetXaxis()->GetBinCenter(bin),/*ite->second[trkVec[it]]._impPos.getY() - ite->second[trkVec[it]]._cryOrigin.getY() - deltaYc );*/ite->second[trkVec[it]]._impPos.getY() - ite->second[trkVec[it]]._cluCog.getY() );
                                        }
                                        _hTHistDeltaZquality->Fill(_hTHistDeltaZquality->GetXaxis()->GetBinCenter(bin),/* ite->second[trkVec[it]]._impPos.getZ() - ite->second[trkVec[it]]._cryOrigin.getZ());*/ite->second[trkVec[it]]._impPos.getZ() - ite->second[trkVec[it]]._cluCog.getZ() );

                                }
                                ++ite;
                        }

                        if(trovato){
                                std::vector<unsigned int>::iterator er = trkVec.begin();
                                er +=canc;
                                trkVec.erase(er);
                                size = trkVec.size();
                                if(size2 == 0) break;
                        }else{
                                ++it;
                        }

                }

                //                cout <<"after.... trkVec.size = "<< trkVec.size() <<endl;
                //                cout<<"...\n...\n "<<endl;

                //_hTHistEff->SetBinContent(bin, NrecTot);

        }//end for(_hTHistEff->GetNbinsX())

        for(unsigned int b=1; b <=(unsigned int) _hTHistEff->GetNbinsX(); ++b){
                double tmpN=0.0;
                double tmpN2 = 0.0;
                for(int v = 0; v < nVane; ++v){
                        tmpN +=vecRec[v][b-1];
                        tmpN2 += vecRecNorm[v][b-1];
                        //cout<< "------> tmpN = "<< tmpN<<endl;
                }
                globalVecRec[b-1] += tmpN;

                //                cout<<"@@@@@@@@@@@@@"<<endl;
                //                cout <<"RecCaloCut[ "<< b-1 << " ]"<< " = "<< tmpN <<", NtrkCut = "<<NtrkCut <<endl;
                //                cout <<"CaloCut[ "<< b-1 << " ]"<< " = "<< tmpN2 <<", NtrkCut = "<<NtrkCut <<endl;
                //                cout << "globalCaloCut[ "<< b-1<<" ] = " <<globalCaloCut[b-1]<< ", globalNtrkCut = "<< globalNtrkCut<<endl;
                //                cout << "globalVecRec[ "<< b-1<<" ] = " <<globalVecRec[b-1]<< ", globalNtrkCut = "<< globalNtrkCut<<endl;

                _hTHistEff->Fill(_hTHistEff->GetXaxis()->GetBinCenter(b), tmpN / (double)NtrkCut );

                _hTHistEffTrk->Fill(_hTHistEffTrk->GetXaxis()->GetBinCenter(b), tmpN2 / (double)NtrkCut );


        }
        cout << "Event "<<evt.id().event()<<" CaloClusterEff done..."<<endl;

}

}


using mu2e::CaloClusterEff;
DEFINE_ART_MODULE(CaloClusterEff);

