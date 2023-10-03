//
// TODO:  This file will analyze the outgoing photons energy, momentum and position by looking at the indicent electron/positron
// Original author S. Middleton
//

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryPoint.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <memory>
#include <utility>
#include <functional>
#include <float.h>
#include <vector>
#include <map>

//ROOT
#include "TStyle.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TTree.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TProfile.h"
using namespace std;
namespace mu2e{

  class PhotonAna : public art::EDAnalyzer {
	  public:
		  struct Config {
		  using Name=fhicl::Name;
		  using Comment=fhicl::Comment;
		  fhicl::Atom<int> diag{Name("diag"), Comment("Create diag histograms"),0};
	    fhicl::Atom<art::InputTag> StrawToken{Name("StrawGasStepCollection"),Comment("tag for straw gas step collection")};
		  fhicl::Atom<art::InputTag> MCToken{Name("MCTrajectoryCollection"),Comment("tag for MC trajectory collection")};
		  fhicl::Atom<art::InputTag> SimToken{Name("SimParticleCollection"),Comment("tag for Sim collection")};
  };
	  typedef art::EDAnalyzer::Table<Config> Parameters;
	  explicit PhotonAna(const Parameters& conf);
	  virtual void beginJob() override;
	  virtual void analyze(const art::Event& e);

  private:
	  Config _conf;
	  int _diagLevel;
	  art::InputTag _StrawToken;
	  art::InputTag _MCToken;
	  art::InputTag _SimToken;
	  const StrawGasStepCollection* _StrawCol;
	  const MCTrajectoryCollection* _TrajCol;
	  const SimParticleCollection* _SimCol;
	  
  //Declaring variables: Daughter particles 
	  Float_t _pdgid;
	  //TTree *_photon_analyzer;
	  TTree *_simpart_analyzer;
	  //StrawId strawid;
	  Int_t nparts;
	  
	  Float_t ElecPairMom;
	  Float_t PosPairMom;
	  Float_t TotalPairMom;
	  
	  Float_t ElecStartMomx;
	  Float_t ElecStartMomy;
	  Float_t ElecStartMomz;
	  
	  Float_t ElecStartPosx;
	  Float_t ElecStartPosy;
	  Float_t ElecStartPosz;

	  Float_t ElecTheta;
    Float_t ElecPhi;
	  Float_t ElecTransMom;
	  Float_t ElecPhotDirMom;
	  
	  Float_t PosStartMomx;
	  Float_t PosStartMomy;
	  Float_t PosStartMomz;
	  
	  Float_t PosTheta;
    Float_t PosPhi;
	  Float_t PosTransMom;
	  Float_t PosPhotDirMom;
	  
	  Float_t PosStartPosx;
	  Float_t PosStartPosy;
	  Float_t PosStartPosz;
	  	  
	  Float_t PhotStartMomx;
	  Float_t PhotStartMomy;
	  Float_t PhotStartMomz;
	  Float_t PhotTotMom;
	  
	  Float_t PhotTheta;
	  Float_t PhotPhi;
	  
	  Int_t nParts;
	  
	  unsigned int nConv = 0;
	  unsigned int nConvParticles = 0;
	  unsigned int nBoth = 0;
	  unsigned int nSameMother = 0;
	  unsigned int nSameMotherPhoton = 0;
	  };

  PhotonAna::PhotonAna(const Parameters& conf) :
  art::EDAnalyzer(conf),
  _diagLevel(conf().diag()),
  _StrawToken(conf().StrawToken()),
  _MCToken(conf().MCToken()),
  _SimToken(conf().SimToken())
  {}

  void PhotonAna::beginJob() { //TODO - can add TTree and THistograms here if required
	  art::ServiceHandle<art::TFileService> tfs;
  //1D Historgrams

	  _simpart_analyzer = tfs->make<TTree>("simpart_analyzer"," Diagnostics for Photon Conversion Track Fitting");
	  _simpart_analyzer->Branch("ElecPairMom",&ElecPairMom,"elecPairMom/F");
	  _simpart_analyzer->Branch("PosPairMom",&PosPairMom,"posPairMom/F");
	  _simpart_analyzer->Branch("TotalPairMom",&TotalPairMom,"totalPairMom/F");
	  _simpart_analyzer->Branch("Particles", &nParts, "nParts/I");
	  
	  _simpart_analyzer->Branch("ElecStartMomx",&ElecStartMomx,"ElecStartMomx/F");
	  _simpart_analyzer->Branch("ElecStartMomy",&ElecStartMomy,"ElecStartMomy/F");
	  _simpart_analyzer->Branch("ElecStartMomz",&ElecStartMomz,"ElecStartMomz/F");
	  
	  _simpart_analyzer->Branch("ElecStartPosx",&ElecStartPosx,"ElecStartPosx/F");
	  _simpart_analyzer->Branch("ElecStartPosy",&ElecStartPosy,"ElecStartPosy/F");
	  _simpart_analyzer->Branch("ElecStartPosz",&ElecStartPosz,"ElecStartPosz/F");
	  
	  _simpart_analyzer->Branch("PosStartMomx",&PosStartMomx,"PosStartMomx/F");
	  _simpart_analyzer->Branch("PosStartMomy",&PosStartMomy,"PosStartMomy/F");
	  _simpart_analyzer->Branch("PosStartMomz",&PosStartMomz,"PosStartMomz/F");
	  
	  _simpart_analyzer->Branch("PosStartPosx",&PosStartPosx,"PosStartPosx/F");
	  _simpart_analyzer->Branch("PosStartPosy",&PosStartPosy,"PosStartPosy/F");
	  _simpart_analyzer->Branch("PosStartPosz",&PosStartPosz,"PosStartPosz/F");
	  
	  _simpart_analyzer->Branch("PhotStartMomx",&PhotStartMomx,"PhotStartMomx/F");
	  _simpart_analyzer->Branch("PhotStartMomy",&PhotStartMomy,"PhotStartMomy/F");
	  _simpart_analyzer->Branch("PhotStartMomz",&PhotStartMomz,"PhotStartMomz/F");
	  _simpart_analyzer->Branch("PhotTotMom",&PhotTotMom,"PhotTotMom/F");
	  
	  _simpart_analyzer->Branch("ElecTheta",&ElecTheta,"ElecTheta/F");
	  _simpart_analyzer->Branch("ElecPhi",&ElecPhi,"ElecPhi/F");
	  _simpart_analyzer->Branch("ElecTransMom",&ElecTransMom,"ElecTransMom/F");	  
	  _simpart_analyzer->Branch("ElecPhotDirMom",&ElecPhotDirMom,"ElecPhotDirMom/F");
	   
		_simpart_analyzer->Branch("PosTheta",&PosTheta,"PosTheta/F");
	  _simpart_analyzer->Branch("PosPhi",&PosPhi,"PosPhi/F");
	  _simpart_analyzer->Branch("PosTransMom",&PosTransMom,"PosTransMom/F");	  
	  _simpart_analyzer->Branch("PosPhotDirMom",&PosPhotDirMom,"PosPhotDirMom/F");
	  
	  _simpart_analyzer->Branch("PhotTheta",&PhotTheta,"PhotTheta/F");
	  _simpart_analyzer->Branch("PhotPhi",&PhotPhi,"PhotPhi/F");
  }
  void PhotonAna::analyze(const art::Event& event) {
    
    int eventid_ = event.id().event(); 
    int runid_ = event.run();
    int subrunid_ = event.subRun();
    if(_diagLevel > 1){
      std::cout<<"+=======================================================+"<<std::endl;
      std::cout<<"event : "<<eventid_<<" subrun : "<<subrunid_<<" run : "<<runid_<<std::endl;
    }
    //------------SimParticles-------------//
    auto sH = event.getValidHandle<mu2e::SimParticleCollection>(_SimToken);
    _SimCol = sH.product();

    //std::cout<<"SimCol Size : "<<_SimCol->size()<<std::endl;
    
    // some counters and bools:
    nConvParticles = 0;
    bool isConv = false;
    bool hasElectron = false;
    bool hasPositron = false;
    bool hasMotherPhoton = false;
    double ElectronMom = 0;
    double PositronMom = 0;
    
    std::vector<double> ElecMotherPos; //TODO update plots to use these positions
    std::vector<double> PosMotherPos;

    bool hasBoth = false;
    //loop over Simparticles:
    for ( SimParticleCollection::const_iterator i=_SimCol->begin(); i!=_SimCol->end(); ++i ){
      SimParticle const& sim = i->second;
      
      //check the particle is e-/e+ and created by conversion
      if(abs(sim.pdgId()) == 11 and sim.creationCode()==13){//(sim.creationCode()==173 or sim.creationCode()==174)
          
          double TotalMomentum = sim.startMomentum().rho();
          
          // fill counters
          if(isConv ==false) nConv ++; // avoids double counting for each of the pair
          isConv = true;
          nConvParticles ++;

          //print information on daughter:
          if(_diagLevel > 1){
            std::cout<<"Daughter PDG ID : "<<sim.pdgId()<<" Creation code "<<sim.creationCode()<<std::endl;
            std::cout<<"Daughter Start Pos. : "<<sim.startPosition().x()<<" , "<<sim.startPosition().y()<<" , "<<sim.startPosition().z()<<std::endl;
            std::cout<<"Daughter End Pos. : "<<sim.endPosition().x()<<" , "<<sim.endPosition().y()<<" , "<<sim.endPosition().z()<<std::endl;
            std::cout<<"Daughter total mom : "<<TotalMomentum<<std::endl;
          }

          // parent information on mother:
          art::Ptr<SimParticle> const& parent = sim.parent();
          
          PhotStartMomx = parent->startMomentum().x();
          PhotStartMomy = parent->startMomentum().y();
          PhotStartMomz = parent->startMomentum().z();
          PhotTotMom = parent->startMomentum().rho();
          PhotTheta = parent->startMomentum().theta();
          PhotPhi = parent->startMomentum().phi() ;
          
          // information on grandmother
          /*art::Ptr<SimParticle> grandparent = parent->parent();
          if(_diagLevel > 1) std::cout<<"ancestor id "<<grandparent->parent()->pdgId()<<std::endl;
          if(!grandparent->isPrimary() and abs(grandparent->pdgId()!=13)){
            int id = grandparent->parent()->pdgId();
            std::vector<int> ancestors;
            ancestors.push_back(id);
            while(abs(id) !=13){
              id = grandparent->parent()->pdgId();
              grandparent =  grandparent->parent();
              if(_diagLevel > 1) std::cout<<"ancestor id "<<id<<std::endl;
            }
          }*/
			  
          
          if(_diagLevel > 1){
            std::cout<<"Mother PDG ID : "<<parent->pdgId()<<std::endl;
            std::cout<<"Mother End Pos. : "<<parent->endPosition().x()<<" , "<<parent->endPosition().y()<<" , "<<parent->endPosition().z()<<std::endl;
          }
          
          //check mother ID is photon
          if(parent->pdgId() == 22 and parent->creationCode()==175) hasMotherPhoton = true; //TODO keep this at ==13?
          
          // check electron and positrons and save parent positions for later checks
          if(sim.pdgId() == 11) {
            hasElectron  = true;
            ElectronMom = TotalMomentum; //TODO - add these to the TTree
            ElecMotherPos.push_back(parent->endPosition().x());
            ElecMotherPos.push_back(parent->endPosition().y());
            ElecMotherPos.push_back(parent->endPosition().z());

            ElecStartMomx  =  sim.startMomentum().x();
            ElecStartMomy =  sim.startMomentum().y();
            ElecStartMomz =  sim.startMomentum().z();

            ElecStartPosx = sim.startPosition().x();
            ElecStartPosy = sim.startPosition().y();
            ElecStartPosz = sim.startPosition().z();

            ElecTheta = sim.startMomentum().theta();
            ElecPhi = sim.startMomentum().phi();
            ElecTransMom = TotalMomentum*(sin(ElecTheta));
           
            ElecPhotDirMom = sim.startMomentum().dot(parent->startMomentum());//(ElecMom[0]*ElecMotherMom[0]+ElecMom[1]*ElecMotherMom[1]+ElecMom[2]*ElecMotherMom[2])/ElecMotherAbsMom;

           }
          if(sim.pdgId() == -11){
            hasPositron  = true;
            PositronMom = TotalMomentum;
            PosMotherPos.push_back(parent->endPosition().x());
            PosMotherPos.push_back(parent->endPosition().y());
            PosMotherPos.push_back(parent->endPosition().z());

            PosStartMomx  =  sim.startMomentum().x();
            PosStartMomy =  sim.startMomentum().y();
            PosStartMomz =  sim.startMomentum().z();

            PosStartPosx = sim.startPosition().x();
            PosStartPosy = sim.startPosition().y();
            PosStartPosz = sim.startPosition().z();

            PosTheta = sim.startMomentum().theta();
            PosPhi = sim.startMomentum().phi();
            PosTransMom = TotalMomentum*(sin(PosTheta));

            PosPhotDirMom = sim.startMomentum().dot(parent->startMomentum());

          }

        }
      }
      
      // fill counters
      if(hasElectron and hasPositron) hasBoth = true;
      if(isConv and hasBoth) nBoth ++;
      if(ElecMotherPos.size() !=0 and PosMotherPos.size() != 0 and ElecMotherPos[0] == PosMotherPos[0] and ElecMotherPos[1] == PosMotherPos[1] and ElecMotherPos[2] == PosMotherPos[2]) nSameMother++;
      if(ElecMotherPos.size() !=0 and PosMotherPos.size() != 0 and ElecMotherPos[0] == PosMotherPos[0] and ElecMotherPos[1] == PosMotherPos[1] and ElecMotherPos[2] == PosMotherPos[2] and hasMotherPhoton) nSameMotherPhoton++;
      
      // only print if its an interesting event but iterate without printing others
      if(_diagLevel > 0 and isConv and hasBoth and ElecMotherPos.size() !=0 and PosMotherPos.size() != 0 and ElecMotherPos[0] == PosMotherPos[0] and ElecMotherPos[1] == PosMotherPos[1] and ElecMotherPos[2] == PosMotherPos[2] and hasMotherPhoton and nConvParticles==2){ //nConv == 2
        if(_diagLevel == 1){ std::cout<<"+=======================================================+"<<std::endl;
        std::cout<<"event : "<<eventid_<<" subrun : "<<subrunid_<<" run : "<<runid_<<std::endl;}
        std::cout<<"number of conversion particles in this event "<<nConvParticles<<std::endl;
        std::cout<<"number of total conversions "<<nConv<<std::endl;
        std::cout<<"number of events with both e+/e-s conversions "<<nBoth<<std::endl;
        std::cout<<"number of events with both e+/e-s conversions from same mother "<<nSameMother<<std::endl;
        std::cout<<"number of events with both e+/e-s conversions from same mother and its a photon "<<nSameMotherPhoton<<std::endl;
        std::cout<<"elec mom "<<ElectronMom<<" positron mom "<<PositronMom<<" total mom "<<ElectronMom+PositronMom<<std::endl;
        ElecPairMom = ElectronMom;
        PosPairMom = PositronMom;
        TotalPairMom = ElectronMom + PositronMom;
        nParts = nConvParticles;
        _simpart_analyzer->Fill();
      }
      
        
    
    //------------ MC Trajectory --------------------//
    /*auto chH = event.getValidHandle<mu2e::MCTrajectoryCollection>(_MCToken);
	  _TrajCol = chH.product();
	  std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
	  for(unsigned int k = 0; k < _TrajCol->size(); k++){ 
		  for(trajectoryIter=_TrajCol->begin(); trajectoryIter!=_TrajCol->end(); trajectoryIter++){
		    _pdgid = trajectoryIter->first->pdgId();
		    art::Ptr<SimParticle> const& _parent = trajectoryIter->first->parent();
	      ProcessCode _creationCode = trajectoryIter->first->creationCode();
 
        //initialsation of variables
		    theta = atan((trajectoryIter->first->startMomentum().y())/trajectoryIter->first->startMomentum().x());
		    nparts = _TrajCol->size();
		    startmomentum = trajectoryIter->first->startMomentum().rho();
		    transmomentum = startmomentum*(sin(theta));
    //Fill for e+e- pairs origination from same parent particle 
		    if (_TrajCol->size()==2 && _creationCode==13 && (_pdgid==11 || _pdgid==-11) && _parent->pdgId()==22 && (transmomentum>5 || transmomentum <-5) ) {
          cout << "********************************" << k << "/ " << _TrajCol->size() << "*************************" << endl;
         
		      startmomentum_D = trajectoryIter->first->startMomentum().rho();
		      transmomentum_D = startmomentum_D*(sin(theta));
		      endmomentum_D = trajectoryIter->first->endMomentum().rho();
		      startposx_D = (trajectoryIter->first->startPosition().x())+3903;
		      startposy_D = trajectoryIter->first->startPosition().y();
		      startposz_D = trajectoryIter->first->startPosition().z();
		      startposr_D = sqrt(((startposx_D)*(startposx_D))+startposy_D*startposy_D);
		      endposx_D = (trajectoryIter->first->endPosition().x())+3903;
		      endposy_D = trajectoryIter->first->endPosition().y();
		      endposz_D = trajectoryIter->first->endPosition().z();
		      endposr_D = sqrt(((endposx_D)*(endposx_D))+endposy_D*endposy_D);
		      time_D = trajectoryIter->first->startGlobalTime();

          
          std::cout << "Code: "<< _creationCode << " ID: " << _pdgid << " Parent: " << _parent->pdgId() <<endl;
          std::cout << "Parent: " << _parent->pdgId() <<"  Start Pos X: "<< startposx_P <<" end Pos X: "<< endposx_P <<endl;
          std::cout << "ID: " << _pdgid << "  Start Pos X: "<< startposx_D << " End Pos X: " << endposx_D <<endl;
          std::cout << "ID sim particle:  " << trajectoryIter->first->id() << endl;
          std::cout << "ID parent particle:  " << trajectoryIter->first->parent()->id() << endl;
          
          //_photon_analyzer->Fill(); 
	      };
	    }
    }*/
    
  }
}//end mu2e namespace
using mu2e::PhotonAna;
DEFINE_ART_MODULE(PhotonAna)
