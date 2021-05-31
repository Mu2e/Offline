//
// An EDAnalyzer module that reads back the hits created by the calorimeter and produces an ntuple
//
//
// Original author Bertrand Echenard
//

#include "CLHEP/Units/SystemOfUnits.h"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/CaloShowerStep.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/GenEventCount.hh"


#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <vector>



namespace mu2e {


  class CaloNeutron : public art::EDAnalyzer {

     public:

       explicit CaloNeutron(fhicl::ParameterSet const& pset);
       virtual ~CaloNeutron() { }

       virtual void beginJob();
       virtual void endJob();
       virtual void endSubRun(const art::SubRun& sr);

       // This is called for each event.
       virtual void analyze(const art::Event& e);





     private:

       typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;

       std::string g4ModuleLabel_;
       std::string calorimeterStepPoints_;
       std::string calorimeterROStepPoints_;
       std::string calorimeterROCardStepPoints_;
       std::string calorimeterCrateStepPoints_;
       std::string virtualDetectorLabel_;
       SimParticleTimeOffset toff_;  
       int diagLevel_;
       int nProcess_;
       int numEvents_;
       TTree* Ntup_;
       TH1D *hEvents_;
       double sumEcrys_;

       int   _evt,_run;

       int   _nVd,_vdId[16384],_vdPdgId[16384],_vdenIdx[16384];
       float _vdTime[16384],_vdPosX[16384],_vdPosY[16384],_vdPosZ[16384],_vdMom[16384],_vdMomX[16384],_vdMomY[16384],_vdMomZ[16384];
       
       int   _nCry,_cryId[16384];
       float _cryEdep[16384],_cryTime[16384],_cryPosX[16384],_cryPosY[16384],_cryPosZ[16384];

       int   _nRO,_ROId[16384];
       float _ROEdep[16384],_ROPosX[16384],_ROPosY[16384],_ROPosZ[16384]; 

       int   _nROCard,_ROCardId[16384];
       float _ROCardEdep[16384],_ROCardPosX[16384],_ROCardPosY[16384],_ROCardPosZ[16384]; 

       int   _nCrate,_CrateId[16384];
       float _CrateEdep[16384],_CratePosX[16384],_CratePosY[16384],_CratePosZ[16384]; 


  };


  CaloNeutron::CaloNeutron(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    g4ModuleLabel_(pset.get<std::string>("g4ModuleLabel")),
    calorimeterStepPoints_( pset.get<std::string>("calorimeterStepPoints") ),
    calorimeterROStepPoints_( pset.get<std::string>("calorimeterROStepPoints") ),
    calorimeterROCardStepPoints_( pset.get<std::string>("calorimeterROCardStepPoints") ),
    calorimeterCrateStepPoints_( pset.get<std::string>("calorimeterCrateStepPoints") ),
    virtualDetectorLabel_(pset.get<std::string>("virtualDetectorName")),
    toff_(pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
    diagLevel_(pset.get<int>("diagLevel",0)),
    nProcess_(0),
    numEvents_(0),
    Ntup_(0),
    hEvents_(0),
    sumEcrys_(0)

  {}

  void CaloNeutron::beginJob(){

       art::ServiceHandle<art::TFileService> tfs;
       
       hEvents_ = tfs->make<TH1D>("numEvents", "numEvents", 1, -0.5, 0.5);
       
       
       Ntup_  = tfs->make<TTree>("Calo", "Calo");

       Ntup_->Branch("evt",          &_evt ,        "evt/I");
       Ntup_->Branch("run",          &_run ,        "run/I");

       Ntup_->Branch("nVd",      &_nVd ,     "nVd/I");
       Ntup_->Branch("vdId",     &_vdId ,    "vdId[nVd]/I");
       Ntup_->Branch("vdPdgId",  &_vdPdgId , "vdPdgId[nVd]/I");
       Ntup_->Branch("vdMom",    &_vdMom ,   "vdMom[nVd]/F");
       Ntup_->Branch("vdMomX",   &_vdMomX ,  "vdMomX[nVd]/F");
       Ntup_->Branch("vdMomY",   &_vdMomY ,  "vdMomY[nVd]/F");
       Ntup_->Branch("vdMomZ",   &_vdMomZ ,  "vdMomZ[nVd]/F");
       Ntup_->Branch("vdPosX",   &_vdPosX ,  "vdPosX[nVd]/F");
       Ntup_->Branch("vdPosY",   &_vdPosY ,  "vdPosY[nVd]/F");
       Ntup_->Branch("vdPosZ",   &_vdPosZ ,  "vdPosZ[nVd]/F");
       Ntup_->Branch("vdTime",   &_vdTime ,  "vdTime[nVd]/F");

       Ntup_->Branch("nCry",     &_nCry ,        "nCry/I");
       Ntup_->Branch("cryId",    &_cryId ,       "cryId[nCry]/I");
       Ntup_->Branch("cryEdep",  &_cryEdep ,     "cryEdep[nCry]/F");
       Ntup_->Branch("cryTime",  &_cryTime ,     "cryTime[nCry]/F");
       Ntup_->Branch("cryPosX",  &_cryPosX ,     "cryPosX[nCry]/F");
       Ntup_->Branch("cryPosY",  &_cryPosY ,     "cryPosY[nCry]/F");
       Ntup_->Branch("cryPosZ",  &_cryPosZ ,     "cryPosZ[nCry]/F");

       Ntup_->Branch("nRO",      &_nRO ,     "nRO/I");
       Ntup_->Branch("ROId",     &_ROId ,    "vdId[nRO]/I");
       Ntup_->Branch("ROEdep",   &_ROEdep ,  "ROEdep[nRO]/F");
       Ntup_->Branch("ROPosX",   &_ROPosX ,  "ROPosX[nRO]/F");
       Ntup_->Branch("ROPosY",   &_ROPosY ,  "ROPosY[nRO]/F");
       Ntup_->Branch("ROPosZ",   &_ROPosZ ,  "ROPosZ[nRO]/F");
       
       Ntup_->Branch("nROCard",      &_nROCard ,     "nROCard/I");
       Ntup_->Branch("ROCardId",     &_ROCardId ,    "vdId[nROCard]/I");
       Ntup_->Branch("ROCardEdep",   &_ROCardEdep ,  "ROCardEdep[nROCard]/F");
       Ntup_->Branch("ROCardPosX",   &_ROCardPosX ,  "ROCardPosX[nROCard]/F");
       Ntup_->Branch("ROCardPosY",   &_ROCardPosY ,  "ROCardPosY[nROCard]/F");
       Ntup_->Branch("ROCardPosZ",   &_ROCardPosZ ,  "ROCardPosZ[nROCard]/F");

       Ntup_->Branch("nCrate",      &_nCrate ,     "nCrate/I");
       Ntup_->Branch("CrateId",     &_CrateId ,    "vdId[nCrate]/I");
       Ntup_->Branch("CrateEdep",   &_CrateEdep ,  "CrateEdep[nCrate]/F");
       Ntup_->Branch("CratePosX",   &_CratePosX ,  "CratePosX[nCrate]/F");
       Ntup_->Branch("CratePosY",   &_CratePosY ,  "CratePosY[nCrate]/F");
       Ntup_->Branch("CratePosZ",   &_CratePosZ ,  "CratePosZ[nCrate]/F");

 
  }



  void CaloNeutron::endJob()
  {
      std::cout<<"We processed "<<numEvents_<<" events"<<std::endl;
      std::cout<<"Porcessed crystal energy sum="<<sumEcrys_<<std::endl;
      hEvents_->Fill(0., numEvents_);
  }

  //================================================================
  void CaloNeutron::endSubRun(const art::SubRun& sr)
  {

    std::vector<art::Handle<GenEventCount> > hh = sr.getMany<GenEventCount>();

    if(hh.size() > 1) 
    {
       std::ostringstream os;
       os<<"GenEventCountReader: multiple GenEventCount objects found in "
         <<sr.id()<<":\n";
       for(const auto& h : hh) {
         os<<"    moduleLabel = "<<h.provenance()->moduleLabel()
           <<", instance = "<<h.provenance()->productInstanceName()
           <<", process = "<<h.provenance()->processName()
           <<"\n";
       }
       os<<"\n";
       throw cet::exception("BADCONFIG")<<os.str();
    }
    else if(hh.empty())
    {
       throw cet::exception("BADCONFIG")
             <<"GenEventCountReader: no GenEventCount record in "<<sr.id()<<"\n";
    }

    mf::LogInfo("INFO")<<"GenEventCount: "
                       <<hh.front()->count()<<" events in "<<sr.id()
                       <<"\n";
    numEvents_ += hh.front()->count();
  }



  void CaloNeutron::analyze(const art::Event& event) {

      ++nProcess_;
      if (nProcess_%10==0 && diagLevel_ > 0) std::cout<<"Processing event from CaloNeutron =  "<<nProcess_ <<std::endl;

      ConditionsHandle<AcceleratorParams> accPar("ignored");
      double _mbtime = accPar->deBuncherPeriod;
      toff_.updateMap(event);

      //Handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());

      //Virtual detector hits
      art::Handle<StepPointMCCollection> vdhits;
      event.getByLabel(g4ModuleLabel_,virtualDetectorLabel_,vdhits);



      // These selectors will select data products with the given instance name, and ignore all other fields of the product ID.
      art::ProductInstanceNameSelector getCrystalSteps(calorimeterStepPoints_);
      art::ProductInstanceNameSelector getROSteps(calorimeterROStepPoints_);
      art::ProductInstanceNameSelector getROCardSteps(calorimeterROCardStepPoints_);
      art::ProductInstanceNameSelector getCrateSteps(calorimeterCrateStepPoints_);


      // Get the StepPointMCs from the event.
      HandleVector crystalStepsHandles, ROStepsHandles,ROCardStepsHandles,CrateStepsHandles;
      crystalStepsHandles = event.getMany<StepPointMCCollection>(getCrystalSteps);
      ROStepsHandles = event.getMany<StepPointMCCollection>(getROSteps);
      ROCardStepsHandles = event.getMany<StepPointMCCollection>(getROCardSteps);
      CrateStepsHandles = event.getMany<StepPointMCCollection>(getCrateSteps);




      _evt = event.id().event();
      _run = event.run();


       _nCry=0;
       for ( HandleVector::const_iterator i=crystalStepsHandles.begin(), e=crystalStepsHandles.end(); i != e; ++i )
       {     
            const art::Handle<StepPointMCCollection>& handle(*i);
            const StepPointMCCollection& steps(*handle);

            for (const auto& step : steps)
            {
                CLHEP::Hep3Vector Pos = cal.geomUtil().mu2eToTracker(step.position());
                
                _cryId[_nCry]   = step.volumeId();
                _cryEdep[_nCry] = step.totalEDep();
                _cryTime[_nCry] = step.time();
                _cryPosX[_nCry] = Pos.x();
                _cryPosY[_nCry] = Pos.y();
                _cryPosZ[_nCry] = Pos.z();
                ++_nCry;
                sumEcrys_ +=step.totalEDep();
            }
       }

       _nRO=0;
       for ( HandleVector::const_iterator i=ROStepsHandles.begin(), e=ROStepsHandles.end(); i != e; ++i )
       {     
            const art::Handle<StepPointMCCollection>& handle(*i);
            const StepPointMCCollection& steps(*handle);

            for (const auto& step : steps)
            {
                CLHEP::Hep3Vector Pos = cal.geomUtil().mu2eToTracker(step.position());
                _ROId[_nRO]=step.volumeId();
                _ROEdep[_nRO]=step.totalEDep();
                _ROPosX[_nRO]=Pos.x();
                _ROPosY[_nRO]=Pos.y();
                _ROPosZ[_nRO]=Pos.z();
                ++_nRO;
            }
       }

       _nROCard=0;
       for ( HandleVector::const_iterator i=ROCardStepsHandles.begin(), e=ROCardStepsHandles.end(); i != e; ++i )
       {     
            const art::Handle<StepPointMCCollection>& handle(*i);
            const StepPointMCCollection& steps(*handle);

            for (const auto& step : steps)
            {
                CLHEP::Hep3Vector Pos = cal.geomUtil().mu2eToTracker(step.position());
                _ROCardId[_nROCard]=step.volumeId();
                _ROCardEdep[_nROCard]=step.totalEDep();
                _ROCardPosX[_nROCard]=Pos.x();
                _ROCardPosY[_nROCard]=Pos.y();
                _ROCardPosZ[_nROCard]=Pos.z();
                ++_nROCard;
            }
       }

       _nCrate=0;
       for ( HandleVector::const_iterator i=CrateStepsHandles.begin(), e=CrateStepsHandles.end(); i != e; ++i )
       {     
            const art::Handle<StepPointMCCollection>& handle(*i);
            const StepPointMCCollection& steps(*handle);

            for (const auto& step : steps)
            {
                CLHEP::Hep3Vector Pos = cal.geomUtil().mu2eToTracker(step.position());
                _CrateId[_nCrate]=step.volumeId();
                _CrateEdep[_nCrate]=step.totalEDep();
                _CratePosX[_nCrate]=Pos.x();
                _CratePosY[_nCrate]=Pos.y();
                _CratePosZ[_nCrate]=Pos.z();
                ++_nCrate;
            }
       }



      //--------------------------  Do virtual detectors --------------------------------
      //73/74/77/78 front back inner outer edges disk 0
      //75/76/79/80 front back inner outer edges disk 1


      _nVd = 0;
      if (vdhits.isValid())
      {
          for (auto iter=vdhits->begin(), ie=vdhits->end(); iter!=ie; ++iter)
            {
              const StepPointMC& hit = *iter;

              //if (hit.volumeId()<VirtualDetectorId::EMC_Disk_0_SurfIn || hit.volumeId()>VirtualDetectorId::EMC_Disk_1_EdgeOut) continue;

              double hitTimeUnfolded = toff_.timeWithOffsetsApplied(hit);
   	      double hitTime         = fmod(hitTimeUnfolded,_mbtime);

              CLHEP::Hep3Vector VDPos = cal.geomUtil().mu2eToTracker(hit.position());
              //CLHEP::Hep3Vector VDPos = cal.toTrackerFrame(hit.position());

              _vdId[_nVd]    = hit.volumeId();
              _vdPdgId[_nVd] = hit.simParticle()->pdgId();
              _vdTime[_nVd]  = hitTime;
              _vdPosX[_nVd]  = VDPos.x(); //tracker frame
              _vdPosY[_nVd]  = VDPos.y();
              _vdPosZ[_nVd]  = VDPos.z();
              _vdMom [_nVd]  = hit.momentum().mag();
              _vdMomX[_nVd]  = hit.momentum().x();
              _vdMomY[_nVd]  = hit.momentum().y();
              _vdMomZ[_nVd]  = hit.momentum().z();
              _vdenIdx[_nVd] = hit.simParticle()->generatorIndex();
              ++_nVd;
            }             
        }


        Ntup_->Fill();

  }



}  

DEFINE_ART_MODULE(mu2e::CaloNeutron);


