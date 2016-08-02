
// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/CaloShowerCollection.hh"
#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"

#include "TH2F.h"
#include "TFile.h"

#include <iostream>
#include <string>


namespace mu2e {



  class OldCaloHitFromShower : public art::EDProducer {

     public:

	 enum processorStrategy {LogNormalFit = 0, RawExtract = 1};

	 explicit OldCaloHitFromShower(fhicl::ParameterSet const& pset) :

	   caloShowerModuleLabel_ (pset.get<std::string>("caloshowerModuleLabel")),
	   diagLevel_             (pset.get<int>        ("diagLevel",0))
	 {
	     produces<CaloRecoDigiCollection>();
	 }

	 virtual ~OldCaloHitFromShower() {}
	 void produce(art::Event& e);



     private:

	std::string  caloShowerModuleLabel_; 
	int          diagLevel_;

	void makeCaloRecoDigis(const CaloShowerCollection& caloShowers,
			       CaloRecoDigiCollection &recoCaloHits);       


  };


  //-------------------------------------------------------
  void OldCaloHitFromShower::produce(art::Event& event) 
  {

      if (diagLevel_ > 2 ) printf("[OldCaloHitFromShower::produce] event %i\n", event.id().event());

      art::Handle<CaloShowerCollection> caloShowerHandle;
      event.getByLabel(caloShowerModuleLabel_, caloShowerHandle);
      const CaloShowerCollection& caloShowers(*caloShowerHandle);


      //Create a new RecoCaloCrystalHit collection and fill it
      std::unique_ptr<CaloRecoDigiCollection> recoCaloDigiColl(new CaloRecoDigiCollection);      

      makeCaloRecoDigis(caloShowers, *recoCaloDigiColl);
      
      event.put(std::move(recoCaloDigiColl));
    
      return;
  }

    
  //--------------------------------------------------------------------------------------
  void OldCaloHitFromShower::makeCaloRecoDigis(const CaloShowerCollection& caloShowers,
				               CaloRecoDigiCollection& recoCaloDigis)						   
  {
      
      Calorimeter const &cal = *(GeomHandle<Calorimeter>());
          
      int nROs = cal.caloGeomInfo().nROPerCrystal();
      art::Ptr<CaloDigi> dummy;     
      
      
      for (const auto& caloShower : caloShowers)
      {
          int ROidBase = cal.ROBaseByCrystal(caloShower.crystalId());
          for (int i=0; i<nROs; ++i) 
	     recoCaloDigis.emplace_back(CaloRecoDigi(ROidBase + i, dummy , caloShower.energy(),0,caloShower.time(),0,-1,-1,0));
      }

  }	  







}

using mu2e::OldCaloHitFromShower;
DEFINE_ART_MODULE(OldCaloHitFromShower);
