//Author: S Middleton
//Date: Nov 2019
//Purpose: To make CaloRecoHits quickly from the Digis using the  data product structure (for Online compatib.)

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/CaloRecoDigi.hh"
#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"

#include "DataProducts/inc/EventWindowMarker.hh"

#include <iostream>
#include <string>


#include <iostream>
#include <string>
   
namespace mu2e {

  class FastCaloRecoDigiFromDigi : public art::EDProducer {

  	public:
		explicit FastCaloRecoDigiFromDigi(fhicl::ParameterSet const& pset) :
		art::EDProducer{pset}, caloDigisToken_{consumes<CaloDigiCollection>(pset.get<std::string>("caloDigiModuleLabel"))},
		digiSampling_       (pset.get<double>("digiSampling")),
		minDigiE_   	    (pset.get<double> ("minDigiE")),
		shiftTime_          (pset.get<double> ("shiftTime")),
		diagLevel_          (pset.get<int>    ("diagLevel",0)),
		ewMarkerTag_(pset.get<art::InputTag>("EventWindowMarkerLabel"))
	     	{
			consumes<EventWindowMarker>(ewMarkerTag_);
	      		produces<CaloRecoDigiCollection>();
		}

    		void produce(art::Event& e) override;

	private:

		art::ProductToken<CaloDigiCollection> const caloDigisToken_;
		double const digiSampling_;
		double 	 minDigiE_ ;
		double 	 shiftTime_ ;
		int      diagLevel_;
   		art::InputTag ewMarkerTag_;// name of the module that makes eventwindowmarkers
   
		void extractRecoDigi(art::ValidHandle<CaloDigiCollection> const& caloDigis, CaloRecoDigiCollection& recoCaloHits, double const& ewOffset);

	};

	void FastCaloRecoDigiFromDigi::produce(art::Event& event)
  	{

		if (diagLevel_ > 0) std::cout<<"[FastRecoDigiFromDigi::produce] begin"<<std::endl;

		auto const& caloDigisH = event.getValidHandle(caloDigisToken_);
		auto recoCaloDigiColl = std::make_unique<CaloRecoDigiCollection>();
		
		double ewmOffset = 0;
		art::Handle<EventWindowMarker> ewMarkerHandle;
		if (event.getByLabel(ewMarkerTag_, ewMarkerHandle)){
			const EventWindowMarker& ewMarker(*ewMarkerHandle);
			ewmOffset = ewMarker.timeOffset();
		}

		extractRecoDigi(caloDigisH, *recoCaloDigiColl, ewmOffset);

		if ( diagLevel_ > 0 )std::cout<<"[FastCaloRecoDigiFromDigi::produce] extracted "<<recoCaloDigiColl->size()<<" RecoDigis"<<std::endl;
		

		event.put(std::move(recoCaloDigiColl));

		if (diagLevel_ > 0) std::cout<<"[FastRecoDigiFromDigi::produce] end"<<std::endl;

		return;
	}

	void FastCaloRecoDigiFromDigi::extractRecoDigi(art::ValidHandle<CaloDigiCollection> const& caloDigisHandle, CaloRecoDigiCollection &recoCaloHits, double const& eventWindowOffset)
  	{
		ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
		auto const& caloDigis = *caloDigisHandle;
		CaloDigi const* base = &caloDigis.front(); 

		for (const auto& caloDigi : caloDigis)
		{
 			int    roId     = caloDigi.roId();
			double t0       = caloDigi.t0();
				// TODO:+ calorimeterCalibrations->timeOffset(roId);
			double adc2MeV  = calorimeterCalibrations->Peak2MeV(roId);
			size_t index = &caloDigi - base;
			art::Ptr<CaloDigi> caloDigiPtr(caloDigisHandle, index);

			double time = t0 + (caloDigi.peakpos()+0.5)*digiSampling_+ eventWindowOffset - shiftTime_; 
			double eDep = (caloDigi.waveform().at(caloDigi.peakpos()))*adc2MeV; 
			if(eDep <  minDigiE_) {continue;}
			double eDepErr =  0*adc2MeV;
			double timeErr = 0;

			if (diagLevel_ > 1)
			  {
			    std::cout<<"[FastRecoDigiFromDigi::extractAmplitude] extracted Digi with roId =  "<<roId<<"  eDep = "<<eDep<<" time = " <<time<<std::endl;
			  }

			recoCaloHits.emplace_back(CaloRecoDigi(roId, caloDigiPtr, eDep,eDepErr,time,timeErr,0,1,false));
			
		}
	}
}

DEFINE_ART_MODULE(mu2e::FastCaloRecoDigiFromDigi);
