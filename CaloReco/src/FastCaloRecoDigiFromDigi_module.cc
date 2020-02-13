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

		if ( diagLevel_ > 0 )
		{
			printf("[FastRecoDigiFromDigi::produce] produced RecoCrystalHits ");
			printf(", recoCaloDigiColl size  = %i \n", int(recoCaloDigiColl->size()));
		}

		event.put(std::move(recoCaloDigiColl));

		if (diagLevel_ > 0) std::cout<<"[FastRecoDigiFromDigi::produce] end"<<std::endl;

		return;
	}

	void FastCaloRecoDigiFromDigi::extractRecoDigi(art::ValidHandle<CaloDigiCollection> const& caloDigisHandle, CaloRecoDigiCollection &recoCaloHits, double const& eventWindowOffset)
  	{

		std::vector<double> x,y;
		ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

		auto const& caloDigis = *caloDigisHandle;
		
		CaloDigi const* base = &caloDigis.front(); 

		for (const auto& caloDigi : caloDigis)
		{
    
			int    roId     = caloDigi.roId();
			double t0       = caloDigi.t0();
			double adc2MeV  = calorimeterCalibrations->Peak2MeV(roId);
			const std::vector<int>& waveform = caloDigi.waveform();

			size_t index = &caloDigi - base;
			art::Ptr<CaloDigi> caloDigiPtr(caloDigisHandle, index);

			x.clear();
			y.clear();
			for (unsigned int i=0;i<waveform.size();++i)
			  {
			   	x.push_back(t0 + (i+0.5)*digiSampling_+eventWindowOffset); 
				y.push_back(waveform.at(i)); 
			  }

			if (diagLevel_ > 3)
			  {
				    std::cout<<"[FastRecoDigiFromDigi::extractRecoDigi] extract amplitude from this set of hits for RoId="<<roId<<" a time "<<t0<<std::endl;
				    for (auto const& val : waveform) {std::cout<< val<<" ";} 
			  }

			
			double eDep= y[caloDigi.peakpos()]*adc2MeV;
			if(eDep <  minDigiE_) {continue;}
			double eDepErr =  0*adc2MeV;
			double time =  x[caloDigi.peakpos()] - shiftTime_;
			double timeErr = 0;

			 if (diagLevel_ > 1)
			  {
			    std::cout<<"[FastRecoDigiFromDigi::extractAmplitude] extract "<<roId<<"  eDep="<<eDep<<" time="<<time<<std::endl;
			  }
			recoCaloHits.emplace_back(CaloRecoDigi(roId, caloDigiPtr, eDep,eDepErr,time,timeErr,0,1,false));
			
		}
	}
}

DEFINE_ART_MODULE(mu2e::FastCaloRecoDigiFromDigi);
