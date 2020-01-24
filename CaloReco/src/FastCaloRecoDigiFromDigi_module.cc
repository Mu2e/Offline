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
		art::EDProducer{pset},
		caloDigisToken_{consumes<CaloDigiCollection>(pset.get<std::string>("caloDigiModuleLabel"))},
		digiSampling_        (pset.get<double>("digiSampling")),
		windowPeak_         (pset.get<int>    ("windowPeak")),
		minPeakAmplitude_   (pset.get<double> ("minPeakAmplitude")),
		shiftTime_          (pset.get<double> ("shiftTime")),
		scaleFactor_        (pset.get<double> ("scaleFactor")),
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
		int 	 windowPeak_ ;
		double 	 minPeakAmplitude_ ;
		double 	 shiftTime_ ;
		double       scaleFactor_;
		int          diagLevel_;
   		art::InputTag ewMarkerTag_;// name of the module that makes eventwindowmarkers
   
		void extractRecoDigi(art::ValidHandle<CaloDigiCollection> const& caloDigis, CaloRecoDigiCollection& recoCaloHits, double const& ewOffset);

	};

	void FastCaloRecoDigiFromDigi::produce(art::Event& event)
  	{

		if (diagLevel_ > 0) std::cout<<"[FastRecoDigiFromDigi::produce] begin"<<std::endl;

		auto const& caloDigisH = event.getValidHandle(caloDigisToken_);
		auto recoCaloDigiColl = std::make_unique<CaloRecoDigiCollection>();
		//Get EW time offset (as per StrawHitReco):
		double ewmOffset = 0;
		art::Handle<EventWindowMarker> ewMarkerHandle;
		if (event.getByLabel(ewMarkerTag_, ewMarkerHandle)){
			const EventWindowMarker& ewMarker(*ewMarkerHandle);
			ewmOffset = ewMarker.timeOffset();
		}

		extractRecoDigi(caloDigisH, *recoCaloDigiColl, ewmOffset);

		if ( diagLevel_ > 3 )
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
		ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");//TODO

		auto const& caloDigis = *caloDigisHandle;
		//if(caloDigis.size() == 0){ continue; } TODO
		CaloDigi const* base = &caloDigis.front(); 

		for (const auto& caloDigi : caloDigis)
		{
    
			int    roId     = caloDigi.roId();
			double t0       = caloDigi.t0();
			//removed crId -> no need as we extract crystal hit in this path
			double adc2MeV  = calorimeterCalibrations->ADC2MeV(roId);
			const std::vector<int>& waveform = caloDigi.waveform();

			size_t index = &caloDigi - base;
			art::Ptr<CaloDigi> caloDigiPtr(caloDigisHandle, index);

			x.clear();
			y.clear();
			for (unsigned int i=0;i<waveform.size();++i)
			  {
			   	x.push_back(t0 + (i+0.5)*digiSampling_+eventWindowOffset); 
				y.push_back(waveform.at(i)); //amplitude
			  }

			if (diagLevel_ > 3)
			  {
				    std::cout<<"[FastRecoDigiFromDigi::extractRecoDigi] extract amplitude from this set of hits for RoId="<<roId<<" a time "<<t0<<std::endl;
				    for (auto const& val : waveform) {std::cout<< val<<" ";} std::cout<<std::endl;
			  }

			unsigned int nPeaks_ = caloDigi.peakpos().size(); //should always be 1 (assume no pile up)
			double chi2   = 0;
			int ndf    = 1;
			bool isPileUp;
			if(nPeaks_ > 1) { isPileUp = true ; } //set pile up flag
			else { isPileUp = false; }
			for (unsigned int i=0;i<nPeaks_;++i)
			{ 
		 
				if(y[caloDigi.peakpos()[i]] <  minPeakAmplitude_) {continue;} // New Feature!
				double eDep= scaleFactor_*y[caloDigi.peakpos()[i]]*adc2MeV;
				double eDepErr =  0*adc2MeV;
				double time =  x[caloDigi.peakpos()[i]] - shiftTime_;
				double timeErr = 0;

				 if (diagLevel_ > 1)
				  {
				    std::cout<<"[FastRecoDigiFromDigi::extractAmplitude] extract "<<roId<<"   i="<<i<<"  eDep="<<eDep<<" time="<<time<<"  chi2="<<chi2<<std::endl;
				  }
				recoCaloHits.emplace_back(CaloRecoDigi(roId, caloDigiPtr, eDep,eDepErr,time,timeErr,chi2,ndf,isPileUp));
			}
		}
	}
}

DEFINE_ART_MODULE(mu2e::FastCaloRecoDigiFromDigi);
