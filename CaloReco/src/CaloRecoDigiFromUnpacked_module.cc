//
// Recall, the caloFigiColl structure is the following:
// nTotWords - nWords_roID - roiID - nWord_roID_ihit - time_roID_ihit - Data_roID_ihit - ...
//

// C++ includes.
#include <iostream>
#include <string>

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
#include "CaloReco/inc/WaveformProcessor.hh"
#include "CaloReco/inc/LogNormalProcessor.hh"
#include "CaloReco/inc/FixedFastProcessor.hh"
#include "CaloReco/inc/RawProcessor.hh"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloDigiPacked.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"



namespace mu2e {



  class CaloRecoDigiFromUnpacked : public art::EDProducer {

    public:

        enum processorStrategy {NoChoice, RawExtract, LogNormalFit, FixedFast};

        explicit CaloRecoDigiFromUnpacked(fhicl::ParameterSet const& pset) :
          caloDigiUnpackModuleLabel_ (pset.get<std::string>("caloDigiUnpackModuleLabel")),
          processorStrategy_         (pset.get<std::string>("processorStrategy")),
          digiSampling_              (pset.get<double>     ("digiSampling")),
          maxChi2Cut_                (pset.get<double>     ("maxChi2Cut")),
          diagLevel_                 (pset.get<int>        ("diagLevel",0)),
	  nplot_(0)
        {

            produces<CaloRecoDigiCollection>();

            std::map<std::string, processorStrategy> spmap;
            spmap["RawExtract"]   = RawExtract;
            spmap["LogNormalFit"] = LogNormalFit;
            spmap["FixedFast"]    = FixedFast;            

            switch (spmap[processorStrategy_])
            {
                case RawExtract:
                {
                    auto const& param = pset.get<fhicl::ParameterSet>("RawProcessor",fhicl::ParameterSet());
                    waveformProcessor_ = std::unique_ptr<WaveformProcessor>(new RawProcessor(param));
                    break;
                }               

                case LogNormalFit: 
                {
                    auto const& param = pset.get<fhicl::ParameterSet>("LogNormalProcessor",fhicl::ParameterSet());
                    waveformProcessor_ = std::unique_ptr<WaveformProcessor>(new LogNormalProcessor(param));           
                    break;
                } 

                case FixedFast: 
                {
                    auto const& param = pset.get<fhicl::ParameterSet>("FixedFastProcessor",fhicl::ParameterSet());
                    waveformProcessor_ = std::unique_ptr<WaveformProcessor>(new FixedFastProcessor(param));           
                    break;
                }

                default:
                {
                    throw cet::exception("CATEGORY")<< "Unrecognized processor in CaloHitsFromDigis module";    
                }
            }
        }

        virtual ~CaloRecoDigiFromUnpacked() {}
        
        virtual void beginRun(art::Run& aRun);
	virtual void produce(art::Event& e);




    private:

       std::string  caloDigiUnpackModuleLabel_; 
       std::string  processorStrategy_;
       double       digiSampling_;
       double       maxChi2Cut_;
       int          diagLevel_;
       int          nplot_;

       std::unique_ptr<WaveformProcessor> waveformProcessor_;


       void extractRecoDigi(const art::Handle<CaloDigiCollection>& caloDigisHandle, 
                            CaloRecoDigiCollection &recoCaloHits);



  };


  //-------------------------------------------------------
  void CaloRecoDigiFromUnpacked::produce(art::Event& event) 
  {

       if (diagLevel_ > 3) printf("[CaloRecoDigiFromUnpacked::produce] event %i\n", event.id().event());

       //Get handles to calorimeter RO (aka APD) collection
       art::Handle<CaloDigiCollection> caloDigisHandle;
       event.getByLabel(caloDigiUnpackModuleLabel_, caloDigisHandle);


       std::unique_ptr<CaloRecoDigiCollection> recoCaloDigiColl(new CaloRecoDigiCollection);
       extractRecoDigi( caloDigisHandle, *recoCaloDigiColl);


       if ( diagLevel_ > 3 )
       {
           printf("[CaloRecoDigiFromUnpacked::produce] produced RecoCrystalHits ");
           printf(", recoCaloDigiColl size  = %i \n", int(recoCaloDigiColl->size()));
       }

       event.put(std::move(recoCaloDigiColl));

       return;
  }
  
  //-----------------------------------------------------------------------------
  void CaloRecoDigiFromUnpacked::beginRun(art::Run& aRun)
  {
      waveformProcessor_->initialize(); 
  }

    
  //--------------------------------------------------------------------------------------
  void CaloRecoDigiFromUnpacked::extractRecoDigi(const art::Handle<CaloDigiCollection>& caloDigisHandle, 
                                            CaloRecoDigiCollection &recoCaloHits)
  {
      
      std::vector<double> x,y;
      
      ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

      const CaloDigiCollection& caloDigis(*caloDigisHandle);
      CaloDigi const* base = &caloDigis.front();

    
      for (const auto& caloDigi : caloDigis)
      {
            int    roId     = caloDigi.roId();
            double t0       = caloDigi.t0();
            double adc2MeV  = calorimeterCalibrations->ADC2MeV(roId);
            const std::vector<int> waveform = caloDigi.waveform();
           
            size_t index = &caloDigi - base; 
            art::Ptr<CaloDigi> caloDigiPtr(caloDigisHandle, index);
              
            x.clear();
            y.clear();
            for (unsigned int i=0;i<waveform.size();++i)
            {
                x.push_back(t0 + (i+0.5)*digiSampling_); // add 0.5 to be in middle of bin
                y.push_back(waveform.at(i));            
            }

            if (diagLevel_ > 3)
            {
                std::cout<<"CaloRecoDigiFromUnpacked::extractAmplitude from this set of hits for RoId="<<roId<<" a time "<<t0<<std::endl;
                for (auto const& val : waveform) std::cout<< val<<" "; std::cout<<std::endl;
            }
            
            waveformProcessor_->reset();
            waveformProcessor_->extract(x,y);
            
            for (int i=0;i<waveformProcessor_->nPeaks();++i)
            {
                 double eDep      = waveformProcessor_->amplitude(i)*adc2MeV;
                 double eDepErr   = waveformProcessor_->amplitudeErr(i)*adc2MeV;
                 double time      = waveformProcessor_->time(i);
                 double timeErr   = waveformProcessor_->timeErr(i);
                 bool   isPileUp  = waveformProcessor_->isPileUp(i);
                 double chi2      = waveformProcessor_->chi2();
                 int    ndf       = waveformProcessor_->ndf();               

                 if (diagLevel_ > 1)
                 {
                     std::cout<<"CaloRecoDigiFromUnpacked::extractAmplitude "<<roId<<"   i="<<i<<"  eDep="<<eDep<<" time="<<time<<"  chi2="<<chi2<<std::endl;                   
                 }           

                 if (chi2/ndf > maxChi2Cut_) continue;
                 recoCaloHits.emplace_back(CaloRecoDigi(roId, caloDigiPtr, eDep,eDepErr,time,timeErr,chi2,ndf,isPileUp));
            }

      }
  
  }          

}

using mu2e::CaloRecoDigiFromUnpacked;
DEFINE_ART_MODULE(CaloRecoDigiFromUnpacked);


 /*
 if (waveformProcessor_->chi2()/waveformProcessor_->ndf() > 4 )
 {                 
     std::string pname = "plots/plot_"+std::to_string(nplot_)+".pdf";
     std::cout<<"Saved in file "<<pname<<std::endl;
     waveformProcessor_->plot(pname);          
     ++nplot_;
 }
 */
