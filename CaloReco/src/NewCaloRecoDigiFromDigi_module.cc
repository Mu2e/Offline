#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

#include "CaloReco/inc/WaveformProcessor.hh"
#include "CaloReco/inc/LogNormalProcessor.hh"
#include "CaloReco/inc/FixedFastProcessor.hh"
#include "CaloReco/inc/RawProcessor.hh"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "RecoDataProducts/inc/NewCaloDigi.hh"
#include "RecoDataProducts/inc/NewCaloDigiCollection.hh"
#include "RecoDataProducts/inc/NewCaloRecoDigi.hh"
#include "RecoDataProducts/inc/NewCaloRecoDigiCollection.hh"

#include <iostream>
#include <string>


namespace mu2e {

  class NewCaloRecoDigiFromDigi : public art::EDProducer {

  public:

    enum processorStrategy {NoChoice, RawExtract, LogNormalFit, FixedFast};
    struct Config{
	      using Name=fhicl::Name;
	      using Comment=fhicl::Comment;
	      fhicl::Atom<art::InputTag> digitag{Name("NewCaloDigiCollection"),Comment("tag for digi collection")};
	      fhicl::Atom<double> digiSampling{Name("digiSampling"),Comment("todo")};
	      fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("set to 1 for info"),0};
	      fhicl::Atom<std::string> processorStrategy{Name("processorStrategy"), Comment("todo")};
	     // fhicl::Table<SomeOtherClass::Config> tfit{Name("CallName"), Comment("fit")};
	      fhicl::Atom<double> maxChi2Cut{Name("maxChi2Cut"), Comment("todo")};
              fhicl::Atom<int> nplot{Name("nplot"), Comment("todo"), 0};
	     
       };
      
      
     typedef art::EDProducer::Table<Config> Parameters;

      explicit NewCaloRecoDigiFromDigi(const Parameters& conf);
      virtual ~NewCaloRecoDigiFromDigi() {};
      virtual void produce(art::Event& e);
      virtual void beginRun(art::Run& aRun) override;
   
  private:
    Config _conf;
    

    art::InputTag digitag_;
    double const digiSampling_;
    int          diagLevel_;
    std::string const processorStrategy_;
    double       maxChi2Cut_;   
    int          nplot_;

    std::unique_ptr<WaveformProcessor> waveformProcessor_;

    void extractRecoDigi(art::ValidHandle<NewCaloDigiCollection> const& caloDigis,
                         NewCaloRecoDigiCollection& recoCaloHits);

  };

NewCaloRecoDigiFromDigi::NewCaloRecoDigiFromDigi(const Parameters& conf) :
	   art::EDProducer(conf),
	   digitag_(conf().digitag()),
	   digiSampling_(conf().digiSampling()),
	   diagLevel_(conf().diagLevel()),
    	   processorStrategy_(conf().processorStrategy()),
           maxChi2Cut_(conf().maxChi2Cut()),
           nplot_(conf().nplot())
    {
		 
	      produces<NewCaloRecoDigiCollection>();
             
	      std::map<std::string, processorStrategy> spmap;
	      spmap["RawExtract"]   = RawExtract;
	      spmap["LogNormalFit"] = LogNormalFit;
	      spmap["FixedFast"]    = FixedFast;

	      switch (spmap[processorStrategy_])
		{
		case RawExtract:
		  {
		    auto const& param = pset.get<fhicl::ParameterSet>("RawProcessor", {});
		    waveformProcessor_ = std::make_unique<RawProcessor>(param);
		    break;
		  }

		case LogNormalFit:
		  {
		    auto const& param = pset.get<fhicl::ParameterSet>("LogNormalProcessor", {});
		    waveformProcessor_ = std::make_unique<LogNormalProcessor>(param);
		    break;
		  }

		case FixedFast:
		  {
		    auto const& param = pset.get<fhicl::ParameterSet>("FixedFastProcessor", {});
		    waveformProcessor_ = std::make_unique<FixedFastProcessor>(param);
		    break;
		  }

		default:
		  {
		    throw cet::exception("CATEGORY")<< "Unrecognized processor in CaloHitsFromDigis module";
		  }
		}
	
  }
  //-------------------------------------------------------
  void NewCaloRecoDigiFromDigi::produce(art::Event& event)
  {

    if (diagLevel_ > 0) std::cout<<"[NewCaloRecoDigiFromDigi::produce] begin"<<std::endl;

    auto const& digiH = event.getValidHandle<NewCaloDigiCollection>(digitag_);
	
    //Get the calorimeter Digis
    //auto const& caloDigisH = event.getValidHandle(caloDigisToken_);

    auto recoCaloDigiColl = std::make_unique<NewCaloRecoDigiCollection>();
    extractRecoDigi(digiH, *recoCaloDigiColl);

    if ( diagLevel_ > 3 )
      {
        printf("[NewCaloRecoDigiFromDigi::produce] produced RecoCrystalHits ");
        printf(", recoCaloDigiColl size  = %i \n", int(recoCaloDigiColl->size()));
      }

    event.put(std::move(recoCaloDigiColl));

    if (diagLevel_ > 0) std::cout<<"[NewCaloRecoDigiFromDigi::produce] end"<<std::endl;

    return;
  }

  //-----------------------------------------------------------------------------
  void NewCaloRecoDigiFromDigi::beginRun(art::Run& aRun)
  {
    waveformProcessor_->initialize();
  }

  //--------------------------------------------------------------------------------------
  void NewCaloRecoDigiFromDigi::extractRecoDigi(art::ValidHandle<NewCaloDigiCollection> const& caloDigisHandle, NewCaloRecoDigiCollection &recoCaloHits)
  {

    std::vector<double> x,y;
    ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

    auto const& caloDigis = *caloDigisHandle;
    NewCaloDigi const* base = &caloDigis.front(); // What if caloDigis is empty?
    
    for (const auto& caloDigi : caloDigis)
      {
       
        int    roId     = caloDigi.roId();
        double t0       = caloDigi.t0();
        double adc2MeV  = calorimeterCalibrations->ADC2MeV(roId);
        const std::vector<int>& waveform = caloDigi.waveform();

        size_t index = &caloDigi - base;
        art::Ptr<NewCaloDigi> caloDigiPtr(caloDigisHandle, index);

        x.clear();
        y.clear();
        for (unsigned int i=0;i<waveform.size();++i)
          {
            x.push_back(t0 + (i+0.5)*digiSampling_); // add 0.5 to be in middle of bin
            y.push_back(waveform.at(i));
          }

        if (diagLevel_ > 3)
          {
            std::cout<<"[NewCaloRecoDigiFromDigi::extractRecoDigi] extract amplitude from this set of hits for RoId="<<roId<<" a time "<<t0<<std::endl;
            for (auto const& val : waveform) {std::cout<< val<<" ";} std::cout<<std::endl;
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
                std::cout<<"[NewCaloRecoDigiFromDigi::extractAmplitude] extract "<<roId<<"   i="<<i<<"  eDep="<<eDep<<" time="<<time<<"  chi2="<<chi2<<std::endl;
              }

            if (chi2/ndf > maxChi2Cut_) continue;

            recoCaloHits.emplace_back(NewCaloRecoDigi(roId, caloDigiPtr, eDep,eDepErr,time,timeErr,chi2,ndf,isPileUp));
          }

      }

  }

}

DEFINE_ART_MODULE(mu2e::NewCaloRecoDigiFromDigi);
