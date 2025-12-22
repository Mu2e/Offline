// Ed Callaghan
// Collate digis from multiple sources, and apply collision resolution
// August 2024

// stl
#include <string>
#include <vector>

// art
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"

// mu2e
#include "Offline/MCDataProducts/inc/CrvDigiMC.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawElectronics.hh"
#include "Offline/TrackerMC/inc/StrawDigiBundle.hh"
#include "Offline/TrackerMC/inc/StrawDigiBundleCollection.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/CaloMC/inc/CaloDigiWrapperCollection.hh"
#include "Offline/CRVConditions/inc/CRVCalib.hh"
#include "Offline/CRVConditions/inc/CRVADCRange.hh"
#include "Offline/CRVResponse/inc/CrvMCHelper.hh"

namespace mu2e{
  class MergeDigis: public art::EDProducer{
    public:
      struct Config{
        fhicl::Sequence<art::InputTag> tracker_digi_tags{
          fhicl::Name("StrawDigiCollections"),
          fhicl::Comment("art::InputTags of source StrawDigi and StrawDigiADCWaveforms")
        };
        fhicl::Atom<bool> tracker_mc{
          fhicl::Name("MergeStrawDigiMCs"),
          fhicl::Comment("True/false to merge tracker MC truth")
        };
        fhicl::Sequence<art::InputTag> calo_digi_tags{
          fhicl::Name("CaloDigiCollections"),
          fhicl::Comment("art::InputTags of source CaloDigis")
        };
        fhicl::Atom<unsigned int> calo_adc_bits{
          fhicl::Name("CalorimeterADCBitDepth"),
          fhicl::Comment("Bit depth of calorimeter adc readings (temporary)")
        };
        fhicl::Sequence<art::InputTag> crv_digi_tags{
          fhicl::Name("CrvDigiCollections"),
          fhicl::Comment("art::InputTags of source CrvDigi")
        };
        fhicl::Sequence<art::InputTag> crv_digimc_tags{   //needs to be empty or index-matched with CrvDigiCollections
          fhicl::Name("CrvDigiMCCollections"),
          fhicl::Comment("art::InputTags of source CrvDigiMC")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      MergeDigis(const Parameters&);

    protected:
      // tracker
      std::vector<art::InputTag> _tracker_digi_tags;
      bool _tracker_mc;
      ProditionsHandle<StrawElectronics> _tracker_conditions_handle;
      // crv
      std::vector<art::InputTag> _crv_digi_tags;
      std::vector<art::InputTag> _crv_digimc_tags;
      ProditionsHandle<CRVCalib>  _crv_calib_handle;

      // calorimeter
      std::vector<art::InputTag> _calo_digi_tags;
      CaloDigiWrapper::sample_t _calo_max_adc;

    private:
      void produce(art::Event&);
      void mergeCrvDigis(art::Event&);
  };

  // constructor
  MergeDigis::MergeDigis(const Parameters& config):
      art::EDProducer(config),
      _tracker_digi_tags(config().tracker_digi_tags()),
      _tracker_mc(config().tracker_mc()),
      _calo_digi_tags(config().calo_digi_tags()),
      _calo_max_adc(1 << config().calo_adc_bits()),
      _crv_digi_tags(config().crv_digi_tags()),
      _crv_digimc_tags(config().crv_digimc_tags()){
    // tracker
    for (const auto& tag: _tracker_digi_tags){
      this->consumes<StrawDigiCollection>(tag);
      this->consumes<StrawDigiADCWaveformCollection>(tag);
    }
    this->produces<StrawDigiCollection>();
    this->produces<StrawDigiADCWaveformCollection>();

    // calorimeter...
    for (const auto& tag: _calo_digi_tags){
      this->consumes<CaloDigiCollection>(tag);
    }
    this->produces<CaloDigiCollection>();

    // crv...
    for (const auto& tag: _crv_digi_tags) this->consumes<CrvDigiCollection>(tag);
    for (const auto& tag: _crv_digimc_tags) this->consumes<CrvDigiMCCollection>(tag);
    this->produces<CrvDigiCollection>();


    // mc truth
    if (_tracker_mc){
      for (const auto& tag: _tracker_digi_tags){
        this->consumes<StrawDigiMCCollection>(tag);
      }
      this->produces<StrawDigiMCCollection>();
    }
    if(!_crv_digimc_tags.empty()){
      this->produces<CrvDigiMCCollection>();
    }
  }

  void MergeDigis::produce(art::Event& event){
    // tracker: two easy steps:
    //   i) read all digis into a StrawDigiBundleCollection
    //  ii) defer collision resolution to that collection
    StrawDigiBundleCollection bundles;
    for (const auto& tag: _tracker_digi_tags){
      auto digi_handle = event.getValidHandle<StrawDigiCollection>(tag);
      auto adcs_handle = event.getValidHandle<StrawDigiADCWaveformCollection>(tag);
      if (_tracker_mc){
        auto dgmc_handle = event.getValidHandle<StrawDigiMCCollection>(tag);
        bundles.Append(*digi_handle, *adcs_handle, *dgmc_handle);
      }
      else{
        bundles.Append(*digi_handle, *adcs_handle);
      }
    }
    const auto& electronics = _tracker_conditions_handle.get(event.id());
    StrawDigiBundleCollection resolved;
    bundles.ResolveCollisions(electronics, resolved);
    auto digis = resolved.GetStrawDigiPtrs();
    auto adcs = resolved.GetStrawDigiADCWaveformPtrs();
    event.put(std::move(digis));
    event.put(std::move(adcs));
    if (_tracker_mc){
      auto dgmc = resolved.GetStrawDigiMCPtrs();
      event.put(std::move(dgmc));
    }

    // calorimeter: two easy steps:
    //   i) read all digis into a CaloDigiWrapperCollection
    //  ii) defer collision resolution to that collection
    CaloDigiWrapperCollection wrappers;
    for (const auto& tag: _calo_digi_tags){
      auto handle = event.getValidHandle<CaloDigiCollection>(tag);
      wrappers.Append(*handle);
    }
    CaloDigiWrapperCollection calo_resolved;
    wrappers.ResolveCollisions(_calo_max_adc, calo_resolved);
    auto calo_digis = calo_resolved.GetDigis();
    event.put(std::move(calo_digis));

    // crv...
    mergeCrvDigis(event);
  }

  void MergeDigis::mergeCrvDigis(art::Event& event)
  {
//assumptions
//(1) the order in which these collections were mixed is identical for CrvDigiMC and CrvDigi
//(2) the order in which these collections are merged is identical for CrvDigiMC and CrvDigi
//(3) ZS and NZS collections need to be kept separate  //TODO //FIXME

    bool hasCrvDigiMCs = !_crv_digimc_tags.empty();

    if(_crv_digi_tags.size()!=_crv_digimc_tags.size() && hasCrvDigiMCs)
      throw cet::exception("MergeDigis::mergeCrvDigis") << "mismatched number of CrvDigi and CrvDigiMC collections" << std::endl;

    size_t nCollections = _crv_digi_tags.size();

//output collections
    std::unique_ptr<mu2e::CrvDigiCollection> crvDigisNew(new mu2e::CrvDigiCollection);
    std::unique_ptr<mu2e::CrvDigiMCCollection> crvDigisMCNew(new mu2e::CrvDigiMCCollection);

//get calibration constants (needed for the pedestal)
    auto const& calib = _crv_calib_handle.get(event.id());

//Sort digis by channel number
//Keep digi and digiMC as a pair to keep both index-matched during ordering
//Or have no digiMCs at all
    typedef std::pair<mu2e::CrvDigi, std::optional<mu2e::CrvDigiMC> > DigiPair;
    std::map<size_t, std::vector<DigiPair> > digiPairMap;
    for(size_t iTag=0; iTag<nCollections; ++iTag)
    {
      const auto crvDigis = event.getValidHandle<CrvDigiCollection>(_crv_digi_tags.at(iTag));
      std::optional<const art::ValidHandle<CrvDigiMCCollection>> crvDigisMC = (hasCrvDigiMCs ? std::optional<const art::ValidHandle<CrvDigiMCCollection>>(event.getValidHandle<CrvDigiMCCollection>(_crv_digimc_tags.at(iTag))) : std::nullopt);

      if(hasCrvDigiMCs)
      {
        if(crvDigis->size()!=crvDigisMC.value()->size())
          throw cet::exception("MergeDigis::mergeCrvDigis") << "mismatched CrvDigi/CrvDigiMC collection sizes" << std::endl;
      }

      for(size_t i=0; i<crvDigis->size(); ++i)
      {
        const auto &digi=crvDigis->at(i);
        const CRSScintillatorBarIndex &barIndex = digi.GetScintillatorBarIndex();
        int SiPM = digi.GetSiPMNumber();
        size_t channel = barIndex.asUint()*CRVId::nChanPerBar + SiPM;

        if(hasCrvDigiMCs) digiPairMap[channel].emplace_back(digi,crvDigisMC.value()->at(i));
        else digiPairMap[channel].emplace_back(digi,std::nullopt);
      }
    }//loop through collections

//Loop through all channels
    for(auto digiPairsIter=digiPairMap.begin(); digiPairsIter!=digiPairMap.end(); ++digiPairsIter)
    {
      size_t channel  = digiPairsIter->first;
      std::vector<DigiPair> &digiPairs = digiPairsIter->second;

//Sort digis of a channel by TDC start time and keep digiMCs index-matched
      std::sort(digiPairs.begin(),digiPairs.end(), [](const DigiPair &d1, const DigiPair &d2) {return d1.first.GetStartTDC()<d2.first.GetStartTDC();});

//Loop through digis of a channel to solve overlap problems that may have occured during merging
      for(size_t i=1; i<digiPairs.size(); ++i)
      {
        const auto &digi1=digiPairs.at(i-1).first;
        const auto &digi2=digiPairs.at(i).first;
        const auto &digiMC1=digiPairs.at(i-1).second;
        const auto &digiMC2=digiPairs.at(i).second;
        size_t start1=digi1.GetStartTDC();
        size_t start2=digi2.GetStartTDC();
        size_t length1=digi1.GetADCs().size();
        size_t length2=digi2.GetADCs().size();
        if(start1+length1>start2)
        {
//found an overlap between digiPair_1 and digiPair_2
          double pedestal = calib.pedestal(channel);
          size_t lengthTotal=std::max(start1+length1,start2+length2)-start1;
          size_t startOverlap=start2-start1;

//combine ADCs and voltages
          std::vector<int16_t> ADCs;
          std::vector<double>  voltages;
          ADCs.resize(lengthTotal);
          voltages.resize(lengthTotal);
          for(size_t j=0; j<lengthTotal; ++j)
          {
            bool hasDigi1 = (j<length1);
            bool hasDigi2 = (startOverlap<=j && j<startOverlap+length2);
            if(hasDigi1 && !hasDigi2)
            {
              ADCs.at(j)=digi1.GetADCs().at(j);
              if(hasCrvDigiMCs) voltages.at(j)=digiMC1.value().GetVoltages().at(j);
            }
            if(!hasDigi1 && hasDigi2)
            {
              ADCs.at(j)=digi2.GetADCs().at(j-startOverlap);
              if(hasCrvDigiMCs) voltages.at(j)=digiMC2.value().GetVoltages().at(j-startOverlap);
            }
            if(hasDigi1 && hasDigi2)
            {
              //pedestal is included in both ADC values, but can be used only once. so it needs to be subtracted.
              int16_t ADC = digi1.GetADCs().at(j) + digi2.GetADCs().at(j-startOverlap) - pedestal;
              ADCs.at(j)=ADC;
              if(hasCrvDigiMCs)
              {
                double voltage = digiMC1.value().GetVoltages().at(j) + digiMC2.value().GetVoltages().at(j-startOverlap);
                voltages.at(j)=voltage;
              }
            }
          }
          for(size_t j=0; j<ADCs.size(); ++j)
          {
            if(ADCs.at(j)>CRVMaxADC) ADCs.at(j)=CRVMaxADC;
          }

//combine other values
          uint16_t startTDC     = start1;
          bool     NZS          = digi1.IsNZS();
          bool     oddTimestamp = digi1.HasOddTimestamp() || digi2.HasOddTimestamp();
          mu2e::CRSScintillatorBarIndex scintillatorBarIndex = digi1.GetScintillatorBarIndex();
          uint8_t  SiPMNumber   = digi1.GetSiPMNumber();
          uint8_t  ROC          = digi1.GetROC();
          uint8_t  FEB          = digi1.GetFEB();
          uint8_t  FEBchannel   = digi1.GetFEBchannel();

          double   startTime    = 0;
          double   TDC0Time     = 0;
          std::vector<art::Ptr<CrvStep> > steps;
          art::Ptr<SimParticle> simParticle;

          if(hasCrvDigiMCs)
          {
            startTime = digiMC1.value().GetStartTime();
            TDC0Time  = digiMC1.value().GetTDC0Time();

            //no duplicate crvSteps expected
            steps.insert(steps.begin(),digiMC1.value().GetCrvSteps().begin(),digiMC1.value().GetCrvSteps().end());
            steps.insert(steps.begin(),digiMC2.value().GetCrvSteps().begin(),digiMC2.value().GetCrvSteps().end());

            //find simparticle with biggest contribution of visible deposited energy
            std::map<art::Ptr<SimParticle>,double> simParticleMap;
            for(auto stepsIter=steps.begin(); stepsIter!=steps.end(); ++stepsIter)
            {
              if((*stepsIter).isNonnull()) simParticleMap[(*stepsIter)->simParticle()]+=(*stepsIter)->visibleEDep();
            }
            double simParticleDepEnergy=0;
            for(auto simParticleIter=simParticleMap.begin(); simParticleIter!=simParticleMap.end(); ++simParticleIter)
            {
              if(simParticleIter->second>simParticleDepEnergy)
              {
                simParticleDepEnergy=simParticleIter->second;
                simParticle=simParticleIter->first;
              }
            }
          }

//erase both digiPairs (i-1 and i)
          digiPairs.erase(digiPairs.begin()+i);
          digiPairs.erase(digiPairs.begin()+i-1);

//insert new digiPair
          if(hasCrvDigiMCs)
          {
            DigiPair newDigiPair(mu2e::CrvDigi(ADCs,startTDC,NZS,oddTimestamp,scintillatorBarIndex,SiPMNumber,ROC,FEB,FEBchannel),
                                 mu2e::CrvDigiMC(voltages,steps,simParticle,startTime,TDC0Time,NZS,scintillatorBarIndex,SiPMNumber));
            digiPairs.insert(digiPairs.begin()+i-1, newDigiPair);
          }
          else
          {
            DigiPair newDigiPair(mu2e::CrvDigi(ADCs,startTDC,NZS,oddTimestamp,scintillatorBarIndex,SiPMNumber,ROC,FEB,FEBchannel),std::nullopt);
            digiPairs.insert(digiPairs.begin()+i-1, newDigiPair);
          }

//the pair i (formerly pair i+1) needs to be checked with the newly created pair i-1
          --i; //for loop will increase i by 1 so that we are again back at the same index.

        }//found an overlap
      }//loop over digi pairs

      for(size_t i=0; i<digiPairs.size(); ++i)
      {
        const auto &digi=digiPairs.at(i).first;
        const auto &digiMC=digiPairs.at(i).second;

        crvDigisNew->emplace_back(digi);
        if(hasCrvDigiMCs) crvDigisMCNew->emplace_back(digiMC.value());
      }
    }//loop over channel

    event.put(std::move(crvDigisNew));
    if(hasCrvDigiMCs) event.put(std::move(crvDigisMCNew));
  }//mergeCrvDigis

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::MergeDigis);
