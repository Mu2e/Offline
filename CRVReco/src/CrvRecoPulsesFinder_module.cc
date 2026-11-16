//
// A module to find the times and PEs of Crv pulses
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CRVReco/inc/MakeCrvRecoPulses.hh"

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"

#include "Offline/CRVConditions/inc/CRVDigitizationPeriod.hh"
#include "Offline/CRVConditions/inc/CRVCalib.hh"
#include "Offline/CRVConditions/inc/CRVStatus.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulseFlags.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <string>
#include <numeric>

#include <TMath.h>

namespace mu2e
{
  class CrvRecoPulsesFinder : public art::EDProducer
  {

    public:
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> crvDigiModuleLabel{Name("crvDigiModuleLabel"), Comment("module label for CrvDigis")};
      fhicl::Atom<bool> NZSdata{Name("NZSdata"), Comment("use digis with the NZS instance label")};  //false

      //for calibration
      fhicl::Atom<float>  minPedestalDifference{Name("minPedestalDifference"), Comment("minimum difference between first ADC sample and DB pedestal"), -10000};
      fhicl::Atom<float>  maxPedestalDifference{Name("maxPedestalDifference"), Comment("maximum difference between first ADC sample and DB pedestal"),  10000};
      fhicl::Atom<float>  minAdjustedPedestalDifference{Name("minAdjustedPedestalDifference"), Comment("minimum difference between adjusted pedestal and DB pedestal"), -10000};
      fhicl::Atom<float>  maxAdjustedPedestalDifference{Name("maxAdjustedPedestalDifference"), Comment("maximum difference between adjusted pedestal and DB pedestal"),  10000};
      fhicl::Atom<size_t> maxPulseLength{Name("maxPulseLength"), Comment("maximum pulse length for waveform to be considered"),  10000};
      fhicl::Atom<bool>   pedestalAdjustment{Name("pedestalAdjustment"), Comment("precision pedestal adjustment for dark counts"),  false};
      fhicl::Atom<size_t> pedestalAdjustmentFirstPeakPos{Name("pedestalAdjustmentFirstPeakPos"), Comment("first allowed peak pos for precision pedestal adjustment for dark counts"), 2};
      fhicl::Atom<size_t> pedestalAdjustmentLastPeakPos{Name("pedestalAdjustmentLastPeakPos"), Comment("last allowed peak pos for precision pedestal adjustment for dark counts"), 3};
      fhicl::Atom<bool>   usePedestalAdjustmentSample10{Name("usePedestalAdjustmentSample10"), Comment("use sample 10 for pedestal adjustment for precision pedestal adjustment for dark counts"), false};
      fhicl::Atom<bool>   usePedestalAdjustmentSample11{Name("usePedestalAdjustmentSample11"), Comment("use sample 11 for pedestal adjustment for precision pedestal adjustment for dark counts"), true};
      fhicl::Atom<float>  pedestalAdjustmentFactorSample10{Name("pedestalAdjustmentFactorSample10"), Comment("pedestal adjustment factor for sample 10 for precision pedestal adjustment for dark counts"), 0.015};
      fhicl::Atom<float>  pedestalAdjustmentFactorSample11{Name("pedestalAdjustmentFactorSample11"), Comment("pedestal adjustment factor for sample 11 for precision pedestal adjustment for dark counts"), 0.023};
      fhicl::Atom<float>  pedestalAdjustmentMaxFluctuation{Name("pedestalAdjustmentMaxFluctuation"), Comment("maximum pedestal fluctuation for precision pedestal adjustment for dark counts"), 5};

      fhicl::Atom<float> pedestalOffset{Name("pedestalOffset"), Comment("offset added to pedestal"), 0};
      fhicl::Atom<float> pedestalUndershootThreshold{Name("pedestalUndershootThreshold"), Comment("if 1st ADC sample is below DB pedestal by this threshold value, use 1st ADC sample as pedestal. E.g. -10000 forces it to always use the 1st ADC sample.")};  //6.0
      fhicl::Atom<float> minADCdifference{Name("minADCdifference"), Comment("minimum ADC difference above pedestal to be considered for reconstruction")};  //40
      fhicl::Atom<float> defaultBeta{Name("defaultBeta"), Comment("initialization value for fit and default value for invalid fits (regular pulses: 19.0ns, dark counts for calibration: 12.0ns)")};
      fhicl::Atom<float> minBeta{Name("minBeta"), Comment("smallest accepted beta for valid fit [ns]")};  //5.0ns
      fhicl::Atom<float> maxBeta{Name("maxBeta"), Comment("largest accepted beta for valid fit [ns]")}; //50.0ns
      fhicl::Atom<float> maxTimeDifference{Name("maxTimeDifference"), Comment("largest accepted difference between time of largest ADC value and fitted peak [ns]")}; //20.0ns
      fhicl::Atom<float> minPulseHeightRatio{Name("minPulseHeightRatio"), Comment("smallest accepted ratio between largest ADC value and fitted peak")}; //0.7
      fhicl::Atom<float> maxPulseHeightRatio{Name("maxPulseHeightRatio"), Comment("largest accepted ratio between largest ADC value and fitted peak")}; //1.5
      fhicl::Atom<float> LEtimeFactor{Name("LEtimeFactor"), Comment("time of leading edge is peakTime-LEtimeFactor*beta (0.985,1.385,1.587 for a leading edge of 0.5,0.2,0.1 pulse height")};
      fhicl::Atom<int>   samplesBefore{Name("samplesBefore"), Comment("number of samples before the peak to be included in fit")};
      fhicl::Atom<int>   samplesAfter{Name("samplesAfter"), Comment("number of samples after the peak to be included in fit")};
      fhicl::Atom<art::InputTag> eventWindowMarkerTag{Name("eventWindowMarkerTag"), Comment("EventWindowMarker producer"),"EWMProducer"};
      fhicl::Atom<art::InputTag> protonBunchTimeTag{Name("protonBunchTimeTag"), Comment("ProtonBunchTime producer"),"EWMProducer"};
      fhicl::Atom<bool> useTimeOffsetDB{Name("useTimeOffsetDB"), Comment("apply time offsets from the DB")}; //true
      fhicl::Atom<float> timeOffsetScale{Name("timeOffsetScale"), Comment("scale factor for time offsets from database (use 1.0, if measured values)")}; //1.0
      fhicl::Atom<bool> ignoreChannels{Name("ignoreChannels"), Comment("ignore channels that have status 2 (bit 1) in CRVstatus DB")}; //true
    };

    typedef art::EDProducer::Table<Config> Parameters;

    explicit CrvRecoPulsesFinder(const Parameters& config);
    void produce(art::Event& e) override;
    void beginJob() override;
    void beginRun(art::Run &run) override;
    void endJob() override;

    private:
    boost::shared_ptr<mu2eCrv::MakeCrvRecoPulses> _makeCrvRecoPulses;

    std::string _crvDigiModuleLabel;
    bool        _NZSdata;
    float       _minPedestalDifference, _maxPedestalDifference;
    float       _minAdjustedPedestalDifference, _maxAdjustedPedestalDifference;
    size_t      _maxPulseLength;
    bool        _pedestalAdjustment;
    size_t      _pedestalAdjustmentFirstPeakPos, _pedestalAdjustmentLastPeakPos;
    bool        _usePedestalAdjustmentSample10, _usePedestalAdjustmentSample11;
    float       _pedestalAdjustmentFactorSample10, _pedestalAdjustmentFactorSample11;
    float       _pedestalAdjustmentMaxFluctuation;
    float       _pedestalOffset;
    float       _pedestalUndershootThreshold;
    art::InputTag _eventWindowMarkerTag;
    art::InputTag _protonBunchTimeTag;

    bool  _useTimeOffsetDB;
    float _timeOffsetScale;

    bool  _ignoreChannels;

    ProditionsHandle<CRVCalib>  _calib;
    ProditionsHandle<CRVStatus> _sipmStatus;
    ProditionsHandle<EventTiming> _eventTiming;
  };


  CrvRecoPulsesFinder::CrvRecoPulsesFinder(const Parameters& conf) :
    art::EDProducer(conf),
    _crvDigiModuleLabel(conf().crvDigiModuleLabel()),
    _NZSdata(conf().NZSdata()),
    _minPedestalDifference(conf().minPedestalDifference()),
    _maxPedestalDifference(conf().maxPedestalDifference()),
    _minAdjustedPedestalDifference(conf().minAdjustedPedestalDifference()),
    _maxAdjustedPedestalDifference(conf().maxAdjustedPedestalDifference()),
    _maxPulseLength(conf().maxPulseLength()),
    _pedestalAdjustment(conf().pedestalAdjustment()),
    _pedestalAdjustmentFirstPeakPos(conf().pedestalAdjustmentFirstPeakPos()),
    _pedestalAdjustmentLastPeakPos(conf().pedestalAdjustmentLastPeakPos()),
    _usePedestalAdjustmentSample10(conf().usePedestalAdjustmentSample10()),
    _usePedestalAdjustmentSample11(conf().usePedestalAdjustmentSample11()),
    _pedestalAdjustmentFactorSample10(conf().pedestalAdjustmentFactorSample10()),
    _pedestalAdjustmentFactorSample11(conf().pedestalAdjustmentFactorSample11()),
    _pedestalAdjustmentMaxFluctuation(conf().pedestalAdjustmentMaxFluctuation()),
    _pedestalOffset(conf().pedestalOffset()),
    _pedestalUndershootThreshold(conf().pedestalUndershootThreshold()),
    _eventWindowMarkerTag(conf().eventWindowMarkerTag()),
    _protonBunchTimeTag(conf().protonBunchTimeTag()),
    _useTimeOffsetDB(conf().useTimeOffsetDB()),
    _timeOffsetScale(conf().timeOffsetScale()),
    _ignoreChannels(conf().ignoreChannels())
  {
    produces<CrvRecoPulseCollection>(_NZSdata?"NZS":"");
    _makeCrvRecoPulses=boost::shared_ptr<mu2eCrv::MakeCrvRecoPulses>(new mu2eCrv::MakeCrvRecoPulses(conf().minADCdifference(),
                                                                                                    conf().defaultBeta(),
                                                                                                    conf().minBeta(),
                                                                                                    conf().maxBeta(),
                                                                                                    conf().maxTimeDifference(),
                                                                                                    conf().minPulseHeightRatio(),
                                                                                                    conf().maxPulseHeightRatio(),
                                                                                                    conf().LEtimeFactor(),
                                                                                                    conf().samplesBefore(),
                                                                                                    conf().samplesAfter()));
  }

  void CrvRecoPulsesFinder::beginJob()
  {
  }

  void CrvRecoPulsesFinder::endJob()
  {
  }

  void CrvRecoPulsesFinder::beginRun(art::Run &run)
  {
  }

  void CrvRecoPulsesFinder::produce(art::Event& event)
  {
    std::unique_ptr<CrvRecoPulseCollection> crvRecoPulseCollection(new CrvRecoPulseCollection);

    double TDC0time = 0;

    art::Handle<EventWindowMarker> eventWindowMarker;
    if(event.getByLabel(_eventWindowMarkerTag,eventWindowMarker))
    {
      EventWindowMarker::SpillType spillType = eventWindowMarker->spillType();
      if(spillType==EventWindowMarker::SpillType::onspill)
      {
        art::Handle<ProtonBunchTime> protonBunchTime;
        event.getByLabel(_protonBunchTimeTag, protonBunchTime);
        if(protonBunchTime.isValid())
        {
          TDC0time = -protonBunchTime->pbtime_; //200ns...225ns (only for onspill)
        }
      }
    }

    art::Handle<CrvDigiCollection> crvDigiCollection;
    event.getByLabel(_crvDigiModuleLabel,(_NZSdata?"NZS":""),crvDigiCollection);

    auto const& calib = _calib.get(event.id());
    auto const& sipmStatus = _sipmStatus.get(event.id());
    auto const& eventTiming = _eventTiming.get(event.id());
    TDC0time += eventTiming.crvTrackerTimeOffset();

    size_t waveformIndex = 0;
    while(waveformIndex<crvDigiCollection->size())
    {
      const CrvDigi &digi = crvDigiCollection->at(waveformIndex);
      const CRSScintillatorBarIndex &barIndex = digi.GetScintillatorBarIndex();
      uint8_t SiPM = digi.GetSiPMNumber();
      uint8_t ROC = digi.GetROC();
      uint8_t FEB = digi.GetFEB();
      uint8_t FEBchannel = digi.GetFEBchannel();
      uint16_t startTDC = digi.GetStartTDC();
      std::vector<int16_t> ADCs=digi.GetADCs();
      std::vector<size_t> waveformIndices;
      waveformIndices.push_back(waveformIndex);

      //checking following digis whether they are a continuation of the current digis
      //if that is the case, append the next digis
      //requires that digis belonging to the same channel are grouped together and ordered by timestamp
      while(++waveformIndex<crvDigiCollection->size())
      {
        const CrvDigi &nextDigi = crvDigiCollection->at(waveformIndex);
        if(barIndex!=nextDigi.GetScintillatorBarIndex()) break;
        if(SiPM!=nextDigi.GetSiPMNumber()) break;
        if(startTDC+ADCs.size()!=nextDigi.GetStartTDC()) break;  //checks that the next waveform begins right after the previous one
        ADCs.insert(ADCs.end(),nextDigi.GetADCs().begin(),nextDigi.GetADCs().end()); //append the waveform
        waveformIndices.push_back(waveformIndex);
      }

      size_t channel = barIndex.asUint()*CRVId::nChanPerBar + SiPM;

      if(_ignoreChannels)
      {
        std::bitset<16> status(sipmStatus.status(channel));
        if(status.test(CRVStatus::Flags::ignoreChannel)) continue; //ignore this channel (bit 1)
      }

      //get default pedestal from DB
      double pedestal = calib.pedestal(channel);
      bool pedestalFromDB = true;

      //for dark counts (for calibration): ignore pulses where the first ADC sample is outside of a range away from the DB pedestal
      if(((float)ADCs.at(0))-pedestal<_minPedestalDifference || ((float)ADCs.at(0))-pedestal>_maxPedestalDifference) continue;

      //for dark counts (for calibration): ignore long pulses (i.e. waveforms with more than one hit), since their pedestals cannot be determined precisely.
      if(ADCs.size()>_maxPulseLength) continue;

      //for dark counts (for calibration): get a more precise pedestal using the 1st, 11th, 12th (indices 0, 10, 11) ADC sample.
      //the last ADC sample is typically in an undershoot which depends on the pulse height.
      //a correction for this undershoot needs to be applies.
      //this correction factor needs to be found independently and requires that the peak is at a certain index.
      //ignore pulse if peak is outside of index range or if the difference between 1st sample and adjusted 12th sample is too large.
      if(_pedestalAdjustment)
      {
        if(ADCs.size()<12 && (_usePedestalAdjustmentSample10 || _usePedestalAdjustmentSample11)) continue;

        auto ADCpeak=std::max_element(ADCs.begin(),ADCs.end());
        size_t  ADCpeakPos=std::distance(ADCs.begin(),ADCpeak);
        if(ADCpeakPos<_pedestalAdjustmentFirstPeakPos || ADCpeakPos>_pedestalAdjustmentLastPeakPos) continue;

        std::vector<float> samples;
        samples.push_back((float)ADCs.at(0));

        float ADCpeakHeight=*ADCpeak-ADCs.at(0);
        if(_usePedestalAdjustmentSample10) samples.push_back(ADCs.at(10)+_pedestalAdjustmentFactorSample10*ADCpeakHeight);
        if(_usePedestalAdjustmentSample11) samples.push_back(ADCs.at(11)+_pedestalAdjustmentFactorSample11*ADCpeakHeight);

        auto [minSample, maxSample] = std::minmax_element(samples.begin(),samples.end());
        if(*maxSample-*minSample>_pedestalAdjustmentMaxFluctuation) continue;
        float adjustedPedestal = std::accumulate(samples.begin(),samples.end(), 0.0)/samples.size();

        //if(fabs((float)ADCs.at(0)-adjustedSample11)>_pedestalAdjustmentMaxFluctuation) continue;
        //float adjustedPedestal=(ADCs.at(0)+adjustedSample11)/2.0;

        if(adjustedPedestal-pedestal<_minAdjustedPedestalDifference || adjustedPedestal-pedestal>_maxAdjustedPedestalDifference) continue;
        pedestal=adjustedPedestal;
        pedestalFromDB=false;
      }

      //for signals: get pedestal from 1st ADC sample instead of DB for pulses that seem to be in an undershoot.
      //adjustment using the last ADC sample as done for dark counts cannot be done, since signal pulses are much longer.
      //due to strong sample-to-sample fluctuations, this is not a very precise pedestal value but it is sufficient for pulses > 5PEs.
      //values for _pedestalUndershootThreshold at 100000 turns it off. a value of -100000 uses the first sample as pedestal for every pulse.
      if(pedestal-((float)ADCs.at(0))>_pedestalUndershootThreshold && !_pedestalAdjustment)  //apply only if not doing _pedestalAdjustement above
      {
        pedestal=ADCs.at(0);
        pedestalFromDB=false;
      }

      //global pedestal adjustment
      pedestal+=_pedestalOffset;
      if(_pedestalOffset!=0) pedestalFromDB=false;

      //get calibration values
      double calibPulseArea = calib.pulseArea(channel);
      double calibPulseHeight = calib.pulseHeight(channel);
      double timeOffset = 0.0;
      if(_useTimeOffsetDB)
      {
        timeOffset = calib.timeOffset(channel);
        timeOffset*=_timeOffsetScale;   //random time offsets can be scaled to a wider or smaller spread
      }

      //do the fit
      _makeCrvRecoPulses->SetWaveform(ADCs, startTDC, CRVDigitizationPeriod, pedestal, calibPulseArea, calibPulseHeight);

      //extract fit results
      size_t n = _makeCrvRecoPulses->GetPEs().size();
      for(size_t j=0; j<n; ++j)
      {
        //the TDC times were recorded with respect to the event window start.
        //need to shift the times back to the original time scale (i.e. microbunch time)
        double pulseTime   = _makeCrvRecoPulses->GetPulseTimes().at(j) + TDC0time + timeOffset;
        double LEtime      = _makeCrvRecoPulses->GetLEtimes().at(j) + TDC0time + timeOffset;
        float  PEs         = _makeCrvRecoPulses->GetPEs().at(j);
        float  PEsPulseHeight = _makeCrvRecoPulses->GetPEsPulseHeight().at(j);
        float  pulseHeight = _makeCrvRecoPulses->GetPulseHeights().at(j);
        float  pulseBeta   = _makeCrvRecoPulses->GetPulseBetas().at(j);
        float  pulseFitChi2= _makeCrvRecoPulses->GetPulseFitChi2s().at(j);

        bool   failedFit              = _makeCrvRecoPulses->GetFailedFits().at(j);
        bool   zeroNdf                = _makeCrvRecoPulses->GetZeroNdfs().at(j);
        CrvRecoPulseFlags flags;
        if(failedFit)              flags.set(CrvRecoPulseFlagEnums::failedFit);
        if(zeroNdf)                flags.set(CrvRecoPulseFlagEnums::zeroNdf);

        if(calibPulseArea<=0) {PEs=0; flags.set(CrvRecoPulseFlagEnums::noCalibConstPulseArea);}
        if(calibPulseHeight<=0) {PEsPulseHeight=0; flags.set(CrvRecoPulseFlagEnums::noCalibConstPulseHeight);}

        crvRecoPulseCollection->emplace_back(PEs, PEsPulseHeight, pulseTime, pulseHeight, pulseBeta, pulseFitChi2, LEtime, j, flags,
                                             waveformIndices, barIndex, SiPM, ROC, FEB, FEBchannel, pedestal, pedestalFromDB);
      }

    }

    event.put(std::move(crvRecoPulseCollection),(_NZSdata?"NZS":""));
  } // end produce

} // end namespace mu2e

using mu2e::CrvRecoPulsesFinder;
DEFINE_ART_MODULE(CrvRecoPulsesFinder)
