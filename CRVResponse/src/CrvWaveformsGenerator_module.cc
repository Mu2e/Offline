//
// A module to create CRV waveforms from CRV SiPMCharges
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CRVResponse/inc/MakeCrvWaveforms.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CRVConditions/inc/CRVCalib.hh"
#include "Offline/CRVConditions/inc/CRVDigitizationPeriod.hh"
#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/MCDataProducts/inc/CrvSiPMCharges.hh"
#include "Offline/MCDataProducts/inc/CrvDigiMC.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/SeedService/inc/SeedService.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Random/Randomize.h"

#include <string>

#include <TMath.h>

namespace mu2e
{
  class CrvWaveformsGenerator : public art::EDProducer
  {

    public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config
    {
      fhicl::Atom<std::string> crvSiPMChargesModuleLabel{Name("crvSiPMChargesModuleLabel")};
      fhicl::Atom<std::string> singlePEWaveformFileName{Name("singlePEWaveformFileName")};
      fhicl::Atom<double> digitizationStart{Name("digitizationStart"), Comment("start of digitization after DAQ event window start")}; //200ns (400ns...425ns after POT)
      fhicl::Atom<double> digitizationEnd{Name("digitizationEnd"), Comment("end of digitization after DAQ event window start")}; //1500ns (1700ns...1725ns after POT)
      fhicl::Atom<art::InputTag> eventWindowMarkerTag{Name("eventWindowMarkerTag"), Comment("EventWindowMarker producer"),"EWMProducer" };
      fhicl::Atom<art::InputTag> protonBunchTimeMCTag{Name("protonBunchTimeMCTag"), Comment("ProtonBunchTimeMC producer"),"EWMProducer" };
      fhicl::Atom<double> minVoltage{Name("minVoltage")};                       //0.022V (corresponds to 3.5PE)
      fhicl::Atom<double> noise{Name("noise")};                                 //4.0e-4V
      fhicl::Atom<double> timeOffsetScale{Name("timeOffsetScale")};             // 1.0ns (scale factor applied to the database values)
      fhicl::Atom<double> timeOffsetCutoffLow{Name("timeOffsetCutoffLow")};     //-3.0ns
      fhicl::Atom<double> timeOffsetCutoffHigh{Name("timeOffsetCutoffHigh")};   // 3.0ns
                                                                                //note: if measured time offsets are used,
                                                                                //the cutoffs should be set to the maximum values
      fhicl::Atom<bool> useTimeOffsetDB{Name("useTimeOffsetDB")};  //false, will be applied at reco
      fhicl::Atom<double> singlePEWaveformMaxTime{Name("singlePEWaveformMaxTime")};        //100ns
      fhicl::Atom<double> singlePEWaveformPrecision{Name("singlePEWaveformPrecision")};    //1.0 ns
      fhicl::Atom<double> singlePEWaveformStretchFactor{Name("singlePEWaveformStretchFactor")};    //1.047
      fhicl::Atom<double> singlePEReferenceCharge{Name("singlePEReferenceCharge")}; //2.652e-13 C (the charge which was used to generate the above 1PE waveform)

      fhicl::Atom<int> numberSamplesZS{Name("numberSamplesZS")}; //12
      fhicl::Atom<int> numberSamplesNZS{Name("numberSamplesNZS")}; //134
      fhicl::Atom<bool> simulateNZS{Name("simulateNZS")}; //false
      fhicl::Atom<int> prescalingFactorNZS{Name("prescalingFactorNZS")}; //10
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit CrvWaveformsGenerator(const Parameters& conf);
    void produce(art::Event& e) override;
    void beginRun(art::Run &run) override;

    private:
    std::string _crvSiPMChargesModuleLabel;
    std::string _singlePEWaveformFileName;
    art::InputTag _eventWindowMarkerTag;
    art::InputTag _protonBunchTimeMCTag;

    boost::shared_ptr<mu2eCrv::MakeCrvWaveforms> _makeCrvWaveforms;

    double                              _digitizationStart;
    double                              _digitizationEnd;
    double                              _minVoltage;
    double                              _noise;
    double                              _timeOffsetScale;
    double                              _timeOffsetCutoffLow;
    double                              _timeOffsetCutoffHigh;
    bool                                _useTimeOffsetDB;
    double                              _singlePEWaveformMaxTime;

    int                                 _numberSamplesZS;
    int                                 _numberSamplesNZS;
    bool                                _simulateNZS;
    int                                 _prescalingFactorNZS;

    CLHEP::HepRandomEngine&             _engine;
    CLHEP::RandFlat                     _randFlat;
    CLHEP::RandGaussQ                   _randGaussQ;

    ProditionsHandle<CRVCalib>         _calib;
    ProditionsHandle<CRVOrdinal>       _crvChannelMap;
    std::vector<double> _digitizationPointShiftFEBs;

    static constexpr int nDigiPeriods = 4; //number of digitization periods to check for charges

    struct ChargeCluster
    {
      std::vector<std::pair<double,double> > timesAndCharges;
    };

    void FindChargeClusters(const std::vector<CrvSiPMCharges::SingleCharge> &timesAndCharges,
                            std::vector<ChargeCluster> &chargeClusters, double timeOffset);
    bool SingleWaveformStart(std::vector<double> &fullWaveform, size_t i);
  };

  CrvWaveformsGenerator::CrvWaveformsGenerator(const Parameters& conf) :
    art::EDProducer{conf},
    _crvSiPMChargesModuleLabel(conf().crvSiPMChargesModuleLabel()),
    _singlePEWaveformFileName(conf().singlePEWaveformFileName()),
    _eventWindowMarkerTag(conf().eventWindowMarkerTag()),
    _protonBunchTimeMCTag(conf().protonBunchTimeMCTag()),
    _digitizationStart(conf().digitizationStart()),
    _digitizationEnd(conf().digitizationEnd()),
    _minVoltage(conf().minVoltage()),
    _noise(conf().noise()),
    _timeOffsetScale(conf().timeOffsetScale()),
    _timeOffsetCutoffLow(conf().timeOffsetCutoffLow()),
    _timeOffsetCutoffHigh(conf().timeOffsetCutoffHigh()),
    _useTimeOffsetDB(conf().useTimeOffsetDB()),
    _singlePEWaveformMaxTime(conf().singlePEWaveformMaxTime()),
    _numberSamplesZS(conf().numberSamplesZS()),
    _numberSamplesNZS(conf().numberSamplesNZS()),
    _simulateNZS(conf().simulateNZS()),
    _prescalingFactorNZS(conf().prescalingFactorNZS()),
    _engine{createEngine(art::ServiceHandle<SeedService>()->getSeed())},
    _randFlat{_engine},
    _randGaussQ{_engine}
  {
    double singlePEWaveformPrecision(conf().singlePEWaveformPrecision());
    double singlePEWaveformStretchFactor(conf().singlePEWaveformStretchFactor());
    double singlePEReferenceCharge(conf().singlePEReferenceCharge());
    ConfigFileLookupPolicy configFile;
    _singlePEWaveformFileName = configFile(_singlePEWaveformFileName);
    _makeCrvWaveforms = boost::shared_ptr<mu2eCrv::MakeCrvWaveforms>(new mu2eCrv::MakeCrvWaveforms());
    _makeCrvWaveforms->LoadSinglePEWaveform(_singlePEWaveformFileName, singlePEWaveformPrecision, singlePEWaveformStretchFactor,
                                            _singlePEWaveformMaxTime, singlePEReferenceCharge);
    produces<CrvDigiMCCollection>();
    if(_simulateNZS) produces<CrvDigiMCCollection>("NZS");
  }

  void CrvWaveformsGenerator::beginRun(art::Run &run)
  {
  }

  void CrvWaveformsGenerator::produce(art::Event& event)
  {
    art::Handle<EventWindowMarker> eventWindowMarker;
    event.getByLabel(_eventWindowMarkerTag,eventWindowMarker);
    EventWindowMarker::SpillType spillType = eventWindowMarker->spillType();
    double eventWindowLength = eventWindowMarker->eventLength(); //onspill: 1675ns/1700ns, offspill: 100000ns

    //offspill
    double eventWindowStart=0;
    double digitizationStart=0;
    double digitizationEnd=eventWindowLength;

    //onspill
    if(spillType==EventWindowMarker::SpillType::onspill)
    {
      art::Handle<ProtonBunchTimeMC> protonBunchTimeMC;
      event.getByLabel(_protonBunchTimeMCTag, protonBunchTimeMC);
      eventWindowStart = -protonBunchTimeMC->pbtime_; //200ns...225ns

      //for ZS data
      digitizationStart=eventWindowStart+_digitizationStart; //400ns...425ns
      digitizationEnd=eventWindowStart+_digitizationEnd; //1700ns...1725ns
    }

    auto const& calib = _calib.get(event.id());
    auto const& crvChannelMap = _crvChannelMap.get(event.id());

    std::unique_ptr<CrvDigiMCCollection> crvDigiMCCollection(new CrvDigiMCCollection);
    std::unique_ptr<CrvDigiMCCollection> crvDigiMCCollectionNZS(new CrvDigiMCCollection);

    art::Handle<CrvSiPMChargesCollection> crvSiPMChargesCollection;
    event.getByLabel(_crvSiPMChargesModuleLabel,"",crvSiPMChargesCollection);

    _digitizationPointShiftFEBs.clear();
    unsigned int nFEBs = CRVId::nFEBPerROC*CRVId::nROC;  //a few more FEBs than we need
    for(unsigned int i=0; i<nFEBs; ++i)
    {
      //the closest digitization point with respect to a certain time
      //can happen anywhere within the digitization period of 12.5ns.
      //this time is different for each FEB.
      _digitizationPointShiftFEBs.emplace_back(_randFlat.fire()*CRVDigitizationPeriod);
    }

    for(CrvSiPMChargesCollection::const_iterator iter=crvSiPMChargesCollection->begin();
        iter!=crvSiPMChargesCollection->end(); iter++)
    {
      int SiPM = iter->GetSiPMNumber();
      CRSScintillatorBarIndex barIndex = iter->GetScintillatorBarIndex();

      //get FEB from database
      uint16_t offlineChannel = barIndex.asUint()*4 + SiPM;
      CRVROC   onlineChannel  = crvChannelMap.online(offlineChannel);
      uint16_t FEB            = onlineChannel.FEB();
      uint16_t FEBchannel     = onlineChannel.FEBchannel();  //will be used for the NZS data below

      double timeOffset=0.0;
      if(_useTimeOffsetDB)
      {
        //the FEBs will be synchronized to account for cable length differences etc.,
        //but there may still be small time differences between the FEBs.
        //get the numbers from the database (either measured values of random values)
        timeOffset = calib.timeOffset(barIndex.asUint()*CRVId::nChanPerBar + SiPM);
        timeOffset*=_timeOffsetScale;   //random time offsets can be scaled to a wider or smaller spread
        if(timeOffset<_timeOffsetCutoffLow)  timeOffset=_timeOffsetCutoffLow;  //random time offsets can be cutoff at some limit
        if(timeOffset>_timeOffsetCutoffHigh) timeOffset=_timeOffsetCutoffHigh;
      }

      //need to find where this FEB's TDC=0 (first point after the event window start, i.e. first clock tick after POT) is located with respect to the global time
      //can be anywhere within the digitization period
      double digitizationPointShiftFEB=_digitizationPointShiftFEBs[FEB];
      double TDC0time=eventWindowStart+digitizationPointShiftFEB;  //that's the time when TDC=0 for this FEB

      //charges and times
      const std::vector<CrvSiPMCharges::SingleCharge> &timesAndCharges = iter->GetCharges();

      //zero suppressed data
      std::vector<ChargeCluster> chargeClusters;
      FindChargeClusters(timesAndCharges, chargeClusters, timeOffset);

      for(size_t iCluster=0; iCluster<chargeClusters.size(); ++iCluster)
      {
        //if the number of charges in this cluster cannot achieve the minimum voltage, skip this cluster
        if(chargeClusters[iCluster].timesAndCharges.size()*_makeCrvWaveforms->GetSinglePEMaxVoltage()<_minVoltage) continue;

        //find the TDC time when the first charge occurs (adjusted for this FEB)
        double firstChargeTime=chargeClusters[iCluster].timesAndCharges.front().first;
        firstChargeTime-=1.0*CRVDigitizationPeriod;  //start somewhere before the first charge
        double TDCstartTime=ceil((firstChargeTime-TDC0time)/CRVDigitizationPeriod) * CRVDigitizationPeriod + TDC0time;

        //first create the full waveform
        std::vector<double> fullWaveform;
        _makeCrvWaveforms->MakeWaveform(chargeClusters[iCluster].timesAndCharges,
                                        fullWaveform, TDCstartTime, CRVDigitizationPeriod);
        _makeCrvWaveforms->AddElectronicNoise(fullWaveform, _noise, _randGaussQ);

        //break the waveform apart into short pieces (_numberSampleZS) and apply the zero suppression
        //don't digitize outside of digitizationStart and digitizationEnd
        for(size_t i=0; i<fullWaveform.size(); ++i)
        {
          if(SingleWaveformStart(fullWaveform, i)) //acts as a zero suppression
          {
            //start new single waveform
            double digiStartTime=TDCstartTime+i*CRVDigitizationPeriod;
            if(digiStartTime<digitizationStart) continue; //digis cannot start before the digitization interval
            if(digiStartTime>digitizationEnd) continue; //digis cannot start after the digitization interval
//            if(digiStartTime+(_numberSamplesZS-1)*CRVDigitizationPeriod>digitizationEnd) continue; //digis cannot end after the digitization interval

            //collect voltages
            std::vector<double> voltages;
            voltages.resize(_numberSamplesZS);
            for(int singleWaveformIndex=0; singleWaveformIndex<_numberSamplesZS; ++i, ++singleWaveformIndex)
            {
//              if(i<fullWaveform.size() && TDCstartTime+i*CRVDigitizationPeriod<=digitizationEnd) voltages[singleWaveformIndex]=fullWaveform[i];  //cuts off pulse in the middle of the hit
              if(i<fullWaveform.size()) voltages[singleWaveformIndex]=fullWaveform[i];
              else voltages[singleWaveformIndex]=0.0;  //so that all unused single waveform samples are set to zero
            }

            //collect CrvSteps and SimParticles responsible for this single waveform
            std::set<art::Ptr<CrvStep> > steps;  //use a set to remove dublicate steppoints
            std::map<art::Ptr<SimParticle>, int> simparticles;
            for(size_t j=0; j<timesAndCharges.size(); ++j)
            {
              if(timesAndCharges[j]._time>=digiStartTime-_singlePEWaveformMaxTime &&
                 timesAndCharges[j]._time<=digiStartTime+_numberSamplesZS*CRVDigitizationPeriod)
              {
                steps.insert(timesAndCharges[j]._step);
                if(timesAndCharges[j]._step.isNonnull()) simparticles[timesAndCharges[j]._step->simParticle()]++;
              }
            }

            //loop through the steps to fill the single waveform
            std::vector<art::Ptr<CrvStep> > stepVector;
            std::set<art::Ptr<CrvStep> >::iterator stepIter;
            for(stepIter=steps.begin(); stepIter!=steps.end(); stepIter++) stepVector.push_back(*stepIter);

            //find the most likely SimParticle
            //if no SimParticle was recorded for this single waveform, then it was caused either by noise hits (if the threshold is low enough),
            //or is the tail end of the peak. in that case, _simparticle will be null (set by the default constructor of art::Ptr)
            art::Ptr<SimParticle> simParticle;
            std::map<art::Ptr<SimParticle>,int >::iterator simparticleIter;
            int simparticleCount=0;
            for(simparticleIter=simparticles.begin(); simparticleIter!=simparticles.end(); ++simparticleIter)
            {
              if(simparticleIter->second>simparticleCount)
              {
                simparticleCount=simparticleIter->second;
                simParticle=simparticleIter->first;
              }
            }

            --i;
            crvDigiMCCollection->emplace_back(voltages, stepVector, simParticle, digiStartTime, TDC0time, false, barIndex, SiPM);
          }
        } //waveform
      } //charge cluster for zero suppressed data

      if(_simulateNZS &&
         event.event()%CRVId::nChanPerFPGA==FEBchannel%CRVId::nChanPerFPGA &&  //only one of the 16 channels of an FPGA of an FEB will record the NZS data
         event.event()%_prescalingFactorNZS==0) //a pre-scaling of the NZS data will be applied
      {
        //non-zero suppressed data
        //use only charges for the first 134 samples (_numberSamplesNZS)
        std::vector<std::pair<double,double> > timesAndChargesNZS;
        for(size_t iCharge=0; iCharge<timesAndCharges.size(); ++iCharge)
        {
          if(timesAndCharges[iCharge]._time>=TDC0time-_singlePEWaveformMaxTime &&
             timesAndCharges[iCharge]._time<=TDC0time+CRVDigitizationPeriod*_numberSamplesNZS)
          {
            timesAndChargesNZS.emplace_back(timesAndCharges[iCharge]._time,timesAndCharges[iCharge]._charge);
          }
        }

        //create NZS waveform
        std::vector<double> fullWaveformNZS;
        _makeCrvWaveforms->MakeWaveform(timesAndChargesNZS, fullWaveformNZS, TDC0time, CRVDigitizationPeriod);
        //record first 134 samples (_numberSamplesNZS)
        //assumed to be always within eventWindowStart and eventWindowEnd
        fullWaveformNZS.resize(_numberSamplesNZS);
        _makeCrvWaveforms->AddElectronicNoise(fullWaveformNZS, _noise, _randGaussQ);

        //use empty step vector and empty simParticle for NZS
        std::vector<art::Ptr<CrvStep> > stepVectorNZS;
        art::Ptr<SimParticle> simParticleNZS;  //will be set to Null automatically

        crvDigiMCCollectionNZS->emplace_back(fullWaveformNZS, stepVectorNZS, simParticleNZS, TDC0time, TDC0time, true, barIndex, SiPM);
      } //NZS
    } //SiPMCharges of one channel

    event.put(std::move(crvDigiMCCollection));
    if(_simulateNZS) event.put(std::move(crvDigiMCCollectionNZS),"NZS");
  } // end produce

  void CrvWaveformsGenerator::FindChargeClusters(const std::vector<CrvSiPMCharges::SingleCharge> &timesAndCharges,
                                                 std::vector<ChargeCluster> &chargeClusters, double timeOffset)
  {
    chargeClusters.reserve(timesAndCharges.size());
    for(size_t i=0; i<timesAndCharges.size(); ++i)
    {
      //No check whether times are within digitizationStart-_digitizationMargin and digitizationEnd
      double timeTmp=timesAndCharges[i]._time+timeOffset;  //apply timeOffset to account for inaccuraries in the FEB time calibration

      if(chargeClusters.empty())
      {
        chargeClusters.resize(1);
        chargeClusters.back().timesAndCharges.reserve(timesAndCharges.size());
      }
      else
      {
        if(timeTmp-chargeClusters.back().timesAndCharges.back().first>_singlePEWaveformMaxTime+nDigiPeriods*CRVDigitizationPeriod)
        //if the difference b/w the time of the next charge and the time of the last charge
        //is greater than a full single PE waveform plus four additional digitization periods
        //-->start a new charge cluster
        {
          chargeClusters.resize(chargeClusters.size()+1);
          chargeClusters.back().timesAndCharges.reserve(timesAndCharges.size());
        }
      }
      chargeClusters.back().timesAndCharges.emplace_back(timeTmp,timesAndCharges[i]._charge);
    }
  }

  bool CrvWaveformsGenerator::SingleWaveformStart(std::vector<double> &fullWaveform, size_t i)
  {
    //test up to 2 points before and after i, if the waveform is/was above _minVoltage
    //to decide whether to start or continue to record digis
    for(size_t j=std::max<size_t>(i,2)-2; j<std::min<size_t>(i+3,fullWaveform.size()); ++j)
    {
      if(fullWaveform[j]>_minVoltage) return true;
    }

    return false;
  }

} // end namespace mu2e

using mu2e::CrvWaveformsGenerator;
DEFINE_ART_MODULE(CrvWaveformsGenerator)
