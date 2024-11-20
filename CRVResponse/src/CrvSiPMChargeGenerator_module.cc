//
// A module to create CRV SiPM charges from CRV photons
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CRVConditions/inc/CRVDigitizationPeriod.hh"
#include "Offline/CRVConditions/inc/CRVStatus.hh"
#include "Offline/CRVResponse/inc/MakeCrvSiPMCharges.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/MCDataProducts/inc/CrvPhotons.hh"
#include "Offline/MCDataProducts/inc/CrvSiPMCharges.hh"
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
#include <bitset>

#include <TMath.h>

namespace mu2e
{
  class CrvSiPMChargeGenerator : public art::EDProducer
  {
    public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config
    {
      fhicl::Atom<std::string> crvPhotonsModuleLabel{Name("crvPhotonsModuleLabel")};
      fhicl::Atom<int> nPixelsX{Name("nPixelsX")};                            //40
      fhicl::Atom<int> nPixelsY{Name("nPixelsY")};                            //40
      fhicl::Atom<double> overvoltage{Name("overvoltage")};                   //3.0V
      fhicl::Atom<double> timeConstant{Name("timeConstant")};                 //12.0ns
      fhicl::Atom<double> capacitance{Name("capacitance")};                   //8.84e-14F (per pixel)
      fhicl::Atom<double> digitizationStart{Name("digitizationStart"),
                                            Comment("start of digitization after DAQ event window start")};
                                            //200ns (400ns...425ns after POT)
      fhicl::Atom<double> digitizationEnd{Name("digitizationEnd"),
                                          Comment("end of digitization after DAQ event window start")};
                                            //1500ns (700ns...1725ns after POT)
      fhicl::Atom<double> digitizationStartMargin{Name("digitizationStartMargin"),
                                                  Comment("time window before digitization starts to account for photon travel time and electronics response.")};
                                                  //50.0ns  start recording earlier to account for electronics response times
      fhicl::Atom<int> numberSamplesNZS{Name("numberSamplesNZS")};           //134
      fhicl::Atom<bool> simulateNZS{Name("simulateNZS")};                    //false
      fhicl::Atom<art::InputTag> eventWindowMarkerTag{Name("eventWindowMarkerTag"), Comment("EventWindowMarker producer"),"EWMProducer" };
      fhicl::Atom<art::InputTag> protonBunchTimeMCTag{Name("protonBunchTimeMCTag"), Comment("ProtonBunchTimeMC producer"),"EWMProducer" };
      fhicl::Atom<bool> useSipmStatusDB{Name("useSipmStatusDB")};             //false (all channels will be simulated. channels with status bit 1 can be ignored in reco)
      fhicl::Sequence<fhicl::Sequence<int,2u> > inactivePixels{Name("inactivePixels")};      //{18,18},....,{21,21}
      fhicl::Atom<std::string> photonMapFileName{Name("photonMapFileName")};
      fhicl::Atom<double> avalancheProbParam1{Name("AvalancheProbParam1")};  //0.65
      fhicl::Atom<double> avalancheProbParam2{Name("AvalancheProbParam2")};  //2.7
      fhicl::Atom<double> trapType0Prob{Name("TrapType0Prob")};              //0
      fhicl::Atom<double> trapType1Prob{Name("TrapType1Prob")};              //0
      fhicl::Atom<double> trapType0Lifetime{Name("TrapType0Lifetime")};      //5ns
      fhicl::Atom<double> trapType1Lifetime{Name("TrapType1Lifetime")};      //50ns
      fhicl::Atom<double> thermalRate{Name("ThermalRate")};                  //1.0e-4 ns^-1   100kHz for entire SiPM
      fhicl::Atom<double> crossTalkProb{Name("CrossTalkProb")};              //0.04
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit CrvSiPMChargeGenerator(const Parameters& conf);
    void produce(art::Event& e);
    void beginRun(art::Run &run);

    private:
    std::string _crvPhotonsModuleLabel;
    int         _nPixelsX;
    int         _nPixelsY;
    double      _overvoltage;
    double      _timeConstant;
    double      _capacitance;
    double      _digitizationStart, _digitizationEnd, _digitizationStartMargin;
    int         _numberSamplesNZS;
    bool        _simulateNZS;
    art::InputTag _eventWindowMarkerTag;
    art::InputTag _protonBunchTimeMCTag;

    mu2e::ProditionsHandle<mu2e::CRVStatus> _sipmStatus;
    bool                                    _useSipmStatusDB;

    mu2eCrv::MakeCrvSiPMCharges::ProbabilitiesStruct _probabilities;
    std::vector<std::pair<int,int> >   _inactivePixels;

    boost::shared_ptr<mu2eCrv::MakeCrvSiPMCharges> _makeCrvSiPMCharges;

    CLHEP::HepRandomEngine& _engine;
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;

    std::string             _photonMapFileName;
    ConfigFileLookupPolicy  _resolveFullPath;
  };

  CrvSiPMChargeGenerator::CrvSiPMChargeGenerator(const Parameters& conf) :
    art::EDProducer{conf},
    _crvPhotonsModuleLabel(conf().crvPhotonsModuleLabel()),
    _nPixelsX(conf().nPixelsX()),
    _nPixelsY(conf().nPixelsY()),
    _overvoltage(conf().overvoltage()),
    _timeConstant(conf().timeConstant()),
    _capacitance(conf().capacitance()),
    _digitizationStart(conf().digitizationStart()),
    _digitizationEnd(conf().digitizationEnd()),
    _digitizationStartMargin(conf().digitizationStartMargin()),
    _numberSamplesNZS(conf().numberSamplesNZS()),
    _simulateNZS(conf().simulateNZS()),
    _eventWindowMarkerTag(conf().eventWindowMarkerTag()),
    _protonBunchTimeMCTag(conf().protonBunchTimeMCTag()),
    _useSipmStatusDB(conf().useSipmStatusDB()),
    _engine{createEngine(art::ServiceHandle<SeedService>()->getSeed())},
    _randFlat{_engine},
    _randPoissonQ{_engine},
    _photonMapFileName(conf().photonMapFileName())
  {
    produces<CrvSiPMChargesCollection>();
    const auto inactivePixelsTmp = conf().inactivePixels();
    _inactivePixels.resize(inactivePixelsTmp.size());
    for(size_t i=0; i<inactivePixelsTmp.size(); ++i)
    {
      _inactivePixels[i]=std::pair<int,int>(inactivePixelsTmp.at(i).at(0),inactivePixelsTmp.at(i).at(1));
    }
    _probabilities._avalancheProbParam1=conf().avalancheProbParam1();  //0.65
    _probabilities._avalancheProbParam2=conf().avalancheProbParam2();  //2.7
    _probabilities._trapType0Prob=conf().trapType0Prob();              //0
    _probabilities._trapType1Prob=conf().trapType1Prob();              //0
    _probabilities._trapType0Lifetime=conf().trapType0Lifetime();      //5.0ns
    _probabilities._trapType1Lifetime=conf().trapType1Lifetime();      //50.0ns
    _probabilities._thermalRate=conf().thermalRate();                  //1.0e-4 ns^-1   100kHz for entire SiPM
    _probabilities._crossTalkProb=conf().crossTalkProb();              //0.05

    std::string fullPhotonMapFileName(_resolveFullPath(_photonMapFileName));
    _makeCrvSiPMCharges = boost::shared_ptr<mu2eCrv::MakeCrvSiPMCharges>(new mu2eCrv::MakeCrvSiPMCharges(_randFlat, _randPoissonQ, fullPhotonMapFileName));
  }

  void CrvSiPMChargeGenerator::beginRun(art::Run &run)
  {
    _makeCrvSiPMCharges->SetSiPMConstants(_nPixelsX, _nPixelsY, _overvoltage, _timeConstant, _capacitance, _probabilities, _inactivePixels);
  }

  void CrvSiPMChargeGenerator::produce(art::Event& event)
  {
    std::unique_ptr<CrvSiPMChargesCollection> crvSiPMChargesCollection(new CrvSiPMChargesCollection);

    art::Handle<CrvPhotonsCollection> crvPhotonsCollection;
    event.getByLabel(_crvPhotonsModuleLabel,crvPhotonsCollection);

    art::Handle<EventWindowMarker> eventWindowMarker;
    event.getByLabel(_eventWindowMarkerTag,eventWindowMarker);
    EventWindowMarker::SpillType spillType = eventWindowMarker->spillType();
    double eventWindowLength = eventWindowMarker->eventLength(); //onspill: 1675ns/1700ns, offspill: 100000ns

    //offspill
    double eventWindowStart=0;
    double startTime=0;
    double endTime=eventWindowLength;

    //onspill
    if(spillType==EventWindowMarker::SpillType::onspill)
    {
      art::Handle<ProtonBunchTimeMC> protonBunchTimeMC;
      event.getByLabel(_protonBunchTimeMCTag, protonBunchTimeMC);
      eventWindowStart = -protonBunchTimeMC->pbtime_; //200ns...225ns

      startTime=eventWindowStart+_digitizationStart-_digitizationStartMargin; //300ns...325ns
      endTime=eventWindowStart+_digitizationEnd; //1700ns...1725ns

      if(_simulateNZS)
      {
        startTime=eventWindowStart-_digitizationStartMargin; //100ns...125ns
        endTime=eventWindowStart+_numberSamplesNZS*CRVDigitizationPeriod; //1875ns...1900s
      }
    }

    auto const& sipmStatus = _sipmStatus.get(event.id());

    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    std::vector<std::shared_ptr<CRSScintillatorBar> >::const_iterator iter;
    for(iter=counters.begin(); iter!=counters.end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = (*iter)->index();

      for(size_t SiPM=0; SiPM<CRVId::nChanPerBar; SiPM++)
      {
        if(!(*iter)->getBarDetail().hasCMB(SiPM%CRVId::nSidesPerBar)) continue;  //no SiPM charges at non-existing SiPMs
                                                               //SiPM%2 returns the side of the CRV counter
                                                               //0 ... negative side
                                                               //1 ... positive side

        if(_useSipmStatusDB)
        {
          size_t channel = barIndex.asUint()*CRVId::nChanPerBar + SiPM;
          std::bitset<16> status(sipmStatus.status(channel));
          if(status.test(CRVStatus::Flags::notConnected) || status.test(CRVStatus::Flags::noData)) continue; //SiPM not connected (bit 0) or no data (bit 2)
        }

        //time wrapping happened in the photon generator
        std::vector<std::pair<double,size_t> > photonTimesNew;   //pair of photon time and index in the original photon vector
        CrvPhotonsCollection::const_iterator crvPhotons;
        for(crvPhotons=crvPhotonsCollection->begin(); crvPhotons!=crvPhotonsCollection->end(); crvPhotons++)
        {
          if(crvPhotons->GetScintillatorBarIndex()==barIndex && crvPhotons->GetSiPMNumber()==(int)SiPM)
          {
            const std::vector<CrvPhotons::SinglePhoton> &photonTimes = crvPhotons->GetPhotons();
            for(size_t iphoton=0; iphoton<photonTimes.size(); iphoton++)
            {
              //No check whether photons are within startTime and endTime
              double time = photonTimes[iphoton]._time;
              photonTimesNew.emplace_back(time,iphoton);
            }
            break;
          }
        }

        std::vector<mu2eCrv::SiPMresponse> SiPMresponseVector;
        _makeCrvSiPMCharges->Simulate(photonTimesNew, SiPMresponseVector, startTime, endTime);

        if(SiPMresponseVector.size()>0)
        {
          crvSiPMChargesCollection->emplace_back(barIndex,SiPM);
          std::vector<CrvSiPMCharges::SingleCharge> &charges = crvSiPMChargesCollection->back().GetCharges();

          std::vector<mu2eCrv::SiPMresponse>::const_iterator responseIter;
          for(responseIter=SiPMresponseVector.begin(); responseIter!=SiPMresponseVector.end(); responseIter++)
          {
            double time=responseIter->_time;
            double charge=responseIter->_charge;
            double chargeInPEs=responseIter->_chargeInPEs;
            int photonIndex=responseIter->_photonIndex;
            bool darkNoise=responseIter->_darkNoise;
            if(!darkNoise)
            {
              const std::vector<CrvPhotons::SinglePhoton> &photonTimes = crvPhotons->GetPhotons();
              charges.emplace_back(time, charge, chargeInPEs, photonTimes[photonIndex]._step);
            }
            else charges.emplace_back(time, charge, chargeInPEs);
          }
        }//non-empty SiPM charges
      }//SiPM
    }//barIndex

    event.put(std::move(crvSiPMChargesCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvSiPMChargeGenerator;
DEFINE_ART_MODULE(CrvSiPMChargeGenerator)
