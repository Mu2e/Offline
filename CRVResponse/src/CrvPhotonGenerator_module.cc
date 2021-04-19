//
// A module to create CRV photons arriving at the SiPMs (using StepPointMCs)
//
//
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvPhotons.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvStep.hh"
#include "MCDataProducts/inc/CrvPhotons.hh"
#include "MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "SeedService/inc/SeedService.hh"

#include "ProditionsService/inc/ProditionsHandle.hh"
#include "DAQConditions/inc/EventTiming.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/CrvParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "DataProducts/inc/EventWindowMarker.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Random/Randomize.h"

#include <string>

#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>


namespace mu2e
{
  class CrvPhotonGenerator : public art::EDProducer
  {

    public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config 
    {
      fhicl::Sequence<std::string> moduleLabels{ Name("crvStepModuleLabels"), Comment("CrvStepModule labels")};
      fhicl::Sequence<std::string> processNames{ Name("crvStepProcessNames"), Comment("process names of CrvSteps")};
      fhicl::Sequence<std::string> CRVSectors{ Name("CRVSectors"), Comment("Crv sectors")};
      fhicl::Sequence<int> reflectors{ Name("reflectors"), Comment("location of reflectors at Crv sectors")};
      fhicl::Sequence<std::string> lookupTableFileNames{ Name("lookupTableFileNames"), Comment("lookup tables for Crv sectors")};
      fhicl::Sequence<double> scintillationYields{ Name("scintillationYields"), Comment("scintillation yields at Crv sectors")};
      fhicl::Atom<double> scintillationYieldScaleFactor{ Name("scintillationYieldScaleFactor"), 
                                                        Comment("scale factor for scintillation yield")};
      fhicl::Atom<double> scintillationYieldVariation{ Name("scintillationYieldVariation"), 
                                                      Comment("sigma of gaussian variation of scintillation yield")};
      fhicl::Atom<double> scintillationYieldVariationCutoffLow{ Name("scintillationYieldVariationCutoffLow"), 
                                                               Comment("lower cutoff at scintillation yield variation")};
      fhicl::Atom<double> scintillationYieldVariationCutoffHigh{ Name("scintillationYieldVariationCutoffHigh"), 
                                                                Comment("upper cutoff at scintillation yield variation")};
      fhicl::Atom<double> digitizationStart{ Name("digitizationStart"), Comment("start of digitization")};
      fhicl::Atom<double> digitizationEnd{ Name("digitizationEnd"), Comment("end of digitization")};
      fhicl::Atom<double> crvStepMargin{ Name("crvStepMargin"), Comment("time window allowed for crvSteps before digitization starts")};
      fhicl::Atom<art::InputTag> eventWindowMarkerTag{ Name("eventWindowMarkerTag"), Comment("EventWindowMarker producer"),"EWMProducer" };
      fhicl::Atom<art::InputTag> protonBunchTimeMCTag{ Name("protonBunchTimeMCTag"), Comment("ProtonBunchTimeMC producer"),"EWMProducer" };
      fhicl::Sequence<art::InputTag> timeOffsets { Name("timeOffsets"), Comment("Sim Particle Time Offset Maps")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit CrvPhotonGenerator(const Parameters& conf);
    void produce(art::Event& e);
    void beginRun(art::Run& r);

    private:
    std::vector<std::string> _moduleLabels;
    std::vector<std::string> _processNames;
    std::vector<std::unique_ptr<art::Selector> > _selectors;

    ConfigFileLookupPolicy                                     _resolveFullPath;
    std::vector<std::string>                                   _CRVSectors;
    std::vector<int>                                           _reflectors;
    std::vector<std::string>                                   _lookupTableFileNames;
    std::vector<double>                                        _scintillationYields;
    std::vector<boost::shared_ptr<mu2eCrv::MakeCrvPhotons> >   _makeCrvPhotons;

    double      _scintillationYieldScaleFactor;
    double      _scintillationYieldVariation;
    double      _scintillationYieldVariationCutoffLow;
    double      _scintillationYieldVariationCutoffHigh;

    //On-spill:
    //Photons will be recorded during the digitization window (400ns...1750ns). The CRV digitization window needs to start
    //about 100ns before the tracker digitization (which start at 500ns) to catch cosmic ray muons which may cause signals
    //in the tracker.
    //CrvSteps need to be recorded starting about 50ns earlier to acount for the travel time of the photons inside the CRV
    //counters. This gives a start time of 350ns. It will avoid recording hits occurring around the time of the beam flash.
    //All photon times need to be folded modulus 1695ns (=microbunch time). Since the the digitization window extends past
    //the microbunch time (for 1750ns-1695ns=55ns), photons that are recorded between 0ns and 55ns will be copied to times
    //after the microbunch ("ghost hits"). Therefore, CrvSteps also need be recorded before 55ns.
    //Summary:
    //-record CrvSteps before t_digitization_end-length_microbunch (i.e. 1750ns-1695ns=55ns)
    //-record CrvSteps after t_digitization_start-length_crvstepmargin (i.e. 400ns-50ns=350ns)
    //
    //Off-spill:
    //The digitization window is equal to the event time window in off-spill. Photons are only recorded during the event 
    //time window, but the CrvSteps need to be recorded starting about 50ns earlier to account for the travel time of the 
    //photons inside the CRV counters. There is no time folding in off-spill.
    //-record CrvSteps after eventwindow_start-length_crvstepmargin until eventwindow_end

    double      _digitizationStart; //400ns
    double      _digitizationEnd;   //1750ns
    double      _crvStepMargin;     //50ns
    art::InputTag _eventWindowMarkerTag;
    art::InputTag _protonBunchTimeMCTag;
    double      _microBunchPeriod;

    SimParticleTimeOffset _timeOffsets;

    CLHEP::HepRandomEngine& _engine;
    CLHEP::RandFlat       _randFlat;
    CLHEP::RandGaussQ     _randGaussQ;
    CLHEP::RandPoissonQ   _randPoissonQ;

    std::map<CRSScintillatorBarIndex,double>  _scintillationYieldsAdjusted;
  };

  CrvPhotonGenerator::CrvPhotonGenerator(const Parameters& conf) :
    art::EDProducer{conf},
    _moduleLabels(conf().moduleLabels()),
    _processNames(conf().processNames()),
    _CRVSectors(conf().CRVSectors()),
    _reflectors(conf().reflectors()),
    _lookupTableFileNames(conf().lookupTableFileNames()),
    _scintillationYields(conf().scintillationYields()),
    _scintillationYieldScaleFactor(conf().scintillationYieldScaleFactor()),
    _scintillationYieldVariation(conf().scintillationYieldVariation()),
    _scintillationYieldVariationCutoffLow(conf().scintillationYieldVariationCutoffLow()),
    _scintillationYieldVariationCutoffHigh(conf().scintillationYieldVariationCutoffHigh()),
    _digitizationStart(conf().digitizationStart()),
    _digitizationEnd(conf().digitizationEnd()),
    _crvStepMargin(conf().crvStepMargin()),
    _eventWindowMarkerTag(conf().eventWindowMarkerTag()),
    _protonBunchTimeMCTag(conf().protonBunchTimeMCTag()),
    _timeOffsets(conf().timeOffsets()),
    _engine{createEngine(art::ServiceHandle<SeedService>()->getSeed())},
    _randFlat(_engine),
    _randGaussQ(_engine),
    _randPoissonQ(_engine)
  {
    if(_moduleLabels.size()==0) throw std::logic_error("ERROR: a list of crvSteps module labels needs to be provided");
    if(_moduleLabels.size()!=_processNames.size()) throw std::logic_error("ERROR: mismatch between specified selectors (crvStepModuleLabels/crvStepProcessNames)");
    for(size_t i=0; i<_moduleLabels.size(); ++i)
    {
      if(_moduleLabels[i]!="*")
        _selectors.push_back(std::unique_ptr<art::Selector>(new art::Selector(art::ProductInstanceNameSelector("") &&
                                                            art::ModuleLabelSelector(_moduleLabels[i]) &&
                                                            art::ProcessNameSelector(_processNames[i]))));
      else
        _selectors.push_back(std::unique_ptr<art::Selector>(new art::Selector(art::ProductInstanceNameSelector("") &&
                                                            art::ProcessNameSelector(_processNames[i]))));
    }

    if(_lookupTableFileNames.size()!=_CRVSectors.size()) throw std::logic_error("ERROR: mismatch between specified CRV sector names and lookup table list");
    if(_reflectors.size()!=_CRVSectors.size()) throw std::logic_error("ERROR: mismatch between specified CRV sector names and reflector list");
    if(_scintillationYields.size()!=_CRVSectors.size()) throw std::logic_error("ERROR: mismatch between specified CRV sector names and scintillation yield list");

    for(size_t i=0; i<_lookupTableFileNames.size(); ++i)
    {
      _scintillationYields[i]*=_scintillationYieldScaleFactor;

      bool tableLoaded=false;
      for(size_t j=0; j<i; ++j)
      {
        if(_lookupTableFileNames[i]==_lookupTableFileNames[j])
        {
           tableLoaded=true;
           _makeCrvPhotons.emplace_back(_makeCrvPhotons[j]);
           std::cout<<"CRV sector "<<i<<" ("<<_CRVSectors[i]<<") uses "<<_makeCrvPhotons.back()->GetFileName()<<std::endl;
           break;
        }
      }
      if(tableLoaded) continue;

      _makeCrvPhotons.emplace_back(boost::shared_ptr<mu2eCrv::MakeCrvPhotons>(new mu2eCrv::MakeCrvPhotons(_randFlat, _randGaussQ, _randPoissonQ)));
      boost::shared_ptr<mu2eCrv::MakeCrvPhotons> &photonMaker=_makeCrvPhotons.back();
      photonMaker->LoadLookupTable(_resolveFullPath(_lookupTableFileNames[i]));
      photonMaker->SetScintillationYield(_scintillationYields[i]);
      std::cout<<"CRV sector "<<i<<" ("<<_CRVSectors[i]<<") uses "<<_makeCrvPhotons.back()->GetFileName()<<" with scintillation yield of "<<_scintillationYields[i]<<" photons/MeV"<<std::endl;
    }

    produces<CrvPhotonsCollection>();
  }

  void CrvPhotonGenerator::beginRun(art::Run& rr)
  {
    GeomHandle<CosmicRayShield> CRS;
    std::vector<CRSScintillatorShield> const &shields = CRS->getCRSScintillatorShields();
    if(shields.size()!=_CRVSectors.size()) throw std::logic_error("ERROR: mismatch between the geometry and the specified lookup table CRVSectors");

    for(size_t i=0; i<shields.size(); ++i)
    {
      if(shields[i].getCRSScintillatorBarDetail().getMaterialName()!="G4_POLYSTYRENE")
        throw std::logic_error("ERROR: scintillator material is not the expected G4_POLYSTYRENE which is used in the look-up tables");
      if(shields[i].getName().substr(4)!=_CRVSectors[i]) throw std::logic_error("ERROR: mismatch between the geometry and the specified lookup table CRVSectors");
                            //substr(4) removes the "CRV_" part of the sector name
    }

    mu2e::ConditionsHandle<mu2e::AcceleratorParams> accPar("ignored");
    _microBunchPeriod = accPar->deBuncherPeriod;
  }

  void CrvPhotonGenerator::produce(art::Event& event)
  {
    _timeOffsets.updateMap(event);

    _scintillationYieldsAdjusted.clear();

    std::unique_ptr<CrvPhotonsCollection> crvPhotonsCollection(new CrvPhotonsCollection);

    std::map<std::pair<mu2e::CRSScintillatorBarIndex,int>,std::vector<CrvPhotons::SinglePhoton> > photonMap;

    GeomHandle<CosmicRayShield> CRS;
    GlobalConstantsHandle<ParticleDataTable> particleDataTable;

    art::Handle<EventWindowMarker> eventWindowMarker;
    event.getByLabel(_eventWindowMarkerTag,eventWindowMarker);
    EventWindowMarker::SpillType spillType = eventWindowMarker->spillType();
    double eventWindowLength = eventWindowMarker->eventLength();

    art::Handle<ProtonBunchTimeMC> protonBunchTimeMC;
    event.getByLabel(_protonBunchTimeMCTag, protonBunchTimeMC);
    double eventWindowStart = -protonBunchTimeMC->pbtime_;
    double eventWindowEnd = eventWindowStart + eventWindowLength;

    ProditionsHandle<EventTiming> eventTimingHandle;
    const EventTiming &eventTiming = eventTimingHandle.get(event.id());
    double jitter = eventWindowStart - eventTiming.timeFromProtonsToDRMarker();

    double digitizationStart=_digitizationStart+jitter;
    double digitizationEnd=_digitizationEnd+jitter;

    for(size_t j=0; j<_selectors.size(); ++j)
    {
      std::vector<art::Handle<CrvStepCollection> > CrvStepsVector = event.getMany<CrvStepCollection>(*(_selectors.at(j)));
      for(size_t i=0; i<CrvStepsVector.size(); ++i)
      {
        const art::Handle<CrvStepCollection> &CrvSteps = CrvStepsVector[i];
        for(size_t istep=0; istep<CrvSteps->size(); ++istep)
        {
          CrvStep const& step(CrvSteps->at(istep));

          double timeOffset = _timeOffsets.totalTimeOffset(step.simParticle());
          double t1 = step.startTime()+timeOffset;
          double t2 = step.endTime()+timeOffset;
          if(isnan(t1) || isnan(t2)) continue;  //This situation was observed once. Not sure how it happened.

          //see explanation above
          //On-spill: No photons in the blind time between before digitizationStart and up to 55ns after 0
          //(which will moved past the microbunch end after time folding) but allow a little bit more for the crvSteps
          //-record CrvSteps before digitizationEnd-microBunchPeriod (i.e. 1750ns-1695ns=55ns)
          //-record CrvSteps after digitizationStart-crvStepMargin (i.e. 400ns-50ns=350ns)
          //Off-spill: Only photons within eventWindow, but allow a little bit more for the crvSteps
          //-record CrvSteps after eventwindowStart-crvStepMargin until eventWindowEnd
          if(spillType==EventWindowMarker::SpillType::onspill)
          {
            if(t1>digitizationEnd-_microBunchPeriod && t2<digitizationStart-_crvStepMargin) continue;
          }
          else
          {
            if(t2<eventWindowStart-_crvStepMargin || t1>eventWindowEnd) continue;
          }

          CLHEP::Hep3Vector pos1 = step.startPosition();  //TODO: Need to convert everything into XYZVec, so that const references can be used
          CLHEP::Hep3Vector pos2 = step.endPosition();

          int PDGcode = step.simParticle()->pdgId();
          ParticleDataTable::maybe_ref particle = particleDataTable->particle(PDGcode);
          if(!particle)
          {
            std::cerr<<"Error in CrvPhotonGenerator: Found a PDG code which is not in the GEANT particle table: ";
            std::cerr<<PDGcode<<std::endl;
            continue;
          }
          double mass = particle.ref().mass();  //MeV/c^2
          double charge = particle.ref().charge(); //in units of elementary charges

          double energy1   = sqrt(step.startMom().mag2() + mass*mass); //MeV
          double energy2   = sqrt(step.endMom()*step.endMom() + mass*mass);
          double avgEnergy = 0.5*(energy1+energy2);
          double avgGamma  = avgEnergy/mass;
          double avgBeta   = sqrt(1.0-1.0/(avgGamma*avgGamma));

          const CRSScintillatorBar &CRSbar = CRS->getBar(step.barIndex());
          const CLHEP::Hep3Vector &pos1Local = CRSbar.toLocal(pos1);
          const CLHEP::Hep3Vector &pos2Local = CRSbar.toLocal(pos2);

          const CRSScintillatorBarId &barId = CRSbar.id();
          int CRVSectorNumber=barId.getShieldNumber();
          if(_scintillationYieldsAdjusted.find(step.barIndex())==_scintillationYieldsAdjusted.end())
          {
            double sectorScintillationYield=_scintillationYields[CRVSectorNumber];
            double adjustedYield=0;
            do
            {
              adjustedYield=_randGaussQ.fire(sectorScintillationYield, sectorScintillationYield*_scintillationYieldVariation);
            } while(adjustedYield<sectorScintillationYield*_scintillationYieldVariationCutoffLow ||
                    adjustedYield>sectorScintillationYield*_scintillationYieldVariationCutoffHigh);

            _scintillationYieldsAdjusted[step.barIndex()] = adjustedYield;
          }
          double currentAdjustedYield = _scintillationYieldsAdjusted[step.barIndex()];

          boost::shared_ptr<mu2eCrv::MakeCrvPhotons> &photonMaker=_makeCrvPhotons.at(CRVSectorNumber);
          photonMaker->SetScintillationYield(currentAdjustedYield);
          photonMaker->MakePhotons(pos1Local, pos2Local, t1, t2,
                                        avgBeta, charge,
                                        step.visibleEDep(),
                                        step.pathLength(),
                                        _reflectors[CRVSectorNumber]);

          art::Ptr<CrvStep> crvStepPtr(CrvSteps,istep);
          for(int SiPM=0; SiPM<4; ++SiPM)
          {
            std::pair<CRSScintillatorBarIndex,int> barIndexSiPMNumber(step.barIndex(),SiPM);
            const std::vector<double> &times=photonMaker->GetArrivalTimes(SiPM);
            if(times.empty()) continue;
            std::vector<CrvPhotons::SinglePhoton> &photons = photonMap[barIndexSiPMNumber];
            for(size_t itime=0; itime<times.size(); ++itime)
            {
              double timeTmp=times[itime];
              if(spillType==EventWindowMarker::SpillType::onspill)
              {
                timeTmp = fmod(timeTmp,_microBunchPeriod);
                //photons before the digitization start get removed except photons 
                //in the first 55ns which get moved to the interval between the end of 
                //the microbunch period and the digitization end
                if(timeTmp<digitizationEnd-_microBunchPeriod) timeTmp+=_microBunchPeriod;  
                if(timeTmp<digitizationStart) continue;
              }
              else
              {              
                //photons outside the eventWindow get removed
                if(timeTmp<eventWindowStart || timeTmp>eventWindowEnd) continue;
              }
              photons.emplace_back(timeTmp,crvStepPtr);
            }
          }

        } //loop over StepPointMCs in the StepPointMC collection
      } //loop over all StepPointMC collections
    } //loop over all module labels / process names from the fcl file

    for(auto p=photonMap.begin(); p!=photonMap.end(); ++p)
    {
      crvPhotonsCollection->emplace_back(p->first.first,p->first.second,p->second);
    }

    event.put(std::move(crvPhotonsCollection));
  }

}

using mu2e::CrvPhotonGenerator;
DEFINE_ART_MODULE(CrvPhotonGenerator)
