//
// A module to create CRV photons arriving at the SiPMs (using StepPointMCs)
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "CRVResponse/inc/MakeCrvPhotonArrivals.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CrvPhotonArrivalsCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "SeedService/inc/SeedService.hh"

#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Random/Randomize.h"

#include <string>

#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>


namespace mu2e 
{
  class CrvPhotonArrivalsGenerator : public art::EDProducer 
  {

    public:
    explicit CrvPhotonArrivalsGenerator(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void endJob();
    void beginRun(art::Run& r);

    private:
    std::vector<std::string> _g4ModuleLabels;
    std::vector<std::string> _processNames;

    std::vector<std::string>                                             _lookupTableFileNames;
    std::vector<double>                                                  _lookupTableCounterLengths;
    std::map<double, boost::shared_ptr<mu2eCrv::MakeCrvPhotonArrivals> > _makeCrvPhotonArrivals;

    double      _scintillationYield;
    double      _scintillationYieldTolerance;
    double      _scintillatorBirksConstant;
    double      _scintillatorRatioFastSlow;
    double      _scintillatorDecayTimeFast;
    double      _scintillatorDecayTimeSlow;
    double      _fiberDecayTime;

    double      _startTime;             //StepPoint times before this time will be ignored to reduce computation times
                                        //(in particular by ignoring hits during the beam flash).
                                        //This time should be at least 100ns before the end of the SiPM's blind time
                                        //to account for the travel time of the photons inside the CRV bar.
                                        //Default is 0.

    SimParticleTimeOffset _timeOffsets;

    std::string                           _backgroundSampleFileName;
    int                                   _countersInBackgroundSample;
    int                                   _backgroundSampleFactor;
    double                                _maxBackgroundTimeShift;

    struct Background
    {
      CLHEP::Hep3Vector _p1Local;
      CLHEP::Hep3Vector _p2Local;
      double            _t1, _t2;
      double            _beta, _charge, _totalE, _nonIonizingE;
      int               _PDGcode;
      Background(double t1, double x1, double y1, double z1,
                 double t2, double x2, double y2, double z2,
                 double beta, double charge, double totalE, double nonIonizingE, int PDGcode) :
                 _p1Local(x1,y1,z1), _p2Local(x2,y2,z2), _t1(t1), _t2(t2), 
                 _beta(beta), _charge(charge), _totalE(totalE), _nonIonizingE(nonIonizingE), _PDGcode(PDGcode) {}
    };

    std::vector<std::vector<std::vector<Background> > > _backgroundEvents;
    std::vector<std::vector<std::vector<Background> > >::const_iterator _currentBackgroundEvent;

    CLHEP::RandFlat       _randFlat;
    CLHEP::RandGaussQ     _randGaussQ;

    std::map<CRSScintillatorBarIndex,double>  _scintillationYieldAdjustments;
  };

  CrvPhotonArrivalsGenerator::CrvPhotonArrivalsGenerator(fhicl::ParameterSet const& pset) :
    _g4ModuleLabels(pset.get<std::vector<std::string> >("g4ModuleLabels")),
    _processNames(pset.get<std::vector<std::string> >("processNames")),
    _lookupTableFileNames(pset.get<std::vector<std::string> >("lookupTableFileNames")),
    _lookupTableCounterLengths(pset.get<std::vector<double> >("lookupTableCounterLengths")),
    _scintillationYield(pset.get<double>("scintillationYield")),    //5000.0 photons per MeV
    _scintillationYieldTolerance(pset.get<double>("scintillationYieldTolerance")),    //0.0%
    _scintillatorBirksConstant(pset.get<double>("scintillatorBirksConstant")), //0.126 mm/MeV
    _scintillatorRatioFastSlow(pset.get<double>("scintillatorRatioFastSlow")), //1.0
    _scintillatorDecayTimeFast(pset.get<double>("scintillatorDecayTimeFast")), //10.0 ns, includes WLS components in the scintillator
    _scintillatorDecayTimeSlow(pset.get<double>("scintillatorDecayTimeSlow")), //100.0 ns, unknown, not used
    _fiberDecayTime(pset.get<double>("fiberDecayTime")),     //7.4 ns
    _startTime(pset.get<double>("startTime")),               //0.0 ns
    _timeOffsets(pset.get<fhicl::ParameterSet>("timeOffsets", fhicl::ParameterSet())),
    _backgroundSampleFileName(pset.get<std::string>("backgroundSampleFileName","")),
    _countersInBackgroundSample(pset.get<int>("countersInBackgroundSample",0)),
    _backgroundSampleFactor(pset.get<int>("backgroundSampleFactor",0)),
    _maxBackgroundTimeShift(pset.get<double>("maxBackgroundTimeShift",0)),
    _randFlat(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    _randGaussQ(art::ServiceHandle<art::RandomNumberGenerator>()->getEngine())
  {
    if(_g4ModuleLabels.size()!=_processNames.size()) throw std::logic_error("ERROR: mismatch between specified selectors (g4ModuleLabels/processNames)");

    if(_lookupTableFileNames.size()!=_lookupTableCounterLengths.size()) throw std::logic_error("ERROR: mismatch between specified lookup tables (lookupTableFileNames/lookupTableCounterLengths)");
    for(unsigned int i=0; i<_lookupTableFileNames.size(); i++)
    {
      double counterLength = _lookupTableCounterLengths[i];
      _makeCrvPhotonArrivals.emplace(counterLength, boost::shared_ptr<mu2eCrv::MakeCrvPhotonArrivals>(new mu2eCrv::MakeCrvPhotonArrivals(_randFlat)));
      std::map<double, boost::shared_ptr<mu2eCrv::MakeCrvPhotonArrivals> >::iterator iterCPA=_makeCrvPhotonArrivals.find(counterLength);
      iterCPA->second->LoadLookupTable(_lookupTableFileNames[i].c_str());
      iterCPA->second->SetScintillationYield(_scintillationYield);
      iterCPA->second->SetScintillatorBirksConstant(_scintillatorBirksConstant);
      iterCPA->second->SetScintillatorRatioFastSlow(_scintillatorRatioFastSlow);
      iterCPA->second->SetScintillatorDecayTimeFast(_scintillatorDecayTimeFast);
      iterCPA->second->SetScintillatorDecayTimeSlow(_scintillatorDecayTimeSlow);
      iterCPA->second->SetFiberDecayTime(_fiberDecayTime);
    }

    if(_backgroundSampleFileName!="" && _countersInBackgroundSample!=0 && _backgroundSampleFactor!=0)
    {
      std::cout<<"CRVResponse uses a background overlay ("<<_backgroundSampleFileName<<")"<<std::endl;

      TDirectory *directory = gDirectory;

      TFile *backgroundFile = TFile::Open(_backgroundSampleFileName.c_str());
      gDirectory->cd("background/background");
      TTree *tree = dynamic_cast<TTree*>(gDirectory->FindObjectAny("background"));

      ULong64_t     eventNumber, volumeId;
      double        t1, x1, y1, z1;
      double        t2, x2, y2, z2;
      double        beta, charge, totalE, nonIonizingE;
      int           PDGcode;
      tree->SetBranchAddress("eventNumber",&eventNumber);
      tree->SetBranchAddress("volumeId",&volumeId);
      tree->SetBranchAddress("t1", &t1);
      tree->SetBranchAddress("x1", &x1);
      tree->SetBranchAddress("y1", &y1);
      tree->SetBranchAddress("z1", &z1);
      tree->SetBranchAddress("t2", &t2);
      tree->SetBranchAddress("x2", &x2);
      tree->SetBranchAddress("y2", &y2);
      tree->SetBranchAddress("z2", &z2);
      tree->SetBranchAddress("beta", &beta);
      tree->SetBranchAddress("charge", &charge);
      tree->SetBranchAddress("totalE", &totalE);
      tree->SetBranchAddress("nonIonizingE", &nonIonizingE);
      tree->SetBranchAddress("PDGcode", &PDGcode);

      int n = tree->GetEntries();

      std::map<unsigned long, std::map<unsigned long, std::vector<Background> > > backgroundMap;
      for(int i=0; i<n; i++)
      {
        tree->GetEntry(i);
        std::map<unsigned long, std::vector<Background> > &backgroundEvent = backgroundMap[eventNumber];
        std::vector<Background> &backgroundVector = backgroundEvent[volumeId];
        backgroundVector.emplace_back(t1, x1, y1, z1, t2, x2, y2, z2, beta, charge, totalE, nonIonizingE, PDGcode);
      }
      backgroundFile->Close();

      //convert maps into vectors to decrease the access times (the actual event number and volume Id is not needed)
      std::map<unsigned long, std::map<unsigned long, std::vector<Background> > >::const_iterator backgroundIter;
      for(backgroundIter=backgroundMap.begin(); backgroundIter!=backgroundMap.end(); backgroundIter++)
      {
        _backgroundEvents.resize(_backgroundEvents.size()+1);

        std::map<unsigned long, std::vector<Background> >::const_iterator backgroundIter2;
        for(backgroundIter2=backgroundIter->second.begin(); backgroundIter2!=backgroundIter->second.end(); backgroundIter2++)
        {
          const std::vector<Background> &backgroundVector = backgroundIter2->second;
          _backgroundEvents.back().push_back(backgroundVector);
        }
      }

      _currentBackgroundEvent = _backgroundEvents.begin();

      directory->cd();
    }

    produces<CrvPhotonArrivalsCollection>();
  }

  void CrvPhotonArrivalsGenerator::beginJob()
  {
  }

  void CrvPhotonArrivalsGenerator::endJob()
  {
  }

  void CrvPhotonArrivalsGenerator::beginRun(art::Run& rr)
  {
    GeomHandle<CosmicRayShield> CRS;
    std::vector<CRSScintillatorShield> const &shields = CRS->getCRSScintillatorShields();
    std::vector<CRSScintillatorShield>::const_iterator ishield;
    for(ishield=shields.begin(); ishield!=shields.end(); ++ishield) 
    {
      if(ishield->getCRSScintillatorBarDetail().getMaterialName()!="G4_POLYSTYRENE")
      {
        throw std::logic_error("scintillator material is not the expected G4_POLYSTYRENE which is used in the look-up tables");
      }
    }
  }

  void CrvPhotonArrivalsGenerator::produce(art::Event& event) 
  {
    _timeOffsets.updateMap(event);
    _scintillationYieldAdjustments.clear();

    std::unique_ptr<CrvPhotonArrivalsCollection> crvPhotonArrivalsCollection(new CrvPhotonArrivalsCollection);

    GeomHandle<CosmicRayShield> CRS;

    std::vector<art::Handle<StepPointMCCollection> > CRVStepsVector;
    std::unique_ptr<art::Selector> selector;
    for(size_t j=0; j<_g4ModuleLabels.size(); j++)
    {
      if(_g4ModuleLabels[j]!="" && _g4ModuleLabels[j]!="*")
        selector = std::unique_ptr<art::Selector>(new art::Selector(art::ProductInstanceNameSelector("CRV") &&
                                                                    art::ModuleLabelSelector(_g4ModuleLabels[j]) && 
                                                                    art::ProcessNameSelector(_processNames[j])));
      else
        selector = std::unique_ptr<art::Selector>(new art::Selector(art::ProductInstanceNameSelector("CRV") &&
                                                                    art::ProcessNameSelector(_processNames[j])));
      //the ProcessNameSelector allows "*" and ""

      event.getMany(*selector, CRVStepsVector);
      for(size_t i=0; i<CRVStepsVector.size(); i++)
      {
        const art::Handle<StepPointMCCollection> &CRVSteps = CRVStepsVector[i];

        for(StepPointMCCollection::const_iterator iter=CRVSteps->begin(); iter!=CRVSteps->end(); iter++)
        {
          StepPointMC const& step(*iter);

          double t1 = _timeOffsets.timeWithOffsetsApplied(step);
          if(t1<_startTime) continue;   //Ignore this StepPoint to reduce computation time.

          const CLHEP::Hep3Vector &p1 = step.position();
          CLHEP::Hep3Vector p2 = p1 + step.momentum().unit()*step.stepLength();
          double energyDepositedTotal= step.totalEDep();
          double energyDepositedNonIonizing = step.nonIonizingEDep();

          GlobalConstantsHandle<ParticleDataTable> particleDataTable;
          int PDGcode = step.simParticle()->pdgId();
          ParticleDataTable::maybe_ref particle = particleDataTable->particle(PDGcode);
          if(!particle) 
          {
            std::cerr<<"Error in CrvPhotonArrivalsGenerator: Found a PDG code which is not in the GEANT particle table: ";
            std::cerr<<PDGcode<<std::endl;
            continue;
          }
          double mass = particle.ref().mass();  //MeV/c^2
          double charge = particle.ref().charge(); //in units of elementary charges 

          double momentum1 = step.momentum().mag(); //MeV/c
          double energy1 = sqrt(momentum1*momentum1 + mass*mass); //MeV
//FIXME: does not take the energy of daughter particles into account
          double energy2 = energy1 - energyDepositedTotal; //MeV  
          if(energy2<mass) energy2=mass;

          double gamma1 = energy1 / mass;
          double gamma2 = energy2 / mass;
          double beta1 = sqrt(1.0-1.0/(gamma1*gamma1));
          double beta2 = sqrt(1.0-1.0/(gamma2*gamma2));
          double beta = (beta1+beta2)/2.0;
          double velocity = beta*CLHEP::c_light;
          double t2 = t1 + step.stepLength()/velocity;

          const CRSScintillatorBar &CRSbar = CRS->getBar(step.barIndex());
          const CLHEP::Hep3Vector &p1Local = CRSbar.toLocal(p1);
          const CLHEP::Hep3Vector &p2Local = CRSbar.toLocal(p2);

          if(_scintillationYieldAdjustments.find(step.barIndex())==_scintillationYieldAdjustments.end())
          {
            double adjustment = _randGaussQ.fire(0, _scintillationYield*_scintillationYieldTolerance);
            _scintillationYieldAdjustments[step.barIndex()] = adjustment;
          }
          double scintillationYieldAdjustment = _scintillationYieldAdjustments[step.barIndex()];

          double counterLength = CRSbar.getHalfLength()*2.0;
          std::map<double, boost::shared_ptr<mu2eCrv::MakeCrvPhotonArrivals> >::iterator iterCPA=_makeCrvPhotonArrivals.find(counterLength);
          if(iterCPA!=_makeCrvPhotonArrivals.end())
          {
            iterCPA->second->MakePhotons(p1Local, p2Local, t1, t2,  
                                        PDGcode, beta, charge,
                                        energyDepositedTotal,
                                        energyDepositedNonIonizing,
                                        scintillationYieldAdjustment);

            CrvPhotonArrivals &crvPhotons = (*crvPhotonArrivalsCollection)[step.barIndex()];
            for(int SiPM=0; SiPM<4; SiPM++)
            {
              const std::vector<double> &times=iterCPA->second->GetArrivalTimes(SiPM);
              crvPhotons.GetPhotonArrivalTimes(SiPM).insert(crvPhotons.GetPhotonArrivalTimes(SiPM).end(),times.begin(),times.end());
            }
          }
          else std::cout<<"ERROR: Did not find the matching lookuptable file for the CRV counter with the length "<<counterLength<<"."<<std::endl;

        } //loop over StepPointMCs in the StepPointMC collection
      } //loop over all StepPointMC collections
    } //loop over all module labels / process names from the fcl file

/* overlay background */
/* should only be used if background samples and actual "detector" use only one counter length */

   
    if(_backgroundSampleFileName!="" && _countersInBackgroundSample!=0 && _backgroundSampleFactor!=0)
    {
      const std::vector<std::vector<Background> > &backgroundSampleCounters = *_currentBackgroundEvent;
      const std::vector<std::shared_ptr<CRSScintillatorBar> > &actualCounters = CRS->getAllCRSScintillatorBars();
      for(int i=0; i<_backgroundSampleFactor; i++)
      {
        for(size_t j=0; j<actualCounters.size(); j++)
        {
          CRSScintillatorBarIndex jj(j);
          if(_scintillationYieldAdjustments.find(jj)==_scintillationYieldAdjustments.end())
          {
            double adjustment = _randGaussQ.fire(0, _scintillationYield*_scintillationYieldTolerance);
            _scintillationYieldAdjustments[jj] = adjustment;
          }
          double scintillationYieldAdjustment = _scintillationYieldAdjustments[jj];

          double counterLength = actualCounters[j]->getHalfLength()*2.0;
          std::map<double, boost::shared_ptr<mu2eCrv::MakeCrvPhotonArrivals> >::iterator iterCPA=_makeCrvPhotonArrivals.find(counterLength);

          if(iterCPA!=_makeCrvPhotonArrivals.end())
          {
            unsigned int backgroundCounter = _randFlat.fire(_countersInBackgroundSample);
            if(backgroundCounter<backgroundSampleCounters.size())
            {
              const std::vector<Background> &backgroundVector = backgroundSampleCounters[backgroundCounter];
              double timeShift = _randFlat.fire(_maxBackgroundTimeShift);
              for(size_t k=0; k<backgroundVector.size(); k++)
              {
                const Background &b=backgroundVector[k];
                iterCPA->second->MakePhotons(b._p1Local, b._p2Local, b._t1-timeShift, b._t2-timeShift,  
                                             b._PDGcode, b._beta, b._charge,
                                             b._totalE, b._nonIonizingE,
                                             scintillationYieldAdjustment);

                CrvPhotonArrivals &crvPhotons = (*crvPhotonArrivalsCollection)[jj];
                for(int SiPM=0; SiPM<4; SiPM++)
                {
                  const std::vector<double> &times=iterCPA->second->GetArrivalTimes(SiPM);
                  crvPhotons.GetPhotonArrivalTimes(SiPM).insert(crvPhotons.GetPhotonArrivalTimes(SiPM).end(),times.begin(),times.end());
                }
              }
            }
          }
        }
      }

      _currentBackgroundEvent++;
      if(_currentBackgroundEvent==_backgroundEvents.end()) _currentBackgroundEvent = _backgroundEvents.begin();
    }

/* photnns into the event */

    event.put(std::move(crvPhotonArrivalsCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvPhotonArrivalsGenerator;
DEFINE_ART_MODULE(CrvPhotonArrivalsGenerator)
