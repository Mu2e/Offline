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
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "SeedService/inc/SeedService.hh"

#include "canvas/Persistency/Common/Ptr.h"
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
  class CrvPhotonGenerator : public art::EDProducer
  {

    public:
    explicit CrvPhotonGenerator(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginRun(art::Run& r);

    private:
    std::vector<std::string> _moduleLabels;
    std::vector<std::string> _processNames;
    std::vector<std::unique_ptr<art::Selector> > _selectors;

    ConfigFileLookupPolicy                                     _resolveFullPath;
    std::vector<std::string>                                   _lookupTableFileNames;
    std::vector<int>                                           _lookupTableReflectors;
    std::vector<std::string>                                   _lookupTableCRVSectors;
    std::vector<boost::shared_ptr<mu2eCrv::MakeCrvPhotons> >   _makeCrvPhotons;
    std::vector<double>                                        _scintillationYields;

    double      _scintillationYieldScaleFactor;
    double      _scintillationYieldVariation;
    double      _scintillationYieldVariationCutoffLow;
    double      _scintillationYieldVariationCutoffHigh;

    double      _startTime;             //StepPoint times before this time will be ignored to reduce computation times
                                        //(in particular by ignoring hits during the beam flash).
                                        //This time should be at least 100ns before the end of the SiPM's blind time
                                        //to account for the travel time of the photons inside the CRV bar.
                                        //Default is 0.

    SimParticleTimeOffset _timeOffsets;

    CLHEP::HepRandomEngine& _engine;
    CLHEP::RandFlat       _randFlat;
    CLHEP::RandGaussQ     _randGaussQ;
    CLHEP::RandPoissonQ   _randPoissonQ;

    std::map<CRSScintillatorBarIndex,double>  _scintillationYieldsAdjusted;
  };

  CrvPhotonGenerator::CrvPhotonGenerator(fhicl::ParameterSet const& pset) :
    EDProducer{pset},
    _moduleLabels(pset.get<std::vector<std::string> >("crvStepModuleLabels")),
    _processNames(pset.get<std::vector<std::string> >("crvStepProcessNames")),
    _lookupTableFileNames(pset.get<std::vector<std::string> >("lookupTableFileNames")),
    _lookupTableReflectors(pset.get<std::vector<int> >("reflectors")),
    _lookupTableCRVSectors(pset.get<std::vector<std::string> >("CRVSectors")),
    _scintillationYields(pset.get<std::vector<double> >("scintillationYields")),    //39400 photons per MeV
    _scintillationYieldScaleFactor(pset.get<double>("scintillationYieldScaleFactor")),    //100%
    _scintillationYieldVariation(pset.get<double>("scintillationYieldVariation")),    //20.0%
    _scintillationYieldVariationCutoffLow(pset.get<double>("scintillationYieldVariationCutoffLow")),    //60.0%
    _scintillationYieldVariationCutoffHigh(pset.get<double>("scintillationYieldVariationCutoffHigh")),    //120.0%
    _startTime(pset.get<double>("startTime")),               //0.0 ns
    _timeOffsets(pset.get<fhicl::ParameterSet>("timeOffsets", fhicl::ParameterSet())),
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

    if(_lookupTableFileNames.size()!=_lookupTableCRVSectors.size()) throw std::logic_error("ERROR: mismatch between specified CRV sector names and lookup table list");
    if(_lookupTableReflectors.size()!=_lookupTableCRVSectors.size()) throw std::logic_error("ERROR: mismatch between specified CRV sector names and reflector list");
    if(_scintillationYields.size()!=_lookupTableCRVSectors.size()) throw std::logic_error("ERROR: mismatch between specified CRV sector names and scintillation yield list");

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
           std::cout<<"CRV sector "<<i<<" ("<<_lookupTableCRVSectors[i]<<") uses "<<_makeCrvPhotons.back()->GetFileName()<<std::endl;
           break;
        }
      }
      if(tableLoaded) continue;

      _makeCrvPhotons.emplace_back(boost::shared_ptr<mu2eCrv::MakeCrvPhotons>(new mu2eCrv::MakeCrvPhotons(_randFlat, _randGaussQ, _randPoissonQ)));
      boost::shared_ptr<mu2eCrv::MakeCrvPhotons> &photonMaker=_makeCrvPhotons.back();
      photonMaker->LoadLookupTable(_resolveFullPath(_lookupTableFileNames[i]));
      photonMaker->SetScintillationYield(_scintillationYields[i]);
      std::cout<<"CRV sector "<<i<<" ("<<_lookupTableCRVSectors[i]<<") uses "<<_makeCrvPhotons.back()->GetFileName()<<" with scintillation yield of "<<_scintillationYields[i]<<" photons/MeV"<<std::endl;
    }

    produces<CrvPhotonsCollection>();
  }

  void CrvPhotonGenerator::beginRun(art::Run& rr)
  {
    GeomHandle<CosmicRayShield> CRS;
    std::vector<CRSScintillatorShield> const &shields = CRS->getCRSScintillatorShields();
    if(shields.size()!=_lookupTableCRVSectors.size()) throw std::logic_error("ERROR: mismatch between the geometry and the specified lookup table CRVSectors");

    for(size_t i=0; i<shields.size(); ++i)
    {
      if(shields[i].getCRSScintillatorBarDetail().getMaterialName()!="G4_POLYSTYRENE")
        throw std::logic_error("ERROR: scintillator material is not the expected G4_POLYSTYRENE which is used in the look-up tables");
      if(shields[i].getName().substr(4)!=_lookupTableCRVSectors[i]) throw std::logic_error("ERROR: mismatch between the geometry and the specified lookup table CRVSectors");
                            //substr(4) removes the "CRV_" part of the sector name
    }
  }

  void CrvPhotonGenerator::produce(art::Event& event)
  {
    _timeOffsets.updateMap(event);

    _scintillationYieldsAdjusted.clear();

    std::unique_ptr<CrvPhotonsCollection> crvPhotonsCollection(new CrvPhotonsCollection);

    std::map<std::pair<mu2e::CRSScintillatorBarIndex,int>,std::vector<CrvPhotons::SinglePhoton> > photonMap;

    GeomHandle<CosmicRayShield> CRS;
    GlobalConstantsHandle<ParticleDataTable> particleDataTable;

    for(size_t j=0; j<_selectors.size(); ++j)
    {
      std::vector<art::Handle<CrvStepCollection> > CrvStepsVector;
      event.getMany(*(_selectors.at(j)), CrvStepsVector);
      for(size_t i=0; i<CrvStepsVector.size(); ++i)
      {
        const art::Handle<CrvStepCollection> &CrvSteps = CrvStepsVector[i];
        for(size_t istep=0; istep<CrvSteps->size(); ++istep)
        {
          CrvStep const& step(CrvSteps->at(istep));

          double timeOffset = _timeOffsets.totalTimeOffset(step.simParticle());
          double t1 = step.startTime()+timeOffset;
          double t2 = step.endTime()+timeOffset;
          if(t1<_startTime) continue;   //Ignore this StepPoint to reduce computation time.

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
                                        _lookupTableReflectors[CRVSectorNumber]);

          art::Ptr<CrvStep> crvStepPtr(CrvSteps,istep);
          for(int SiPM=0; SiPM<4; ++SiPM)
          {
            std::pair<CRSScintillatorBarIndex,int> barIndexSiPMNumber(step.barIndex(),SiPM);
            const std::vector<double> &times=photonMaker->GetArrivalTimes(SiPM);
            if(times.empty()) continue;
            std::vector<CrvPhotons::SinglePhoton> &photons = photonMap[barIndexSiPMNumber];
            for(size_t itime=0; itime<times.size(); ++itime) photons.emplace_back(times[itime],crvStepPtr);
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
