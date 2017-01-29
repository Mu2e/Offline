//
// A module to collect background StepPoints and put their information into a TTree
// The StepPoint positions are stored relative to the counter center
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include <string>

#include <TTree.h>
#include <TVector2.h>


namespace mu2e 
{
  class CrvBackgroundCollector : public art::EDAnalyzer 
  {

    public:
    explicit CrvBackgroundCollector(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& e);

    private:
    std::vector<std::string> _g4ModuleLabels;
    std::vector<std::string> _processNames;

    size_t                   _volumeIdStart, _volumeIdEnd;
    double                   _startTime;

    TTree                   *_tree;
    ULong64_t                _eventId;
    //StepPointMC information
    ULong64_t                _volumeId;
    double                   _totalEDep;
    double                   _nonIonizingEDep;
    double                   _time;
    double                   _proper;
    double                   _positionX, _positionY, _positionZ;   //these are positions relative to the center of the counter
    double                   _momentumX, _momentumY, _momentumZ;
    double                   _stepLength;
    int                      _endProcessCode;
    //relevant SimParticle information
    ULong64_t                _id;
    int                      _pdgId;
  };

  CrvBackgroundCollector::CrvBackgroundCollector(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _g4ModuleLabels(pset.get<std::vector<std::string> >("g4ModuleLabels")),
    _processNames(pset.get<std::vector<std::string> >("processNames")),
    _volumeIdStart(pset.get<size_t>("volumeIdStart")),
    _volumeIdEnd(pset.get<size_t>("volumeIdEnd")),
    _startTime(pset.get<double>("startTime"))
  {
    if(_g4ModuleLabels.size()!=_processNames.size()) throw std::logic_error("ERROR: mismatch between specified selectors (g4ModuleLabels/processNames)");

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("background");
    _tree = tfdir.make<TTree>("background","background");
    _tree->Branch("eventId", &_eventId);
    _tree->Branch("volumeId", &_volumeId);
    _tree->Branch("totalEDep", &_totalEDep);
    _tree->Branch("nonIonizingEDep", &_nonIonizingEDep);
    _tree->Branch("time", &_time);
    _tree->Branch("proper", &_proper);
    _tree->Branch("positionX", &_positionX);
    _tree->Branch("positionY", &_positionY);
    _tree->Branch("positionZ", &_positionZ);
    _tree->Branch("momentumX", &_momentumX);
    _tree->Branch("momentumY", &_momentumY);
    _tree->Branch("momentumZ", &_momentumZ);
    _tree->Branch("stepLength", &_stepLength);
    _tree->Branch("endProcessCode", &_endProcessCode);
    _tree->Branch("id", &_id);
    _tree->Branch("pdgId", &_pdgId);

    _tree->GetUserInfo()->Add(new TVector2(_volumeIdStart,_volumeIdEnd));
  }

  void CrvBackgroundCollector::analyze(const art::Event& event) 
  {
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

          _time            = step.time();
          _volumeId        = step.barIndex().asInt();
          if(_time<_startTime) continue;
          if(_volumeId<_volumeIdStart || _volumeId>_volumeIdEnd) continue;

          _eventId         = event.event();
          _totalEDep       = step.totalEDep();
          _nonIonizingEDep = step.nonIonizingEDep();
          _proper          = step.properTime();
          
          const CLHEP::Hep3Vector &position = step.position();
          const CRSScintillatorBar &CRSbar = CRS->getBar(CRSScintillatorBarIndex(_volumeId));
          const CLHEP::Hep3Vector &positionLocal = CRSbar.toLocal(position);

          _positionX       = positionLocal.x();
          _positionY       = positionLocal.y();
          _positionZ       = positionLocal.z();
          _momentumX       = step.momentum().x();
          _momentumY       = step.momentum().y();
          _momentumZ       = step.momentum().z();
          _stepLength      = step.stepLength();
          _endProcessCode  = step.endProcessCode();
          _id    = step.simParticle()->id().asInt();
          _pdgId = step.simParticle()->pdgId();

          _tree->Fill();
        }
      }
    }
  }

} // end namespace mu2e

using mu2e::CrvBackgroundCollector;
DEFINE_ART_MODULE(CrvBackgroundCollector)
