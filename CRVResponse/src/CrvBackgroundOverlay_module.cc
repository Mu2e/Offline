//
// A module to read the background overlay StepPointMC information from a TTree, and recreates a new StepPointMC collection
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
//
// Original Author: Ralf Ehrlich

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include <string>

#include <TDirectory.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector2.h>

namespace mu2e
{
  class CrvBackgroundOverlay : public art::EDProducer
  {

    public:
    explicit CrvBackgroundOverlay(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);

    private:
    TFile                   *_file;
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

    std::string              _backgroundFile;
    int                      _overlayFactor;
    size_t                   _volumeIdStart, _volumeIdEnd;
    size_t                   _currentEntry;
    size_t                   _nEntries;
  };

  CrvBackgroundOverlay::CrvBackgroundOverlay(fhicl::ParameterSet const& pset) :
    _backgroundFile(pset.get<std::string>("backgroundFile")),
    _overlayFactor(pset.get<int>("overlayFactor"))
  {
    ConfigFileLookupPolicy configFile;
    _backgroundFile = configFile(_backgroundFile);

    TDirectory *directory = gDirectory;
    _file = TFile::Open(_backgroundFile.c_str());
    gDirectory->cd("background/background");
    _tree = dynamic_cast<TTree*>(gDirectory->FindObjectAny("background"));

    _tree->SetBranchAddress("eventId",&_eventId);
    _tree->SetBranchAddress("volumeId",&_volumeId);
    _tree->SetBranchAddress("totalEDep",&_totalEDep);
    _tree->SetBranchAddress("nonIonizingEDep",&_nonIonizingEDep);
    _tree->SetBranchAddress("time",&_time);
    _tree->SetBranchAddress("proper",&_proper);
    _tree->SetBranchAddress("positionX",&_positionX);
    _tree->SetBranchAddress("positionY",&_positionY);
    _tree->SetBranchAddress("positionZ",&_positionZ);
    _tree->SetBranchAddress("momentumX",&_momentumX);
    _tree->SetBranchAddress("momentumY",&_momentumY);
    _tree->SetBranchAddress("momentumZ",&_momentumZ);
    _tree->SetBranchAddress("stepLength",&_stepLength);
    _tree->SetBranchAddress("endProcessCode",&_endProcessCode);
    _tree->SetBranchAddress("id",&_id);
    _tree->SetBranchAddress("pdgId",&_pdgId);

    _volumeIdStart = (size_t)((TVector2*)_tree->GetUserInfo()->At(0))->X()+0.5;
    _volumeIdEnd = (size_t)((TVector2*)_tree->GetUserInfo()->At(0))->Y()+0.5;

    std::cout<<"CRVResponse uses a background overlay ("<<_backgroundFile<<")"<<std::endl;
    std::cout<<"from counters "<<_volumeIdStart<<" to "<<_volumeIdEnd<<std::endl;
    std::cout<<"with a factor of "<<_overlayFactor<<"."<<std::endl;

    _currentEntry = 0;
    _nEntries = _tree->GetEntries();

    directory->cd();

    produces<StepPointMCCollection>("CRV");
    produces<SimParticleCollection>();

  }

  void CrvBackgroundOverlay::produce(art::Event& event)
  {
    GeomHandle<CosmicRayShield> CRS;

    art::ProductID simPartId(getProductID<SimParticleCollection>());
    SimParticleHelper spHelper(0, simPartId, event);

    std::unique_ptr<StepPointMCCollection> stepPointMCs(new StepPointMCCollection);
    std::unique_ptr<SimParticleCollection> simParticles(new SimParticleCollection);

    TDirectory *directory = gDirectory;
    _file->cd();

    std::map<size_t,int> particleMap;  //can't insert the SimParticles directly into the SimParticleCollection, since it expects the track _ids to be sorted

    int    overlayCount=0;
    size_t currentEventId;
    while(_currentEntry<=_nEntries)
    {
      _tree->GetEntry(_currentEntry);
      if(overlayCount==0 || currentEventId!=_eventId)
      {
        currentEventId=_eventId;
        overlayCount++;
        if(overlayCount>_overlayFactor) break;  //only store StepPointMCs and SimParticles of the number of overlay events indicated by overlayFactor
      }

      _id+=1e10*overlayCount;     //ensures that the same track ID isn't used again for overlays from different overlay events

      particleMap[_id]=_pdgId;    //stores every _id only once, and orderes the _ids

      art::Ptr<SimParticle> particle(spHelper.productID(), art::Ptr<SimParticle>::key_type(_id), spHelper.productGetter());

      size_t nVolumeIds=_volumeIdEnd-_volumeIdStart+1;
      size_t nNewVolumeIds=CRS->getAllCRSScintillatorBars().size();
      for(size_t newVolumeId=_volumeId-_volumeIdStart; newVolumeId<nNewVolumeIds; newVolumeId+=nVolumeIds)
      {
        const CRSScintillatorBar &CRSbar = CRS->getBar(CRSScintillatorBarIndex(newVolumeId));
        const CLHEP::Hep3Vector localPosition(_positionX,_positionY,_positionZ);
        const CLHEP::Hep3Vector &position = CRSbar.toWorld(localPosition);

        CLHEP::Hep3Vector momentum(_momentumX,_momentumY,_momentumZ);

        stepPointMCs->push_back(StepPointMC(particle,
                                            newVolumeId,
                                            _totalEDep,
                                            _nonIonizingEDep,
                                            _time,
                                            _proper,
                                            position,
                                            momentum,
                                            _stepLength,
                                            (mu2e::ProcessCode::enum_type)_endProcessCode));
      }

      _currentEntry++;
      if(_currentEntry==_nEntries)
      {
        _currentEntry=0;
        break;
      }
    }


    std::map<size_t,int>::const_iterator particleIter;
    for(particleIter=particleMap.begin(); particleIter!=particleMap.end(); particleIter++)
    {
      size_t id=particleIter->first;
      size_t idTmp=id - simParticles->delta();  //to counteract the weird way how the map_vector::push_back function changes the first component of the std::pair
      size_t pdgId=particleIter->second;

      simParticles->push_back(std::pair<cet::map_vector_key,mu2e::SimParticle>(cet::map_vector_key(idTmp), SimParticle(cet::map_vector_key(id),    //track ID
                                                                                                                       art::Ptr<SimParticle>(),    //parent SimParticle
                                                                                                                       (mu2e::PDGCode::type)pdgId, //PDG ID
                                                                                                                       art::Ptr<GenParticle>(),    //GenParticle
                                                                                                                       CLHEP::Hep3Vector(),        //position
                                                                                                                       CLHEP::HepLorentzVector(),  //momentum,
                                                                                                                       NAN,                        //startGlobalTime
                                                                                                                       NAN,                        //startProperTime
                                                                                                                       0,                          //startVolumeIndex
                                                                                                                       0,                          //startG4Status
                                                                                                                       mu2e::ProcessCode::unknown, //creationCode
                                                                                                                       1.0)));                     //weight

    }

    directory->cd();

    event.put(std::move(stepPointMCs),"CRV");
    event.put(std::move(simParticles));
  } // end produce

} // end namespace mu2e

using mu2e::CrvBackgroundOverlay;
DEFINE_ART_MODULE(CrvBackgroundOverlay)
