// Mu2e includes.
#include "ConfigTools/inc/SimpleConfig.hh"
#include "ConfigTools/inc/requireUniqueKey.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"
#include "SeedService/inc/SeedService.hh"

// Includes from art and its toolchain.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandFlat.h"

#include "Mu2eUtilities/inc/VectorVolume.hh"


using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;


namespace mu2e {

  class CorsikaEventGenerator : public art::EDProducer {
    public:

      struct Config
      {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<std::string> corsikaModuleLabel{Name("corsikaModuleLabel"), Comment("Reference point in the coordinate system"), "FromCorsikaBinary"};
        fhicl::Atom<std::string> refPointChoice{Name("refPointChoice"), Comment("Reference point in the coordinate system"), "UNDEFINED"};
        fhicl::Atom<bool> projectToTargetBox{Name("projectToTargetBox"), Comment("Store only events that cross the target box"), false};
        fhicl::Atom<float> targetBoxXmin{Name("targetBoxXmin"), Comment("Extension of the generation plane"), -5000};
        fhicl::Atom<float> targetBoxXmax{Name("targetBoxXmax"), Comment("Extension of the generation plane"), 5000};
        fhicl::Atom<float> targetBoxYmin{Name("targetBoxYmin"), Comment("Extension of the generation plane"), -5000};
        fhicl::Atom<float> targetBoxYmax{Name("targetBoxYmax"), Comment("Extension of the generation plane"), 5000};
        fhicl::Atom<float> targetBoxZmin{Name("targetBoxZmin"), Comment("Extension of the generation plane"), -5000};
        fhicl::Atom<float> targetBoxZmax{Name("targetBoxZmax"), Comment("Extension of the generation plane"), 5000};
      };
      typedef art::EDProducer::Table<Config> Parameters;

      explicit CorsikaEventGenerator(const Parameters &conf);
      // Accept compiler written d'tor.  Modules are never moved or copied.
      virtual void produce (art::Event& e);
      virtual void beginSubRun(art::SubRun &sr);

    private:

      std::vector<CLHEP::Hep3Vector> _targetBoxIntersections;
      std::vector<CLHEP::Hep3Vector> _worldIntersections;

      bool  _geomInfoObtained = false;
      Config _conf;
      std::string _corsikaModuleLabel =  "FromCorsikaBinary";

      float _worldXmin = 0;
      float _worldXmax = 0;
      float _worldYmin = 0;
      float _worldYmax = 0;
      float _worldZmin = 0;
      float _worldZmax = 0;
      float _targetBoxXmin = 0;
      float _targetBoxXmax = 0;
      float _targetBoxYmin = 0;
      float _targetBoxYmax = 0;
      float _targetBoxZmin = 0;
      float _targetBoxZmax = 0;
      bool _projectToTargetBox = false;
      Hep3Vector _cosmicReferencePointInMu2e;
      std::string _refPointChoice;

  };

  CorsikaEventGenerator::CorsikaEventGenerator(const Parameters &conf) : EDProducer{conf},
                                                                         _conf(conf()),
                                                                         _corsikaModuleLabel(_conf.corsikaModuleLabel()),
                                                                         _targetBoxXmin(_conf.targetBoxXmin()),
                                                                         _targetBoxXmax(_conf.targetBoxXmax()),
                                                                         _targetBoxYmin(_conf.targetBoxYmin()),
                                                                         _targetBoxYmax(_conf.targetBoxYmax()),
                                                                         _targetBoxZmin(_conf.targetBoxZmin()),
                                                                         _targetBoxZmax(_conf.targetBoxZmax()),
                                                                         _projectToTargetBox(_conf.projectToTargetBox()),
                                                                         _refPointChoice(_conf.refPointChoice())
  {
    produces<GenParticleCollection>();
  }


  void CorsikaEventGenerator::beginSubRun( art::SubRun &subrun){

  }

  void CorsikaEventGenerator::produce(art::Event &evt)
  {
    if (!_geomInfoObtained) {
      GeomHandle<Mu2eEnvelope> env;
      GeomHandle<WorldG4> worldGeom;
      GeomHandle<DetectorSystem> detsys;

      const float deltaX = 1; // mm
      _worldXmin = worldGeom->mu2eOriginInWorld().x() - worldGeom->halfLengths()[0] + deltaX;
      _worldXmax = worldGeom->mu2eOriginInWorld().x() + worldGeom->halfLengths()[0] - deltaX;
      _worldYmin = env->ymin() + deltaX;
      _worldYmax = env->ymax() - deltaX;
      _worldZmin = worldGeom->mu2eOriginInWorld().z() - worldGeom->halfLengths()[2] + deltaX;
      _worldZmax = worldGeom->mu2eOriginInWorld().z() + worldGeom->halfLengths()[2] - deltaX;

      if (_refPointChoice == "TRACKER")
      {
        _cosmicReferencePointInMu2e = Hep3Vector(detsys->getOrigin().x(),
                                                 0, detsys->getOrigin().z());
      }
      else if (_refPointChoice == "EXTMONFNAL")
      {
        GeomHandle<ExtMonFNAL::ExtMon> extMonFNAL;
        _cosmicReferencePointInMu2e =
            Hep3Vector(extMonFNAL->detectorCenterInMu2e().x(),
                       0, extMonFNAL->detectorCenterInMu2e().z());
      }
      else if (_refPointChoice == "CALO")
      {
        GeomHandle<Calorimeter> calorimeter;
        _cosmicReferencePointInMu2e = Hep3Vector(detsys->getOrigin().x(),
                                                 0, calorimeter->disk(0).geomInfo().origin().z());
      }
      else if (_refPointChoice == "UNDEFINED")
        _cosmicReferencePointInMu2e = Hep3Vector(0., 0, 0.);

      _geomInfoObtained = true;
    }

    art::Handle<mu2e::GenParticleCollection> corsikaParticles;
    evt.getByLabel(_corsikaModuleLabel, corsikaParticles);
    const GenParticleCollection &particles(*corsikaParticles);
    std::unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);

    for (GenParticleCollection::const_iterator i = particles.begin(); i != particles.end(); ++i)
    {
      GenParticle const &particle = *i;
      const HepLorentzVector mom4 = particle.momentum();

      const Hep3Vector position(
          particle.position().x() + _cosmicReferencePointInMu2e.x(),
          particle.position().y(),
          particle.position().z() + _cosmicReferencePointInMu2e.z());

      if (_projectToTargetBox)
      {
        _worldIntersections.clear();
        VectorVolume particleWorld(particle.position(), particle.momentum().vect(),
                                   _worldXmin, _worldXmax, _worldYmin, _worldYmax, _worldZmin, _worldZmax);
        particleWorld.calIntersections(_worldIntersections);

        if (_worldIntersections.size() > 0)
        {
          // Being inside the world volume, the intersection can be only one.
          const Hep3Vector projectedPos = _worldIntersections.at(0);
          genParticles->push_back(GenParticle(static_cast<PDGCode::type>(particle.pdgId()),
                                              GenId::cosmicCORSIKA, projectedPos, mom4,
                                              particle.time()));
        }
      } else {
        genParticles->push_back(GenParticle(static_cast<PDGCode::type>(particle.pdgId()),
                  GenId::cosmicCORSIKA, position, mom4,
                  particle.time()));
      }
    }

    evt.put(std::move(genParticles));
  }

}


using mu2e::CorsikaEventGenerator;
DEFINE_ART_MODULE(CorsikaEventGenerator);
