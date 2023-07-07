// Mu2e includes.
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/ConfigTools/inc/requireUniqueKey.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/G4BeamlineInfo.hh"
#include "Offline/SeedService/inc/SeedService.hh"

// Includes from art and its toolchain.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/Mu2eEnvelope.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/CosmicLivetime.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimTimeOffset.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "Offline/GeneralUtilities/inc/safeSqrt.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandFlat.h"

#include "Offline/Mu2eUtilities/inc/VectorVolume.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;


namespace mu2e {

  class CorsikaEventGenerator : public art::EDProducer {
    public:

      struct Config
      {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<std::string> corsikaModuleLabel{Name("corsikaModuleLabel"), Comment("Name of CORSIKA module label"), "FromCorsikaBinary"};
        fhicl::Atom<std::string> refPointChoice{Name("refPointChoice"), Comment("Reference point in the coordinate system"), "UNDEFINED"};
        fhicl::Atom<bool> projectToTargetBox{Name("projectToTargetBox"), Comment("Store only events that cross the target box"), false};
        fhicl::Atom<float> targetBoxYmax{Name("targetBoxYmax"), Comment("Top coordinate of the target box on the Y axis")};
        fhicl::Atom<float> intDist{Name("intDist"), Comment("Maximum distance of closest approach for backtracked trajectories")};
        fhicl::OptionalAtom<art::InputTag> simTimeOffset{Name("simTimeOffset"), Comment("Simulation time offset to apply (optional)")};
      };
      typedef art::EDProducer::Table<Config> Parameters;

      explicit CorsikaEventGenerator(const Parameters &conf);
      // Accept compiler written d'tor.  Modules are never moved or copied.
      virtual void produce (art::Event& e);
      virtual void endSubRun(art::SubRun &sr);
      virtual void beginSubRun(art::SubRun &sr);

    private:

      std::vector<CLHEP::Hep3Vector> _worldIntersections;

      bool  _geomInfoObtained = false;
      std::string _corsikaModuleLabel =  "FromCorsikaBinary";

      float _worldXmin = 0;
      float _worldXmax = 0;
      float _worldYmin = 0;
      float _worldYmax = 0;
      float _worldZmin = 0;
      float _worldZmax = 0;

      float _targetBoxYmax;

      bool _projectToTargetBox = false;
      Hep3Vector _cosmicReferencePointInMu2e;
      std::string _refPointChoice;

      float _intDist = 0;

      unsigned long long _primaries = 0;
      float _area = 0;
      float _lowE = 0;
      float _highE = 0;
      float _fluxConstant = 0;

      bool _applyTimeOffset;
      art::InputTag _timeOffsetTag;
      SimTimeOffset _stoff; // time offset for SimParticles and StepPointMCs

  };

  CorsikaEventGenerator::CorsikaEventGenerator(const Parameters &conf) : EDProducer{conf},
                                                                         _corsikaModuleLabel(conf().corsikaModuleLabel()),
                                                                         _targetBoxYmax(conf().targetBoxYmax()),
                                                                         _projectToTargetBox(conf().projectToTargetBox()),
                                                                         _refPointChoice(conf().refPointChoice()),
                                                                         _intDist(conf().intDist()),
                                                                         _applyTimeOffset(conf().simTimeOffset.hasValue()),
                                                                         _stoff(0.0)
  {
    produces<GenParticleCollection>();
    produces<CosmicLivetime,art::InSubRun>();
    if (_applyTimeOffset) {
      _timeOffsetTag = conf().simTimeOffset().value();
    }
  }

  void CorsikaEventGenerator::beginSubRun(art::SubRun &subrun)
  {
    _primaries = 0;
  }

  void CorsikaEventGenerator::endSubRun(art::SubRun &subrun)
  {
    std::unique_ptr<CosmicLivetime> livetime(new CosmicLivetime(_primaries, _area, _lowE, _highE, _fluxConstant));
    std::cout << *livetime << std::endl;
    subrun.put(std::move(livetime), art::fullSubRun());
  }

  void CorsikaEventGenerator::produce(art::Event &evt)
  {

    if (_applyTimeOffset) {
      // find the time offset in the event, and copy it locally
      const auto& stoH = evt.getValidHandle<SimTimeOffset>(_timeOffsetTag);
      _stoff = *stoH;
    }

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

    art::Handle<mu2e::CosmicLivetime> livetime;
    evt.getByLabel(_corsikaModuleLabel, livetime);
    _area = (*livetime).area();
    _lowE = (*livetime).lowE();
    _highE = (*livetime).highE();
    _fluxConstant = (*livetime).fluxConstant();
    _primaries += (*livetime).primaries();

    /*
     * In this block we check if there are two or more particles that intersect
     * before the roof of the Mu2e world. If this happens, the particles are assumed to come
     * from a common mother particle. This is to avoid simulating unphysical processes
     * that can happen between the roof and the interaction point.
     */
    std::vector<CLHEP::Hep3Vector> endPoints(particles.size());
    std::vector<bool> doNotBacktrack(particles.size(), false);

    if (_intDist >= 0) {
      for (GenParticleCollection::const_iterator i = particles.begin(); i != particles.end(); ++i)
      {
        for (GenParticleCollection::const_iterator j = i+1; j != particles.end(); ++j)
        {

          GenParticle const &particle1 = *i;
          GenParticle const &particle2 = *j;

          TwoLinePCA twoLine(particle1.position(), -particle1.momentum().vect(), particle2.position(), -particle2.momentum().vect());
          const CLHEP::Hep3Vector point1 = twoLine.point1();
          const CLHEP::Hep3Vector point2 = twoLine.point2();

          /*
          * If their distance of closest approach is smaller than _intDist
          * and the point of closest approach is between the roof and the
          * top of the target box then we stop the backtracking at interaction
          * point
          */
          if (twoLine.dca() < _intDist &&
              point1.y() > _targetBoxYmax &&
              point1.y() < _worldYmax &&
              point2.y() > _targetBoxYmax &&
              point2.y() < _worldYmax)
          {
            doNotBacktrack[i - particles.begin()] = true;
            doNotBacktrack[j - particles.begin()] = true;
            endPoints[i - particles.begin()] = point1;
            endPoints[j - particles.begin()] = point2;
          }
        }
      }
    }

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
        // Check if the particle needs to be backtracked to the roof or not
        if (doNotBacktrack[i - particles.begin()]) {
          Hep3Vector projectedPos = endPoints[i - particles.begin()];

          projectedPos.setX(projectedPos.x() + _cosmicReferencePointInMu2e.x());
          projectedPos.setZ(projectedPos.z() + _cosmicReferencePointInMu2e.z());

          genParticles->push_back(GenParticle(static_cast<PDGCode::type>(particle.pdgId()),
                                              GenId::cosmicCORSIKA, projectedPos, mom4,
                                              particle.time() + _stoff.timeOffset_));
        } else {
          _worldIntersections.clear();
          VectorVolume particleWorld(position, particle.momentum().vect(),
                                    _worldXmin, _worldXmax, _worldYmin, _worldYmax, _worldZmin, _worldZmax);
          particleWorld.calIntersections(_worldIntersections);

          if (_worldIntersections.size() > 0)
          {
            // Being inside the world volume, the intersection can be only one.
            const Hep3Vector projectedPos = _worldIntersections.at(0);
            genParticles->push_back(GenParticle(static_cast<PDGCode::type>(particle.pdgId()),
                                                GenId::cosmicCORSIKA, projectedPos, mom4,
                                                particle.time() + _stoff.timeOffset_));
          }
        }
      } else {
        genParticles->push_back(GenParticle(static_cast<PDGCode::type>(particle.pdgId()),
                  GenId::cosmicCORSIKA, position, mom4,
                  particle.time() + _stoff.timeOffset_));
      }
    }

    evt.put(std::move(genParticles));
  }

}


using mu2e::CorsikaEventGenerator;
DEFINE_ART_MODULE(CorsikaEventGenerator)
