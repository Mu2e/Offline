// Shifts and scales the StepPointMCs associated with a circular virtual detector and moves them to a new location
// The new location is defined as (OutputX, OutputY, OutputZ) and has radius OutputRadius
//
// Original author: Pawel Plesniak

// stdlib includes
#include <iostream>
#include <string>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"

// Offline includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"


typedef unsigned long VolumeId_type;

namespace mu2e {
  class ShiftVirtualDetectorStepPointMCs : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config
      {
        fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Input tag of StepPointMCs you want to shift")};
        fhicl::Atom<VolumeId_type> VirtualDetectorID{Name("VirtualDetectorID"), Comment("ID of the virtual detector to shift")};
        fhicl::Atom<float> InputX{Name("InputX"), Comment("Input x location of the StepPointMC distribution")};
        fhicl::Atom<float> InputY{Name("InputY"), Comment("Input y location of the StepPointMC distribution")};
        fhicl::Atom<float> InputZ{Name("InputZ"), Comment("Input z location of the StepPointMC distribution")};
        fhicl::Atom<float> OutputX{Name("OutputX"), Comment("Target x location of the StepPointMC distribution")};
        fhicl::Atom<float> OutputY{Name("OutputY"), Comment("Target y location of the StepPointMC distribution")};
        fhicl::Atom<float> OutputZ{Name("OutputZ"), Comment("Target z location of the StepPointMC distribution")};
        fhicl::Atom<float> InputRadius{Name("InputRadius"), Comment("Initial radius of the virutal detector")};
        fhicl::Atom<float> OutputRadius{Name("OutputRadius"), Comment("Target virtual detector radius")};
        fhicl::Atom<int> pdgID{Name("pdgID"), Comment("pdgID. If set to 0, includes all particles")};
      };
      using Parameters=art::EDProducer::Table<Config>;
      explicit ShiftVirtualDetectorStepPointMCs(const Parameters& pset);
      virtual void produce(art::Event& event) override;
      void makeNewStepPointMC(const StepPointMC &step, art::Ptr<SimParticle> const &particle, StepPointMC &newStepPointer, CLHEP::Hep3Vector newPosition, CLHEP::Hep3Vector newPostPosition, CLHEP::Hep3Vector newMomentum, CLHEP::Hep3Vector newPostMomentum);
    private:
      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      VolumeId_type VirtualDetectorID = 0; // Filter out all the StepPointMCs from VD101 for resampling
      const double xHPGeAbsorber = -3944.6, HPGeAbsorberHalfWidth = 25, SSCApertureSpacing = 81.2;
      double inputX = 0.0,  inputY = 0.0,  inputZ = 0.0; // Old position coordinates
      double outputX = 0.0, outputY = 0.0, outputZ = 0.0; // New position coordinates
      double x = 0.0, y = 0.0, z = 0.0; // Buffer variables
      double mass = 0.0, preE = 0.0, preE2 = 0.0, postE = 0.0, postE2 = 0.0, newPrePz = 0.0, newPrePz2 = 0.0, newPostPz = 0.0, newPostPz2 = 0.0;
      double inputRadius = 0.0, outputRadius = 0.0; // Radii of VD101 and SSC apertures
      double scaleFactor = 0.0; // Scale the distance from the centre of VD101 to the hit position by this factor, ratio of VD101 to SSCAperture radii.
      int pdgID = 0;
      CLHEP::Hep3Vector oldPosition, oldPostPosition, oldMomentum, oldPostMomentum, newPosition, newPostPosition, newMomentum, newPostMomentum, oldCentre, distanceFromOldCentre, newCentre, distanceFromNewCentre;
      StepPointMC newStep;
      GlobalConstantsHandle<ParticleDataList> pdt;
  };

  ShiftVirtualDetectorStepPointMCs::ShiftVirtualDetectorStepPointMCs(const Parameters& conf) :
    art::EDProducer{conf},
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    VirtualDetectorID(conf().VirtualDetectorID()),
    inputX(conf().InputX()),
    inputY(conf().InputY()),
    inputZ(conf().InputZ()),
    outputX(conf().OutputX()),
    outputY(conf().OutputY()),
    outputZ(conf().OutputZ()),
    inputRadius(conf().InputRadius()),
    outputRadius(conf().OutputRadius()),
    pdgID(conf().pdgID()) {
      produces<StepPointMCCollection>();
      scaleFactor = outputRadius/inputRadius;
      oldCentre.setX(inputX);
      oldCentre.setY(inputY);
      oldCentre.setZ(inputZ);
      newCentre.setX(outputX);
      newCentre.setY(outputY);
      newCentre.setZ(outputZ);
    };

  void ShiftVirtualDetectorStepPointMCs::produce(art::Event& event) {
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    if (StepPointMCs.empty())
      throw cet::exception("DataError", "Requested data product not found\n");

    // Define the StepPointMCCollection to be added to the event
    std::unique_ptr<StepPointMCCollection> _outputStepPointMCs(new StepPointMCCollection);

    // Check if the event has a hit in VirtualDetectorFilterID. If so add it to the collection
    for (const StepPointMC& step : StepPointMCs) {
      if (step.volumeId() != VirtualDetectorID)
        continue;
      art::Ptr<SimParticle> const &particle = step.simParticle();
      if ((particle->pdgId() != pdgID) && (pdgID != 0))
        continue;

      mass = pdt->particle(particle->pdgId()).mass();
      oldPosition = step.position();
      oldPostPosition = step.postPosition();
      oldMomentum = step.momentum();
      oldPostMomentum = step.postMomentum();

      // If the step is in the area of the HPGe absorber, skip these particles to avoid double counting the absorber particles and replace them with those symmetrically across in front of the LaBr aperture
      if (std::abs(oldPosition.x()-xHPGeAbsorber) < HPGeAbsorberHalfWidth && std::abs(oldPosition.y()) < HPGeAbsorberHalfWidth)
        continue;

      distanceFromOldCentre = oldPosition - oldCentre; // Calculate the vector from the hit position to the center of VD101
      distanceFromNewCentre = distanceFromOldCentre * scaleFactor; // Scale it to determine the new offset vector
      newPosition = newCentre + distanceFromNewCentre; // Calculate the effective same position in the new distribution
      newPosition.setZ(newCentre.z());

      newPostPosition = newPosition + (oldPostPosition - oldPosition) * scaleFactor;
      newPostPosition.setZ(newCentre.z() + (oldPostPosition - oldPosition).z());

      preE = std::sqrt(oldMomentum.mag2() + mass * mass);
      newMomentum.setX(oldMomentum.x() * scaleFactor);
      newMomentum.setY(oldMomentum.y() * scaleFactor);
      newMomentum.setZ(preE - oldMomentum.perp() * scaleFactor);
      if (newMomentum.mag() - oldMomentum.mag() > std::numeric_limits<double>::epsilon())
        throw cet::exception("LogicError", "Difference in momentum magnitude is large");

      postE = std::sqrt(oldPostMomentum.mag2()+mass*mass);
      newPostMomentum.setX(oldPostMomentum.x() * scaleFactor);
      newPostMomentum.setY(oldPostMomentum.y() * scaleFactor);
      newPostMomentum.setZ(postE - oldPostMomentum.perp() * scaleFactor);
      if (newPostMomentum.mag() - oldPostMomentum.mag() > std::numeric_limits<double>::epsilon())
        throw cet::exception("LogicError", "Difference in momentum magnitude is large");

      makeNewStepPointMC(step, particle, newStep, newPosition, newPostPosition, newMomentum, newPostMomentum); // Convert it into a StepPointMC
      _outputStepPointMCs->emplace_back(newStep); // Add it to the collection

      // If the step is in the area of the HPGe absorber but placed symmetrically across in the equivalent space upstream of the LaBr detector, double count it to simulate the effect of removing the absorber
      // First move the point and shift its momentum
      if (std::abs(oldPosition.x() - xHPGeAbsorber - SSCApertureSpacing) < HPGeAbsorberHalfWidth && std::abs(oldPosition.y()) < HPGeAbsorberHalfWidth)
        {
          oldPosition.setX(oldCentre.x()*2-oldPosition.x());
          oldPostPosition.setX(oldCentre.x()*2-oldPostPosition.x());

          distanceFromOldCentre = oldPosition - oldCentre; // Calculate the vector from the hit position to the center of VD101
          distanceFromNewCentre = distanceFromOldCentre * scaleFactor; // Scale it to determine the new offset vector
          newPosition = newCentre + distanceFromNewCentre; // Calculate the effective same position in the new distribution
          newPosition.setZ(newCentre.z());

          newPostPosition = newPosition + (oldPostPosition-oldPosition)*scaleFactor;
          newPostPosition.setZ(newCentre.z()+(oldPostPosition-oldPosition).z());

          newMomentum.setX(-1*scaleFactor*oldMomentum.x());
          newMomentum.setY(scaleFactor*oldMomentum.y());
          newMomentum.setZ(newPrePz);

          newPostMomentum.setX(-1*scaleFactor*oldPostMomentum.x());
          newPostMomentum.setY(scaleFactor*oldPostMomentum.y());
          newPostMomentum.setZ(newPostPz);

          makeNewStepPointMC(step, particle, newStep, newPosition, newPostPosition, newMomentum, newPostMomentum); // Convert it into a StepPointMC
          _outputStepPointMCs->emplace_back(newStep); // Add it to the collection
        }; // end if position is in the absorber shadow
      }; // end for step : StepPointMCs

    event.put(std::move(_outputStepPointMCs));
    return;
  };

  void ShiftVirtualDetectorStepPointMCs::makeNewStepPointMC(const StepPointMC &step, art::Ptr<SimParticle> const &particle, StepPointMC &newStepPointer, CLHEP::Hep3Vector newPosition, CLHEP::Hep3Vector newPostPosition, CLHEP::Hep3Vector newMomentum, CLHEP::Hep3Vector newPostMomentum) {
    StepPointMC newStep(
                        particle,
                        step.volumeId(),
                        step.totalEDep(),
                        step.nonIonizingEDep(),
                        step.visibleEDep(),
                        step.time(),
                        step.properTime(),
                        newPosition,
                        newPostPosition,
                        newMomentum,
                        newPostMomentum,
                        step.stepLength(),
                        step.endProcessCode()
    );
    newStepPointer = newStep;
    return;
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ShiftVirtualDetectorStepPointMCs)
