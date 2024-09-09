// Shifts and scales the StepPointMCs associated with a circular virtual detector and moves them to a new location
// The new location is defined as (OutputX, OutputY, OutputZ) and has radius OutputRadius
// STMStudy is a parameter for the STM rate investigation. It removes the effect of the HPGe absorber
//
// Pawel Plesniak

// stdlib includes
#include <iostream>
#include <string>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "canvas/Utilities/InputTag.h"

// Offline includes
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;

namespace mu2e{
  class ShiftVirtualDetectorStepPointMCs : public art::EDProducer
  {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    typedef unsigned long VolumeId_type;
    struct Config
    {
      fhicl::Atom<std::string> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Input tag of StepPointMCs you want to shift")};
      fhicl::Atom<VolumeId_type> VirtualDetectorID{Name("VirtualDetectorID"), Comment("ID of the virtual detector to shift")};
      fhicl::Atom<float> InputX{Name("InputX"), Comment("Input x location of the StepPointMC distribution")};
      fhicl::Atom<float> InputY{Name("InputY"), Comment("Input y location of the StepPointMC distribution")};
      fhicl::Atom<float> InputZ{Name("InputZ"), Comment("Input z location of the StepPointMC distribution")};
      fhicl::Atom<float> OutputX{Name("OutputX"), Comment("Target x location of the StepPointMC distribution")};
      fhicl::Atom<float> OutputY{Name("OutputY"), Comment("Target y location of the StepPointMC distribution")};
      fhicl::Atom<float> OutputZ{Name("OutputZ"), Comment("Target z location of the StepPointMC distribution")};
      fhicl::Atom<float> InputRadius{Name("InputRadius"), Comment("Initial radius of the virutal detector")};
      fhicl::Atom<float> OutputRadius{Name("OutputRadius"), Comment("Target virtual detector radius")};
      fhicl::Atom<bool> STMStudy{Name("STMStudy"), Comment("If this is true, will reflect the StepPointMCs associated with the upstream of LaBr to mitigate the effect of the HPGe absorber")};
      fhicl::OptionalAtom<bool> verbose{Name("verbose"), Comment("Prints summary")};
    };

    using Parameters=art::EDProducer::Table<Config>;

    explicit ShiftVirtualDetectorStepPointMCs(const Parameters& pset);
    virtual void produce(art::Event& event) override;
    virtual void beginJob() override;
    void makeNewStepPointMC(const StepPointMC &step, art::Ptr<SimParticle> const& particle, StepPointMC &newStepPointer, CLHEP::Hep3Vector newPosition, CLHEP::Hep3Vector newPostPosition, CLHEP::Hep3Vector newMomentum, CLHEP::Hep3Vector newPostMomentum);

  private:
    art::ProductToken<StepPointMCCollection> StepPointMCsToken;
    VolumeId_type VirtualDetectorID = 0; // Filter out all the StepPointMCs from VD101 for resampling
    const float xHPGeAbsorber = -3944.6, HPGeAbsorberHalfWidth = 25, SSCApertureSpacing = 81.2;
    float InputX = 0.0, InputY = 0.0, InputZ = 0.0; // Old position coordinates
    float OutputX = 0.0, OutputY = 0.0, OutputZ = 0.0; // New position coordinates
    float x = 0.0, y = 0.0, z = 0.0; // Buffer variables
    float InputRadius = 0.0, OutputRadius = 0.0; // Radii of VD101 and SSC apertures
    float scale_factor = 0.0; // Scale the distance from the centre of VD101 to the hit position by this factor, ratio of VD101 to SSCAperture radii.
    bool STMStudy = true, verbose = false;
    CLHEP::Hep3Vector oldPosition, oldPostPosition, oldMomentum, oldPostMomentum, newPosition, newPostPosition, newMomentum, newPostMomentum, oldCentre, distanceFromOldCentre, newCentre, distanceFromNewCentre;
    StepPointMC newStep;
  };
  // ===================================================
  ShiftVirtualDetectorStepPointMCs::ShiftVirtualDetectorStepPointMCs(const Parameters& conf) :
    art::EDProducer{conf},
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    VirtualDetectorID(conf().VirtualDetectorID()),
    InputX(conf().InputX()),
    InputY(conf().InputY()),
    InputZ(conf().InputZ()),
    OutputX(conf().OutputX()),
    OutputY(conf().OutputY()),
    OutputZ(conf().OutputZ()),
    InputRadius(conf().InputRadius()),
    OutputRadius(conf().OutputRadius()),
    STMStudy(conf().STMStudy())
    {
      produces<StepPointMCCollection>();
      scale_factor = OutputRadius/InputRadius;
      oldCentre.setX(InputX);
      oldCentre.setY(InputY);
      oldCentre.setZ(InputZ);
      newCentre.setX(OutputX);
      newCentre.setY(OutputY);
      newCentre.setZ(OutputZ);

      auto _verbose = conf().verbose();
      if(_verbose)verbose = *_verbose;
    };
  // ===================================================
  void ShiftVirtualDetectorStepPointMCs::beginJob()
  {
    if (verbose == true)
      std::cout << "Moving StepPointMCs in VD" << VirtualDetectorID << " from " << oldCentre << " to " << newCentre << std::endl;
  };
  // ===================================================
  void ShiftVirtualDetectorStepPointMCs::produce(art::Event& event)
  {
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);

    // Define the StepPointMCCollection to be added to the event
    std::unique_ptr<StepPointMCCollection> _outputStepPointMCs(new StepPointMCCollection);

    // Check if the event has a hit in VirtualDetectorFilterID. If so add it to the collection
    for (const StepPointMC& step : StepPointMCs)
      {
        if (step.volumeId() == VirtualDetectorID) // This could be virtualDetectorId from StepPointMC
          {
            const art::Ptr<SimParticle> particle = step.simParticle();
            oldPosition = step.position();
            oldPostPosition = step.postPosition();
            oldMomentum = step.momentum();
            oldPostMomentum = step.postMomentum();

            // If the step is in the area of the HPGe absorber, skip these particles to avoid double counting the absorber particles and replace them with those symmetrically across in front of the LaBr aperture
            if (STMStudy && (std::abs(oldPosition.x()-xHPGeAbsorber) < HPGeAbsorberHalfWidth) && (std::abs(oldPosition.y()) < HPGeAbsorberHalfWidth))
              continue;

            distanceFromOldCentre = oldPosition - oldCentre; // Calculate the vector from the hit position to the center of VD101
            distanceFromNewCentre = distanceFromOldCentre * scale_factor; // Scale it to determine the new offset vector
            newPosition = newCentre + distanceFromNewCentre; // Calculate the effective same position in the new distribution
            newPosition.setZ(newCentre.z());
            newPostPosition = newPosition + (oldPostPosition-oldPosition)*scale_factor;
            newPostPosition.setZ(newCentre.z());

            newMomentum = oldMomentum;
            newMomentum.setX(scale_factor*newMomentum.x());
            newMomentum.setY(scale_factor*newMomentum.y());

            newPostMomentum = oldPostMomentum;
            newPostMomentum.setX(scale_factor*newPostMomentum.x());
            newPostMomentum.setY(scale_factor*newPostMomentum.y());

            makeNewStepPointMC(step, particle, newStep, newPosition, newPostPosition, newMomentum, newPostMomentum); // Convert it into a StepPointMC
            _outputStepPointMCs->emplace_back(newStep); // Add it to the collection

            // If the step is in the area of the HPGe absorber but placed symmetrically across in the equivalent space upstream of the LaBr detector, double count it to simulate the effect of removing the absorber
            // First move the point and shift its momentum
            if (STMStudy && (std::abs(oldPosition.x()-xHPGeAbsorber-SSCApertureSpacing) < HPGeAbsorberHalfWidth) && (std::abs(oldPosition.y()) < HPGeAbsorberHalfWidth))
              {
                oldPosition.setX(oldCentre.x()*2-oldPosition.x());
                distanceFromOldCentre = oldPosition - oldCentre; // Calculate the vector from the hit position to the center of VD101
                distanceFromNewCentre = distanceFromOldCentre * scale_factor; // Scale it to determine the new offset vector
                newPosition = newCentre + distanceFromNewCentre; // Calculate the effective same position in the new distribution
                newPosition.setZ(newCentre.z());
                newPostPosition = newPosition + (oldPostPosition-oldPosition)*scale_factor;
                newPostPosition.setZ(newCentre.z());

                newMomentum = oldMomentum;
                newMomentum.setX(-1*scale_factor*newMomentum.x());
                newMomentum.setY(scale_factor*newMomentum.y());

                newPostMomentum = oldPostMomentum;
                newPostMomentum.setX(-1*scale_factor*newPostMomentum.x());
                newPostMomentum.setY(scale_factor*newPostMomentum.y());

                makeNewStepPointMC(step, particle, newStep, newPosition, newPostPosition, newMomentum, newPostMomentum); // Convert it into a StepPointMC
                _outputStepPointMCs->emplace_back(newStep); // Add it to the collection
            };
          };
      };
    event.put(std::move(_outputStepPointMCs));
    return;
  };
  // ===================================================
  void ShiftVirtualDetectorStepPointMCs::makeNewStepPointMC(const StepPointMC &step, art::Ptr<SimParticle> const& particle, StepPointMC &newStepPointer, CLHEP::Hep3Vector newPosition, CLHEP::Hep3Vector newPostPosition, CLHEP::Hep3Vector newMomentum, CLHEP::Hep3Vector newPostMomentum)
  {
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
  }
}

DEFINE_ART_MODULE(mu2e::ShiftVirtualDetectorStepPointMCs)
