// Cuts that can be used in Mu2eG4SteppingAction.
//
// Andrei Gaponenko, 2015

#ifndef Mu2eG4_Mu2eG4Cuts_hh
#define Mu2eG4_Mu2eG4Cuts_hh

#include "Mu2eG4/inc/IMu2eG4Cut.hh"

#include <string>
#include <memory>
#include <array>

#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "MCDataProducts/inc/StepPointMCCollection.hh"

namespace mu2e {
  namespace SteppingCuts {

    //================================================================
    // A common implementation for some of the required IMu2eG4Cut methods
    class IOHelper: virtual public IMu2eG4Cut {
    public:
      virtual void declareProducts(art::EDProducer *parent) override;
      virtual void finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) override;
      virtual void beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) override;
      virtual void put(art::Event& event) override;

    protected:
      explicit IOHelper(const std::string& outputName);
      std::string outputName_;
      std::unique_ptr<StepPointMCCollection> output_;
      CLHEP::Hep3Vector mu2eOrigin_;
      const SimParticleHelper *spHelper_;

      void addHit(const G4Step *aStep, ProcessCode endCode);
    };

    //================================================================
    class Union: virtual public IMu2eG4Cut,
                 public IOHelper
    {
    public:
      virtual bool evaluate(const G4Step *aStep);
      static std::unique_ptr<Union> maybe_instance(const fhicl::ParameterSet& pset);

      // Sequences need a different implementation
      virtual void declareProducts(art::EDProducer *parent) override;
      virtual void beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) override;
      virtual void put(art::Event& event) override;
      virtual void finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) override;

    private:
      explicit Union(const fhicl::ParameterSet& pset);
      std::vector<std::unique_ptr<IMu2eG4Cut> > cuts_;
    };

    //================================================================
    class Intersection: virtual public IMu2eG4Cut,
                        public IOHelper
    {
    public:
      virtual bool evaluate(const G4Step *aStep);
      static std::unique_ptr<Intersection> maybe_instance(const fhicl::ParameterSet& pset);

      // Sequences need a different implementation
      virtual void declareProducts(art::EDProducer *parent) override;
      virtual void beginEvent(const art::Event& evt, const SimParticleHelper& spHelper) override;
      virtual void put(art::Event& event) override;
      virtual void finishConstruction(const CLHEP::Hep3Vector& mu2eOriginInWorld) override;

    private:
      explicit Intersection(const fhicl::ParameterSet& pset);
      std::vector<std::unique_ptr<IMu2eG4Cut> > cuts_;
    };

    //================================================================
    class Plane: virtual public IMu2eG4Cut,
                 public IOHelper
    {
    public:
      virtual bool evaluate(const G4Step *aStep);
      static std::unique_ptr<Plane> maybe_instance(const fhicl::ParameterSet& pset);
    private:
      explicit Plane(const fhicl::ParameterSet& pset);
      std::array<double,3> normal_;
      double offset_;
    };


    //================================================================
    class Constant: virtual public IMu2eG4Cut,
                    public IOHelper
    {
    public:
      virtual bool evaluate(const G4Step *aStep);
      static std::unique_ptr<Constant> maybe_instance(const fhicl::ParameterSet& pset);
      explicit Constant(bool val);
    private:
      explicit Constant(const fhicl::ParameterSet& pset);
      bool value_;
    };

    //================================================================
    std::unique_ptr<IMu2eG4Cut> createCuts(const fhicl::ParameterSet& pset);

  }

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4Cuts_hh */
