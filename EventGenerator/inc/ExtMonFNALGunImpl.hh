#ifndef ExtMonFNALGunImpl_hh
#define ExtMonFNALGunImpl_hh

// A generator to conveniently generate particles in the acceptance of
// the ExtMonFNAl filter.  The particles are generated at the entrance
// to the filter channel, in a cone around the axis of the first
// collimator.
//
//
// Original author Andrei Gaponenko, 2012
//
// The position is given in the Mu2e coordinate system.
//

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/TupleAs.h"

#include "EventGenerator/inc/ParticleGunImpl.hh"
#include "MCDataProducts/inc/GenParticle.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

// Forward references.
namespace art{ class Run; }

namespace mu2e {

  class ExtMonFNALBuilding;

  class ExtMonFNALGunImpl {

  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      using Hep3Vector_t = CLHEP::Hep3Vector(double, double, double);

      fhicl::Atom<double> multiplicity {
        Name("multiplicity"),
          Comment("The number of generated particles is sampled from the Poisson\n"
                  "distribution with mean==multiplicity for mulitplicity>0,\n"
                  "or is fixed to |multiplicity| for negative integer inputs."
                  ),
          -1
          };

      fhicl::Atom<int> pdgId {
        Name("pdgId"),
          Comment("PDG ID of the produced GenParticles"),
          2212 // proton
          };

      fhicl::OptionalAtom<double> pmin {
        Name("pmin"),
          Comment("Override generated particle min momentum.  The default is the nominal ExtMon signal momentum.")
          };

      fhicl::OptionalAtom<double> pmax {
        Name("pmax"),
          Comment("Override generated particle max momentum.  The default is the nominal ExtMon signal momentum.")
          };

      fhicl::Atom<double> coneAngleMin {
        Name("coneAngleMin"),
          Comment("Angle, in radians, w.r.t. the reference axis.\n"
                  "Particle momenta point between cones with min and max given angles."),
          0.
          };

      fhicl::Atom<double> coneAngleMax {
        Name("coneAngleMax"),
          Comment("Angle, in radians, w.r.t. the reference axis.\n"
                  "Particle momenta point between cones with min and max given angles."),
          0.
          };

      fhicl::Atom<double> tmin { Name("tmin"), Comment("Earliest generated time."), 0. };
      fhicl::Atom<double> tmax { Name("tmax"), Comment("Latest generated time."), 0. };

      fhicl::TupleAs<Hep3Vector_t> offset {
        Name("offset"),
          Comment("The 3D offset of the center of generated particle origin distribution from the reference point."),
          CLHEP::Hep3Vector(0.,0.,0.)
          };

      fhicl::TupleAs<Hep3Vector_t> halfSize {
        Name("halfSize"),
          Comment("The origin of generated particles is sampled uniformly inside a box of the give half size."),
          CLHEP::Hep3Vector(0.,0.,0.)
          };

      fhicl::Atom<std::string> histDirName {
        Name("histDirName"),
          Comment("A non-empty string will activate histogramming and write out results into the named ROOT subdir."),
          ""
          };

      fhicl::Atom<std::string> reference {
        Name("reference"),
          Comment("Particle origin reference point.  Allowed values are:\n"
                  "filter, detector, productionTargetEntrance, productionTargetExit, productionTargetCenter."
                  )
          };

      fhicl::Atom<int> verbosityLevel {
        Name("verbosityLevel"),
          Comment("A number from 0 to 3, larger values mean more printouts."),
          0
          };

    };

    explicit ExtMonFNALGunImpl(CLHEP::HepRandomEngine& engine, const Config& conf);

    // adds generated particles to the collection
    virtual void generate(GenParticleCollection& out);

  private:
    ParticleGunImpl m_gun;

    CLHEP::HepRotation m_rotation;
    CLHEP::Hep3Vector  m_translation;

    void initGeom(const std::string& refPoint);

    // From local coordinates used by our m_gun to mu2e frame
    GenParticle transform(const GenParticle& in) const;

    static double getpmin(const Config& conf, const ExtMonFNALBuilding& emb);
    static double getpmax(const Config& conf, const ExtMonFNALBuilding& emb);
  };

} // end namespace mu2e,

#endif /* ExtMonFNALGunImpl_hh */
