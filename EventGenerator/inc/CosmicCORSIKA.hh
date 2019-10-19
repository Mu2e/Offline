#ifndef CORSIKAEventGenerator_CosmicCORSIKA_hh
#define CORSIKAEventGenerator_CosmicCORSIKA_hh

// Cosmic rays generator using CORSIKA

#include <vector>


#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

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

namespace art
{
  class Run;
}


class GenParticleCollection;

namespace mu2e {

  class CosmicCORSIKA{

    public:
      CosmicCORSIKA(art::Run &run, const SimpleConfig &, CLHEP::HepRandomEngine &engine);
      ~CosmicCORSIKA();
      const double getLiveTime();
      const unsigned int getNumShowers();
      static int corsikaToHepevtID(const int corsikaID);
      virtual void generate(GenParticleCollection &);

    private:

      void genEvent(GenParticleCollection &genParts);
      static double wrapvar( const double var, const double low, const double high);

      std::vector<CLHEP::Hep3Vector> _targetBoxIntersections;
      std::vector<CLHEP::Hep3Vector> _worldIntersections;

      static void calIntersections(CLHEP::Hep3Vector orig, CLHEP::Hep3Vector dir,
                            std::vector<CLHEP::Hep3Vector> &intersections,
                            double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
      static bool pointInBox(double x, double y, double x0, double y0, double x1, double z1);
      static double distance(const CLHEP::Hep3Vector &u, const CLHEP::Hep3Vector &v);

      const float _GeV2MeV = CLHEP::GeV / CLHEP::MeV;

      int _showerInputs; ///< Number of shower inputs to process from

      CLHEP::Hep3Vector _cosmicReferencePointInMu2e;

      double _refY0; ///< Height to which particles will be projected [cm]
      bool _projectToTargetBox;
      std::vector< std::string > _showerInputFiles; ///< Set of CORSIKA shower data files to use
      double _tOffset; ///< Time offset of sample, defaults to zero (no offset) [s]
      double _buffer; ///< Buffer box extensions to cryostat in each direction (6 of them: x_lo,x_hi,y_lo,y_hi,z_lo,z_hi) [mm]
      double _showerAreaExtension; ///< Extend distribution of corsika particles in x,z by this much (e.g. 1000 will extend 10 m in -x, +x, -z, and +z) [mm]

      std::string _refPointChoice;
      bool _geomInfoObtained;
      double _targetBoxXmin;
      double _targetBoxXmax;
      double _targetBoxYmin;
      double _targetBoxYmax;
      double _targetBoxZmin;
      double _targetBoxZmax;
      double _worldXmin;
      double _worldXmax;
      double _worldYmin;
      double _worldYmax;
      double _worldZmin;
      double _worldZmax;
      double _randomXZshift;

      int _fileIndex;

      FILE *in;
      int _loops = 0; // number of loops, necessary to skip the garbage data between blocks
      float _garbage;

      unsigned int _primaries = 0;
      CLHEP::RandFlat _flat;
      double _fluxConstant;
  };  // CosmicCORSIKA

}
#endif
