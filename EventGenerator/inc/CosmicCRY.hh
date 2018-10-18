#ifndef CRYEventGenerator_CosmicCRY_hh
#define CRYEventGenerator_CosmicCRY_hh

// Cosmic rays generator using CRY

#include <vector>

#include "CLHEP/Random/RandEngine.h"

namespace art{
  class Run;
}

class CRYSetup;
class CRYGenerator;
class GenParticleCollection;

namespace mu2e {

  class CosmicCRY{

    public:
      CosmicCRY(art::Run& run, const SimpleConfig& config,
          CLHEP::HepRandomEngine& engine);

      virtual void generate( GenParticleCollection&  );

    private:
      int  _verbose;
      // CRY output options
      bool _returnMuons;
      bool _returnNeutrons;
      bool _returnProtons;
      bool _returnGammas;
      bool _returnElectrons;
      bool _returnPions;
      bool _returnKaons;

      // CRY input parameters:
      // - the date that we want to simulate flux, month-day-year
      // - latitude of the detector site, 41.8 for Fermilab
      // - altitude in meter, CRY accepts 3 values: 0, 2100, 11300, default to 0
      // - sub box length
      int _month;
      int _day;
      int _year;
      double _latitude;
      int _altitude;
      double _subboxLength;

      std::string _setupString;
      std::string _cryDataPath;

      CRYSetup * _crySetup;
      std::shared_ptr<CRYGenerator> _cryGen;

      double _refY0;
      std::string _refPointChoice;
      std::string _directionChoice;
      CLHEP::Hep3Vector _cosmicReferencePointInMu2e;
      bool _vertical;

      bool _projectToEnvelope;

      bool _geomInfoObtained;
      double _envXmin;
      double _envXmax;
      double _envYmin;
      double _envYmax;
      double _envZmin;
      double _envZmax;

      double _GeV2MeV;
      double _m2mm;

      std::vector<CLHEP::Hep3Vector> _envIntersections;
      void calIntersections(CLHEP::Hep3Vector orig, CLHEP::Hep3Vector dir);
      bool pointInBox(double x, double y, double x0, double y0, double x1, double z1);

      void createSetupString();
  };  // CosmicCRY

}
#endif /* end of include guard: CRYEventGenerator_CosmicCRY_hh */
