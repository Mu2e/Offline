#ifndef TrackerConditions_StrawPhysics_hh
#define TrackerConditions_StrawPhysics_hh
//
// StrawPhysics collects the electronics response behavior of a Mu2e straw in
// several functions and parameters
//

// C++ includes
#include <iostream>
#include <array>
#include <vector>
#include <utility>
#include <algorithm>
// CLHEP
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"
// Mu2e includes
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"

#include "Offline/TrackerConditions/inc/StrawDrift.hh"

namespace mu2e {
  class StrawPhysics : virtual public ProditionsEntity {
    public:
      typedef std::shared_ptr<StrawPhysics> ptr_t;
      typedef std::shared_ptr<const StrawPhysics> cptr_t;
      constexpr static const char* cxname = {"StrawPhysics"};

      StrawPhysics(
          std::vector<double> EIonize,
          double meanpath, double eKin, double Qe,
          std::vector<double> intNProb,
          double gasgain, double polyaA,
          double gslope, unsigned nggauss, double vprop,
          double NAverage,
          double EAverage,
          std::vector<double> cdpoly,
          std::vector<double> dtvar,
          bool nonlindrift, double bz,
          StrawDrift::cptr_t strawDrift) :
        ProditionsEntity(cxname),
        _EIonize(EIonize),
        _meanpath(meanpath),  _eKin(eKin), _Qe(Qe),
        _intNProb(intNProb),
        _gasgain(gasgain), _polyaA(polyaA),
        _gslope(gslope), _nggauss(nggauss), _vprop(vprop),
        _NAverage(NAverage),
        _EAverage(EAverage),
        _cdpoly(cdpoly),
        _dtvar(dtvar),
        _nonlindrift(nonlindrift), _bz(bz),
        _strawDrift(strawDrift) {}

      virtual ~StrawPhysics() = default;

      // models.  Note these are different from the corresponding
      // functions used in reconstruction, as those can be wire-
      // dependent and have emergent properties
      unsigned nePerIon(double urand) const; // number of electrons for a given ionization, given a flat random number 0 < x < 1
      unsigned nePerEIon(double Eion) const; // number of electrons for a given ionization energy (approximate)
      double strawGain() const { return _gasgain; } // nominal gain
      double clusterGain(CLHEP::RandGaussQ& rgauss, CLHEP::RandFlat& rflat, unsigned nele) const;
      double driftDistanceToTime(double ddist, double phi) const;  // single cluster!
      double testStrawDrift(double ddist, double phi);
      double driftTimeSpread(double ddist) const; // single electron sigma
      double propagationTime(double wdist) const;
      double meanFreePath() const { return _meanpath; }
      double ionizationEnergy(unsigned nele=1) const;
      double ionizationEnergy(double q) const { return _EAverage*q/_Qe; }  // energy to produce a given charge.  This assumes the internal distribution of the number of electrons/ionization
      double ionizationCharge(double eion) const { return _Qe*eion/_EAverage; }
      double ionizationCharge(unsigned nele=1) const { return _Qe*nele; }
      double meanElectronEnergy() const { return _EAverage; }
      double meanElectronCount() const { return _NAverage; }
      void print(std::ostream& os) const;

    private:

      std::vector<double> _EIonize; // cumulative energy to create N ionization electrons (MeV)
      double _meanpath; // mean free path (mm)
      double _eKin; // average kinetic energy of ioniztion electrons (MeV)
      double _Qe; // charge of a single ionization electron (=e, pC)
      std::vector<double> _intNProb; // integrated probability distribution of the number of e produced per cluster
      double _gasgain; // nominal (average) avalanche gain
      double _polyaA; // 'A' parameter of Polya function used in gain fluctuation
      double _gslope; // slope of gain relative RMS vs 1/sqrt(n)
      unsigned _nggauss; // number of electrons/cluster to switch to a Gaussian model
      // attenuation length of charge down the wire; note
      // there is a short and a long component, each with it's own amplitude
      double _vprop; // (average) propagation velocity
      double _NAverage; // Average number of ionization electrons/ionization
      double _EAverage; // Average energy of a single ionization electron (MeV)
      // parameters describing cluster DtoT
      std::vector<double> _cdpoly;
      std::vector<double> _dtvar; // variance depencence on drift time
      bool _nonlindrift;
      double _bz;

      StrawDrift::cptr_t _strawDrift;
  };
}
#endif

