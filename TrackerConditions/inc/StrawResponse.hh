#ifndef TrackerConditions_StrawResponse_hh
#define TrackerConditions_StrawResponse_hh
//
// StrawResponse collects the net response features of straws used in reconstruction 
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>
#include <vector>
// Mu2e includes
#include "TrackerConditions/inc/Types.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
#include "TrackerConditions/inc/StrawPhysics.hh"
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {
  class Straw;
  class StrawDrift;
  class StrawId;
  class StrawResponse : virtual public ConditionsEntity {
    public:
      // construct from parameters
      explicit StrawResponse(fhicl::ParameterSet const& pset);
      virtual ~StrawResponse();
      bool wireDistance(Straw const& straw, float edep, float dt, float& wdist, float& wderr) const;
      bool useDriftError() const { return _usederr; } 
      bool useNonLinearDrift() const { return _usenonlindrift; }
      double Mint0doca() const { return _mint0doca;}
      double strawGain() const { return _gasGain;}

      float halfPropV(StrawId strawId, float kedep) const;

      double driftDistanceToTime(StrawId strawId, double ddist, double phi) const;
      double driftTimeToDistance(StrawId strawId, double dtime, double phi) const;
      double driftInstantSpeed(StrawId strawId, double ddist, double phi) const;
      double driftConstantSpeed() const {return _lindriftvel;} // constant value used for annealing errors, should be close to average velocity
      double driftDistanceError(StrawId strawId, double ddist, double phi, float DOCA) const;
      double driftDistanceOffset(StrawId strawId, double ddist, double phi, float DOCA) const;

      double peakMinusPedestalEnergyScale() const { return _pmpEnergyScaleAvg; }
      double peakMinusPedestalEnergyScale(StrawId sid) const { return _pmpEnergyScale[sid.getStraw()]; }
      double analogNoise(StrawElectronics::Path ipath) const { return _analognoise[ipath]; }  // incoherent noise
      double fallTime(StrawElectronics::Path ipath) const { return 22.;} //FIXME
      double currentToVoltage(StrawElectronics::Path ipath) const { return _dVdI[ipath]; }
      double saturatedResponse(double vlin) const;
      double ADCPedestal() const { return _ADCped; };
      double t0shift() const { return _t0shift; }

      // converts times from TDC times to time relative to Event Window
      // removes channel to channel delays and overall electronics time delay
      void calibrateTimes(TrkTypes::TDCValues const& tdc, TrkTypes::TDCTimes &times, const StrawId &id) const;
      // approximate drift distatnce from ToT value
      double driftTime(Straw const& straw, float tot) const;
      double pathLength(Straw const& straw, float tot) const;
//      double pathLength(StrawHit const& strawhit, double theta) const;

      void print(std::ostream& os) const;

      // StrawElectronics functions we are allowed to use
      inline size_t nADCPreSamples() const { return _strawele->nADCPreSamples(); }
      inline double adcLSB() const { return _strawele->adcLSB(); }
      inline double totLSB() const { return _strawele->totLSB(); }
      inline double adcPeriod() const { return _strawele->adcPeriod(); }
      inline uint16_t maxADC() const { return _strawele->maxADC(); }
      // StrawPhysics functions we are allowed to use
      inline double ionizationEnergy(double q) const { return _strawphys->ionizationEnergy(q); }

    private:
// helper functions
      float wpRes(float kedep, float wdist) const;
      static float PieceLine(std::vector<float> const& xvals, std::vector<float> const& yvals, float xval);
      void initializeStrawDrift() const;

      StrawDrift *_strawDrift;

      // parametric data for calibration functions
      // TD reconstruction uses 1/2 the propagation velocity and depends on the
      // Dependence on position and straw length still needed FIXME!
      // (reconstructed) energy deposit
      std::vector<float> _edep; // energy deposit boundaries
      std::vector<float> _halfvp; // effective 1/2 propagation velocity by edep
      float _central; // max wire distance for central wire region
      std::vector<float> _centres; // wire center resolution by edep
      std::vector<float> _resslope; // resolution slope vs position by edep
      bool _usederr; // flag to use the doca-dependent calibration of the drift error
      std::vector<float> _derr; // parameters describing the drift error function
      float _wbuf; // buffer at the edge of the straws, in terms of sigma
      float _slfac; // factor of straw length to set 'missing cluster' hits
      float _errfac; // error inflation for 'missing cluster' hits

      std::string _driftFile;
      float _wirevoltage;
      int _phiBins;
      int _dIntegrationBins;
      bool _usenonlindrift;
      double _lindriftvel;
      double _rres_min;
      double _rres_max;
      double _rres_rad;
      double _mint0doca;  // minimum doca for t0 calculation.  Note this is a SIGNED QUANTITITY

      double _TOTIntercept, _TOTSlope, _TOTmin, _TOTmax;
      double _t0shift;
      double _gasGain;
      std::vector<double> _pmpEnergyScale;
      double _pmpEnergyScaleAvg;
      double _analognoise[StrawElectronics::npaths]; //noise (mVolt) from the straw itself
      double _dVdI[StrawElectronics::npaths]; // scale factor between charge and voltage (milliVolts/picoCoulombs)
      double _vsat;
      double _ADCped;


      double _electronicsTimeDelay;
      std::vector<double> _timeOffsetPanel;
      std::vector<double> _timeOffsetStrawHV;
      std::vector<double> _timeOffsetStrawCal;
      
      ConditionsHandle<StrawElectronics> _strawele;
      ConditionsHandle<StrawPhysics> _strawphys;

  };
}
#endif

