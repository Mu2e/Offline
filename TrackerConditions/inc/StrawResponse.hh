#ifndef TrackerConditions_StrawResponse_hh
#define TrackerConditions_StrawResponse_hh
//
// StrawResponse collects the net response features of straws 
// used in reconstruction 
//

#include <iostream>
#include <vector>
#include <array>
#include "TrackerGeom/inc/Straw.hh"
#include "DataProducts/inc/TrkTypes.hh"
#include "TrackerConditions/inc/StrawDrift.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
#include "TrackerConditions/inc/StrawPhysics.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"


namespace mu2e {
  class Straw;
  class StrawDrift;
  class StrawId;

  class StrawResponse : virtual public ProditionsEntity {
    public:
    typedef std::shared_ptr<StrawResponse> ptr_t;
    typedef std::shared_ptr<const StrawResponse> cptr_t;

    explicit StrawResponse() : _name("StrawResponse") {}
    explicit StrawResponse( StrawDrift::cptr_t strawDrift,
			    StrawElectronics::cptr_t strawElectronics,
			    StrawPhysics::cptr_t strawPhysics,
        std::vector<double> edep, std::vector<double> halfvp, 
	double central, std::vector<double> centres, 
	std::vector<double> resslope, std::vector<double> totdtime, 
	bool usederr, std::vector<double> derr, 
        double wbuf, double slfac, double errfac, bool usenonlindrift, 
        double lindriftvel, double rres_min, double rres_max, 
        double rres_rad, double mint0doca, double t0shift, 
        std::vector<double> pmpEnergyScale, 
        std::vector<double> timeOffsetPanel, 
        std::vector<double> timeOffsetStrawHV, 
        std::vector<double> timeOffsetStrawCal, 
	double electronicsTimeDelay, double gasGain, 
        std::array<double,StrawElectronics::npaths> analognoise, 
	std::array<double,StrawElectronics::npaths> dVdI, 
        double vsat, double ADCped, double pmpEnergyScaleAvg) : 
      _name("StrawResponse"), 
      _strawDrift(strawDrift),
      _strawElectronics(strawElectronics),
      _strawPhysics(strawPhysics),
      _edep(edep), _halfvp(halfvp), _central(central), _centres(centres), 
      _resslope(resslope), _totdtime(totdtime), _usederr(usederr), 
      _derr(derr), _wbuf(wbuf), _slfac(slfac), _errfac(errfac), 
      _usenonlindrift(usenonlindrift), _lindriftvel(lindriftvel), 
      _rres_min(rres_min), _rres_max(rres_max), _rres_rad(rres_rad), 
      _mint0doca(mint0doca), _t0shift(t0shift), 
      _pmpEnergyScale(pmpEnergyScale), _timeOffsetPanel(timeOffsetPanel), 
      _timeOffsetStrawHV(timeOffsetStrawHV), 
      _timeOffsetStrawCal(timeOffsetStrawCal), 
      _electronicsTimeDelay(electronicsTimeDelay), 
      _gasGain(gasGain), _analognoise(analognoise), 
      _dVdI(dVdI), _vsat(vsat), _ADCped(ADCped), 
      _pmpEnergyScaleAvg(pmpEnergyScaleAvg)  {}

    virtual ~StrawResponse() {}

    bool wireDistance(Straw const& straw, double edep, double dt, 
		    double& wdist, double& wderr, double& halfpv) const;
    bool useDriftError() const { return _usederr; } 
    bool useNonLinearDrift() const { return _usenonlindrift; }
    double Mint0doca() const { return _mint0doca;}
    double strawGain() const { return _strawPhysics->strawGain();}

    double halfPropV(StrawId strawId, double kedep) const;

    double driftDistanceToTime(StrawId strawId, double ddist, double phi) const;
    double driftTimeToDistance(StrawId strawId, double dtime, double phi) const;
    double driftInstantSpeed(StrawId strawId, double ddist, double phi) const;
    double driftConstantSpeed() const {return _lindriftvel;} // constant value used for annealing errors, should be close to average velocity
    double driftDistanceError(StrawId strawId, double ddist, double phi, double DOCA) const;
    double driftDistanceOffset(StrawId strawId, double ddist, double phi, double DOCA) const;

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
    double driftTime(Straw const& straw, double tot, double edep) const;
    double pathLength(Straw const& straw, double tot) const;
    //      double pathLength(StrawHit const& strawhit, double theta) const;

    void print(std::ostream& os) const;
    void printVector(std::ostream& os, std::string const& name, 
		    std::vector<double> const& a) const;

    // StrawElectronics functions we are allowed to use
    inline size_t nADCPreSamples() const { return _strawElectronics->nADCPreSamples(); }
    inline double adcLSB() const { return _strawElectronics->adcLSB(); }
    inline double totLSB() const { return _strawElectronics->totLSB(); }
    inline double adcPeriod() const { return _strawElectronics->adcPeriod(); }
    inline uint16_t maxADC() const { return _strawElectronics->maxADC(); }
    // StrawPhysics functions we are allowed to use
    inline double ionizationEnergy(double q) const { return _strawPhysics->ionizationEnergy(q); }
    
    const std::string& name() const { return _name; }
  private:

    std::string _name;

    // helper functions
    double wpRes(double kedep, double wdist) const;
    static double PieceLine(std::vector<double> const& xvals, 
			    std::vector<double> const& yvals, double xval);

    StrawDrift::cptr_t _strawDrift;
    StrawElectronics::cptr_t _strawElectronics;
    StrawPhysics::cptr_t _strawPhysics;

    // parametric data for calibration functions
    // TD reconstruction uses 1/2 the propagation velocity and depends on the
    // Dependence on position and straw length still needed FIXME!
    // (reconstructed) energy deposit
    std::vector<double> _edep; // energy deposit boundaries
    std::vector<double> _halfvp; // effective 1/2 propagation velocity by edep
    double _central; // max wire distance for central wire region
    std::vector<double> _centres; // wire center resolution by edep
    std::vector<double> _resslope; // resolution slope vs position by edep
    std::vector<double> _totdtime;
    bool _usederr; // flag to use the doca-dependent calibration of the drift error
    std::vector<double> _derr; // parameters describing the drift error function
    double _wbuf; // buffer at the edge of the straws, in terms of sigma
    double _slfac; // factor of straw length to set 'missing cluster' hits
    double _errfac; // error inflation for 'missing cluster' hits
    
    bool _usenonlindrift;
    double _lindriftvel;
    double _rres_min;
    double _rres_max;
    double _rres_rad;
    double _mint0doca;  // minimum doca for t0 calculation.  Note this is a SIGNED QUANTITITY
    
    double _t0shift;
    std::vector<double> _pmpEnergyScale;
    std::vector<double> _timeOffsetPanel;
    std::vector<double> _timeOffsetStrawHV;
    std::vector<double> _timeOffsetStrawCal;
    double _electronicsTimeDelay;
    double _gasGain;
    std::array<double,StrawElectronics::npaths> _analognoise; //noise (mVolt) from the straw itself
    std::array<double,StrawElectronics::npaths> _dVdI; // scale factor between charge and voltage (milliVolts/picoCoulombs)
    double _vsat;
    double _ADCped;
    double _pmpEnergyScaleAvg;


  };
}
#endif

