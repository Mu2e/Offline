
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include <cmath>
#include <algorithm>

using namespace std;

namespace mu2e {


  // simple line interpolation, this should be a utility function, FIXME!
  double StrawResponse::PieceLine(std::vector<double> const& xvals, std::vector<double> const& yvals, double xval){
    int imax = int(xvals.size()-2);
    double xbin = (xvals.back()-xvals.front())/(xvals.size()-1);
    int ibin = min(imax,max(0,int(floor((xval-xvals.front())/xbin))));
    double yval = yvals[ibin];
    int jbin = ibin+1;
    double slope = (yvals[jbin]-yvals[ibin])/(xvals[jbin]-xvals[ibin]);
    yval += (xval-xvals[ibin])*slope;
    return yval;
  }

  double StrawResponse::driftDistanceToTime(StrawId strawId,
      double ddist, double phi) const {
    if(_usenonlindrift){
      return _strawDrift->D2T(ddist,phi);
    }
    else{
      return ddist/_lindriftvel; //or return t assuming a constant drift speed of 0.06 mm/ns (for diagnosis)
    }
  }

  double StrawResponse::driftTimeToDistance(StrawId strawId,
      double dtime, double phi) const {
    if(_usenonlindrift){
      return _strawDrift->T2D(dtime,phi);
    }
    else{
      return dtime*_lindriftvel; //or return t assuming a constant drift speed of 0.06 mm/ns (for diagnosis)
    }
  }

  double StrawResponse::driftInstantSpeed(StrawId strawId,
      double doca, double phi) const {
    if(_usenonlindrift){
      return _strawDrift->GetInstantSpeedFromD(doca);
    }else{
      return _lindriftvel;
    }
  }

  double StrawResponse::driftDistanceError(StrawId strawId,
      double ddist, double phi, double DOCA) const {
    if(useDriftError()){
      // maximum drift is the straw radius.  should come from conditions FIXME!
      static double rstraw(2.5);
      double doca = std::min(fabs(DOCA),rstraw);
      size_t idoca = std::min(_derr.size()-1,size_t(floor(_derr.size()*(doca/rstraw))));
      return _derr[idoca];
    }else{
      double rres = _rres_min;
      if(ddist < _rres_rad){
        rres = _rres_min+_rres_max*(_rres_rad-ddist)/_rres_rad;
      }
      return rres;
    }
  }

  double StrawResponse::driftDistanceOffset(StrawId strawId, double ddist, double phi, double DOCA) const {
    return 0;
  }

  double StrawResponse::driftTimeError(StrawId strawId,
      double ddist, double phi, double DOCA) const {
    if (useParameterizedDriftError()){
      if (DOCA > 2.5)
        DOCA = 2.5;
      return PieceLine(_parDriftDocas, _parDriftRes, DOCA);
    }else{
      return driftDistanceError(strawId, ddist, phi, DOCA) / _lindriftvel;
    }
  }

  double StrawResponse::driftTimeOffset(StrawId strawId, double ddist, double phi, double DOCA) const {
    return PieceLine(_parDriftDocas, _parDriftOffsets, DOCA);
  }



  bool StrawResponse::wireDistance(Straw const& straw, double edep,
      double dt, double& wdist, double& wderr, double &halfpv) const {
    bool retval(true);
    double slen = straw.halfLength();
    // convert edep from Mev to KeV (should be standardized, FIXME!)
    double kedep = 1000.0*edep;
    halfpv = halfPropV(straw.id(),kedep);
    wdist = halfpv*(dt);
    wderr = wpRes(kedep,fabs(wdist));
    // truncate positions that exceed the length of the straw (with a buffer): these come from missing a cluster on one end
    if(fabs(wdist) > slen+_wbuf*wderr){
      // move the position to the correct half of the straw
      wdist = copysign(slen*_slfac,wdist);
      // inflat the error
      wderr *= _errfac;
      retval = false;
    }
    return retval;
  }

  double StrawResponse::halfPropV(StrawId strawId, double kedep) const {
    return PieceLine(_edep,_halfvp,kedep);
  }

  double StrawResponse::wpRes(double kedep,double wlen) const {
    // central resolution depends on edep
    double tdres = PieceLine(_edep,_centres,kedep);
    if( wlen > _central){
      // outside the central region the resolution depends linearly on the distance
      // along the wire.  The slope of that also depends on edep
      double wslope = PieceLine(_edep,_resslope,kedep);
      tdres += (wlen-_central)*wslope;
    }
    return tdres;
  }

  void StrawResponse::calibrateTimes(TrkTypes::TDCValues const& tdc,
      TrkTypes::TDCTimes &times, const StrawId &id) const {
    double electronicsTimeDelay = _strawElectronics->electronicsTimeDelay();
    times[StrawEnd::hv] = tdc[StrawEnd::hv]*_strawElectronics->tdcLSB()
      - electronicsTimeDelay + _timeOffsetPanel[id.getPanel()]
      + _timeOffsetStrawHV[id.uniqueStraw()];
    times[StrawEnd::cal] = tdc[StrawEnd::cal]*_strawElectronics->tdcLSB()
      - electronicsTimeDelay + _timeOffsetPanel[id.getPanel()]
      + _timeOffsetStrawCal[id.uniqueStraw()];
  }


  double StrawResponse::saturatedResponse(double vlin) const {
    if (vlin < _vsat)
      return vlin;
    else
      return _vsat;
  }

  double StrawResponse::driftTime(Straw const& straw,
      double tot, double edep) const {
    // straw is present in case of eventual calibration
    size_t totbin = min(_totTBins-1,static_cast<size_t>(tot/_totTBinWidth));
    size_t ebin = min(_totEBins-1,static_cast<size_t>(edep/_totEBinWidth));
    return _totdtime[totbin*_totEBins+ebin];
  }

  double StrawResponse::pathLength(Straw const& straw, double tot) const {
    // needs to be implemented, FIXME!!
    return 5.0;
  }

  void StrawResponse::print(std::ostream& os) const {
    os << endl << "StrawResponse parameters: "  << std::endl;

    printVector(os,"edep",_edep);
    printVector(os,"halfvp",_halfvp);
    os << "central = " << _central << endl;
    printVector(os,"centres",_centres);
    printVector(os,"resslope",_resslope);
    printVector(os,"totdtime",_totdtime);
    os << "usederr = " << _usederr << endl;
    printVector(os,"derr",_derr);
    os << "wbuf = " << _wbuf << endl;
    os << "slfac = " << _slfac << endl;
    os << "errfac = " << _errfac << endl;
    os << "usenonlindrift = " << _usenonlindrift << endl;
    os << "lindriftvel = " << _lindriftvel << endl;
    os << "rres_min = " << _rres_min << endl;
    os << "rres_max = " << _rres_max << endl;
    os << "mint0doca = " << _mint0doca << endl;
    os << "t0shift = " << _t0shift << endl;
    printArray(os,"pmpEnergyScale",_pmpEnergyScale);
    printArray(os,"timeOffsetPanel",_timeOffsetPanel);
    printArray(os,"timeOffsetStrawHV",_timeOffsetStrawHV);
    printArray(os,"timeOffsetStrawCal",_timeOffsetStrawCal);
    os << "electronicsTimeDelay = " << _electronicsTimeDelay << endl;
    os << "gasGain = " << _gasGain << endl;

    os << "analognoise["<<_analognoise.size()<<"] = ";
    for(auto x: _analognoise) os << x << " " ;
    os << endl;

    os << "dVdI["<<_dVdI.size()<<"] = ";
    for(auto x: _dVdI) os << x << " " ;
    os << endl;

    os << "vsat = " << _vsat << endl;
    os << "ADCped = " << _ADCped << endl;
    os << "pmpEnergyScaleAvg = " << _pmpEnergyScaleAvg << endl;

  }

  void StrawResponse::printVector(std::ostream& os, std::string const& name,
      std::vector<double> const& a) const {
    size_t n = a.size();
    if(n<=4) {
      os << name << " ("<<n<<") = ";
      for(auto x : a) os << x << " ";
      os << endl;
    } else {
      os << name <<" ("<<n<<") = "
        << a[0] << " " << a[1] << " ... "
        << a[n-2] << " " << a[n-1] << endl;
    }
  }
  template<typename T, size_t SIZE>
    void StrawResponse::printArray(std::ostream& os, std::string const& name,
        std::array<T,SIZE> const& a) const {
      size_t n = a.size();
      if(n<=4) {
        os << name << " ("<<n<<") = ";
        for(auto x : a) os << x << " ";
        os << endl;
      } else {
        os << name <<" ("<<n<<") = "
          << a[0] << " " << a[1] << " ... "
          << a[n-2] << " " << a[n-1] << endl;
      }

    }

}
