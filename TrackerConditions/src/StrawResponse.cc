
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "TMath.h"
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

  double ConstrainAngle(double phi) {
    if (phi < 0) {
      phi = -1.0*phi;
    }
    phi = fmod(phi,TMath::Pi());
    if (phi > TMath::Pi()/2.0) {
      phi = phi - 2.0*fmod(phi,TMath::Pi()/2.0);
    }
    return phi;
  }

  double StrawResponse::calibrateDriftDistanceToT2D(double ddist) const {
    if (ddist < 0)
      return _dc[0] + ddist*_dc[2]/(_dc[2]+_dc[0]);

    return sqrt(pow(_dc[1]*ddist+_dc[2],2.0)+_dc[0]*(_dc[0]+2*_dc[2]))-_dc[2];
  }

  double StrawResponse::calibrateT2DToDriftDistance(double t2d) const {
    if (t2d <= _dc[0]){
      return (t2d-_dc[0])*(_dc[2]+_dc[0])/_dc[2];
    }
    return (sqrt(4*pow(_dc[1]*_dc[2],2.0)-4*pow(_dc[1],2.0)*(pow(_dc[0],2.0)+2*_dc[0]*_dc[2]-2*_dc[2]*t2d-pow(t2d,2.0)))-2*_dc[1]*_dc[2])/(2*pow(_dc[1],2.0));
  }

  double StrawResponse::calibrateDriftDistanceToT2DDerivative(double ddist) const {
    if (ddist < 0)
      return _dc[2]/(_dc[2]+_dc[0]);
    return _dc[1]*(_dc[1]*ddist+_dc[2])/sqrt(pow(_dc[0],2)+2*_dc[0]*_dc[2]+pow(_dc[1]*ddist+_dc[2],2));
  }

  double StrawResponse::calibrateT2DToDriftDistanceDerivative(double t2d) const {
    if (t2d <= _dc[0]){
      return (_dc[2]+_dc[0])/_dc[2];
    }
    return (_dc[2] + t2d)/sqrt(pow(_dc[1],2)-(-1*pow(_dc[0],2)-2*_dc[0]*_dc[2]+pow(_dc[2]+t2d,2)));
  }

  StrawResponse::DriftInfo StrawResponse::driftInfoAtDistance(StrawId strawId,
      double ddist, double dtime, double phi, bool forceOld) const {
    StrawResponse::DriftInfo info;
    if (!_useOldDrift && !forceOld){
      //FIXME
      info.distance = driftTimeToDistance(strawId,dtime,phi,forceOld);
      info.speed = driftInstantSpeed(strawId,ddist,phi,forceOld);
      info.variance = pow(driftDistanceError(strawId,ddist,phi,forceOld),2);
    }else{
      //FIXME to be deprecated
      info.distance = driftTimeToDistance(strawId,dtime,phi,true);
      info.speed = driftInstantSpeed(strawId,ddist,phi,true);
      info.variance = pow(driftDistanceError(strawId,ddist,phi,true),2);
    }
    return info;
  }

  double StrawResponse::driftDistanceToTime(StrawId strawId,
      double ddist, double phi, bool forceOld) const {
    if (!_useOldDrift && !forceOld){
      if (_driftIgnorePhi)
        phi = 0;
      double calibrated_ddist = calibrateDriftDistanceToT2D(ddist);
      double dtime = _strawDrift->D2T(calibrated_ddist,phi);
      return dtime;
    }else{
      //FIXME to be deprecated
      if(_usenonlindrift){
        return _strawDrift->D2T(ddist,phi);
      }
      else{
        return ddist/_lindriftvel; //or return t assuming a constant drift speed of 0.06 mm/ns (for diagnosis)
      }
    }
  }

  double StrawResponse::driftTimeError(StrawId strawId,
      double ddist, double phi, bool forceOld) const {
    if (!_useOldDrift && !forceOld){
      if (_driftResIsTime){
        if (ddist > 2.5)
          ddist = 2.5;
        return PieceLine(_parDriftDocas, _parDriftRes, ddist);
      }else{
        double distance_error = driftDistanceError(strawId,ddist,phi,forceOld);
        double speed_at_ddist = driftInstantSpeed(strawId,ddist,phi,forceOld);
        return distance_error/speed_at_ddist;
      }
    }else{
      //FIXME to be deprecated
      if (useParameterizedDriftError()){
        if (ddist > 2.5)
          ddist = 2.5;
        return PieceLine(_parDriftDocas, _parDriftRes, ddist);
      }else{
        return driftDistanceError(strawId, ddist, phi, forceOld) / _lindriftvel;
      }
    }
  }

  double StrawResponse::driftInstantSpeed(StrawId strawId,
      double ddist, double phi, bool forceOld) const {
    if (!_useOldDrift && !forceOld){
      if (_driftIgnorePhi)
        phi = 0;
      double calibrated_ddist = calibrateDriftDistanceToT2D(ddist);
      return _strawDrift->GetInstantSpeedFromD(calibrated_ddist)/calibrateDriftDistanceToT2DDerivative(ddist);
    }else{
      // FIXME to be deprecated
      if(_usenonlindrift){
        return _strawDrift->GetInstantSpeedFromD(ddist);
      }else{
        return _lindriftvel;
      }
    }
  }

  double StrawResponse::driftTimeToDistance(StrawId strawId,
      double dtime, double phi, bool forceOld) const {
    if (!_useOldDrift && !forceOld){
      if (_driftIgnorePhi)
        phi = 0;
      double t2d_driftonly = _strawDrift->T2D(dtime,phi,false);
      double ddist = calibrateT2DToDriftDistance(t2d_driftonly);
      return ddist;
    }else{
      if(_usenonlindrift){
        return _strawDrift->T2D(dtime,phi);
      }
      else{
        return dtime*_lindriftvel; //or return t assuming a constant drift speed of 0.06 mm/ns (for diagnosis)
      }
    }
  }

  double StrawResponse::driftDistanceError(StrawId strawId,
      double ddist, double phi, bool forceOld) const {
    if (!_useOldDrift && !forceOld){
      if (_driftIgnorePhi)
        phi = 0;
      if (_driftResIsTime){
        double time_error = driftTimeError(strawId,ddist,phi,forceOld);
        double speed_at_ddist = driftInstantSpeed(strawId,ddist,phi,forceOld);
        return time_error*speed_at_ddist;
      }else{
        if (ddist > 2.5)
          ddist = 2.5;
        return PieceLine(_parDriftDocas, _parDriftRes, ddist);
      }
    }else{
      // maximum drift is the straw radius.  should come from conditions FIXME!
      static double rstraw(2.5);
      ddist = std::min(fabs(ddist),rstraw);
      size_t idoca = std::min(_derr.size()-1,size_t(floor(_derr.size()*(ddist/rstraw))));
      return _derr[idoca];
    }
  }

  // FIXME to be deprecated
  double StrawResponse::driftDistanceOffset(StrawId strawId, double ddist, double phi, double DOCA) const {
    return 0;
  }

  // FIXME to be deprecated
  double StrawResponse::driftTimeOffset(StrawId strawId, double ddist, double phi, double DOCA) const {
    if (!_useOldDrift){
      return 0;
    }else{
      return PieceLine(_parDriftDocas, _parDriftOffsets, DOCA);
    }
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
    if(fabs(wdist) > slen+_wbuf*wderr && _truncateLongitudinal){
      // move the position to the correct half of the straw
      wdist = copysign(slen*_slfac,wdist);
      // inflat the error
      wderr *= _errfac;
      retval = false;
    }
    return retval;
  }

  double StrawResponse::halfPropV(StrawId strawId, double kedep) const {
    double mean_prop_v = _strawHalfvp[strawId.uniqueStraw()];
    return PieceLine(_edep,_halfvpscale,kedep)*mean_prop_v;
  }

  double StrawResponse::wpRes(double kedep,double wlen) const {
    // central resolution depends on edep
    double tdres = PieceLine(_edep,_centres,kedep);
    if (_rmsLongErrors){
      if( wlen > _central){
        // outside the central region the resolution depends linearly on the distance
        // along the wire.  The slope of that also depends on edep
        double wslope = PieceLine(_edep,_resslope,kedep);
        tdres += (wlen-_central)*wslope;
      }
      return tdres;
    }else{
      double wslope = PieceLine(_edep,_resslope,kedep);

      return tdres + wslope*wlen*wlen;
    }
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

  double StrawResponse::TOTdriftTime(Straw const& straw,
      double tot, double edep) const {
    // straw is present in case of eventual calibration
    size_t totbin = min(_totTBins-1,static_cast<size_t>(tot/_totTBinWidth));
    size_t ebin = min(_totEBins-1,static_cast<size_t>(edep/_totEBinWidth));
    return _totdtime[totbin*_totEBins+ebin];
  }

  double StrawResponse::TOTdriftTimeError(Straw const& straw,
      double tot, double edep) const {
    size_t totbin = min(_totTBins-1,static_cast<size_t>(tot/_totTBinWidth));
    size_t ebin = min(_totEBins-1,static_cast<size_t>(edep/_totEBinWidth));
    return _totderror[totbin*_totEBins+ebin];
  }

  double StrawResponse::pathLength(Straw const& straw, double tot) const {
    // needs to be implemented, FIXME!!
    return 5.0;
  }

  void StrawResponse::print(std::ostream& os) const {
    os << endl << "StrawResponse parameters: "  << std::endl;

    printVector(os,"edep",_edep);
    printVector(os,"ehalfvpscale",_halfvpscale);
    os << "central = " << _central << endl;
    printArray(os,"halfvp",_strawHalfvp);
    printVector(os,"centres",_centres);
    printVector(os,"resslope",_resslope);
    printVector(os,"totdtime",_totdtime);
    printVector(os,"derr",_derr);
    os << "wbuf = " << _wbuf << endl;
    os << "slfac = " << _slfac << endl;
    os << "errfac = " << _errfac << endl;
    os << "usenonlindrift = " << _usenonlindrift << endl;
    os << "lindriftvel = " << _lindriftvel << endl;
    os << "mint0doca = " << _mint0doca << endl;
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
