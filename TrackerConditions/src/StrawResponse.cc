#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "cetlib_except/exception.h"
#include "TMath.h"
#include <cmath>
#include <algorithm>
#include "gsl/gsl_fit.h"

using namespace std;

namespace mu2e {
  double StrawResponse::rstraw_(2.5);  // should come from geometry, TODO

  // simple line interpolation, this should be a utility function, TODO
  // This only works if the bins are uniform, should be rewritten TODO
  double StrawResponse::PieceLine(std::vector<double> const& xvals, std::vector<double> const& yvals, double xval){
    int imax = int(xvals.size()-2);
    double xbin = (xvals.back()-xvals.front())/(xvals.size()-1);
    int ibin = min(imax,max(0,int(floor((xval-xvals.front())/xbin))));
    double yval = yvals[ibin];
    auto jbin = ibin+1;
    double slope = (yvals[jbin]-yvals[ibin])/xbin;
    yval += (xval-xvals[ibin])*slope;
    return yval;
  }

  double StrawResponse::PieceLineDrift(std::vector<double> const& bins,std::vector<double> const& yvals, double xval){
      int imax = yvals.size()-1;
    double xbin = (bins[1]-bins[0])/yvals.size();
    int ibin = min(imax,max(0,int(floor((xval-bins[0])/xbin))));
    double yval = yvals[ibin];
    double slope;
    if(ibin < imax) {
      auto jbin = ibin+1;
      slope = (yvals[jbin]-yvals[ibin])/xbin;
    } else {
      auto jbin = ibin-1;
      slope = (yvals[ibin]-yvals[jbin])/xbin;
    }
    yval += (xval-(ibin+0.5)*xbin)*slope;
    return yval;
  }

  void StrawResponse::interpolateCalib(std::vector<double> const& bins,std::vector<double> const& yvals, double xval,
      int halfrange, double& value, double& slope) {
    int maxindex = yvals.size()-1;
    double xbin = (bins[1]-bins[0])/yvals.size();
    int ibin = min(maxindex,max(0,int(floor((xval-bins[0])/xbin))));
// find a range of N bins about the central bin
    int imin = std::max(0,ibin-halfrange);
    int imax = std::min(maxindex,ibin+halfrange);
    std::vector<double> xfit, yfit;
    for(int jbin=imin;jbin<=imax;++jbin){
      xfit.push_back(bins[0]+xbin*(jbin+0.5)); // y values are at bin center, x vaues are bin edges
      yfit.push_back(yvals[jbin]);
    }
    double c0,c1;
    double cov00, cov01, cov11, sumsq;
    auto fitok = gsl_fit_linear(xfit.data(),1,yfit.data(),1,xfit.size(),
        &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
// test fitok
    if(fitok != 0)throw cet::exception("RECO")<<"mu2e::StrawResponse: calibration interpolation faiulre" << endl;
    value = c0 + c1*xval;
    slope = c1;
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
    if (ddist <= 0)
      return _dc[0] + ddist*calibrateDriftDistanceToT2DDerivative(ddist);

    return sqrt(pow(_dc[1]*ddist+_dc[2],2.0)+_dc[0]*(_dc[0]+2*_dc[2]))-_dc[2];
  }

  double StrawResponse::calibrateT2DToDriftDistance(double t2d) const {
    if (t2d <= _dc[0])
      return (t2d-_dc[0])*calibrateT2DToDriftDistanceDerivative(t2d);

    return (sqrt( pow(t2d+_dc[2],2) - _dc[0]*(_dc[0]+2*_dc[2]))-_dc[2])/_dc[1];
  }

  double StrawResponse::calibrateDriftDistanceToT2DDerivative(double ddist) const {
    if (ddist <= 0){
      if (_dc[2] + _dc[0] == 0)
        return _dc[1];
      else
        return _dc[1]*_dc[2]/(_dc[2]+_dc[0]);
    }

    return _dc[1]*(_dc[1]*ddist+_dc[2])/sqrt(pow(_dc[1]*ddist+_dc[2],2.0)+_dc[0]*(_dc[0]+2*_dc[2]));
  }

  double StrawResponse::calibrateT2DToDriftDistanceDerivative(double t2d) const {
    if (t2d <= _dc[0]){
      if (_dc[2] == 0 && _dc[0] == 0)
        return 1.0/_dc[1];
      else if (_dc[2]== 0)
        return 1.0e8;
      else
        return (_dc[2]+_dc[0])/(_dc[1]*_dc[2]);
    }
    return (_dc[2] + t2d)/sqrt(pow(_dc[1],2)*(-1*pow(_dc[0],2)-2*_dc[0]*_dc[2]+pow(_dc[2]+t2d,2)));
  }

  DriftInfo StrawResponse::driftInfo(StrawId strawId, double dtime, double phi, bool calibrated) const {
    if (_driftIgnorePhi) phi = 0;
    DriftInfo dinfo;
    dinfo.LorentzAngle_ = phi;
    double cdist = _strawDrift->T2D(dtime,phi,false); // allow values outside the physical range at this point
    double derr, derrslope;
    int halfrange(2); // should be a parameter TODO
    interpolateCalib(_driftResBins,_driftResRMS, cdist, halfrange, derr, derrslope);
    dinfo.driftDistanceError_ = derr;
    if(calibrated){
      double dcorr, dcorrslope;
      interpolateCalib(_driftResBins,_driftResOffset, cdist, halfrange, dcorr, dcorrslope);
      dinfo.driftDistance_ = cdist - dcorr;
      dinfo.driftVelocity_ = _strawDrift->GetInstantSpeedFromD(dinfo.driftDistance_);
      //      std::cout << "Drift time " << dtime << " Cluster distance " << cdist
      //         << " offset " << dcorr << " slope " << dcorrslope
      //         << " ddist " << dinfo.driftDistance_ << " dvel " << dinfo.driftVelocity_ << std::endl;
    } else {
      dinfo.driftDistance_ = cdist;
      dinfo.driftVelocity_ = _strawDrift->GetInstantSpeedFromD(dinfo.driftDistance_);
    }
    // force into range
//    dinfo.driftDistance_ = std::min(rstraw_,std::max(0.0,dinfo.driftDistance_));
//    dinfo.driftVelocity_ = std::max(0.0,dinfo.driftVelocity_);
    return dinfo;
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
        ddist = std::max(0.0,std::min(rstraw_,ddist));
        return PieceLineDrift(_driftResBins, _driftResRMS, ddist);
      }else{
        double distance_error = driftDistanceError(strawId,ddist,phi,forceOld);
        double speed_at_ddist = driftInstantSpeed(strawId,ddist,phi,forceOld);
        return distance_error/speed_at_ddist;
      }
    }else{
      //FIXME to be deprecated
      if (useParameterizedDriftError()){
        ddist = std::max(0.0,std::min(rstraw_,ddist));
        return PieceLineDrift(_driftResBins, _driftResRMS, ddist);
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
//      double calibrated_ddist = calibrateDriftDistanceToT2D(ddist);
      return _strawDrift->GetInstantSpeedFromD(ddist)/calibrateDriftDistanceToT2DDerivative(ddist);
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
      double ddist = _strawDrift->T2D(dtime,phi,false);
//      double ddist = calibrateT2DToDriftDistance(t2d_driftonly);
      ddist -= driftDistanceOffset(ddist);
      ddist = std::min(std::max(ddist,0.0),rstraw_); // truncate
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
      //  ddist = std::max(0.0,std::min(rstraw_,ddist));
        return PieceLineDrift(_driftResBins,_driftResRMS, ddist);
      }
    }else{
      // maximum drift is the straw radius.  should come from conditions FIXME!
      ddist = std::min(fabs(ddist),rstraw_);
      size_t idoca = std::min(_derr.size()-1,size_t(floor(_derr.size()*(ddist/rstraw_))));
      return _derr[idoca];
    }
  }

  double StrawResponse::driftDistanceOffset(double ddist) const {
    if(_driftResIsTime)
      return 0;
    else{
      ddist = std::max(0.0,std::min(rstraw_,ddist));
      return PieceLineDrift(_driftResBins,_driftResOffset, ddist);
    }
  }


  // FIXME to be deprecated
  double StrawResponse::driftTimeOffset(StrawId strawId, double ddist, double phi, double DOCA) const {
    return 0;
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
    }else{
      double wslope = PieceLine(_edep,_resslope,kedep);
      tdres += wslope*wlen*wlen;
    }
    // insure a minimum value
    tdres = std::max(30.0,tdres);
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
