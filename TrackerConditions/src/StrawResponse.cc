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
      int imax = yvals.size()-2;
    double xbin = (bins[1]-bins[0])/yvals.size();
    int ibin = min(imax,max(0,int(floor((xval-bins[0])/xbin))));
    double yval = yvals[ibin];
    auto jbin = ibin+1;
    double slope = (yvals[jbin]-yvals[ibin])/xbin;
    yval += (xval-(bins[0]+xbin*ibin))*slope;
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

  DriftInfo StrawResponse::driftInfo(StrawId strawId, double dtime, double phi) const {
    if (_driftIgnorePhi) phi = 0;
    DriftInfo dinfo;
    dinfo.LorentzAngle_ = phi;
    dinfo.cDrift_ = _strawDrift->T2D(dtime,phi,false); // allow values outside the physical range at this point
    int halfrange(2); // should be a parameter TODO
    double dcorr, dcorrslope;
    interpolateCalib(_driftOffBins,_driftOffset, dinfo.cDrift_, halfrange, dcorr, dcorrslope);
    dinfo.rDrift_ = dinfo.cDrift_ -dcorr;
    // note 'Velocity' is really dR/dt (change in calibrated drift distance WRT measured time), not a true physical velocity
    dinfo.driftVelocity_ = _strawDrift->GetInstantSpeedFromD(dinfo.cDrift_)*(1.0 - dcorrslope)*_dRdTScale;
    double serrslope,uerrslope;
    interpolateCalib(_driftRMSBins,_signedDriftRMS, dinfo.rDrift_, halfrange, dinfo.signedDriftError_, serrslope);
    interpolateCalib(_driftRMSBins,_unsignedDriftRMS, dinfo.rDrift_, halfrange, dinfo.unsignedDriftError_ , uerrslope);
    return dinfo;
  }

  double StrawResponse::driftDistanceToTime(StrawId strawId, double ddist, double phi) const {
    if (_driftIgnorePhi)
      phi = 0;
    if(_usenonlindrift){
      return  _strawDrift->D2T(ddist,phi);
    }else{
      return ddist/_lindriftvel; //or return t assuming a constant drift speed of 0.06 mm/ns (for diagnosis)
    }
  }

  double StrawResponse::driftTimeOffset(StrawId strawId, double ddist, double phi) const {
    return PieceLineDrift(_llDriftTimeOffBins,_llDriftTimeOffset, ddist);
  }

  double StrawResponse::driftTimeError(StrawId strawId, double ddist, double phi) const {
    ddist = std::max(0.0,std::min(rstraw_,ddist));
    return PieceLineDrift(_llDriftTimeRMSBins, _llDriftTimeRMS, ddist);
  }

  double StrawResponse::driftInstantSpeed(StrawId strawId, double ddist, double) const {
    if(_usenonlindrift){
      return _strawDrift->GetInstantSpeedFromD(ddist);
    }else{
      return _lindriftvel;
    }
  }

  double StrawResponse::driftTimeToDistance(StrawId strawId, double dtime, double phi) const {
    if (_driftIgnorePhi)
      phi = 0;
    if(_usenonlindrift){
      return _strawDrift->T2D(dtime,phi,false);
    }
    else{
      return dtime*_lindriftvel; //or return t assuming a constant drift speed of 0.06 mm/ns (for diagnosis)
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
      - electronicsTimeDelay + _timeOffsetPanel[id.uniquePanel()]
      + _timeOffsetStrawHV[id.uniqueStraw()];
    times[StrawEnd::cal] = tdc[StrawEnd::cal]*_strawElectronics->tdcLSB()
      - electronicsTimeDelay + _timeOffsetPanel[id.uniquePanel()]
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
    return _totdtime[totbin*_totEBins+ebin] - 0.1; // temporary kludge
  }

  double StrawResponse::TOTdriftTimeError(Straw const& straw,
      double tot, double edep) const {
    size_t totbin = min(_totTBins-1,static_cast<size_t>(tot/_totTBinWidth));
    size_t ebin = min(_totEBins-1,static_cast<size_t>(edep/_totEBinWidth));
    return _totderror[totbin*_totEBins+ebin];
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

   double StrawResponse::BTrk_driftDistanceToTime(StrawId strawId, double ddist, double phi) const {
     return _strawDrift->D2T(ddist,phi);
  }

   double StrawResponse::BTrk_driftInstantSpeed(StrawId strawId, double ddist, double phi) const {
     return _strawDrift->GetInstantSpeedFromD(ddist);
   }

   double StrawResponse::BTrk_driftDistanceError(StrawId strawId, double ddist, double phi) const {
     static const std::vector<double> derr = {0.363559, 0.362685, 0.359959, 0.349385,
       0.336731, 0.321784, 0.302363, 0.282691, 0.268223, 0.252673, 0.238557,
       0.229172, 0.2224, 0.219224, 0.217334, 0.212797, 0.210303, 0.209876,
       0.208739, 0.207411, 0.208738, 0.209646, 0.210073, 0.207101, 0.20431,
       0.203994, 0.202931, 0.19953, 0.196999, 0.194559, 0.191766, 0.187725,
       0.185959, 0.181423, 0.17848, 0.171357, 0.171519, 0.168422, 0.161338,
       0.156641, 0.151196, 0.146546, 0.144069, 0.139858, 0.135838, 0.13319,
       0.132159, 0.130062, 0.123545, 0.120212 };
     static double rstraw(2.5);
     double doca = std::min(fabs(ddist),rstraw);
     size_t idoca = std::min(derr.size()-1,size_t(floor(derr.size()*(doca/rstraw))));
     return derr[idoca];
   }

   double StrawResponse::BTrk_driftTimeToDistance(StrawId strawId, double dtime, double phi) const {
     return _strawDrift->T2D(dtime,phi);
   }

}
