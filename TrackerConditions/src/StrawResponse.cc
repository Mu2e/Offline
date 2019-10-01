
#include "TrackerConditions/inc/StrawResponse.hh"
#include <math.h>
#include <algorithm>

using namespace std;

namespace mu2e {


  // simple line interpolation, this should be a utility function, FIXME!
  double StrawResponse::PieceLine(std::vector<double> const& xvals, std::vector<double> const& yvals, double xval){
    double yval;
    if(xvals.size() != yvals.size() || xvals.size() < 2)
      std::cout << "size error " << std::endl;
    int imax = int(xvals.size()-1);
    // approximate constant binning to get initial guess
    double xbin = (xvals.back()-xvals.front())/(xvals.size()-1);
    int ibin = min(imax,max(0,int(floor((xval-xvals.front())/xbin))));
    // scan to exact bin
    while(ibin > 0 && xval < xvals[ibin])
      --ibin;
    while(ibin < imax && xval > xvals[ibin])
      ++ibin;
    // interpolate
    double slope(0.0);
    if(ibin >= 0 && ibin < imax){
      yval = yvals[ibin];
      int jbin = ibin+1;
      slope = (yvals[jbin]-yvals[ibin])/(xvals[jbin]-xvals[ibin]);
      yval += (xval-xvals[ibin])*slope;
    } else if(ibin >= imax){ 
      yval = yvals[imax];
      slope = (yvals[imax]-yvals[imax-1])/(xvals[imax]-xvals[imax-1]);
      yval += (xval-xvals[imax])*slope;
    } else {
      yval = yvals[0];
      slope = (yvals[1]-yvals[0])/(xvals[1]-xvals[0]);
      yval += (xval-xvals[0])*slope;
    }
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
    // fixme
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
      + _timeOffsetStrawHV[id.getStraw()];
    times[StrawEnd::cal] = tdc[StrawEnd::cal]*_strawElectronics->tdcLSB() 
      - electronicsTimeDelay + _timeOffsetPanel[id.getPanel()] 
      + _timeOffsetStrawCal[id.getStraw()];
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
    size_t totbin = (size_t) (tot/4.);
    size_t ebin = (size_t) (edep/0.00025);
    if (totbin > 15)
      totbin = 15;
    if (ebin > 9)
      ebin = 9;
    return _totdtime[totbin*10+ebin];
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
    printVector(os,"pmpEnergyScale",_pmpEnergyScale);
    printVector(os,"timeOffsetPanel",_timeOffsetPanel);
    printVector(os,"timeOffsetStrawHV",_timeOffsetStrawHV);
    printVector(os,"timeOffsetStrawCal",_timeOffsetStrawCal);
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

}
