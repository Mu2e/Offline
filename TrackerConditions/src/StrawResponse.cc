//
// Reconstruction functions for straws
// Original author David Brown, LBNL
//
#include "TrackerConditions/inc/StrawResponse.hh"
// data products
#include "RecoDataProducts/inc/StrawHit.hh"
#include <math.h>
#include <algorithm>

#include "TrackerConditions/inc/StrawDrift.hh"

#include "BFieldGeom/inc/BFieldManager.hh"
#include "BTrk/BField/BField.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CLHEP/Matrix/Vector.h"


using namespace std;
namespace mu2e {
  StrawResponse::StrawResponse(fhicl::ParameterSet const& pset) :
    _strawDrift(new StrawDrift()),
    _edep(pset.get<vector<float> >("EDep",
	  vector<float>{0.214044 , 0.311385 , 0.405932 , 0.503041 , 0.601217 , 0.700171 , 0.799435 , 0.898859 , 0.998941 , 1.09876 , 1.19853 , 1.29892 , 1.39884 , 1.49879 , 1.59901 , 1.69917 , 1.799 , 1.8992 , 1.99935 , 2.09929 , 2.19977 , 2.29964 , 2.39963 , 2.49954 , 2.59964 , 2.70024 , 2.79984 , 2.89971 , 3.00015 , 3.10025 , 3.20018 , 3.3004 , 3.40047 , 3.50056 , 3.60085 , 3.70057 , 3.80023 , 3.89995 , 3.99972 , 4.1001 , 4.19994 , 4.29969 , 4.39939 , 4.49958 , 4.59949 , 4.70012 , 4.79972 , 4.89955 , 5.00006 , 5.09999 , 5.19956 , 5.29984 , 5.39941 , 5.49948 , 5.59953 , 5.69931 , 5.79903 , 5.89868 , 5.9986 , 6.09881 , 6.19851 , 6.29715 , 6.39709 , 6.49732 , 6.5966 , 6.69613 , 6.79413 , 7.11646 })), //KeV
    _halfvp(pset.get<vector<float> >("HalfPropVelocity",
	  vector<float>{78.2426 , 77.959 , 79.6768 , 81.9859 , 84.0565 , 85.7715 , 87.4105 , 88.9116 , 90.3429 , 91.6604 , 92.789 , 93.8731 , 94.9901 , 95.8371 , 96.7258 , 97.7119 , 98.4066 , 99.1643 , 100.082 , 100.903 , 101.591 , 102.627 , 103.56 , 104.628 , 105.924 , 106.967 , 108.189 , 108.931 , 109.836 , 110.911 , 112.013 , 113.362 , 114.558 , 116.156 , 117.553 , 119.312 , 120.478 , 121.325 , 122.056 , 122.498 , 122.872 , 123.221 , 123.242 , 123.056 , 122.908 , 122.972 , 122.926 , 123.092 , 122.91 , 123.313 , 123.347 , 123.571 , 123.815 , 124.01 , 124.2 , 124.678 , 124.802 , 124.926 , 125.056 , 125.349 , 125.687 , 126.092 , 125.699 , 125.636 , 125.562 , 125.487 , 123.53 , 118.477 })), // mm/ns 
    _central(pset.get<float>("CentralWirePos",65.0)), // mm
    _centres(pset.get<vector<float> >("TDCentralRes",
	  vector<float>{76.2103 , 66.2778 , 58.3431 , 53.2878 , 49.9309 , 47.1836 , 45.13 , 42.9479 , 41.208 , 39.6033 , 38.1616 , 37.1039 , 35.9969 , 35.0815 , 34.0823 , 33.4877 , 32.9431 , 32.3805 , 32.071 , 31.4614 , 31.2776 , 30.635 , 30.1217 , 29.7487 , 29.1309 , 28.4732 , 28.1865 , 27.6331 , 27.5013 , 27.0485 , 26.4937 , 25.9912 , 25.5595 , 25.0649 , 24.0109 , 23.1914 , 22.2386 , 21.7544 , 21.3721 , 21.0055 , 20.5813 , 20.1229 , 20.2912 , 20.615 , 20.4152 , 20.6086 , 20.6452 , 20.6372 , 20.6297 , 20.4731 , 20.401 , 20.174 , 20.2679 , 20.3807 , 20.1217 , 19.7474 , 19.8524 , 19.8093 , 19.8447 , 19.8788 , 19.5326 , 19.1781 , 19.7947 , 19.3696 , 18.6408 , 18.0923 , 18.1123 , 21.371 })),
    _resslope(pset.get<vector<float> >("TDResSlope",
	  vector<float>{0.0160651 , 0.023741 , 0.0421304 , 0.0537863 , 0.0633849 , 0.0697326 , 0.0725653 , 0.0736391 , 0.0746531 , 0.0745937 , 0.0750631 , 0.0738193 , 0.0740473 , 0.0761049 , 0.0785487 , 0.0783369 , 0.0802156 , 0.081579 , 0.0815734 , 0.0830998 , 0.0824664 , 0.0847393 , 0.0857842 , 0.0846082 , 0.0866494 , 0.0879027 , 0.0883831 , 0.0917368 , 0.0909706 , 0.0899623 , 0.0899677 , 0.0890875 , 0.0888003 , 0.0865446 , 0.0847627 , 0.0794976 , 0.0802735 , 0.0797023 , 0.0798181 , 0.0750883 , 0.073853 , 0.0751729 , 0.0737468 , 0.071133 , 0.0698957 , 0.0661819 , 0.0642328 , 0.0635607 , 0.0623173 , 0.0606496 , 0.0574217 , 0.0562718 , 0.0532419 , 0.0510767 , 0.0493828 , 0.0472257 , 0.0447443 , 0.0425958 , 0.0419509 , 0.0391222 , 0.0372135 , 0.0399352 , 0.0340117 , 0.0375563 , 0.0368355 , 0.052548 , 0.0649657 , 0.168875 })),
    _usederr(pset.get<bool>("UseDriftErrorCalibration",false)),
    _derr(pset.get<vector<float> >("DriftErrorParameters",vector<float>{0.2631,-0.05315, 2.55, 0.1978, 0.5964})),
    _wbuf(pset.get<float>("WireLengthBuffer",2.0)), //sigma
    _slfac(pset.get<float>("StrawLengthFactor",0.9)),
    _errfac(pset.get<float>("ErrorFactor",1.0)),
    _driftFile(pset.get<string>("DriftFile","TrackerConditions/data/E2v.tbl")),
    _wirevoltage(pset.get<double>("WireVoltage",1400)),
    _phiBins(pset.get<int>("DriftPhiBins",20)),
    _dIntegrationBins(pset.get<int>("DriftIntegrationBins",50)),
    _usenonlindrift(pset.get<bool>("UseNonLinearDrift",false)),
    _lindriftvel(pset.get<double>("LinearDriftVelocity",0.0625)), // mm/ns, only used if nonlindrift==0
    _rres_min(pset.get<double>("MinDriftRadiusResolution",0.2)), //mm
    _rres_max(pset.get<double>("MaxDriftRadiusResolution",0.2)), //mm
    _rres_rad(pset.get<double>("DriftRadiusResolutionRadius",-1)), //mm
    _mint0doca(pset.get<double>("minT0DOCA", -0.2)), //FIXME should be moved to a reconstruction configuration 
    _pmpEnergyScale(pset.get<vector<double> >("peakMinusPedestalEnergyScale",vector<double>(96,0.0042))), // fudge factor for peak minus pedestal energy method
    _timeOffsetPanel(pset.get<vector<double> >("TimeOffsetPanel",vector<double>(240,0))),
    _timeOffsetStrawHV(pset.get<vector<double> >("TimeOffsetStrawHV",vector<double>(96,0))),
    _timeOffsetStrawCal(pset.get<vector<double> >("TimeOffsetStrawHV",vector<double>(96,0)))
    {
      _strawele = ConditionsHandle<StrawElectronics>("ignored");
      _strawphys = ConditionsHandle<StrawPhysics>("ignored");

      _timeOffsetBeam = pset.get<double>("TimeOffsetBeam",_strawele->clockStart());
      _gasGain = pset.get<double>("GasGain",_strawphys->strawGain());
      _analognoise[TrkTypes::thresh] = pset.get<double>("thresholdAnalogNoise",_strawele->analogNoise(TrkTypes::thresh));
      _analognoise[TrkTypes::adc] = pset.get<double>("adcAnalogNoise",_strawele->analogNoise(TrkTypes::adc));
      _dVdI[TrkTypes::thresh] = pset.get<double>("thresholddVdI",_strawele->currentToVoltage(StrawId(0,0,0),TrkTypes::thresh));
      _dVdI[TrkTypes::adc] = pset.get<double>("adcdVdI",_strawele->currentToVoltage(StrawId(0,0,0),TrkTypes::adc));
      _vsat = pset.get<double>("SaturationVoltage",_strawele->saturationVoltage()); // mVolt
      _ADCped = pset.get<unsigned>("ADCPedestal",_strawele->ADCPedestal(StrawId(0,0,0)));


    _pmpEnergyScaleAvg = 0;
      for (size_t i=0;i<_pmpEnergyScale.size();i++)
        _pmpEnergyScaleAvg += _pmpEnergyScale[i];
      _pmpEnergyScaleAvg /= (double) _pmpEnergyScale.size();

    }

  StrawResponse::~StrawResponse(){}

  void StrawResponse::initializeStrawDrift() const
  {
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0,0.0,0.0));
    CLHEP::Hep3Vector b0 = bfmgr->getBField(vpoint_mu2e);
    float Bz = b0.z();
    _strawDrift->Initialize(_driftFile, _wirevoltage, _phiBins, _dIntegrationBins, Bz);
  }


  // simple line interpolation, this should be a utility function, FIXME!
  float StrawResponse::PieceLine(std::vector<float> const& xvals, std::vector<float> const& yvals, float xval){
    float yval;
    if(xvals.size() != yvals.size() || xvals.size() < 2)
      std::cout << "size error " << std::endl;
    int imax = int(xvals.size()-1);
    // approximate constant binning to get initial guess
    float xbin = (xvals.back()-xvals.front())/(xvals.size()-1);
    int ibin = min(imax,max(0,int(floor((xval-xvals.front())/xbin))));
    // scan to exact bin
    while(ibin > 0 && xval < xvals[ibin])
      --ibin;
    while(ibin < imax && xval > xvals[ibin])
      ++ibin;
    // interpolate
    float slope(0.0);
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

  double StrawResponse::driftDistanceToTime(StrawId strawId, double ddist, double phi) const {
    if(_usenonlindrift){
      if (!_strawDrift->isInitialized())
        initializeStrawDrift();
      return _strawDrift->D2T(ddist,phi);
    }
    else{
      return ddist/_lindriftvel; //or return t assuming a constant drift speed of 0.06 mm/ns (for diagnosis)
    }
  }

  double StrawResponse::driftTimeToDistance(StrawId strawId, double dtime, double phi) const {
    if(_usenonlindrift){
      if (!_strawDrift->isInitialized())
        initializeStrawDrift();
      return _strawDrift->T2D(dtime,phi);
    }
    else{
      return dtime*_lindriftvel; //or return t assuming a constant drift speed of 0.06 mm/ns (for diagnosis)
    }
  }

  double StrawResponse::driftInstantSpeed(StrawId strawId, double doca, double phi) const {
    if(_usenonlindrift){
      if (!_strawDrift->isInitialized())
        initializeStrawDrift();
      return _strawDrift->GetInstantSpeedFromD(doca);
    }else{
      return _lindriftvel;
    }
  }

  double StrawResponse::driftDistanceError(StrawId strawId, double ddist, double phi, float DOCA) const {
    if(useDriftError()){
      // maximum drift is the straw radius.  should come from conditions FIXME!
      static float rstraw(2.5);
      DOCA = std::min(DOCA,rstraw);
      // drift errors are modeled as a truncated line + exponential
      return std::min(_derr[4],_derr[0]+_derr[1]*DOCA + _derr[2]*exp(-DOCA/_derr[3]));
    }else{
      double rres = _rres_min;
      if(ddist < _rres_rad){
        rres = _rres_min+_rres_max*(_rres_rad-ddist)/_rres_rad;
      }
      return rres;
    }
  }

  double StrawResponse::driftDistanceOffset(StrawId strawId, double ddist, double phi, float DOCA) const {
    // fixme
    return 0;
  }

  bool StrawResponse::wireDistance(StrawHit const& strawhit, float slen, float& wdist, float& wderr) const {
    bool retval(true);
    // convert edep from Mev to KeV (should be standardized, FIXME!)
    float kedep = 1000.0*strawhit.energyDep();
    wdist = halfPropV(strawhit.strawId(),kedep)*(strawhit.dt());
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

  float StrawResponse::halfPropV(StrawId strawId, float kedep) const {
    return PieceLine(_edep,_halfvp,kedep);
  }

  float StrawResponse::wpRes(float kedep,float wlen) const {
  // central resolution depends on edep
    float tdres = PieceLine(_edep,_centres,kedep);
    if( wlen > _central){
    // outside the central region the resolution depends linearly on the distance
    // along the wire.  The slope of that also depends on edep
      float wslope = PieceLine(_edep,_resslope,kedep);
      tdres += (wlen-_central)*wslope;
    }
    return tdres;
  }

  void StrawResponse::calibrateTimes(TrkTypes::TDCValues const& tdc, TrkTypes::TDCTimes &times, const StrawId &id) const {
    
    times[TrkTypes::hv] = tdc[TrkTypes::hv]*_strawele->tdcLSB() +  _timeOffsetBeam + _timeOffsetPanel[id.getPanel()] + _timeOffsetStrawHV[id.getStraw()];
    times[TrkTypes::cal] = tdc[TrkTypes::cal]*_strawele->tdcLSB() + _timeOffsetBeam + _timeOffsetPanel[id.getPanel()] + _timeOffsetStrawCal[id.getStraw()];
  }
 

  double StrawResponse::saturatedResponse(double vlin) const {
    if (vlin < _vsat)
      return vlin;
    else
      return _vsat;
  }
 
}
