
#include "Offline/TrackerConditions/inc/StrawResponseMaker.hh"
// data products
#include <cmath>
#include <algorithm>
#include <TMath.h>
#include "cetlib_except/exception.h"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/TrackerConditions/inc/StrawDrift.hh"
#include "Offline/GeneralUtilities/inc/SplineInterpolation.hh"

#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "BTrk/BField/BField.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "CLHEP/Matrix/Vector.h"


using namespace std;
namespace mu2e {
  StrawResponse::ptr_t StrawResponseMaker::fromFcl(
      StrawDrift::cptr_t strawDrift,
      StrawElectronics::cptr_t strawElectronics,
      StrawPhysics::cptr_t strawPhysics)     {
    auto thresh = StrawElectronics::thresh;
    auto adc = StrawElectronics::adc;

    // if these value are not defined in fcl, take them
    // from StrawElectronics and StrawPhysics
    double x;
    double electronicsTimeDelay = strawElectronics->electronicsTimeDelay();
    if(_config.electronicsTimeDelay(x)) electronicsTimeDelay = x;
    double gasGain = strawPhysics->strawGain();
    if(_config.gasGain(x)) gasGain = x;
    std::array<double,StrawElectronics::npaths> analognoise;
    analognoise[thresh] = strawElectronics->analogNoise(thresh);
    if(_config.thresholdAnalogNoise(x)) analognoise[thresh] = x;
    analognoise[adc] = strawElectronics->analogNoise(adc);
    if(_config.adcAnalogNoise(x)) analognoise[adc] = x;
    std::array<double,StrawElectronics::npaths> dVdI;
    dVdI[thresh] = strawElectronics->currentToVoltage(StrawId(0,0,0),thresh);
    if(_config.defaultThresholddVdI(x)) dVdI[thresh] = x;
    dVdI[adc] = strawElectronics->currentToVoltage(StrawId(0,0,0),adc);
    if(_config.defaultAdcdVdI(x)) dVdI[adc] = x;
    double vsat = strawElectronics->saturationVoltage();
    if(_config.saturationVoltage(x)) vsat = x;
    double ADCped = strawElectronics->ADCPedestal(StrawId(0,0,0));
    if(_config.ADCPedestal(x)) ADCped = x;

    double pmpEnergyScaleAvg = 0;
    std::array<double, StrawId::_nustraws> pmpEnergyScale;
    if (_config.peakMinusPedestalEnergyScale().size() == 0){
      pmpEnergyScale.fill(_config.defaultPeakMinusPedestalEnergyScale());
    }else{
      for (size_t i=0;i<pmpEnergyScale.size();i++) {
        pmpEnergyScale[i] = _config.peakMinusPedestalEnergyScale()[i];
      }
    }
    for (size_t i=0;i<pmpEnergyScale.size();i++) {
      pmpEnergyScaleAvg += pmpEnergyScale[i];
    }
    pmpEnergyScaleAvg /= (double) pmpEnergyScale.size();

    std::vector<double> _parDriftDocas, _parDriftOffsets, _parDriftRes;

    double sigma = _config.parameterizedDriftSigma();
    double tau = _config.parameterizedDriftTau();
    int parameterizedDriftBins = _config.parameterizedDriftBins();

    _parDriftDocas.reserve(parameterizedDriftBins);
    _parDriftOffsets.reserve(parameterizedDriftBins);
    _parDriftRes.reserve(parameterizedDriftBins);

    for (int i=0;i<parameterizedDriftBins;i++){
      double doca = i*2.5/parameterizedDriftBins;
      double hypotenuse = sqrt(pow(doca,2) + pow(tau*_config.linearDriftVelocity(),2));
      double tau_eff = hypotenuse/_config.linearDriftVelocity() - doca/_config.linearDriftVelocity();

      double sumw = 0;
      double sumwx = 0;
      double sumwx2 = 0;
      double tresid = -20.005;
      for (int it=0;it<10000;it++){
        double weight = exp(sigma*sigma/(2*tau_eff*tau_eff)-tresid/tau_eff)*(1-TMath::Erf((sigma*sigma-tau_eff*tresid)/(sqrt(2)*sigma*tau_eff)));
        sumw += weight;
        sumwx += weight*tresid;
        sumwx2 += weight*tresid*tresid;
        tresid += 0.01;
      }
      double mean = sumwx/sumw;
      double stddev = sqrt(sumwx2/sumw-mean*mean);

      _parDriftDocas.push_back(doca);
      _parDriftOffsets.push_back(mean);
      _parDriftRes.push_back(stddev);
    }


    if (_config.useDriftSplines()){
      if (_config.driftSplineA().size() != (_config.driftSplineBins().size()-1) ||
          _config.driftSplineB().size() != (_config.driftSplineBins().size()-1) ||
          _config.driftSplineC().size() != (_config.driftSplineBins().size()-1) ||
          _config.driftSplineD().size() != (_config.driftSplineBins().size()-1) ||
          _config.driftResSplineA().size() != (_config.driftResSplineBins().size()-1) ||
          _config.driftResSplineB().size() != (_config.driftResSplineBins().size()-1) ||
          _config.driftResSplineC().size() != (_config.driftResSplineBins().size()-1) ||
          _config.driftResSplineD().size() != (_config.driftResSplineBins().size()-1)) {
        throw cet::exception("BADCONFIG")
          << "StrawResponse drift spline vector lengths incorrect" << "\n";
      }
    }
    SplineInterpolation driftspline(_config.driftSplineBins(),_config.driftSplineA(),_config.driftSplineB(),_config.driftSplineC(),_config.driftSplineD(),true);
    SplineInterpolation resspline(_config.driftResSplineBins(),_config.driftResSplineA(),_config.driftResSplineB(),_config.driftResSplineC(),_config.driftResSplineD(),false);
    std::vector<double> d2t_x;
    std::vector<double> d2t_y;
    std::vector<std::vector<double> > d2t_2d(_config.driftSplineNumPhiBins(), std::vector<double> ());
    std::vector<double> t2d_x;
    std::vector<std::vector<double> > t2d_2d(_config.driftSplineNumPhiBins(), std::vector<double> ());
    size_t pieceLineBins = 5001.;
    double tolerance = 1e-5;
    double startingVelocity = 1.0;
    if (_config.splineIsD2T()){
      double deltad = (_config.driftSplineBins().back()-_config.driftSplineBins().front())/(pieceLineBins - 1);
      d2t_x.push_back(-1*deltad);
      d2t_y.push_back(0);
      for (size_t i=0;i<pieceLineBins;i++){
        double dist = deltad*i;
        d2t_x.push_back(dist);
        d2t_y.push_back(driftspline.interpolate(dist));
      }
      startingVelocity = 1.0/driftspline.derivative(_config.driftSplineBins().front());
    }else{
      double deltad = (driftspline.interpolate(_config.driftSplineBins().back())-driftspline.interpolate(_config.driftSplineBins().front()))/(pieceLineBins - 1);
      d2t_x.push_back(-1*deltad);
      d2t_y.push_back(0);
      for (size_t i=0;i<pieceLineBins;i++){
        double dist = deltad*i;
        double t0 = 0;
        double t1 = 100;
        while (fabs(t0-t1) > tolerance){
          double t2 = (t0+t1)/2.;
          double d2 = driftspline.interpolate(t2);
          if (d2 > dist){
            t1 = t2;
          }else{
            t0 =t2;
          }
        }
        d2t_x.push_back(dist);
        d2t_y.push_back((t0+t1)/2.);
      }
      startingVelocity = driftspline.derivative(_config.driftSplineBins().front());
    }

    double deltaPhi = (TMath::Pi()/2.0)/(_config.driftSplineNumPhiBins()-1);
    double min_t = 999;
    double max_t = -999;
    for (int i=0;i<_config.driftSplineNumPhiBins();i++){
      d2t_2d[i].push_back(0);
      double phi = deltaPhi*i;
      for (size_t j=1;j<d2t_x.size();j++){
        d2t_2d[i].push_back(d2t_y[j] + (strawDrift->D2T(d2t_x[j],phi)-strawDrift->D2T(d2t_x[j],0))*_config.driftSplinePhiScaling());
        if (d2t_2d[i][j] > max_t){
          max_t = d2t_2d[i][j];
        }
        if (d2t_2d[i][j] < min_t){
          min_t = d2t_2d[i][j];
        }
      }
    }
    double deltat = (max_t-min_t)/(pieceLineBins-1);
    t2d_x.push_back(min_t - deltat);
    for (size_t i=0;i<pieceLineBins;i++){
      double time = min_t + deltat*i;
      t2d_x.push_back(time);
    }
    for (int i=0;i<_config.driftSplineNumPhiBins();i++){
      t2d_2d[i].push_back(0);
      size_t current_index = 2;
      for (size_t j=0;j<pieceLineBins;j++){
        double time = t2d_x[j+1];
        while (current_index < pieceLineBins-1 && d2t_2d[i][current_index] < time){
          current_index++;
        }
        double d0 = d2t_x[current_index-1];
        double d1 = d2t_x[current_index];
        double t0 = d2t_2d[i][current_index-1];
        double t1 = d2t_2d[i][current_index];
        t2d_2d[i].push_back(d0 + (d1-d0)*(time-t0)/(t1-t0));
      }
    }

    std::vector<std::vector<double> > d2t_2dv(_config.driftSplineNumPhiBins(), std::vector<double> ());
    for (int i=0;i<_config.driftSplineNumPhiBins();i++){
      d2t_2dv[i].push_back(startingVelocity);
      for (size_t j=0;j<pieceLineBins;j++){
        if (j < pieceLineBins-1){
          d2t_2dv[i].push_back((d2t_x[j+1]-d2t_x[j])/(d2t_2d[i][j+1]-d2t_2d[i][j]));
        }else{
          d2t_2dv[i].push_back(d2t_2dv[i][d2t_2dv[i].size()-1]);
        }
      }
    }

    for (int i=0;i<_config.driftSplineNumPhiBins();i++){
      t2d_2d[i][0] = t2d_2d[i][1] - startingVelocity*deltat;
      d2t_2d[i][0] = d2t_2d[i][1] - 1/startingVelocity*(d2t_x[1]-d2t_x[0]);
    }
    std::vector<double> edep;
    for (int i=0;i<_config.eBins();i++)
      edep.push_back(_config.eBinWidth()*i);

    if ((int) _config.ehalfPVScale().size() != _config.eBins() ||
        (int) _config.tdCentralRes().size() != _config.eBins() ||
        (int) _config.tdResSlope().size() != _config.eBins() ||
        (int) _config.totDriftTime().size() != _config.totTBins()*_config.totEBins()){
      throw cet::exception("BADCONFIG")
        << "StrawResponse calibration vector lengths incorrect" << "\n";
    }

    std::array<double, StrawId::_nustraws> strawHalfPropVelocity;
    if (_config.strawHalfPropVelocity().size() == 0){
      strawHalfPropVelocity.fill(_config.defaultHalfPropVelocity());
    }else if ((_config.strawHalfPropVelocity().size() % StrawId::_nstraws) == 0){
      for (size_t i=0;i<StrawId::_nustraws;i++) {
        size_t index = i%_config.strawHalfPropVelocity().size();
        strawHalfPropVelocity[i] = _config.strawHalfPropVelocity()[index];
      }
    }else if (_config.strawHalfPropVelocity().size() == StrawId::_nustraws){
      for (size_t i=0;i<StrawId::_nustraws;i++) {
        strawHalfPropVelocity[i] = _config.strawHalfPropVelocity()[i];
      }
    }else{
      throw cet::exception("BADCONFIG")
        << "StrawResponse calibration vector lengths incorrect" << "\n";
    }



    auto ptr = std::make_shared<StrawResponse>(
        strawDrift,strawElectronics,strawPhysics,
        _config.eBins(), _config.eBinWidth(), edep, _config.ehalfPVScale(),
        _config.centralWirePos(), _config.tdCentralRes(),
        _config.tdResSlope(), _config.truncateLongitudinal(),
        _config.rmsLongErrors(), _config.totTBins(), _config.totTBinWidth(),
        _config.totEBins(), _config.totEBinWidth(), _config.totDriftTime(),
        _config.driftErrorParameters(),
        _config.useParameterizedDriftErrors(), _parDriftDocas, _parDriftOffsets, _parDriftRes,
        _config.wireLengthBuffer(), _config.strawLengthFactor(),
        _config.errorFactor(), _config.useNonLinearDrift(),
        _config.linearDriftVelocity(),
        _config.minT0DOCA(),
        pmpEnergyScale,
        electronicsTimeDelay,
        gasGain, analognoise, dVdI, vsat, ADCped,
        pmpEnergyScaleAvg, strawHalfPropVelocity,
        _config.useDriftSplines(), _config.driftIgnorePhi(),
        resspline,
        deltaPhi,
        d2t_x, d2t_2d, d2t_2dv, t2d_x, t2d_2d);

    std::array<double, StrawId::_nupanels> timeOffsetPanel;
    std::array<double, StrawId::_nustraws> timeOffsetStrawHV, timeOffsetStrawCal;
    if (_config.timeOffsetPanel().size() > 0){
      for (size_t i=0;i<timeOffsetPanel.size();i++)
        timeOffsetPanel[i] = _config.timeOffsetPanel()[i];
    }else{
      for (size_t i=0;i<timeOffsetPanel.size();i++)
        timeOffsetPanel[i] = strawElectronics->getTimeOffsetPanel(i);
    }
    if (_config.timeOffsetStrawHV().size() > 0){
      for(size_t i=0;i<timeOffsetStrawHV.size();i++){
        timeOffsetStrawHV[i] = _config.timeOffsetStrawHV()[i];
        timeOffsetStrawCal[i] = _config.timeOffsetStrawCal()[i];
      }
    }else{
      for (size_t i=0;i<timeOffsetStrawHV.size();i++){
        timeOffsetStrawHV[i] = strawElectronics->getTimeOffsetStrawHV(i);
        timeOffsetStrawCal[i] = strawElectronics->getTimeOffsetStrawCal(i);
      }
    }
    ptr->setOffsets( timeOffsetPanel,
        timeOffsetStrawHV,
        timeOffsetStrawCal );

    return ptr;
  }

}
