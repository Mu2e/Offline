
#include "Offline/TrackerConditions/inc/StrawResponseMaker.hh"
// data products
#include <cmath>
#include <algorithm>
#include <TMath.h>
#include "cetlib_except/exception.h"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/TrackerConditions/inc/StrawDrift.hh"
#include "Offline/GeneralUtilities/inc/TwoDimSpline.hh"

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

    double deltaPhi = 0;
    std::vector<TwoDimSpline> timesplines;
    std::vector<TwoDimSpline> ressplines;

    if (_config.useDriftSplines()){
      std::vector<double> driftSplineDoca;
      for (int i=0;i<_config.driftSplineDocaBins();i++){
        driftSplineDoca.push_back(i*2.5/(_config.driftSplineDocaBins()-1));
      }
      deltaPhi = (TMath::Pi()/2.0)/(_config.driftSplinePhiBins()-1);
      if ((int) _config.driftSplineA().size() != (_config.driftSplineDocaBins()-1)*_config.driftSplinePhiBins() ||
          (int) _config.driftSplineB().size() != (_config.driftSplineDocaBins()-1)*_config.driftSplinePhiBins() ||
          (int) _config.driftSplineC().size() != (_config.driftSplineDocaBins()-1)*_config.driftSplinePhiBins() ||
          (int) _config.driftSplineD().size() != (_config.driftSplineDocaBins()-1)*_config.driftSplinePhiBins() ||
          (int) _config.driftResSplineA().size() != (_config.driftSplineDocaBins()-1)*_config.driftSplinePhiBins() ||
          (int) _config.driftResSplineB().size() != (_config.driftSplineDocaBins()-1)*_config.driftSplinePhiBins() ||
          (int) _config.driftResSplineC().size() != (_config.driftSplineDocaBins()-1)*_config.driftSplinePhiBins() ||
          (int) _config.driftResSplineD().size() != (_config.driftSplineDocaBins()-1)*_config.driftSplinePhiBins()){
        throw cet::exception("BADCONFIG")
          << "StrawResponse drift spline vector lengths incorrect" << "\n";
      }

      for (int i=0;i<_config.driftSplinePhiBins();i++){
        std::vector<double> driftSplineA, driftSplineB, driftSplineC, driftSplineD;
        std::vector<double> driftResSplineA, driftResSplineB, driftResSplineC, driftResSplineD;
        for (int j=0;j<_config.driftSplineDocaBins()-1;j++){
          driftSplineA.push_back(_config.driftSplineA()[i*(_config.driftSplineDocaBins()-1) + j]);
          driftSplineB.push_back(_config.driftSplineB()[i*(_config.driftSplineDocaBins()-1) + j]);
          driftSplineC.push_back(_config.driftSplineC()[i*(_config.driftSplineDocaBins()-1) + j]);
          driftSplineD.push_back(_config.driftSplineD()[i*(_config.driftSplineDocaBins()-1) + j]);
          driftResSplineA.push_back(_config.driftResSplineA()[i*(_config.driftSplineDocaBins()-1) + j]);
          driftResSplineB.push_back(_config.driftResSplineB()[i*(_config.driftSplineDocaBins()-1) + j]);
          driftResSplineC.push_back(_config.driftResSplineC()[i*(_config.driftSplineDocaBins()-1) + j]);
          driftResSplineD.push_back(_config.driftResSplineD()[i*(_config.driftSplineDocaBins()-1) + j]);
        }
        timesplines.push_back(TwoDimSpline(driftSplineDoca,driftSplineA,driftSplineB,driftSplineC,driftSplineD,true));
        ressplines.push_back(TwoDimSpline(driftSplineDoca,driftResSplineA,driftResSplineB,driftResSplineC,driftResSplineD,false));
      }
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
        deltaPhi,timesplines,ressplines);

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
