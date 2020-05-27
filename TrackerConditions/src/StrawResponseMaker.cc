
#include "TrackerConditions/inc/StrawResponseMaker.hh"
// data products
#include <cmath>
#include <algorithm>
#include <TH1D.h>
#include <TMath.h>
#include "DataProducts/inc/StrawId.hh"
#include "TrackerConditions/inc/StrawDrift.hh"

#include "BFieldGeom/inc/BFieldManager.hh"
#include "BTrk/BField/BField.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
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
    auto pmpEnergyScale = _config.peakMinusPedestalEnergyScale();
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
    
    TH1D h("","",10000,-20,80);
    for (int i=0;i<parameterizedDriftBins;i++){
      double doca = i*2.5/parameterizedDriftBins;
      double hypotenuse = sqrt(pow(doca,2) + pow(tau*_config.linearDriftVelocity(),2));
      double tau_eff = hypotenuse/_config.linearDriftVelocity() - doca/_config.linearDriftVelocity();

      for (int it=0;it<h.GetNbinsX();it++){
        double tresid = h.GetBinCenter(it+1);
        h.SetBinContent(it+1,exp(sigma*sigma/(2*tau_eff*tau_eff)-tresid/tau_eff)*(1-TMath::Erf((sigma*sigma-tau_eff*tresid)/(sqrt(2)*sigma*tau_eff))));
      }
      h.Scale(1.0/h.Integral());

      _parDriftDocas.push_back(doca);
      _parDriftOffsets.push_back(h.GetMean());
      _parDriftRes.push_back(h.GetRMS());
      h.Reset();
    }
    
    auto ptr = std::make_shared<StrawResponse>(
	 strawDrift,strawElectronics,strawPhysics,
	 _config.eDep(), _config.halfPropVelocity(), 
	 _config.centralWirePos(), _config.tdCentralRes(), 
	 _config.tdResSlope(), _config.totDriftTime(), 
	 _config.useDriftErrorCalibration(), _config.driftErrorParameters(), 
         _config.useParameterizedDriftErrors(), _parDriftDocas, _parDriftOffsets, _parDriftRes,
	 _config.wireLengthBuffer(), _config.strawLengthFactor(), 
	 _config.errorFactor(), _config.useNonLinearDrift(), 
	 _config.linearDriftVelocity(), _config.minDriftRadiusResolution(), 
	 _config.maxDriftRadiusResolution(), 
	 _config.driftRadiusResolutionRadius(), _config.minT0DOCA(), 
	 _config.t0shift(), _config.peakMinusPedestalEnergyScale(), 
	 _config.timeOffsetPanel(), _config.timeOffsetStrawHV(), 
	 _config.timeOffsetStrawCal(), electronicsTimeDelay, 
	  gasGain, analognoise, dVdI, vsat, ADCped,
	 pmpEnergyScaleAvg );

    return ptr;
  }
     
}
