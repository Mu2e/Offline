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
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "CLHEP/Matrix/Vector.h"


using namespace std;
namespace mu2e {

  namespace {
    // a per-channel fcl list is either empty (no override) or covers every channel
    void checkChannelListLength(char const* name, size_t size, size_t nchannels) {
      if (size != 0 && size != nchannels)
        throw cet::exception("BADCONFIG")
          << "StrawResponse fcl parameter " << name << " has " << size
          << " entries; it must be empty or have " << nchannels << "\n";
    }
  }

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
    auto const peakMinusPedestalEnergyScale = _config.peakMinusPedestalEnergyScale();
    checkChannelListLength("peakMinusPedestalEnergyScale", peakMinusPedestalEnergyScale.size(), StrawId::_nustraws);
    std::array<double, StrawId::_nustraws> pmpEnergyScale;
    if (peakMinusPedestalEnergyScale.size() == 0){
      pmpEnergyScale.fill(_config.defaultPeakMinusPedestalEnergyScale());
    }else{
      for (size_t i=0;i<pmpEnergyScale.size();i++) {
        pmpEnergyScale[i] = peakMinusPedestalEnergyScale[i];
      }
    }
    for (size_t i=0;i<pmpEnergyScale.size();i++) {
      pmpEnergyScaleAvg += pmpEnergyScale[i];
    }
    pmpEnergyScaleAvg /= (double) pmpEnergyScale.size();

    if ( _config.unsignedDriftRMS().size() != _config.signedDriftRMS().size()
        || _config.driftOffBins().size() != 2
        || _config.driftRMSBins().size() != 2
        || _config.llDriftTimeOffBins().size() != 2
        || _config.llDriftTimeRMSBins().size() != 2
        || _config.llDriftTimeOffset().size() < 2
        || _config.llDriftTimeRMS().size() < 2){
      throw cet::exception("BADCONFIG")
        << "StrawResponse drift res vector lengths incorrect" << "\n";
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
        _config.totDriftError(),
        _config.llDriftTimeOffBins(),_config.llDriftTimeOffset(),
        _config.llDriftTimeRMSBins(),_config.llDriftTimeRMS(),
        _config.driftOffBins(),_config.driftOffset(),
        _config.driftRMSBins(),_config.signedDriftRMS(),
        _config.unsignedDriftRMS(),_config.dRdTScale(),
        _config.wireLengthBuffer(), _config.strawLengthFactor(),
        _config.errorFactor(), _config.useNonLinearDrift(),
        _config.linearDriftVelocity(),
        pmpEnergyScale,
        electronicsTimeDelay,
        gasGain, analognoise, dVdI, vsat, ADCped,
        pmpEnergyScaleAvg, strawHalfPropVelocity,
        _config.driftIgnorePhi());

    auto const timeOffsetPanelFcl = _config.timeOffsetPanel();
    auto const timeOffsetStrawHVFcl = _config.timeOffsetStrawHV();
    auto const timeOffsetStrawCalFcl = _config.timeOffsetStrawCal();
    checkChannelListLength("timeOffsetPanel", timeOffsetPanelFcl.size(), StrawId::_nupanels);
    checkChannelListLength("timeOffsetStrawHV", timeOffsetStrawHVFcl.size(), StrawId::_nustraws);
    if (timeOffsetStrawCalFcl.size() != timeOffsetStrawHVFcl.size())
      throw cet::exception("BADCONFIG")
        << "StrawResponse fcl parameters timeOffsetStrawHV and timeOffsetStrawCal must have the same length, not "
        << timeOffsetStrawHVFcl.size() << " and " << timeOffsetStrawCalFcl.size() << "\n";
    std::array<double, StrawId::_nupanels> timeOffsetPanel;
    std::array<double, StrawId::_nustraws> timeOffsetStrawHV, timeOffsetStrawCal;
    if (timeOffsetPanelFcl.size() > 0){
      for (size_t i=0;i<timeOffsetPanel.size();i++)
        timeOffsetPanel[i] = timeOffsetPanelFcl[i];
    }else{
      for (size_t i=0;i<timeOffsetPanel.size();i++)
        timeOffsetPanel[i] = strawElectronics->getTimeOffsetPanel(i);
    }
    if (timeOffsetStrawHVFcl.size() > 0){
      for(size_t i=0;i<timeOffsetStrawHV.size();i++){
        timeOffsetStrawHV[i] = timeOffsetStrawHVFcl[i];
        timeOffsetStrawCal[i] = timeOffsetStrawCalFcl[i];
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
