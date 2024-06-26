// Bundle together references to various StrawDigi products, to circumvent
// triplicated-loops which pop up when iterating over digis
//
// Ed Callaghan, 2023

#include <Offline/DataProducts/inc/TrkTypes.hh>
#include <Offline/MCDataProducts/inc/DigiProvenance.hh>
#include <Offline/TrackerMC/inc/StrawDigiBundle.hh>

namespace mu2e{
  // dummy MC which signals that a given digi was not produced in simulation
  auto empty_mc = StrawDigiMC(StrawDigiMC(), DigiProvenance::External);
}

namespace mu2e{
  StrawDigiBundle::StrawDigiBundle(const StrawDigi digi,
                                   const StrawDigiADCWaveform adcs,
                                   const StrawDigiMC mc):
                                    _digi(digi),
                                    _adcs(adcs),
                                      _mc(mc){
    /**/
  }

  StrawDigiBundle::StrawDigiBundle(const StrawDigi digi,
                                   const StrawDigiADCWaveform adcs):
                                    _digi(digi),
                                    _adcs(adcs),
                                      _mc(empty_mc){
    /**/
  }

  StrawDigiBundle::StrawDigiBundle(const StrawDigiBundle& bundle):
                                    _digi(bundle.GetStrawDigi()),
                                    _adcs(bundle.GetStrawDigiADCWaveform()),
                                      _mc(bundle.GetStrawDigiMC()){
    /**/
  }

  const StrawDigi& StrawDigiBundle::GetStrawDigi() const{
    const StrawDigi& rv = _digi;
    return rv;
  }

  const StrawDigiADCWaveform& StrawDigiBundle::GetStrawDigiADCWaveform() const{
    const StrawDigiADCWaveform& rv = _adcs;
    return rv;
  }

  const StrawDigiMC& StrawDigiBundle::GetStrawDigiMC() const{
    const StrawDigiMC& rv = _mc;
    return rv;
  }

  // interface for sorting into buckets of overlapping digitization windows
  const double StrawDigiBundle::time() const{
    const auto tdcs = _digi.TDC();
    const auto ptr = std::min_element(tdcs.begin(), tdcs.end());
    const auto first = *ptr;
    const double rv = static_cast<double>(first);
    return rv;
  }
}
