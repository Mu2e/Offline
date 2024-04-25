// Bundle together references to various StrawDigi products, to circumvent
// triplicated-loops which pop up when iterating over digis
//
// Ed Callaghan, 2023

#include <Offline/DataProducts/inc/TrkTypes.hh>
#include <Offline/TrackerMC/inc/StrawDigiBundle.hh>

namespace mu2e{
  StrawDigiMC empty_mc;
}

namespace mu2e{
  StrawDigiBundle::StrawDigiBundle(const StrawDigi digi,
                                   const StrawDigiADCWaveform adcs,
                                   const StrawDigiMC mc):
                                    digi(digi),
                                    adcs(adcs),
                                      mc(mc){
    /**/
  }

  StrawDigiBundle::StrawDigiBundle(const StrawDigi digi,
                                   const StrawDigiADCWaveform adcs):
                                    digi(digi),
                                    adcs(adcs),
                                      mc(empty_mc){
    /**/
  }

  StrawDigiBundle::StrawDigiBundle(const StrawDigiBundle& bundle):
                                    digi(bundle.GetStrawDigi()),
                                    adcs(bundle.GetStrawDigiADCWaveform()),
                                      mc(bundle.GetStrawDigiMC()){
    /**/
  }

  const StrawDigi StrawDigiBundle::GetStrawDigi() const{
    const StrawDigi rv = this->digi;
    return rv;
  }

  const StrawDigiADCWaveform StrawDigiBundle::GetStrawDigiADCWaveform() const{
    const StrawDigiADCWaveform rv = this->adcs;
    return rv;
  }

  const StrawDigiMC StrawDigiBundle::GetStrawDigiMC() const{
    const StrawDigiMC rv = this->mc;
    return rv;
  }

  // interface for sorting into buckets of overlapping digitization windows
  const double StrawDigiBundle::time() const{
    const auto tdcs = this->digi.TDC();
    const auto ptr = std::min_element(tdcs.begin(), tdcs.end());
    const auto first = *ptr;
    const double rv = static_cast<double>(first);
    return rv;
  }
}
