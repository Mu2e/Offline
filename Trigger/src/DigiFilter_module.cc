//
//  Filter for selecting good time cluster: this is part of the track trigger
//  Original author: Dave Brown (LBNL) 3/1/2017
//
// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// mu2e
// data
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
// #include "RecoDataProducts/inc/TriggerInfo.hh"
// c++
#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace std;

namespace {
  //Sum over every CaloDigi of the pedestal-subtracted peak sample. The baseline recipe is
  //the one CaloHitMakerFast uses (CaloReco/src/CaloHitMakerFast_module.cc:144-149), so the
  //two modules agree on what a digi's amplitude is. CaloDigi carries no calibrated energy,
  //only the raw waveform, so this sum is in ADC counts rather than MeV.
  float totalCaloSignal(const mu2e::CaloDigiCollection& digis) {
    float total(0.);
    for(const auto& digi : digis) {
      const std::vector<int>& waveform = digi.waveform();
      const int peak = digi.peakpos();
      if(peak < 0 || size_t(peak) >= waveform.size()) continue; //no usable peak sample
      const int nSamPed = std::min<int>(peak > 3 ? 4 : std::max(peak-1, 1), waveform.size());
      float baseline(0.);
      for(int i = 0; i < nSamPed; ++i) baseline += waveform.at(i);
      baseline /= float(nSamPed);
      total += float(waveform.at(peak)) - baseline;
    }
    return total;
  }
}

namespace mu2e
{
  class DigiFilter : public art::EDFilter
  {
  public:
    explicit DigiFilter(fhicl::ParameterSet const& pset);
    virtual bool filter(art::Event& event) override;
    virtual bool endRun( art::Run& ) override;

  private:
    art::InputTag   _sdTag;
    art::InputTag   _cdTag;
    bool            _useSD;   //flag for using the StrawDigi
    bool            _useCD;   //flag for using the CaloDigi

    //list of the parameters used to perform the filtering
    int             _minnsd;  //minimum number of StrawDigi required
    int             _maxnsd;  //maximum number of StrawDigi required
    int             _minncd;  //minimum number of CaloDigi required
    int             _maxncd;  //maximum number of CaloDigi required
    float           _maxcaloE;//maximum summed CaloDigi amplitude; negative disables the cut

    int             _debug;
    // counters
    unsigned _nevt, _npass;
  };

  DigiFilter::DigiFilter(fhicl::ParameterSet const& pset) :
    art::EDFilter{pset},
    _sdTag    (pset.get<art::InputTag>("strawDigiCollection")),
    _cdTag    (pset.get<art::InputTag>("caloDigiCollection")),
    _useSD    (pset.get<bool>("useStrawDigi")),
    _useCD    (pset.get<bool>("useCaloDigi")),
    _minnsd   (pset.get<int>("minNStrawDigi")),
    _maxnsd   (pset.get<int>("maxNStrawDigi")),
    _minncd   (pset.get<int>("minNCaloDigi")),
    _maxncd   (pset.get<int>("maxNCaloDigi")),
    _maxcaloE (pset.get<float>("maxCaloEnergy")),
    _debug    (pset.get<int>("debugLevel",0)),
    _nevt(0), _npass(0)
  {
    // With neither source selected every event is rejected, silently: retval below stays
    // false. Make that a configuration error rather than an empty trigger path.
    if(!_useSD && !_useCD) {
      throw cet::exception("BADCONFIG") << "DigiFilter: neither useStrawDigi nor useCaloDigi is set, "
                                        << "which rejects every event\n";
    }
  }

  bool DigiFilter::filter(art::Event& event){
    ++_nevt;
    bool retval(false), retvalSD(false), retvalCD(false); // preset to fail
    // find the collection

    int         nsd(0), ncd(0);
    float       caloSignal(0.);

    if (_useSD){
      auto sdH = event.getValidHandle<StrawDigiCollection>(_sdTag);
      nsd   = (int)sdH.product()->size();
    }

    if (_useCD){
      auto cdH = event.getValidHandle<CaloDigiCollection>(_cdTag);
      ncd   = (int)cdH.product()->size();
      //only pay for the waveform loop when the cut is enabled
      if (_maxcaloE >= 0.) caloSignal = totalCaloSignal(*cdH.product());
    }

    if (_useSD) {
      if ( (nsd >= _minnsd) &&
           (nsd <= _maxnsd) ){
        retvalSD = true;
      }
    }

    if (_useCD) {
      //a negative maxCaloEnergy means the amplitude cut is switched off, which is what every
      //configuration in mu2e-trig-config currently sets
      if ( (ncd >= _minncd) &&
           (ncd <= _maxncd) &&
           (_maxcaloE < 0. || caloSignal <= _maxcaloE) ){
        retvalCD = true;
      }
    }

    if (_useSD && _useCD) {
      retval = retvalSD && retvalCD;
    }else if (_useSD){
      retval = retvalSD;
    }else if (_useCD){
      retval = retvalCD;
    }

    if (retval){
      ++_npass;

      if(_debug > 1){
        mf::LogInfo("DigiFilter") << moduleDescription().moduleLabel() << " passed event " << event.id();
      }
    }
    return retval;
  }

  bool DigiFilter::endRun( art::Run& ) {
    if(_debug > 0 && _nevt > 0){
      mf::LogInfo("DigiFilter") << moduleDescription().moduleLabel() << " passed " << _npass << " events out of "
                                << _nevt << " for a ratio of " << float(_npass)/float(_nevt);
    }
    _nevt = 0; _npass = 0; //reset for the next run, as MSDHitFilter does
    return true;
  }
}
using mu2e::DigiFilter;
DEFINE_ART_MODULE(DigiFilter)
