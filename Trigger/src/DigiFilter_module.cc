//
//  Filter for selecting good time cluster: this is part of the track trigger
//  Original author: Dave Brown (LBNL) 3/1/2017
//
// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
// mu2e
// data
#include "RecoDataProducts/inc/StrawDigi.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
// #include "RecoDataProducts/inc/TriggerInfo.hh"
// c++
#include <iostream>
#include <memory>
#include <string> 

using namespace std;

namespace mu2e
{
  class DigiFilter : public art::EDFilter
  {
  public:
    explicit DigiFilter(fhicl::ParameterSet const& pset);
    virtual bool filter(art::Event& event) override;
    virtual bool endRun( art::Run& run ) override;

  private:
    art::InputTag   _sdTag;
    art::InputTag   _cdTag;
    bool            _useSD;   //flag for using the StrawDigi
    bool            _useCD;   //flag for using the CaloDigi
    // std::string     _trigPath;

    //list of the parameters used to perform the filtering
    int             _minnsd;  //minimum number of StrawDigi required
    int             _maxnsd;  //maximum number of StrawDigi required
    int             _minncd;  //minimum number of CaloDigi required
    int             _maxncd;  //maximum number of CaloDigi required
    float           _maxcaloE;//maximum energy, from the sum of all caloDigi

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
    // _trigPath (pset.get<std::string>("triggerPath")),
    _minnsd   (pset.get<int>("minNStrawDigi")),
    _maxnsd   (pset.get<int>("maxNStrawDigi")),
    _minncd   (pset.get<int>("minNCaloDigi")),
    _maxncd   (pset.get<int>("maxNCaloDigi")),
    _maxcaloE (pset.get<float>("maxCaloEnergy")),
    _debug    (pset.get<int>("debugLevel",0)),
    _nevt(0), _npass(0)
  {
    // produces<TriggerInfo>();
  }

  bool DigiFilter::filter(art::Event& event){
    // create output
    // unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool retval(false), retvalSD(false), retvalCD(false); // preset to fail
    // find the collection

    const StrawDigiCollection* sdcol(0);
    const CaloDigiCollection*  cdcol(0);

    int         nsd(0), ncd(0);

    if (_useSD){
      auto sdH = event.getValidHandle<StrawDigiCollection>(_sdTag);
      sdcol = sdH.product();
      nsd   = (int)sdcol->size();
    }
    
    if (_useCD){
      auto cdH = event.getValidHandle<CaloDigiCollection>(_cdTag);
      cdcol = cdH.product();    
      ncd   = (int)cdcol->size();
    }
    
    if (_useSD) {
      if ( (nsd >= _minnsd) && 
	   (nsd <= _maxnsd) ){
	retvalSD = true;
      }
    }

    if (_useCD) {
      if ( (ncd >= _minncd) && 
	   (ncd <= _maxncd) ){
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
      
      // if (retvalSD) triginfo->_triggerBits.merge(TriggerFlag::strawDigis);
      // if (retvalCD) triginfo->_triggerBits.merge(TriggerFlag::caloDigis );
      // triginfo->_triggerPath = _trigPath;

      if(_debug > 1){
	cout << moduleDescription().moduleLabel() << " passed event " << event.id() << endl;
      }
    }
    
    // event.put(std::move(triginfo));
    return retval;
  }

  bool DigiFilter::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      cout << moduleDescription().moduleLabel() << " passed " << _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << endl;
    }
    return true;
  }
}
using mu2e::DigiFilter;
DEFINE_ART_MODULE(DigiFilter);
