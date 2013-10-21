// An EDAnalyzer template to print out a data collection
//
// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Analyses_GenericCollectionPrinter_hh
#define ExtinctionMonitorFNAL_Analyses_GenericCollectionPrinter_hh

#include <iostream>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"

namespace mu2e {
  //================================================================
  template<class Coll>
  class GenericCollectionPrinter : public art::EDAnalyzer {
  protected:
    std::string _inModuleLabel;
    std::string _inInstanceName;

  public:
    explicit GenericCollectionPrinter(const fhicl::ParameterSet& pset);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  template<class Coll>
  GenericCollectionPrinter<Coll>::GenericCollectionPrinter(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , _inModuleLabel(pset.get<std::string>("inputModuleLabel"))
    , _inInstanceName(pset.get<std::string>("inputInstanceName"))
  {}

  //================================================================
  template<class Coll>
  void GenericCollectionPrinter<Coll>::analyze(const art::Event& event) {

    art::Handle<Coll> ih;
    event.getByLabel(_inModuleLabel, _inInstanceName, ih);

    const Coll& inputs(*ih);

    std::cout<<"GenericCollectionPrinter: inModuleLabel = "<<_inModuleLabel
             <<", inInstanceName = "<<_inInstanceName<<std::endl;

    for(typename Coll::const_iterator i=inputs.begin(); i!=inputs.end(); ++i) {
      std::cout<<"event "<<event.id()<<", entry "<<*i<<std::endl;
    }
  }

  //================================================================
} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Analyses_GenericCollectionPrinter_hh*/
