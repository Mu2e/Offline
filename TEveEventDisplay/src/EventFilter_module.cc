#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"

#include<TObject.h>
#include <iostream>

using namespace std;
namespace mu2e {

	class EventFilter : public art::EDFilter {
	  public:
	    struct Config {
	      using Name=fhicl::Name;
	      using Comment=fhicl::Comment;
	      fhicl::Atom<int> event{Name("Event"),Comment("event number"),0};
	      fhicl::Atom<int> run{Name("Run"),Comment("run number"),0};
	      fhicl::Atom<bool> Sequential{Name("Sequential"),Comment("sequential events"),false};
	    };
	    
	    typedef art::EDFilter::Table<Config> Parameters;
	    
	    explicit EventFilter(const Parameters& conf);
	    virtual bool filter (art::Event& event) override;
	    virtual ~EventFilter() {}
	  
	  private:
	    Config _conf;
	    int _event, eventid;
	    int _run, runid;
	    bool _seq;
	};
	
	EventFilter::EventFilter(const Parameters& conf):
    art::EDFilter{conf},
    _event(conf().event()),
    _run(conf().run()),
    _seq(conf().Sequential()){}

	
	bool EventFilter::filter(art::Event& event) {
	  if(_seq) return true;
	  eventid = event.id().event();
	  runid = event.run();
	  std::cout<<"Enter Event Number Selection "<<std::endl;
	  std::cin>>_event;
	  std::cout<<"Enter Run Number Selection "<<std::endl;
	  std::cin>>_run;
	  if(eventid == _event and runid==_run) {
	    return true;
	  }
	    return false;
  }
}
using mu2e::EventFilter;
DEFINE_ART_MODULE(EventFilter);
