#ifndef Mu2eUtilities_TriggerResultsNavigator_hh
#define Mu2eUtilities_TriggerResultsNavigator_hh
//
//
// $Id: $
// $Author: $
// $Date: $
//
// Original author G. Pezzullo
//

#include "canvas/Persistency/Common/TriggerResults.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include <map>
#include <string>
#include <vector>
#include <ostream>

namespace mu2e {

  class TriggerResultsNavigator{

  public:
    TriggerResultsNavigator(const art::TriggerResults* trigResults) :
      _trigResults(trigResults){
      auto const  id   = trigResults->parameterSetID();
      auto const& pset = fhicl::ParameterSetRegistry::get(id);
      //set the vector<string> with the names of the tirgger_paths
      _trigPathsNames  = pset.get<std::vector<std::string>>("trigger_paths");
      
      //loop over trigResults to fill the map <string, unsigned int) 
      for (unsigned int i=0; i< _trigPathsNames.size(); ++i){
    	_trigMap.insert(std::pair<std::string, unsigned int>(_trigPathsNames[i], i));
      }
    }

    
    // Trigger path information for the current process
    size_t  size() const
    {
      return _trigPathsNames.size();
    }
    
    std::vector<std::string> const& getTrigPaths() const  { return _trigPathsNames; }
    std::string              const& getTrigPath(unsigned int const i) const { return _trigPathsNames.at(i); }

    size_t    findTrigPath(std::string const& name) const;
    size_t    find(std::map<std::string, unsigned int> const& posmap, std::string const& name) const;

    // Has ith path accepted the event?
    bool      accept(std::string const& name) const;
    void      print() const;

  private:
    const art::TriggerResults*           _trigResults;
    std::vector<std::string>             _trigPathsNames;
    std::map<std::string, unsigned int>  _trigMap;
  };
  
  
  

} // namespace mu2e

#endif /* Mu2eUtilities_TriggerResultsNavigator_hh */
