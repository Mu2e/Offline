#ifndef Mu2eUtilities_TriggerResultsNavigator_hh
#define Mu2eUtilities_TriggerResultsNavigator_hh
//
// Original author G. Pezzullo
//

#include "canvas/Persistency/Common/TriggerResults.h"

#include <map>
#include <string>
#include <vector>
#include <ostream>

namespace mu2e {

  class TriggerResultsNavigator{

  public:
    TriggerResultsNavigator(const art::TriggerResults* trigResults);
    
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
    bool      accepted(std::string const& name) const;
    art::hlt::HLTState state(std::string const& name) const;
    void      print() const;

  private:
    const art::TriggerResults*           _trigResults;
    std::vector<std::string>             _trigPathsNames;
    std::map<std::string, unsigned int>  _trigMap;
  };
  
  
  

} // namespace mu2e

#endif /* Mu2eUtilities_TriggerResultsNavigator_hh */
