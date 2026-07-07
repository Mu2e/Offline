#ifndef Mu2eUtilities_TriggerResultsNavigator_hh
#define Mu2eUtilities_TriggerResultsNavigator_hh
//
// Original author G. Pezzullo
//

#include "canvas/Persistency/Common/TriggerResults.h"

#include <map>
#include <unordered_map>
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

    std::vector<std::string> const& getTrigPaths() const { return _trigPathsNames; }
    std::string              const& getTrigPathByIndex(unsigned int const index) const { return _trigPathsNames.at(index); }
    std::string              const& getTrigPathNameByIndex(unsigned int const i) const;
    std::string              const& getTrigPathNameByBit(unsigned int const bit) const;
    size_t                          getTrigBitByIndex(unsigned int const pathID) const;
    size_t                          getTrigBitByName(const std::string& name) const;
    size_t                          getTrigPathIndex(std::string const& name) const;
    size_t                          getTrigPathIndex(const size_t bit) const;

    bool      validPath(std::string const& name) const {
      return _trigMap.contains(name);
    }
    bool      validPath(const unsigned int bit) {
      return _bitToPathName.contains(bit);
    }

    // Did the path pass the event
    bool      accepted(std::string const& name) const;
    bool      accepted(unsigned int const bit) const;

    // Was the path run in the event
    bool      wasrun(std::string const& name) const;

    //NOTE: the following three functions can be used only within the same job that runs the
    // trigger paths, otherwise they will fail
    std::vector<std::string>   triggerModules (std::string const& name) const;
    unsigned                   indexLastModule(std::string const& name) const;
    std::string                nameLastModule (std::string const& name) const;

    art::hlt::HLTState state(std::string const& name) const;
    void      print() const;

  private:
    const art::TriggerResults*                     _trigResults;     // trigger results info
    std::vector<std::string>                       _trigPathsNames;  // vector of trigger path names
    std::unordered_map<std::string, unsigned int>  _trigMap;         // map of trigger path name to index in the art::TriggerResults
    std::unordered_map<std::string, unsigned int>  _trigPathMap;     // map of trigger path name to path ID (bit)
    std::unordered_map<unsigned int, std::string>  _indexToPathName; // trigger path index -> name map
    std::unordered_map<unsigned int, std::string>  _bitToPathName;   // trigger path bit -> name map
    std::unordered_map<unsigned int, unsigned int> _bitToIndex;      // trigger path bit -> index map

    static constexpr size_t NOTFOUND = -1;                           // special bit/index value for no trigger was found
  };




} // namespace mu2e

#endif /* Mu2eUtilities_TriggerResultsNavigator_hh */
