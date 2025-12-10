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

    std::vector<std::string> const& getTrigPaths   () const  { return _trigPathsNames; }
    std::string              const& getTrigPath    (unsigned int const i) const { return _trigPathsNames.at(i); }
    std::string              const  getTrigPathName(unsigned int const i) const;
    size_t                          getTrigBit     (unsigned int const pathID) const;
    size_t    findTrigPath(std::string const& name) const;
    size_t    find(std::map<std::string, unsigned int> const& posmap, std::string const& name) const;
    size_t    findTrigPathID(std::string const& name) const;

    bool      validPath(std::string const& name) const {
      size_t path_index = findTrigPath(name);
      return path_index < _trigPathsNames.size();
    }

    size_t    getTrigBit(std::string const& name) const {
      if(!validPath(name))
        throw cet::exception("TRIGGER") << "TriggerResultsNavigator: Path name " <<  name << " not found";
      return findTrigPathID(name);
    }

    std::string const& getTrigNameByBit(size_t const bit) const {
      for (const auto& pair : _trigPathMap) {
        if (pair.second == bit) {
          return pair.first;
        }
      }
      throw cet::exception("TRIGGER") << "TriggerResultsNavigator: Bit " <<  bit << " not found";
    }

    // Has ith path accepted the event?
    bool      accepted(std::string const& name) const;
    bool      accepted(unsigned int const bit) const {
      auto name = getTrigNameByBit(bit);
      return accepted(name);
    }

    bool      wasrun(std::string const& name) const;

    //NOTE: the following three functions can be used only within the same job that runs the
    // trigger paths, otherwise they will fail
    std::vector<std::string>   triggerModules (std::string const& name) const;
    unsigned                   indexLastModule(std::string const& name) const;
    std::string                nameLastModule (std::string const& name) const;

    art::hlt::HLTState state(std::string const& name) const;
    void      print() const;

  private:
    const art::TriggerResults*           _trigResults;
    std::vector<std::string>             _trigPathsNames;  // vector of trigger path names
    std::map<std::string, unsigned int>  _trigMap;         // map of trigger path name to index in TriggerResults
    std::map<std::string, unsigned int>  _trigPathMap;     // map of trigger path name to path ID (bit)
  };




} // namespace mu2e

#endif /* Mu2eUtilities_TriggerResultsNavigator_hh */
