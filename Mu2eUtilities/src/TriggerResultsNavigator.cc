#include <exception>
#include <stddef.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "canvas/Persistency/Common/HLTenums.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "fhiclcpp/coding.h"
#include "fhiclcpp/exception.h"

#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"

namespace mu2e {

  TriggerResultsNavigator::TriggerResultsNavigator(const art::TriggerResults* trigResults):
    _trigResults(trigResults){
    auto const  id   = trigResults->parameterSetID();
    // Simplest thing: ParameterSetRegistry has correct ParameterSet
    fhicl::ParameterSet pset;
    fhicl::ParameterSetRegistry::get(id, pset);

    // Go through ParameterSetRegistry and check for corrupted IDs
    auto const& psets = fhicl::ParameterSetRegistry::get();
    for(auto const& pset_pair : psets) {
      auto nid = pset_pair.second.id();
      assert(nid == pset_pair.first);

      if(nid == id) {
        pset = pset_pair.second;
        if (pset.has_key("trigger_paths")){
          _trigPathsNames = pset.get<std::vector<std::string>>("trigger_paths",std::vector<std::string>());
        }
      }
    }

    //loop over trigResults to fill the map <string, unsigned int)
    std::string   delimeter=":";
    for (unsigned int i=0; i< _trigPathsNames.size(); ++i){
      size_t       pos      = _trigPathsNames[i].find(delimeter);
      unsigned int bit      = std::stoi(_trigPathsNames[i].substr(0, pos));
      std::string  pathName = _trigPathsNames[i].substr(pos+1, _trigPathsNames[i].length());
      _trigMap.insert(std::pair<std::string, unsigned int>(pathName, i));
      _trigPathMap.insert(std::pair<std::string, unsigned int>(pathName, bit));
    }
  }


  std::string const
  TriggerResultsNavigator::getTrigPathName(unsigned int const i) const
  {
    std::string   delimeter =":";
    size_t        pos       = _trigPathsNames[i].find(delimeter);
    if (pos > _trigPathsNames[i].length()) return "TRIG PATH NOT FOUND";
    return _trigPathsNames[i].substr(pos+1, _trigPathsNames[i].length());
  }

  size_t
  TriggerResultsNavigator::getTrigBit(unsigned int const i) const
  {
    if (i>_trigPathsNames.size()) {
      throw cet::exception("TRIGGER") << "TRIG PATHID " << i << " NOT FOUND";
      return 0;
    }
    std::string   delimeter =":";
    size_t        pos       = _trigPathsNames[i].find(delimeter);
    unsigned int  bit       = std::stoi(_trigPathsNames[i].substr(0, pos));
    return bit;
  }

  size_t
  TriggerResultsNavigator::findTrigPath(std::string const& name) const
  {
    return find(_trigMap, name);
  }

  size_t
  TriggerResultsNavigator::find(std::map<std::string, unsigned int> const& posmap, std::string const& name) const
  {
    auto const pos = posmap.find(name);
    if (pos == posmap.cend()) {
      return posmap.size();
    } else {
      return pos->second;
    }
  }

  size_t
  TriggerResultsNavigator::findTrigPathID(std::string const& name) const
  {
    return find(_trigPathMap, name);
  }

  // Has ith path accepted the event?
  bool
  TriggerResultsNavigator::accepted(std::string const& name) const
  {
    size_t index = findTrigPath(name);
    //    return _trigResults->accept(index);
    if (index == _trigResults->size()) return false;
    else                             return _trigResults->accept(index);
  }

  bool
  TriggerResultsNavigator::wasrun(std::string const& name) const
  {
    size_t index = findTrigPath(name);
    return _trigResults->wasrun(index);
  }

  std::vector<std::string>
  TriggerResultsNavigator::triggerModules(std::string const& name) const{
    std::vector<std::string>     modules;

    for ( auto const& i : fhicl::ParameterSetRegistry::get() ){
      auto const  id = i.first;
      if (i.second.has_key(name)){
        auto const &pset = fhicl::ParameterSetRegistry::get(id);
        modules = pset.get<std::vector<std::string>>(name);
        break;
      }
    }
    return modules;
  }

  unsigned
  TriggerResultsNavigator::indexLastModule(std::string const& name) const{
    size_t index = findTrigPath(name);
    return _trigResults->index(index);
   }

  std::string
  TriggerResultsNavigator::nameLastModule (std::string const& name) const{
    unsigned                    indexLast  = indexLastModule(name);
    std::vector<std::string>    modulesVec = triggerModules(name);

    if ( modulesVec.size() == 0) {
      std::string nn = "PATH "+name+" NOT FOUND";
      std::cout << "[TriggerResultsNavigator::nameLastModule] " << nn << std::endl;
      return nn;
    }else {
      return modulesVec[indexLast];
    }
  }

  art::hlt::HLTState
  TriggerResultsNavigator::state(std::string const& name) const{
    size_t index = findTrigPath(name);
    return _trigResults->state(index);
  }

  void
  TriggerResultsNavigator::print() const {
    std::cout << "TriggerResultsNavigator Map" << std::endl;
    std::cout << "//------------------------------------------------//" << std::endl;
    std::cout << "//      trig_pathName           id      accepted  //" << std::endl;
    std::cout << "//------------------------------------------------//" << std::endl;

    for (unsigned i=0; i< getTrigPaths().size(); ++i) {
      const std::string path = getTrigPathName(i);
      const int bit = findTrigPathID(path);
      const bool good = accepted(path);
      std::cout << std::right;
      std::cout <<"//"<<std::setw(40) << path << std::setw(5) << bit << " " << good << " //"<< std::endl;
    }

  }

}
