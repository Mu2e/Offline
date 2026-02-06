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

    // Simplest thing: ParameterSetRegistry has correct ParameterSet
    auto const  id   = trigResults->parameterSetID();
    fhicl::ParameterSet pset;
    fhicl::ParameterSetRegistry::get(id, pset);

    // Go through ParameterSetRegistry and check for corrupted IDs
    auto const& psets = fhicl::ParameterSetRegistry::get();
    for(auto const& pset_pair : psets) {
      auto nid = pset_pair.second.id();
      assert(nid == pset_pair.first);

      // if this is the trigger results parameter set ID, get the trigger path list
      if(nid == id) {
        pset = pset_pair.second;
        if(pset.has_key("trigger_paths")) {
          _trigPathsNames = pset.get<std::vector<std::string>>("trigger_paths");
        }
      }
    }

    //loop over trigResults to fill the map <string, unsigned int>
    std::string   delimeter=":";
    for (size_t index = 0; index < _trigPathsNames.size(); ++index){
      const std::string& path = _trigPathsNames.at(index);
      size_t       pos        = path.find(delimeter);
      if(pos == std::string::npos)
        throw cet::exception("TRIGGER") << "No path bit found for  path " << _trigPathsNames.at(index);
      unsigned int bit        = std::stoi(path.substr(0, pos));
      std::string  pathName   = path.substr(pos+1, path.length());
      _trigMap        .insert(std::pair<std::string, unsigned int>(pathName, index));
      _trigPathMap    .insert(std::pair<std::string, unsigned int>(pathName, bit));
      _indexToPathName.insert(std::pair<unsigned int, std::string>(index, pathName));
      _bitToPathName  .insert(std::pair<unsigned int, std::string>(bit, pathName));
      _bitToIndex     .insert(std::pair<unsigned int, unsigned int>(bit, index));
    }
  }


  std::string const&
  TriggerResultsNavigator::getTrigPathNameByIndex(unsigned int const index) const
  {
    if(!_indexToPathName.contains(index)) {
      throw cet::exception("TRIGGER") << "TRIG PATH INDEX " << index << " NOT FOUND";
    }
    return _indexToPathName.at(index);
  }

  std::string const&
  TriggerResultsNavigator::getTrigPathNameByBit(unsigned int const bit) const
  {
    if(!_bitToPathName.contains(bit)) {
      throw cet::exception("TRIGGER") << "TRIG PATH BIT " << bit << " NOT FOUND";
    }
    return _bitToPathName.at(bit);
  }

  size_t
  TriggerResultsNavigator::getTrigBitByIndex(unsigned int const index) const
  {
    if(!_indexToPathName.contains(index)) {
      throw cet::exception("TRIGGER") << "TRIG PATH INDEX " << index << " NOT FOUND";
    }
    return _trigPathMap.at(_indexToPathName.at(index));
  }

  size_t
  TriggerResultsNavigator::getTrigBitByName(const std::string& name) const
  {
    if(!_trigPathMap.contains(name)) {
      throw cet::exception("TRIGGER") << "TRIG PATH NAME " << name << " NOT FOUND";
    }
    return _trigPathMap.at(name);
  }

  // Return the trigger path index, if found
  size_t
  TriggerResultsNavigator::getTrigPathIndex(std::string const& name) const
  {
    if(!_trigMap.count(name)) return NOTFOUND;
    return _trigMap.at(name);
  }
  size_t
  TriggerResultsNavigator::getTrigPathIndex(const size_t bit) const
  {
    if(!_bitToPathName.count(bit)) return NOTFOUND;
    return _trigMap.at(_bitToPathName.at(bit));
  }

  // Has ith path accepted the event?
  bool
  TriggerResultsNavigator::accepted(std::string const& name) const
  {
    if(!_trigPathMap.count(name)) return false;
    return _trigResults->accept(getTrigPathIndex(name));
  }

  bool
  TriggerResultsNavigator::accepted(unsigned int const bit) const {
    if(!_bitToIndex.count(bit)) return false;
    return _trigResults->accept(getTrigPathIndex(bit));
  }

  bool
  TriggerResultsNavigator::wasrun(std::string const& name) const
  {
    if(!_trigPathMap.count(name)) return false;
    return _trigResults->wasrun(getTrigPathIndex(name));
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
    if(!_trigPathMap.count(name)) return -1;
    return _trigResults->index(getTrigPathIndex(name));
   }

  std::string
  TriggerResultsNavigator::nameLastModule (std::string const& name) const{
    unsigned                    indexLast  = indexLastModule(name);
    std::vector<std::string>    modulesVec;
    if(indexLast != unsigned(-1)) modulesVec = triggerModules(name);

    if(modulesVec.empty()) {
      std::string nn = "PATH "+name+" NOT FOUND";
      std::cout << "[TriggerResultsNavigator::nameLastModule] " << nn << std::endl;
      return nn;
    } else {
      return modulesVec[indexLast];
    }
  }

  art::hlt::HLTState
  TriggerResultsNavigator::state(std::string const& name) const{
    if(!_trigPathMap.count(name))
      throw cet::exception("TRIGGER") << "Path " << name << " not found!";
    return _trigResults->state(getTrigPathIndex(name));
  }

  void
  TriggerResultsNavigator::print() const {
    std::cout << "TriggerResultsNavigator Map" << std::endl;
    std::cout << "//------------------------------------------------//" << std::endl;
    std::cout << "//      trig_pathName           id      accepted  //" << std::endl;
    std::cout << "//------------------------------------------------//" << std::endl;

    const size_t npaths = getTrigPaths().size();
    for (size_t i = 0; i < npaths; ++i) {
      const std::string path = getTrigPathNameByIndex(i);
      const int bit = getTrigBitByName(path);
      const bool good = accepted(path);
      std::cout << std::right;
      std::cout <<"//"<<std::setw(40) << path << std::setw(5) << bit << " " << good << " //"<< std::endl;
    }

  }

}
