//
// Mu2eG4Helper plugin.
//
//
// Original author Rob Kutschke
//

// Framework includes

// Mu2e includes
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

using namespace std;

namespace mu2e {

  Mu2eG4Helper::Mu2eG4Helper(fhicl::ParameterSet const& iPS,
                     art::ActivityRegistry&iRegistry){
  }

  Mu2eG4Helper::~Mu2eG4Helper(){
  }

  // Return the volume info mapped to the given key, throw if the key does not exist.
  VolumeInfo& Mu2eG4Helper::locateVolInfo( const std::string key){
    std::map<std::string,VolumeInfo>::iterator i = _volumeInfoList.find(key);
    if ( i == _volumeInfoList.end() ){
      throw cet::exception("GEOM")
        << "Mu2eG4Helper::locateVolInfo cannot find the volume named: "
        << key
        << "\n";
    }
    return i->second;
  } // end of Mu2eG4Helper::locateVolInfo by key

  // If the key already exists, throw. Otherwise add the (key, value) pair
  // to the map.
  void Mu2eG4Helper::addVolInfo( const VolumeInfo& info ){
    std::map<std::string,VolumeInfo>::iterator i = _volumeInfoList.find(info.name);
    if ( i != _volumeInfoList.end() ){
      throw cet::exception("GEOM")
        << "locateVolInfo already has the key: "
        << info.name
        << "\n";
    }
    _volumeInfoList[info.name] = info;
  } // end of Mu2eG4Helper::addVolInfo

  // Find all info objects whose key matches the supplied regular expression.
  std::vector<VolumeInfo const*> Mu2eG4Helper::locateVolInfo( boost::regex const& expression ) const{

    // Default, empty return value
    std::vector<VolumeInfo const*> infos;

    for ( auto const& i : _volumeInfoList ){

      std::string const& name{i.first};
      VolumeInfo const&  info{i.second};

      if ( boost::regex_match( name.c_str(), expression ) ){
        infos.push_back( &info );
      }
    }

    return infos;

  } // end of Mu2eG4Helper::addlocateVolInfo by regex on key


} // end namespace mu2e

using mu2e::Mu2eG4Helper;
DEFINE_ART_SERVICE(Mu2eG4Helper);
