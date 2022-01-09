
#include <fstream>
#include <iomanip>
#include <iostream>

#include "cetlib_except/exception.h"

#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"


namespace mu2e {

  // *************************************************************

  ParticleDataList::ParticleDataList( SimpleConfig const& config ) {

    ConfigFileLookupPolicy findConfig;
    
    std::string tableFilename = findConfig(
                    config.getString("particleDataList.filename") );
    if( tableFilename.empty() ) {
      throw cet::exception("PARTICLELIST_NO_FILE") 
        << "ParticleDataList: file name not found\n";
    }

    std::ifstream in(tableFilename.c_str());
    if ( !in ) {
      throw cet::exception("PARTICLELIST_OPEN_FAIL")
          << "Unable to open particle data file: "
          << tableFilename << "\n";
    }

    // input text file words on a row
    constexpr size_t nWords = 7;

    std::string line;
    std::string str;
    std::vector<std::string> words;
    words.reserve(nWords);

    while ( std::getline(in,line) ) {
      std::istringstream iss(line);
      while (iss >> str) {
        words.emplace_back(str);
      }
      if ( words.size() != nWords ) {
        throw cet::exception("PARTICLELIST_BAD_LINE")
          << "ParticleList: Not " << nWords << " words on line: "
          << line << "\n";
      }
      int id = std::stoi(words[0]);
      _list.try_emplace(id, id, words[1], words[2], 
         std::stod(words[4]), std::stod(words[5]), std::stod(words[6]) );
      _names.try_emplace(words[1],id);
      if( words[2] != "none" ) _names.try_emplace(words[2],id);
      if( words[3] != "none" ) _names.try_emplace(words[3],id);
      words.clear();
    }

  }

  // *************************************************************

  const ParticleData& ParticleDataList::particle( int id ) const {

    auto it = _list.find(id);
    if( it != _list.end() ) {
      return it->second;
    }

    if ( id <= PDGCode::G4Threshold ){    
      throw cet::exception("PARTICLELIST_BAD_ID")
        << "ParticleList: Failed to look up id: "
        << id << "\n";
    }

    // this is a unknown geant nucleus, add it to the list

    int pA = (std::abs(id)/10)%1000;
    int pZ = (std::abs(id)/10000)%1000;

    // the nuclear excitation level
    int exc = id%10;
    double protonMass = particle(PDGCode::p_plus).mass();
    double mass = (pA+exc)*protonMass;

    std::ostringstream pName;
    pName << "Mu2e_" << id;
    std::string name = pName.str();

    _list.try_emplace(id,id, name, name, double(pZ), mass, 0.0);

    return _list.at(id);
    
  }

  // *************************************************************

  const ParticleData& ParticleDataList::particle( 
                                    const std::string& name ) const {
    auto it = _names.find(name);
    if( it == _names.end() ) {
      throw cet::exception("PARTICLELIST_BAD_NAME")
        << "ParticleList: Failed to look up name: "
        << name << "\n";
    }
    return particle(it->second);
    
  }

  // *************************************************************

  void ParticleDataList::printTable( std::ostream& ostr) const {

    for(auto const& pp: _list) {
      auto const& p = pp.second;
      p.print(ostr);
    }
  }

}  // end namespace mu2e
