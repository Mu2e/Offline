// Andrei Gaponenko, 2020

#include "Mu2eG4/inc/Mu2eG4PrimaryType.hh"

#include <utility>
#include <algorithm>

#include "cetlib_except/exception.h"

namespace mu2e {

  namespace{

    typedef std::pair<Mu2eG4PrimaryType::enum_type,std::string> PT;

    // a plain array to keep memory nice and tight
    const PT typelist[] = {
#define X(x) { Mu2eG4PrimaryType::x, #x},
      MU2EG4_PRIMARY_TYPES
#undef X
    };

    const auto typelist_end = typelist + sizeof(typelist)/sizeof(typelist[0]);
  }

  //================================================================
  Mu2eG4PrimaryType::Mu2eG4PrimaryType(enum_type id) : id_(id) {
    auto res = std::find_if(typelist, typelist_end,
                            [id](const PT& a){ return a.first == id; }
                            );
    if(res == typelist_end) {
      throw cet::exception("CONFIG")
        << "Error: unknown Mu2eG4PrimaryType enum value "<<id<<std::endl;
    }

    name_ = res->second;
  }

  //================================================================
  Mu2eG4PrimaryType::Mu2eG4PrimaryType(std::string s) : name_(s) {
    auto res = std::find_if(typelist, typelist_end,
                            [s](const PT& a){ return a.second == s; }
                            );
    if(res == typelist_end) {
      throw cet::exception("CONFIG")
        << "Error: unknown Mu2eG4PrimaryType name "<<s
        <<".  Recognized values are: "<<
#define X(x) #x " "
      MU2EG4_PRIMARY_TYPES
#undef X
        <<std::endl;
    }

    id_ = res->first;
  }

  //================================================================
  const std::vector<std::string>& Mu2eG4PrimaryType::all_names() {
    static std::vector<std::string> names{
#define X(x) #x,
      MU2EG4_PRIMARY_TYPES
#undef X
    };
    return names;
  }
}
