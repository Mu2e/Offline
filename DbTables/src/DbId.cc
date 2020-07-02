#include <iostream>
#include "cetlib_except/exception.h"
#include "DbTables/inc/DbId.hh"

void mu2e::DbId::setDb(std::string const& name) {
  _name = name;
  if(_name=="mu2e_conditions_prd") {
    _host = "ifdb05";
    _port = "5448";
    _url =        "https://dbdata0vm.fnal.gov:8444/QE/mu2e/prod/app/SQ/query?";
    _urlNoCache = "https://dbdata0vm.fnal.gov:9443/QE/mu2e/prod/app/SQ/query?";
  } else if(_name=="mu2e_conditions_dev") {
    _host = "ifdb04";
    _port = "5444";
    _url =        "https://dbdata0vm.fnal.gov:8444/QE/mu2e/dev/app/SQ/query?";
    _urlNoCache = "https://dbdata0vm.fnal.gov:9443/QE/mu2e/dev/app/SQ/query?";
  } else {
    throw cet::exception("UNKNOWN_DB_ID") 
      << "attempted to create DbId for unknown database "<<_name<<"\n";
  }
}

