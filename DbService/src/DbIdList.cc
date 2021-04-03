#include <iostream>
#include <fstream>
#include "cetlib_except/exception.h"
#include "DbService/inc/DbIdList.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>


mu2e::DbIdList::DbIdList() {

  // read a file listing the connections so jobs can be steered 
  // to new connections in almost real time

   ConfigFileLookupPolicy configFile;
   std::string fn = configFile("Database/connections.txt");

   std::ifstream in( fn );
    if ( !in ) {
      throw cet::exception("DBID_NO_CONNECTION")
        << "DbIdList::DbIdList could not open connections file " << fn << std::endl;
    }

    std::string line;
    std::vector<std::string> words;
    std::string name,host,port,url,urlNoCache;
    while ( in ){
      std::getline(in,line);
      line = line.substr(0, line.find_first_of("#"));      
      boost::trim(line);
      boost::split(words,line, boost::is_any_of(" \t"), boost::token_compress_on);
      // if eof, or starting a new database stanza with a previous database name filled
      if(  in.eof() || ( words.size()==2 && words[0] == "database" && !name.empty()) ) {
	// have read a stanza, check and store it
	bool fail = host.empty() || port.empty() || url.empty() || urlNoCache.empty();
	if ( fail ) {
	  throw cet::exception("DBIDFILE_DATABASE_NOT_LOADED")
	    << "DbIdList::DbIdList could not load database " << name 
	    << " in connections file" << std::endl;
	}
	_ids.emplace_back(name,host,port,url,urlNoCache);
	name=host=port=url=urlNoCache="";
      } // if database
      if(words.size()==2) {
	if(words[0]=="database") {
	  name = words[1];
	} else if(words[0]=="host") {
	  host = words[1];
	} else if(words[0]=="port") {
	  port = words[1];
	} else if(words[0]=="cache") {
	  url = words[1];
	} else if(words[0]=="nocache") {
	  urlNoCache = words[1];
	} else {
	  throw cet::exception("DBIDFILE_BAD_KEYWORD")
	    << "DbIdList::DbIdList unknown keyword: " << words[0] << std::endl;
	}
      } else if (!line.empty()) {
	  throw cet::exception("DBIDFILE_BAD_LINE")
	    << "DbIdList::DbIdList poorly formed line: " << line << std::endl;
      }
    }

    in.close();

}

mu2e::DbId mu2e::DbIdList::getDbId(const std::string& name) {
  for(auto const& id : _ids) {
    if(id.name()==name) {
      return id;
    }
  }
  throw cet::exception("DBIDFILE_DATABASE_NOT_FOUND")
    << "DbIdList::getDbId could not load database " << name << std::endl;

}
