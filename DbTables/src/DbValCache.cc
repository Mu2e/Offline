#include <iostream>
#include <algorithm>
#include "cetlib_except/exception.h"
#include "DbTables/inc/DbValCache.hh"

mu2e::DbTable const& mu2e::DbValCache::asTable(std::string const& tableName) const {
  if(tableName=="ValTables") return _tables;
  if(tableName=="ValCalibrations") return _calibrations;
  if(tableName=="ValIovs") return _iovs;
  if(tableName=="ValGroups") return _groups;
  if(tableName=="ValGroupLists") return _grouplists;
  if(tableName=="ValPurposes") return _purposes;
  if(tableName=="ValLists") return _lists;
  if(tableName=="ValTableLists") return _tablelists;
  if(tableName=="ValVersions") return _versions;
  if(tableName=="ValExtensions") return _extensions;
  if(tableName=="ValExtensionLists") return _extensionlists;
  throw cet::exception("DBVALCACHE_UNKNOWN_TABLE")  
    << "DbValCache::asTable does not recognize table name "
    << tableName << "\n";
}

void mu2e::DbValCache::clear() {
  _tables.clear();
  _calibrations.clear();
  _iovs.clear();
  _groups.clear();
  _grouplists.clear();
  _purposes.clear();
  _lists.clear();
  _tablelists.clear();
  _versions.clear();
  _extensions.clear();
  _extensionlists.clear();
}

size_t mu2e::DbValCache::size() {
  int sum = 0;
  sum += _tables.size();
  sum += _calibrations.size();
  sum += _iovs.size();
  sum += _groups.size();
  sum += _grouplists.size();
  sum += _purposes.size();
  sum += _lists.size();
  sum += _tablelists.size();
  sum += _versions.size();
  sum += _extensions.size();
  sum += _extensionlists.size();
  return sum;
}

void mu2e::DbValCache::print() {

  std::cout << "       Table        rows      size " << std::endl;
  std::cout << std::setw(17) << valTables().name() 
	    << std::setw(6) << valTables().nrow() 
	    << std::setw(12) << valTables().size() << std::endl;
  std::cout << std::setw(17) << valCalibrations().name() 
	    << std::setw(6) << valCalibrations().nrow() 
	    << std::setw(12) << valCalibrations().size() << std::endl;
  std::cout << std::setw(17) << valIovs().name() 
	    << std::setw(6) << valIovs().nrow() 
	    << std::setw(12) << valIovs().size() << std::endl;
  std::cout << std::setw(17) << valGroups().name() 
	    << std::setw(6) << valGroups().nrow() 
	    << std::setw(12) << valGroups().size() << std::endl;
  std::cout << std::setw(17) << valGroupLists().name() 
	    << std::setw(6) << valGroupLists().nrow() 
	    << std::setw(12) << valGroupLists().size() << std::endl;
  std::cout << std::setw(17) << valPurposes().name() 
	    << std::setw(6) << valPurposes().nrow() 
	    << std::setw(12) << valPurposes().size() << std::endl;
  std::cout << std::setw(17) << valLists().name() 
	    << std::setw(6) << valLists().nrow() 
	    << std::setw(12) << valLists().size() << std::endl;
  std::cout << std::setw(17) << valTableLists().name() 
	    << std::setw(6) << valTableLists().nrow() 
	    << std::setw(12) << valTableLists().size() << std::endl;
  std::cout << std::setw(17) << valVersions().name() 
	    << std::setw(6) << valVersions().nrow() 
	    << std::setw(12) << valVersions().size() << std::endl;
  std::cout << std::setw(17) << valExtensions().name() 
	    << std::setw(6) << valExtensions().nrow() 
	    << std::setw(12) << valExtensions().size() << std::endl;
  std::cout << std::setw(17) << valExtensionLists().name() 
	    << std::setw(6) << valExtensionLists().nrow() 
	    << std::setw(12) << valExtensionLists().size() << std::endl;

}
