#include <iostream>
#include "DbTables/inc/DbTable.hh"
#include "DbTables/inc/DbUtil.hh"
#include "cetlib_except/exception.h"
#include <boost/algorithm/string/split.hpp>

int mu2e::DbTable::fill(const std::string& csv, bool saveCsv) {
  std::vector<std::string> lines = DbUtil::splitCsvLines(csv);
  std::vector<std::string> columns; // columns of a row as text
  size_t ncol = 0; // check they all have the same columns
  for(auto const& line: lines) {
    columns = DbUtil::splitCsv(line);
    if(ncol==0) ncol = columns.size();
    if(columns.size()!=ncol) {
      throw cet::exception("DBTABLE_BAD_COLUMN_COUNT")  
	      << "DbTable::fill found "<<columns.size()<<" columns "
	      <<" when "<<ncol<< " was seen in previous rows. Text:"
	      << line <<" \n";
    }
    addRow(columns);
  }

  // if this table has a fixed number of rows, check that
  if(nrowFix()>0 && nrow()!=nrowFix()) {
    throw cet::exception("DBTABLE_BAD_ROW_COUNT") 
      << "DbTable::fill csv line counts is "
      << std::to_string(nrow()) << " but "
      << std::to_string(nrowFix()) << " is required while filling "
      << name();
  }

  // save the plain text
  if(saveCsv) {
    _csv = csv;
  } else {
    _csv.clear();
  }

  return 0;
}

int mu2e::DbTable::toCsv() {
  if(!_csv.empty()) return 0;
  std::ostringstream ss;
  for(std::size_t i=0; i< nrow(); i++) rowToCsv(ss,i);
  _csv = ss.str();
  return 0;
}

void mu2e::DbTable::addRow(const std::vector<std::string>& columns) {
  throw cet::exception("DBTABLE_FUNCTION_NOT_IMPLEMENTED") 
    << "DbTable::addRow must be overridden ";
}

void mu2e::DbTable::rowToCsv(std::ostringstream& stream, size_t irow) const {
  throw cet::exception("DBTABLE_FUNCTION_NOT_IMPLEMENTED") 
    << "DbTable::rowToCsv must be overridden ";
}
