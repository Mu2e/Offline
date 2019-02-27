#ifndef DbTables_TstCalib1_hh
#define DbTables_TstCalib1_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "cetlib_except/exception.h"
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class TstCalib1 : public DbTable {
  public:

    class Row {
    public:
      Row(int channel, int flag, float dtoe):
	_channel(channel),_flag(flag),_dtoe(dtoe) {}
      int  channel() const { return _channel;}
      int     flag() const { return _flag; }
      float   dToE() const { return _dtoe; }
    private:
      int _channel;
      int _flag;
      float _dtoe;
    };


    TstCalib1():DbTable("TstCalib1","tst.calib1","channel,flag,dtoe") {}
    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int channel) const { 
                return _rows.at(_chanIndex.at(channel)); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    //this table should always be 3 rows
    virtual std::size_t nrowFix() const { return 3; }; 
    size_t size() const { return _csv.capacity() + 
	+ nrow()*nrow()/2 + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) {
      int channel = std::stoi(columns[0]);
      // enforce a strict sequential order - optional
      if(channel!=int(_rows.size())) {
	throw cet::exception("TSTCALIB1_BAD_INDEX") 
	  << "TstCalib1::addRow found index out of order: " 
	  <<channel << " != " << _rows.back().channel()+1 <<"\n";
      }
      _rows.emplace_back(channel,
			 std::stoi(columns[1]),
			 std::stof(columns[2]) );
      // add this channel to the map index - optional
      _chanIndex[_rows.back().channel()] = _rows.size()-1;
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.channel()<<",";
      sstream << r.flag()<<",";
      sstream << std::fixed << std::setprecision(3) << r.dToE();
    }

    virtual void clear() { _csv.clear(); _rows.clear(); _chanIndex.clear();}

  private:
    std::vector<Row> _rows;
    std::map<int,std::size_t> _chanIndex;
  };
  
};
#endif
