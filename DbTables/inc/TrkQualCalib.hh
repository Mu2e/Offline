#ifndef DbTables_TrkQualCalib_hh
#define DbTables_TrkQualCalib_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "cetlib_except/exception.h"
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class TrkQualCalib : public DbTable {
  public:

    class Row {
    public:
      Row(int channel, float eff, float cut):
	_channel(channel),_eff(eff),_cut(cut) {}
      int    channel() const { return _channel; }
      float  eff() const { return _eff;}
      float  cut() const { return _cut; }
    private:
      int _channel;
      float _eff;
      float _cut;
    };


    TrkQualCalib():DbTable("TrkQualCalib","trkqual.calib","channel,eff,cut") {}
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
	  << "TrkQualCalib::addRow found index out of order: " 
	  <<channel << " != " << _rows.back().channel()+1 <<"\n";
      }
      _rows.emplace_back(channel,
			 std::stof(columns[1]),
			 std::stof(columns[2]) );
      // add this channel to the map index - optional
      _chanIndex[_rows.back().channel()] = _rows.size()-1;
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.channel()<<",";
      sstream << std::fixed << std::setprecision(3) << r.eff()<<",";
      sstream << std::fixed << std::setprecision(3) << r.cut();
    }

    virtual void clear() { _csv.clear(); _rows.clear(); _chanIndex.clear();}

  private:
    std::vector<Row> _rows;
    std::map<int,std::size_t> _chanIndex;
  };
  
};
#endif
