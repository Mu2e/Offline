#ifndef DbTables_TstCalib2_hh
#define DbTables_TstCalib2_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class TstCalib2 : public DbTable {
  public:

    class Row {
    public:
      Row(int channel, std::string status):
	 _channel(channel),_status(status) {}
      int  channel() const { return _channel;}
      std::string const&  status() const { return _status; }
    private:
      int _channel;
      std::string _status;
    };


    TstCalib2():DbTable("TstCalib2","tst.calib2","channel,status") {}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int channel) const { 
                return _rows.at(_chanIndex.at(channel)); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    // this is not quite right with strings...
    size_t size() const { 
      size_t b = _csv.capacity() + nrow()*nrow()/2 + nrow()*sizeof(Row); 
      for(auto const& r : _rows) b += r.status().capacity(); 
      return b;
    };

    void addRow(const std::vector<std::string>& columns) {
      _rows.emplace_back(std::stoul(columns[0]),
			 columns[1] );
      _chanIndex[_rows.back().channel()] = _rows.size()-1;
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.channel()<<",";
      sstream << r.status();
    }

    virtual void clear() { _csv.clear(); _rows.clear(); _chanIndex.clear();}

  private:
    std::vector<Row> _rows;
    std::map<int,std::size_t> _chanIndex;
  };
  
};
#endif
