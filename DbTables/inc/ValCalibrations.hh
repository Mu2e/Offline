#ifndef DbTables_ValCalibrations_hh
#define DbTables_ValCalibrations_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class ValCalibrations : public DbTable {
  public:

    class Row {
    public:
      Row(int cid, int tid, std::string create_time, std::string create_user):
	_cid(cid),_tid(tid),_create_time(create_time),
	_create_user(create_user) {}
      int  tid() const { return _tid;}
      int  cid() const { return _cid;}
      std::string const& create_time() const { return _create_time; }
      std::string const& create_user() const { return _create_user; }
    private:
      int _cid;
      int _tid;
      std::string _create_time;
      std::string _create_user;
    };


    ValCalibrations():DbTable("ValCalibrations","val.calibrations",
			      "cid,tid,create_time,create_user") {}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int cid) const { return _rows.at(_index.at(cid)); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    size_t size() const { return sizeof(this) + _csv.capacity() 
	+ nrow()*nrow()/2 + nrow()*48; };

    void addRow(const std::vector<std::string>& columns) {
      _rows.emplace_back(std::stoi(columns[0]),std::stoi(columns[1]),
			 columns[2],columns[3]);
      _index[_rows.back().cid()] = _rows.size()-1;
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.cid()<<",";
      sstream << r.tid()<<",";
      sstream << r.create_time()<<",";
      sstream << r.create_user();
    }

    virtual void clear() { _csv.clear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
    std::map<int,std::size_t> _index;
  };
  
};
#endif
