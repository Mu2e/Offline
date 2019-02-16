#ifndef DbTables_ValIovs_hh
#define DbTables_ValIovs_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"
#include "DbTables/inc/DbIoV.hh"

namespace mu2e {

  class ValIovs : public DbTable {
  public:

    class Row {
    public:
      Row(int iid, int cid, 
	  int start_run, int start_subrun, int end_run, int end_subrun,
	  std::string create_time, std::string create_user):
	_iid(iid),_cid(cid),
	_start_run(start_run),_start_subrun(start_subrun),
	_end_run(end_run),_end_subrun(end_subrun),
	_create_time(create_time),_create_user(create_user) {}
      int  iid() const { return _iid;}
      int  cid() const { return _cid;}
      int start_run() const { return _start_run; }
      int start_subrun() const { return _start_subrun; }
      int end_run() const { return _end_run; }
      int end_subrun() const { return _end_subrun; }
      DbIoV iov() const {
	return DbIoV(start_run(),start_subrun(),end_run(),end_subrun()); }
      std::string const& create_time() const { return _create_time; }
      std::string const& create_user() const { return _create_user; }
    private:
      int _iid;
      int _cid;
      int _start_run;
      int _start_subrun;
      int _end_run;
      int _end_subrun;
      std::string _create_time;
      std::string _create_user;
    };


    ValIovs():DbTable("ValIovs","val.iovs",
	  "iid,cid,start_run,start_subrun,end_run,end_subrun,create_time,create_user") {}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int iid) const { return _rows.at(_index.at(iid)); }
    std::size_t nrow() const { return _rows.size(); };
    size_t size() const { return sizeof(this) + _csv.capacity() 
	+ nrow()*nrow()/2 + nrow()*64; };

    void addRow(const std::vector<std::string>& columns) {
      _rows.emplace_back(std::stoi(columns[0]),std::stoi(columns[1]),
			 std::stoi(columns[2]),std::stoi(columns[3]),
			 std::stoi(columns[4]),std::stoi(columns[5]),
			 columns[6],columns[7]);
      _index[_rows.back().iid()] = _rows.size()-1;
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.iid()<<",";
      sstream << r.cid()<<",";
      sstream << r.start_run()<<",";
      sstream << r.start_subrun()<<",";
      sstream << r.end_run()<<",";
      sstream << r.end_subrun()<<",";
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
