#ifndef DbTables_ValGroups_hh
#define DbTables_ValGroups_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class ValGroups : public DbTable {
  public:

    class Row {
    public:
      Row(int gid, std::string create_time, std::string create_user):
	_gid(gid),_create_time(create_time),_create_user(create_user) {}
      int  gid() const { return _gid;}
      std::string const& create_time() const { return _create_time; }
      std::string const& create_user() const { return _create_user; }
    private:
      int _gid;
      std::string _create_time;
      std::string _create_user;
    };


    ValGroups():DbTable("ValGroups","val.groups",
			"gid,create_time,create_user") {}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int gid) const { return _rows.at(_index.at(gid)); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    size_t size() const { return sizeof(this) + _csv.capacity() 
	+ nrow()*44; };

    void addRow(const std::vector<std::string>& columns) {
      _rows.emplace_back(std::stoi(columns[0]),
			 columns[1],columns[2]);
      _index[_rows.back().gid()] = _rows.size()-1;
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.gid()<<",";
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
