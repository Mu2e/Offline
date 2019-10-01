#ifndef DbTables_ValPurposes_hh
#define DbTables_ValPurposes_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class ValPurposes : public DbTable {
  public:

    class Row {
    public:
      Row(int pid, std::string name, std::string comment,
	  std::string create_time, std::string create_user):
	_pid(pid),_name(name),_comment(comment),
	_create_time(create_time),_create_user(create_user) {}
      int  pid() const { return _pid;}
      std::string const& name() const { return _name;}
      std::string const& comment() const { return _comment;}
      std::string const& create_time() const { return _create_time; }
      std::string const& create_user() const { return _create_user; }
    private:
      int _pid;
      std::string _name;
      std::string _comment;
      std::string _create_time;
      std::string _create_user;
    };


    ValPurposes():DbTable("ValPurposes","val.purposes",
			  "pid,name,comment,create_time,create_user") {}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int pid) const { return _rows.at(_index.at(pid)); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    size_t size() const {
      size_t b = sizeof(this) + _csv.capacity() + nrow()*44 + nrow()*nrow()/2;
      for (auto const& r : _rows) b += r.name().capacity() + r.comment().capacity();
      return b;
    };

    void addRow(const std::vector<std::string>& columns) {
      _rows.emplace_back(std::stoi(columns[0]),columns[1],columns[2],
			 columns[3],columns[4]);
      _index[_rows.back().pid()] = _rows.size()-1;
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.pid()<<",";
      sstream << r.name()<<",";
      sstream << r.comment()<<",";
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
