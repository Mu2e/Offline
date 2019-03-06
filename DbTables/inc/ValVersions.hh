#ifndef DbTables_ValVersions_hh
#define DbTables_ValVersions_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class ValVersions : public DbTable {
  public:

    class Row {
    public:
      Row(int vid, int pid, int lid, int major, int minor,
	  std::string comment,
	  std::string create_time, std::string create_user):
	_vid(vid),_pid(pid),_lid(lid),_major(major),_minor(minor),
	_comment(comment),
	_create_time(create_time),_create_user(create_user) {}
      int  vid() const { return _vid;}
      int  pid() const { return _pid;}
      int  lid() const { return _lid;}
      int  major() const { return _major;}
      int  minor() const { return _minor;}
      std::string const& comment() const { return _comment; }
      std::string const& create_time() const { return _create_time; }
      std::string const& create_user() const { return _create_user; }
    private:
      int _vid;
      int _pid;
      int _lid;
      int _major;
      int _minor;
      std::string _comment;
      std::string _create_time;
      std::string _create_user;
    };


    ValVersions():DbTable("ValVersions","val.versions",
      	  "vid,pid,lid,major,minor,comment,create_time,create_user") {}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    size_t size() const {
      size_t b = sizeof(this) + _csv.capacity() + nrow()*60;
      for (auto const& r : _rows) b += r.comment().capacity();
      return b;
    };

    void addRow(const std::vector<std::string>& columns) {
      _rows.emplace_back(std::stoi(columns[0]),std::stoi(columns[1]),
			 std::stoi(columns[2]),std::stoi(columns[3]),
			 std::stoi(columns[4]),columns[5],
			 columns[6],columns[7]);
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.vid()<<",";
      sstream << r.pid()<<",";
      sstream << r.lid()<<",";
      sstream << r.major()<<",";
      sstream << r.minor()<<",";
      sstream << r.comment()<<",";
      sstream << r.create_time()<<",";
      sstream << r.create_user();
    }

    virtual void clear() { _csv.clear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
  };
  
};
#endif
