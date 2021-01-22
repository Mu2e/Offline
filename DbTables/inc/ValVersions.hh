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

    constexpr static const char* cxname = "ValVersions";

    ValVersions():DbTable(cxname,"val.versions",
      	  "vid,pid,lid,major,minor,comment,create_time,create_user") {}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int vid) const { return _rows.at(_index.at(vid)); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override {
      size_t b = baseSize() + sizeof(this) + nrow()*60;
      for (auto const& r : _rows) b += r.comment().capacity();
      return b;
    };

    void addRow(const std::vector<std::string>& columns) override {
      _rows.emplace_back(std::stoi(columns[0]),std::stoi(columns[1]),
			 std::stoi(columns[2]),std::stoi(columns[3]),
			 std::stoi(columns[4]),columns[5],
			 columns[6],columns[7]);
      _index[_rows.back().vid()] = _rows.size()-1;
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
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

    virtual void clear() override { baseClear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
    std::map<int,std::size_t> _index;
 };
  
};
#endif
