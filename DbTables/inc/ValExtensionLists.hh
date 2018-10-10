#ifndef DbTables_ValExtensionLists_hh
#define DbTables_ValExtensionLists_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class ValExtensionLists : public DbTable {
  public:

    class Row {
    public:
      Row(int eid, int gid):
	_eid(eid),_gid(gid) {}
      int  eid() const { return _eid;}
      int  gid() const { return _gid;}
    private:
      int _eid;
      int _gid;
    };


    ValExtensionLists():DbTable("ValExtensionLists","val.extensionlists",
			    "eid,gid") {}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    size_t size() const { return sizeof(this) + _csv.capacity() + nrow()*8; };

    void addRow(const std::vector<std::string>& columns) {
      _rows.emplace_back(std::stoi(columns[0]),std::stoi(columns[1]));
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.eid()<<",";
      sstream << r.gid();
    }

    virtual void clear() { _csv.clear(); _rows.clear(); }

  private:
    std::vector<Row> _rows;
  };
  
};
#endif
