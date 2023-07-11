#ifndef DbTables_ValExtensionLists_hh
#define DbTables_ValExtensionLists_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class ValExtensionLists : public DbTable {
 public:
  class Row {
   public:
    Row(int eid, int gid) : _eid(eid), _gid(gid) {}
    int eid() const { return _eid; }
    int gid() const { return _gid; }

   private:
    int _eid;
    int _gid;
  };

  constexpr static const char* cxname = "ValExtensionLists";

  ValExtensionLists() : DbTable(cxname, "val.extensionlists", "eid,gid") {}

  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  size_t size() const override {
    return baseSize() + sizeof(this) + nrow() * 8;
  };
  const std::string orderBy() const { return std::string("eid,gid"); }

  void addRow(const std::vector<std::string>& columns) override {
    _rows.emplace_back(std::stoi(columns[0]), std::stoi(columns[1]));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.eid() << ",";
    sstream << r.gid();
  }

  virtual void clear() override {
    baseClear();
    _rows.clear();
  }

 private:
  std::vector<Row> _rows;
};

}  // namespace mu2e
#endif
