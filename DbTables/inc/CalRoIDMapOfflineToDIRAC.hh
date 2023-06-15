#ifndef DbTables_CalRoIDMapOfflineToDIRAC_hh
#define DbTables_CalRoIDMapOfflineToDIRAC_hh

#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class CalRoIDMapOfflineToDIRAC : public DbTable {
 public:
  typedef std::shared_ptr<CalRoIDMapOfflineToDIRAC> ptr_t;
  typedef std::shared_ptr<const CalRoIDMapOfflineToDIRAC> cptr_t;
  constexpr static const char* cxname = {"CalRoIDMapOfflineToDIRAC"};

  class Row {
   public:
    Row(int offlineID, uint16_t diracID) :
        _offlineID(offlineID), _diracID(diracID) {}
    int offlineID() const { return _offlineID; }
    uint16_t diracID() const { return _diracID; }

   private:
    int _offlineID;
    uint16_t _diracID;
  };

  CalRoIDMapOfflineToDIRAC() :
      DbTable(cxname, "cal.roidmapofflinetodirac", "offlineid,diracid") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const { return _rows.size(); };
  virtual std::size_t nrowFix() const { return CaloConst::_nRawChannel; };
  size_t size() const { return baseSize() + nrow() * sizeof(Row); };

  void addRow(const std::vector<std::string>& columns) override {
    _rows.emplace_back(std::stoi(columns[0]), std::stoi(columns[1]));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.offlineID() << "," << r.diracID();
  }

  virtual void clear() override {
    baseClear();
    _rows.clear();
  }

 private:
  std::vector<Row> _rows;
};

};  // namespace mu2e
#endif
