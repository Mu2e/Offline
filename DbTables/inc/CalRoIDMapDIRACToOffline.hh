#ifndef DbTables_CalRoIDMapDIRACToOffline_hh
#define DbTables_CalRoIDMapDIRACToOffline_hh

#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class CalRoIDMapDIRACToOffline : public DbTable {
 public:
  typedef std::shared_ptr<CalRoIDMapDIRACToOffline> ptr_t;
  typedef std::shared_ptr<const CalRoIDMapDIRACToOffline> cptr_t;
  constexpr static const char* cxname = {"CalRoIDMapDIRACToOffline"};

  class Row {
   public:
    Row(int diracID, uint16_t caloRoID) :
        _diracID(diracID), _caloRoID(caloRoID) {}
    int diracID() const { return _diracID; }
    uint16_t caloRoID() const { return _caloRoID; }

   private:
    int _diracID;
    uint16_t _caloRoID;
  };

  CalRoIDMapDIRACToOffline() :
      DbTable(cxname, "cal.roidmapdiractooffline", "diracid,caloroid") {}
  const Row& rowAt(const std::size_t diracID) const {
    return _rows.at(diracID);
  }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const { return _rows.size(); };
  virtual std::size_t nrowFix() const { return CaloConst::_nChannel; };
  size_t size() const { return baseSize() + nrow() * sizeof(Row); };

  void addRow(const std::vector<std::string>& columns) override {
    int channel = std::stoi(columns[0]);
    // enforce a strict sequential order - optional
    if (channel != int(_rows.size())) {
      throw cet::exception("CalRoIDMapDIRACToOffline_BAD_INDEX")
          << "CalRoIDMapDIRACToOffline::addRow found index out of order: "
          << channel << " != " << _rows.back().caloRoID() + 1 << "\n";
    }
    _rows.emplace_back(channel, std::stoi(columns[1]));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.diracID() << "," << r.caloRoID();
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
