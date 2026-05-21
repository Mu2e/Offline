#ifndef DbTables_CalEnergyCalibInfo_hh
#define DbTables_CalEnergyCalibInfo_hh

//
// Ad-hoc table linking a CalEnergyCalib CID to the matching
// CalCombinedEnergyCalib CID produced by the same combination pass.
//

#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <sstream>
#include <string>

namespace mu2e {

class CalEnergyCalibInfo : public DbTable {
 public:
  class Row {
   public:
    Row(int energyCalibCid, int combinedEnergyCalibCid,
        const std::string& comment) :
        _energyCalibCid(energyCalibCid),
        _combinedEnergyCalibCid(combinedEnergyCalibCid),
        _comment(comment) {}

    int energyCalibCid() const { return _energyCalibCid; }
    int combinedEnergyCalibCid() const { return _combinedEnergyCalibCid; }
    std::string comment() const { return _comment; }

   private:
    int _energyCalibCid;
    int _combinedEnergyCalibCid;
    std::string _comment;
  };

  constexpr static const char* cxname = "CalEnergyCalibInfo";

  CalEnergyCalibInfo() :
      DbTable(cxname, "cal.energycalibinfo",
              "energycalibcid,combinedenergycalibcid,comment") {}

  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); }
  size_t size() const override { return baseSize() + nrow() * sizeof(Row); }
  tableType type() const override { return Adhoc; }

  void addRow(const std::vector<std::string>& columns) override {
    if (columns.size() < 3) {
      throw cet::exception("CALENERGYCALIBINFO_BAD_COLUMNS")
          << "CalEnergyCalibInfo::addRow expected at least 3 columns, got "
          << columns.size() << "\n";
    }
    _rows.emplace_back(std::stoi(columns[0]), std::stoi(columns[1]),
                       columns[2]);
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.energyCalibCid() << ",";
    sstream << r.combinedEnergyCalibCid() << ",";
    sstream << "\"" << r.comment() << "\"";
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
