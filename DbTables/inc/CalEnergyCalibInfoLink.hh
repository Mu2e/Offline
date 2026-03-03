#ifndef DbTables_CalEnergyCalibInfoLink_hh
#define DbTables_CalEnergyCalibInfoLink_hh

//
// Ad-hoc table linking a CalEnergyCalib CID to the matching
// CalEnergyCalibInfo CID produced by the same combination pass.
//

#include "Offline/DbTables/inc/DbTable.hh"
#include <sstream>
#include <string>
#include <utility>

namespace mu2e {

class CalEnergyCalibInfoLink : public DbTable {
 public:
  class Row {
   public:
    Row(int energyCalibCid, int energyCalibInfoCid, std::string comment) :
        _energyCalibCid(energyCalibCid),
        _energyCalibInfoCid(energyCalibInfoCid),
        _comment(std::move(comment)) {}

    int energyCalibCid() const { return _energyCalibCid; }
    int energyCalibInfoCid() const { return _energyCalibInfoCid; }
    std::string comment() const { return _comment; }

   private:
    int _energyCalibCid;
    int _energyCalibInfoCid;
    std::string _comment;
  };

  constexpr static const char* cxname = "CalEnergyCalibInfoLink";

  CalEnergyCalibInfoLink() :
      DbTable(cxname, "cal.energycalibinfolink",
              "energycalibcid,energycalibinfocid,comment") {}

  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); }
  size_t size() const override { return baseSize() + nrow() * sizeof(Row); }
  tableType type() const override { return Adhoc; }

  void addRow(const std::vector<std::string>& columns) override {
    _rows.emplace_back(std::stoi(columns[0]), std::stoi(columns[1]), columns[2]);
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.energyCalibCid() << ",";
    sstream << r.energyCalibInfoCid() << ",";
    sstream << r.comment();
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
