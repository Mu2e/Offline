#ifndef DbTables_CalCosmicEnergyCalibInfo_hh
#define DbTables_CalCosmicEnergyCalibInfo_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class CalCosmicEnergyCalibInfo : public DbTable {
 public:
  class Row {
   public:
    Row(int FirstCalibRun, int LastCalibRun, std::string EnergyMethod, std::string FitMethod, std::string Comment) :
        _FirstCalibRun(FirstCalibRun),
        _LastCalibRun(LastCalibRun),
        _EnergyMethod(EnergyMethod),
        _FitMethod(FitMethod),
        _Comment(Comment) {}
    int FirstCalibRun() const { return _FirstCalibRun; }
    int LastCalibRun() const { return _LastCalibRun; }
    std::string EnergyMethod() const { return _EnergyMethod; }
    std::string FitMethod() const { return _FitMethod; }
    std::string Comment() const { return _Comment; }

   private:
    int _FirstCalibRun; //First run the calibration was extracted from
    int _LastCalibRun;  //Last run the calibration was extracted from
    std::string _EnergyMethod; // "vmax", "peak", "integral"
    std::string _FitMethod;    // "langaus", "gaus", "..."
    std::string _Comment;
  };

  constexpr static const char* cxname = "CalCosmicEnergyCalibInfo";

  CalCosmicEnergyCalibInfo() : DbTable(cxname, "cal.cosmicenergycalibinfo", "firstcalibrun,lastcalibrun,energymethod,fitmethod,comment") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  size_t size() const override {
    return baseSize() + nrow() * sizeof(Row);
  };
  tableType type() const override { return Adhoc; }

  void addRow(const std::vector<std::string>& columns) override {
    _rows.emplace_back(std::stoi(columns[0]),std::stoi(columns[1]),columns[2],columns[3],columns[4]);
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.FirstCalibRun() << ",";
    sstream << r.LastCalibRun() << ",";
    sstream << r.EnergyMethod() << ",";
    sstream << r.FitMethod() << ",";
    sstream << r.Comment();
  }

  virtual void clear() override { baseClear(); _rows.clear();}

 private:
  std::vector<Row> _rows;
};

}  // namespace mu2e
#endif
