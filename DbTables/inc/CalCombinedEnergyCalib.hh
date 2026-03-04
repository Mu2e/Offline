#ifndef DbTables_CalCombinedEnergyCalib_hh
#define DbTables_CalCombinedEnergyCalib_hh

/*
  Per-SiPM combined calibration metadata table.
  Stores the combined ADC/MeV constants, uncertainties, and status flags
  from the cosmic + source calibration combination algorithm.

  Status-code values and descriptions are defined in:
    Offline/CaloConditions/inc/CalCombinedEnergyCalibStatus.hh

  Author: W. Zhou 2025
*/

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include "cetlib_except/exception.h"
#include "Offline/CaloConditions/inc/CalCombinedEnergyCalibStatus.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"

namespace mu2e {

class CalCombinedEnergyCalib : public DbTable {
 public:
  typedef std::shared_ptr<CalCombinedEnergyCalib> ptr_t;
  typedef std::shared_ptr<const CalCombinedEnergyCalib> cptr_t;

  class Row {
   public:
    Row(CaloSiPMId roid, float ADC2MeV, float ADC2MeV_err,
        int status_code, const std::string& status_message)
        : _roid(roid),
          _ADC2MeV(ADC2MeV),
          _ADC2MeV_err(ADC2MeV_err),
          _status_code(status_code),
          _status_message(status_message) {}

    CaloSiPMId roid() const { return _roid; }
    float ADC2MeV() const { return _ADC2MeV; }
    float ADC2MeV_err() const { return _ADC2MeV_err; }
    int status_code() const { return _status_code; }
    std::string status_message() const { return _status_message; }

   private:
    CaloSiPMId _roid;
    float _ADC2MeV;
    float _ADC2MeV_err;
    int _status_code;
    std::string _status_message;
  };

  constexpr static const char* cxname = "CalCombinedEnergyCalib";

  CalCombinedEnergyCalib()
      : DbTable(cxname, "cal.combinedenergycalib",
                "roid,adc2mev,adc2mev_err,status_code,status_message") {}

  const Row& row(CaloSiPMId roid) const { return _rows.at(roid.id()); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); }
  size_t size() const override { return baseSize() + nrow() * sizeof(Row); }
  virtual std::size_t nrowFix() const override { return CaloConst::_nChannelDB; }
  const std::string orderBy() const { return std::string("roid"); }

  void addRow(const std::vector<std::string>& columns) override {
    std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index >= CaloConst::_nChannelDB || index != _rows.size()) {
      throw cet::exception("CALCOMBINEDENERGYCALIB_BAD_INDEX")
          << "CalCombinedEnergyCalib::addRow found index out of order: "
          << index
          << " != " << _rows.size() << "\n";
    }
    _rows.emplace_back(CaloSiPMId(index), std::stof(columns[1]),
                       std::stof(columns[2]), std::stoi(columns[3]),
                       columns[4]);
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.roid() << ",";
    sstream << std::fixed << std::setprecision(5);
    sstream << r.ADC2MeV() << ",";
    sstream << r.ADC2MeV_err() << ",";
    sstream << r.status_code() << ",";
    sstream << "\"" << r.status_message() << "\"";
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
