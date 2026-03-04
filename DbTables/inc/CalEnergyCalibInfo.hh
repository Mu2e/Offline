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
    std::string comment = columns[2];
    for (std::size_t i = 3; i < columns.size(); ++i) {
      comment += ",";
      comment += columns[i];
    }
    _rows.emplace_back(std::stoi(columns[0]), std::stoi(columns[1]),
                       fromCsvText(comment));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.energyCalibCid() << ",";
    sstream << r.combinedEnergyCalibCid() << ",";
    sstream << toCsvText(r.comment());
  }

  virtual void clear() override {
    baseClear();
    _rows.clear();
  }

 private:
  static std::string toCsvText(const std::string& text) {
    if (text.find_first_of(",\"") == std::string::npos) return text;
    std::string out;
    out.reserve(text.size() + 2);
    out.push_back('"');
    for (char c : text) {
      if (c == '"') out += "\"\"";
      else out.push_back(c);
    }
    out.push_back('"');
    return out;
  }

  static std::string fromCsvText(const std::string& text) {
    if (text.size() < 2 || text.front() != '"' || text.back() != '"') {
      return text;
    }
    std::string out;
    out.reserve(text.size() - 2);
    for (std::size_t i = 1; i + 1 < text.size(); ++i) {
      if (text[i] == '"' && i + 2 < text.size() && text[i + 1] == '"') {
        out.push_back('"');
        ++i;
      } else if (text[i] == '\\' && i + 2 < text.size() &&
                 text[i + 1] == '"') {
        out.push_back('"');
        ++i;
      } else {
        out.push_back(text[i]);
      }
    }
    return out;
  }

  std::vector<Row> _rows;
};

}  // namespace mu2e
#endif
