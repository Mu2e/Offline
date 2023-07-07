#ifndef DbTables_SimEfficiencies2_hh
#define DbTables_SimEfficiencies2_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class SimEfficiencies2 : public DbTable {
 public:
  typedef std::shared_ptr<SimEfficiencies2> ptr_t;
  typedef std::shared_ptr<const SimEfficiencies2> cptr_t;

  class Row {
   public:
    Row(std::string tag, unsigned long numerator, unsigned long denominator,
        double eff) :
        _tag(tag),
        _numerator(numerator), _denominator(denominator), _eff(eff) {}

    std::string tag() const { return _tag; }
    unsigned long long numerator() const { return _numerator; }
    unsigned long long denominator() const { return _denominator; }
    double eff() const { return _eff; }

   private:
    std::string _tag;
    unsigned long long _numerator;
    unsigned long long _denominator;
    double _eff;
  };

  constexpr static const char* cxname = "SimEfficiencies2";

  SimEfficiencies2() :
      DbTable(cxname, "sim.efficiencies2", "tag,numerator,denominator,eff") {}

  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  const Row& row(const int idx) const { return _rows.at(idx); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  //    virtual std::size_t nrowFix() const { return 3; };
  size_t size() const override { return baseSize() + nrow() * sizeof(Row); };

  void addRow(const std::vector<std::string>& columns) override {
    //      int idx = std::stoi(columns[0]);
    _rows.emplace_back(columns[0], std::stoull(columns[1]),
                       std::stoull(columns[2]), std::stof(columns[3]));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.tag() << ",";
    sstream << r.numerator() << ",";
    sstream << r.denominator() << ",";
    sstream << std::fixed << std::setprecision(6) << r.eff();
  }

  void findEff(std::string name, double& eff) const {
    for (const auto& i_row : _rows) {
      if (i_row.tag() == name) {
        eff = i_row.eff();
        return;
      }
    }
    throw cet::exception("SIMEFFICIENCIES2_BAD_TAG")
        << "Efficiency with tag " << name << " not found in database"
        << std::endl;
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
