#ifndef DbTables_CRVScint_hh
#define DbTables_CRVScint_hh

#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace mu2e {

class CRVScint : public DbTable {
 public:
  typedef std::shared_ptr<CRVScint> ptr_t;
  typedef std::shared_ptr<const CRVScint> cptr_t;

  class Row {
   public:
    Row(std::uint16_t channel, float scintYieldDeviation) :
        _channel(channel),
        _scintYieldDeviation(scintYieldDeviation) {}
    std::uint16_t channel() const { return _channel; }
    float scintYieldDeviation() const { return _scintYieldDeviation; }

   private:
    std::uint16_t _channel;
    float _scintYieldDeviation;
  };

  constexpr static const char* cxname = "CRVScint";

  CRVScint() :
      DbTable(cxname, "crv.scint", "channel,scintYieldDeviation") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  const Row& row(std::uint16_t channel) const { return _rows.at(channel); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  std::size_t size() const override {
    return baseSize() + nrow() * sizeof(Row);
  };
  const std::string orderBy() const override { return std::string("channel"); }

  void addRow(const std::vector<std::string>& columns) override {
    std::uint16_t channel = std::stoul(columns[0]);
    // enforce order, so channel can be looked up by index
    if (channel >= CRVId::nBars || channel != _rows.size()) {
      throw cet::exception("CRVSCINT_BAD_CHANNEL")
          << "CRVScint::addRow bad channel, saw " << columns[0] << ", expected "
          << _rows.size() << "\n";
    }
    _rows.emplace_back(std::stoi(columns[0]), std::stof(columns[1]));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.channel() << ",";
    sstream << std::fixed << std::setprecision(3);
    sstream << r.scintYieldDeviation();
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
