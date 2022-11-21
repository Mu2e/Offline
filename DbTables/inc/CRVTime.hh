#ifndef DbTables_CRVTime_hh
#define DbTables_CRVTime_hh

#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class CRVTime : public DbTable {
 public:
  typedef std::shared_ptr<CRVTime> ptr_t;
  typedef std::shared_ptr<const CRVTime> cptr_t;

  class Row {
   public:
    Row(std::size_t channel, float timeOffset) :
        _channel(channel), _timeOffset(timeOffset) {}
    std::size_t channel() const { return _channel; }
    float timeOffset() const { return _timeOffset; }

   private:
    std::size_t _channel;
    float _timeOffset;
  };

  constexpr static const char* cxname = "CRVTime";

  CRVTime() : DbTable(cxname, "crv.time", "channel,timeOffset") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  const Row& row(std::size_t channel) const { return _rows.at(channel); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  virtual std::size_t nrowFix() const override { return CRVId::nChannels; };
  std::size_t size() const override {
    return baseSize() + nrow() * sizeof(Row);
  };
  const std::string orderBy() const override { return std::string("channel"); }

  void addRow(const std::vector<std::string>& columns) override {
    std::size_t channel = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (channel >= CRVId::nChannels || channel != _rows.size()) {
      throw cet::exception("CRVTIME_BAD_CHANNEL")
          << "CRVTime::addRow bad channel, saw " << columns[0] << ", expected "
          << _rows.size() << "\n";
    }
    _rows.emplace_back(channel, std::stof(columns[1]));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.channel() << ",";
    sstream << std::fixed << std::setprecision(3);
    sstream << r.timeOffset();
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
