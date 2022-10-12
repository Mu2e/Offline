#ifndef DbTables_STMPedestals_hh
#define DbTables_STMPedestals_hh

#include "Offline/DataProducts/inc/STMChannel.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class STMPedestals : public DbTable {
 public:
  typedef std::shared_ptr<STMPedestals> ptr_t;
  typedef std::shared_ptr<const STMPedestals> cptr_t;

  class Row {
   public:
    Row(STMChannel const& channel, float pedestal) :
        _channel(channel), _pedestal(pedestal) {}
    STMChannel const& channel() const { return _channel; }
    float pedestal() const { return _pedestal; }

   private:
    STMChannel _channel;
    float _pedestal;
  };

  constexpr static const char* cxname = "STMPedestals";

  STMPedestals() : DbTable(cxname, "stm.pedestals", "channel,pedestal") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  const Row& row(STMChannel channel) const {
    return _rows.at(_chanIndex.at(channel));
  }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  // this table should always be 2 rows
  virtual std::size_t nrowFix() const override { return 2; };
  size_t size() const override {
    return baseSize() + nrow() * nrow() / 2 + nrow() * sizeof(Row);
  };

  void addRow(const std::vector<std::string>& columns) override {
    // enforce channel names
    auto channel = STMChannel::findByName(columns[0]);
    if (!channel.isValid()) {
      throw cet::exception("STMPEDESTALS_BAD_CHANNEL")
          << "STMPedestals::addRow called with bad channel name: " << columns[0]
          << "\n";
    }

    _rows.emplace_back(channel, std::stof(columns[1]));
    // add this channel to the map index
    if (_chanIndex.find(channel) != _chanIndex.end()) {
      throw cet::exception("STMPEDESTALS_DUP_CHANNEL")
          << "STMPedestals::addRow called with duplicate channel name: "
          << columns[0] << "\n";
    }
    _chanIndex[channel] = _rows.size() - 1;
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.channel().name() << ",";
    sstream << std::fixed << std::setprecision(5);
    sstream << r.pedestal();
  }

  virtual void clear() override {
    baseClear();
    _rows.clear();
    _chanIndex.clear();
  }

 private:
  std::vector<Row> _rows;
  std::map<STMChannel, std::size_t> _chanIndex;
};

}  // namespace mu2e
#endif
