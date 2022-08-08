#ifndef DbTables_STMEnergyPar_hh
#define DbTables_STMEnergyPar_hh

#include "Offline/DataProducts/inc/STMChannel.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class STMEnergyPar : public DbTable {
 public:
  typedef std::shared_ptr<STMEnergyPar> ptr_t;
  typedef std::shared_ptr<const STMEnergyPar> cptr_t;

  class Row {
   public:
    Row(STMChannel const& channel, float p0, float p1, float p2) :
        _channel(channel), _p0(p0), _p1(p1), _p2(p2) {}
    STMChannel const& channel() const { return _channel; }
    float p0() const { return _p0; }
    float p1() const { return _p1; }
    float p2() const { return _p2; }

   private:
    STMChannel _channel;
    float _p0;
    float _p1;
    float _p2;
  };

  constexpr static const char* cxname = "STMEnergyPar";

  STMEnergyPar() : DbTable(cxname, "stm.energypar", "channel,p0,p1,p2") {}
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
      throw cet::exception("STMENERGYPAR_BAD_CHANNEL")
          << "STMEnergyPar::addRow called with bad channel name: " << columns[0]
          << "\n";
    }

    _rows.emplace_back(channel, std::stof(columns[1]), std::stof(columns[2]),
                       std::stof(columns[3]));
    // add this channel to the map index
    _chanIndex[channel] = _rows.size() - 1;
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.channel().name() << ",";
    sstream << std::fixed << std::setprecision(5);
    sstream << r.p0() << ",";
    sstream << r.p1() << ",";
    sstream << r.p2();
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
