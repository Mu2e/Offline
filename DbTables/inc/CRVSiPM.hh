#ifndef DbTables_CRVSiPM_hh
#define DbTables_CRVSiPM_hh

#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class CRVSiPM : public DbTable {
 public:
  typedef std::shared_ptr<CRVSiPM> ptr_t;
  typedef std::shared_ptr<const CRVSiPM> cptr_t;

  class Row {
   public:
    Row(std::size_t channel, float pedestal, float height, float area) :
        _channel(channel), _pedestal(pedestal), _height(height), _area(area) {}
    std::size_t channel() const { return _channel; }
    float pedestal() const { return _pedestal; }
    float height() const { return _height; }
    float area() const { return _area; }

   private:
    std::size_t _channel;
    float _pedestal;
    float _height;
    float _area;
  };

  constexpr static const char* cxname = "CRVSiPM";

  CRVSiPM() : DbTable(cxname, "crv.sipm", "channel,pedestal,height,area") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  const Row& row(std::size_t channel) const { return _rows.at(channel); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  size_t size() const override { return baseSize() + nrow() * sizeof(Row); };
  const std::string orderBy() const override { return std::string("channel"); }

  void addRow(const std::vector<std::string>& columns) override {
    std::size_t channel = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (channel >= CRVId::nChannels || channel != _rows.size()) {
      throw cet::exception("CRVSIPM_BAD_CHANNEL")
          << "CRVSiPM::addRow bad channel, saw " << columns[0] << ", expected "
          << _rows.size() << "\n";
    }
    _rows.emplace_back(std::stoi(columns[0]), std::stof(columns[1]),
                       std::stof(columns[2]), std::stof(columns[3]));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.channel() << ",";
    sstream << std::fixed << std::setprecision(3);
    sstream << r.pedestal();
    sstream << r.height();
    sstream << r.area();
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
