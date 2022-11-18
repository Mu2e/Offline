#ifndef DbTables_CRVBadChan_hh
#define DbTables_CRVBadChan_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class CRVBadChan : public DbTable {
 public:
  typedef std::shared_ptr<CRVBadChan> ptr_t;
  typedef std::shared_ptr<const CRVBadChan> cptr_t;

  class Row {
   public:
    Row(std::size_t channel, int status) : _channel(channel), _status(status) {}
    std::size_t channel() const { return _channel; }
    int status() const { return _status; }

   private:
    std::size_t _channel;
    int _status;
  };

  constexpr static const char* cxname = "CRVBadChan";

  CRVBadChan() : DbTable(cxname, "crv.badchan", "channel,status") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  std::size_t size() const override {
    return baseSize() + nrow() * sizeof(Row);
  };
  const std::string orderBy() const override { return std::string("channel"); }

  void addRow(const std::vector<std::string>& columns) override {
    _rows.emplace_back(std::stoul(columns[0]), std::stoi(columns[1]));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.channel() << ",";
    sstream << r.status();
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
