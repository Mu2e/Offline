#ifndef DbTables_TrkDelayPanel_hh
#define DbTables_TrkDelayPanel_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class TrkDelayPanel : public DbTable {
 public:
  typedef std::shared_ptr<TrkDelayPanel> ptr_t;
  typedef std::shared_ptr<const TrkDelayPanel> cptr_t;

  class Row {
   public:
    Row(int index, float delay) : _index(index), _delay(delay) {}
    int index() const { return _index; }
    float delay() const { return _delay; }

   private:
    int _index;
    float _delay;
  };

  constexpr static const char* cxname = "TrkDelayPanel";

  TrkDelayPanel() : DbTable(cxname, "trk.delaypanel", "index,delay") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  virtual std::size_t nrowFix() const override { return 216; };
  size_t size() const override { return baseSize() + nrow() * sizeof(Row); };
  const std::string orderBy() const { return std::string("index"); }

  void addRow(const std::vector<std::string>& columns) override {
    int index = std::stoi(columns[0]);
    // enforce a strict sequential order
    if (index != int(_rows.size())) {
      throw cet::exception("TRKDELAYPANEL_BAD_INDEX")
          << "TrkDelayPanel::addRow found index out of order: " << index
          << " != " << _rows.size() << "\n";
    }
    _rows.emplace_back(index, std::stof(columns[1]));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.index() << ",";
    sstream << std::fixed << std::setprecision(1);
    sstream << r.delay();
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
