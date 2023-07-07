#ifndef DbTables_TrkDelayRStraw_hh
#define DbTables_TrkDelayRStraw_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class TrkDelayRStraw : public DbTable {
 public:
  typedef std::shared_ptr<TrkDelayRStraw> ptr_t;
  typedef std::shared_ptr<const TrkDelayRStraw> cptr_t;

  class Row {
   public:
    Row(int straw, float delay_hv, float delay_cal) :
        _straw(straw), _delay_hv(delay_hv), _delay_cal(delay_cal) {}
    int straw() const { return _straw; }
    float delayHv() const { return _delay_hv; }
    float delayCal() const { return _delay_cal; }

   private:
    int _straw;
    float _delay_hv;
    float _delay_cal;
  };

  constexpr static const char* cxname = "TrkDelayRStraw";

  TrkDelayRStraw() :
      DbTable(cxname, "trk.delayrstraw", "straw,delay_hv,delay_cal") {}
  const Row& rowAt(const std::size_t straw) const { return _rows.at(straw); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  virtual std::size_t nrowFix() const override { return 96; };
  size_t size() const override { return baseSize() + nrow() * sizeof(Row); };
  const std::string orderBy() const { return std::string("straw"); }

  void addRow(const std::vector<std::string>& columns) override {
    int straw = std::stoi(columns[0]);
    // enforce a strict sequential order
    if (straw != int(_rows.size())) {
      throw cet::exception("TRKDELAYRSTRAW_BAD_INDEX")
          << "TrkDelayRStraw::addRow found straw out of order: " << straw
          << " != " << _rows.back().straw() + 1 << "\n";
    }

    _rows.emplace_back(std::stoi(columns[0]), std::stof(columns[1]),
                       std::stof(columns[2]));
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.straw() << ",";
    sstream << std::fixed << std::setprecision(1);
    sstream << r.delayHv() << ",";
    sstream << r.delayCal();
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
