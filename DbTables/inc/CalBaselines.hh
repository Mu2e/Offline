#ifndef DbTables_CalBaselines_hh
#define DbTables_CalBaselines_hh

#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "cetlib_except/exception.h"
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"

namespace mu2e {

  class CalBaselines : public DbTable {
  public:
  typedef std::shared_ptr<CalBaselines> ptr_t;
  typedef std::shared_ptr<const CalBaselines> cptr_t;

    class Row {
    public:
      Row(CaloSiPMId  roid, float baseline, float threshold):_roid(roid),_baseline(baseline),_threshold(threshold) {}
      CaloSiPMId   roid() const { return _roid;}
      float baseline() const { return _baseline; }
      float threshold() const { return _threshold; }

    private:
      CaloSiPMId _roid;
      float _baseline;
      float _threshold;
    };

    constexpr static const char* cxname = "CalBaselines";

    CalBaselines():DbTable(cxname,"cal.baselines","roid,baseline,threshold"){}

    const Row& row(CaloSiPMId id) const {
                return _rows.at(rawid.id());
    }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*sizeof(Row); };
    virtual std::size_t nrowFix() const override { return CaloConst::_nChannelDB; };
    const std::string orderBy() const { return std::string("roid"); }

    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index >= CaloConst::_nChannelDB || index != _rows.size()) {
        throw cet::exception("CALOBASELINES_BAD_INDEX")<<"CalBaselines::addRow found index out of order:"<<index << " != " << int(_rows.size()) <<"\n";
      }
       _rows.emplace_back(CaloSiPMId(index),
       std::stof(columns[1]),
       std::stof(columns[2]));
    }


    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid()<<",";
      sstream << r.baseline()<<",";
      sstream << r.threshold();
    }

    virtual void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;

  };

}
#endif
