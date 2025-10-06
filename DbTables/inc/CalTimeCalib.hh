#ifndef DbTables_CalTimeCalib_hh
#define DbTables_CalTimeCalib_hh


/*
per SiPM time calibration constants reco table -
S Middleton 2023

*/

#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "cetlib_except/exception.h"
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"

namespace mu2e {

  class CalTimeCalib : public DbTable {
  public:
  typedef std::shared_ptr<CalTimeCalib> ptr_t;
  typedef std::shared_ptr<const CalTimeCalib> cptr_t;

    class Row {
    public:
      Row(CaloSiPMId  roid, float tcorr):_roid(roid),_tcorr(tcorr) {}
      CaloSiPMId   roid() const { return _roid;}
      float tcorr() const { return _tcorr; } // correction in ns

    private:
      CaloSiPMId   _roid;
      float _tcorr;
    };

    constexpr static const char* cxname = "CalTimeCalib";

    CalTimeCalib():DbTable(cxname,"cal.timecalib","roid,tcorr"){}

    const Row& row(CaloSiPMId id) const {
                return _rows[id.id()];
    }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*sizeof(Row); };
    virtual std::size_t nrowFix() const override { return CaloConst::_nChannelDB; };
    const std::string orderBy() const { return std::string("roid"); }

    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index!=int(_rows.size())) {
        throw cet::exception("CALOTIMECALIB_BAD_INDEX")<<"CalTimeCalib::addRow found index out of order:"<<index << " != " << int(_rows.size()) <<"\n";
      }
       _rows.emplace_back(CaloSiPMId(index),std::stof(columns[1]));

    }


    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid()<<",";
      sstream << r.tcorr();
    }

    virtual void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;

  };

}
#endif
