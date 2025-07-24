#ifndef DbTables_CalEnergyCalib_hh
#define DbTables_CalEnergyCalib_hh


/*
per SiPM calibration constants reco table -
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

  class CalEnergyCalib : public DbTable {
  public:
  typedef std::shared_ptr<CalEnergyCalib> ptr_t;
  typedef std::shared_ptr<const CalEnergyCalib> cptr_t;

    class Row {
    public:
      Row(CaloSiPMId  roid, float ADC2MeV):_roid(roid),_ADC2MeV(ADC2MeV) {}
      CaloSiPMId   roid() const { return _roid;}
      float ADC2MeV() const { return _ADC2MeV; }

    private:
      CaloSiPMId   _roid;
      float _ADC2MeV;
    };

    constexpr static const char* cxname = "CalEnergyCalib";

    CalEnergyCalib():DbTable(cxname,"cal.energycalib","roid,adc2mev"){}

    const Row& row(CaloSiPMId id) const {
                return _rows[id.id()];
    }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*sizeof(Row); };
    virtual std::size_t nrowFix() const override { return CaloConst::_nChannel; };
    const std::string orderBy() const { return std::string("roid"); }

    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index!=int(_rows.size())) {
        throw cet::exception("CALOENERGYCALIB_BAD_INDEX")<<"CalEnergyCalib::addRow found index out of order:"<<index << " != " << int(_rows.size()) <<"\n";
      }
       _rows.emplace_back(CaloSiPMId(index),std::stof(columns[1]));

    }


    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid()<<",";
      sstream << r.ADC2MeV();
    }

    virtual void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;

  };

}
#endif
