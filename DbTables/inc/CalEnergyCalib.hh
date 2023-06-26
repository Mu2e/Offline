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
      Row(CaloSiPMId  roid, float ADC2MeV, float ErrADC2MeV, int algName):_roid(roid),_ADC2MeV(ADC2MeV), _ErrADC2MeV(ErrADC2MeV),_algName(algName) {}
      CaloSiPMId   roid() const { return _roid;}
      float ADC2MeV() const { return _ADC2MeV; }
      float ErrADC2MeV() const { return _ErrADC2MeV; }
      int algName() const { return _algName; }

    private:
      CaloSiPMId   _roid;
      float _ADC2MeV;
      float _ErrADC2MeV;
      int _algName;
    };

    constexpr static const char* cxname = "CalEnergyCalib";

    CalEnergyCalib():DbTable(cxname,"cal.energycalib","roid,ADC2MeV,ErrADC2MeV,algName"){}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int roid) const { 
                return _rows.at(roid); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*nrow()/2 + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index!=int(_rows.size())) {
        throw cet::exception("CALOENERGYCALIB_BAD_INDEX")<<"CalEnergyCalib::addRow found index out of order:"<<index << " != " << _rows.back().roid()+1 <<"\n";
      }
       _rows.emplace_back(CaloSiPMId (index),std::stoi(columns[1]),std::stof(columns[2]),std::stof(columns[3]));

    }


    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid()<<",";
      sstream << r.ADC2MeV()<<",";
      sstream << r.ErrADC2MeV()<<",";
      sstream << r.algName();
    }

    virtual void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;

  };
  
};
#endif
