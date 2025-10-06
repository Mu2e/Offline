#ifndef DbTables_CalLaserEnergyCalib_hh
#define DbTables_CalLaserEnergyCalib_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"

namespace mu2e {

  class CalLaserEnergyCalib : public DbTable {
  public:

    class Row {
    public:
      Row(CaloSiPMId  roid, float LAS, float ErrLAS, float chisq):_roid(roid),_LAS(LAS), _ErrLAS(ErrLAS), _chisq(chisq){}
      CaloSiPMId   roid() const { return _roid;}
      float LAS() const { return _LAS; }
      float ErrLAS() const { return _ErrLAS; }
      float chisq() const { return _chisq; }

    private:
      CaloSiPMId   _roid;
      float _LAS;
      float _ErrLAS;
      float _chisq;
    };

    constexpr static const char* cxname = "CalLaserEnergyCalib";

    CalLaserEnergyCalib():DbTable(cxname,"cal.laserenergycalib","roid,las,errlas,chisq"){}

    const Row& row(CaloSiPMId  roid) const {
                return _rows.at(roid.id()); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize()  + nrow()*sizeof(Row); };
    virtual std::size_t nrowFix() const override { return CaloConst::_nChannelDB; };
    const std::string orderBy() const { return std::string("roid"); }
    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index >= CaloConst::_nChannelDB  || index != _rows.size()) {
        throw cet::exception("CALOLaserEnergyCALIB_BAD_INDEX")
        << "CalLaserEnergyCalib::addRow found index out of order: "
        <<index<< " != " <<  _rows.size() <<"\n";
      }
       _rows.emplace_back(CaloSiPMId(index),std::stof(columns[1]),std::stof(columns[2]),std::stof(columns[3]));

    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid()<<",";
      sstream << r.LAS()<<",";
      sstream << r.ErrLAS()<<",";
      sstream << r.chisq();
    }

    virtual void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;
    //std::map<int,std::size_t> _chanIndex;
  };
}
#endif
