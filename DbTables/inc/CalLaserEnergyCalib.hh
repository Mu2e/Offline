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
      Row(uint16_t roid, float LAS, float ErrLAS, float chisq):_roid(roid),_LAS(LAS), _ErrLAS(ErrLAS), _chisq(chisq){}
      uint16_t  roid() const { return _roid;}
      float LAS() const { return _LAS; }
      float ErrLAS() const { return _ErrLAS; }
      float chisq() const { return _chisq; }

    private:
      uint16_t  _roid;
      float _LAS;
      float _ErrLAS;
      float _chisq;
    };

    constexpr static const char* cxname = "CalLaserEnergyCalib";

    CalLaserEnergyCalib():DbTable(cxname,"cal.laserenergycalib","roid,LAS,ErrLAS,chisq"){}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int roid) const { return _rows.at(roid); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*nrow()/2 + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index >= CaloConst::_nChannel  || index != _rows.size()) {
        throw cet::exception("CALOLaserEnergyCALIB_BAD_INDEX")
        << "CalLaserEnergyCalib::addRow found index out of order: "
        <<index<< " != " << _rows.back().roid()+1 <<"\n";
      }
       _rows.emplace_back(index,std::stoi(columns[1]),std::stof(columns[2]),std::stof(columns[3]));
      // add this channel to the map index - optional
      //_chanIndex[_rows.back().roid()] = _rows.size()-1;
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
};
#endif
