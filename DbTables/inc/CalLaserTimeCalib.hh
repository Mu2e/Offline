#ifndef DbTables_CalLaserTimeCalib_hh
#define DbTables_CalLaserTimeCalib_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"

namespace mu2e {

  class CalLaserTimeCalib : public DbTable {
    public:

      class Row {
        public:
        Row(CaloSiPMId  roid, double T0, double ErrT0, double chisq):
        _roid(roid),_T0(T0),_ErrT0(ErrT0),_chisq(chisq) {}
        CaloSiPMId       roid()     const { return _roid;} // Offline ID
        float     T0()     const { return _T0; }
        float     ErrT0()  const { return _ErrT0; }
        float     chisq()     const { return _chisq; }

      private:
        CaloSiPMId _roid;
        float _T0;
        float _ErrT0;
        float _chisq;
    };

    constexpr static const char* cxname = "CalLaserTimeCalib";

    CalLaserTimeCalib():DbTable(cxname,"calolasertimecalib",
    "roid,T0,ErrT0,chisq") {}

    const Row& row(CaloSiPMId  _roid) const {
                return _rows.at(_roid.id()); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*sizeof(Row); };
    virtual std::size_t nrowFix() const override { return CaloConst::_nChannel; };
    const std::string orderBy() const { return std::string("roid"); }
    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index >= CaloConst::_nChannel  || index != _rows.size()) {
        throw cet::exception("CalLaserTimeCalib_BAD_INDEX")
        << "CalLaserTimeTable::addRow found index out of order: "
        <<index<< " != " << _rows.size() <<"\n";
      }
      _rows.emplace_back(CaloSiPMId(index),
      std::stof(columns[1]),
      std::stof(columns[2]),
      std::stof(columns[3]));
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid()<<",";
      sstream << r.T0()<<",";
      sstream << r.ErrT0()<<",";
      sstream << r.chisq();
    }

    virtual void clear() override { baseClear(); _rows.clear();}

    private:
      std::vector<Row> _rows;
      //std::map<int,std::size_t> _chanIndex;
  };

}
#endif
