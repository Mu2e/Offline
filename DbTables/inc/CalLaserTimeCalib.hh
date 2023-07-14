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
        Row(CaloSiPMId  roid, double Peak, double ErrPeak, double Width, double ErrWidth, double Chi2):
        _roid(roid),_Peak(Peak),_ErrPeak(ErrPeak),_Width(Width),_ErrWidth(ErrWidth),_Chi2(Chi2) {}
        CaloSiPMId       roid()     const { return _roid;} // Offline ID
        float     Peak()     const { return _Peak; }
        float     ErrPeak()  const { return _ErrPeak; }
        float     Width()    const { return _Width; }
        float     ErrWidth() const { return _ErrWidth; }
        float     Chi2()     const { return _Chi2; }

      private:
        CaloSiPMId    _roid;
        float _Peak;
        float _ErrPeak;
        float _Width;
        float _ErrWidth;
        float _Chi2;
    };

    constexpr static const char* cxname = "CalLaserTimeCalib";

    CalLaserTimeCalib():DbTable(cxname,"calolasertimecalib",
    "roid,Peak,ErrPeak,Width,ErrWidth,Chi2") {}

    const Row& row(const std::uint16_t roid) const { return _rows.at(roid); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*sizeof(Row); };
    virtual std::size_t nrowFix() const override { return CaloConst::_nChannel; };
    const std::string orderBy() const { return std::string("roid"); }
    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index >= CaloConst::_nChannel  || index != _rows.size()) {
        throw cet::exception("CaloLaserTimeCalib_BAD_INDEX")
        << "CaloLaserTimeTable::addRow found index out of order: "
        <<index<< " != " << _rows.size() <<"\n";
      }
      _rows.emplace_back(CaloSiPMId(index),
      std::stof(columns[1]),
      std::stof(columns[2]),
      std::stoi(columns[3]),
      std::stof(columns[4]),
      std::stof(columns[5]) );
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid()<<",";
      sstream << r.Peak()<<",";
      sstream << r.ErrPeak()<<",";
      sstream << r.Width()<<",";
      sstream << r.ErrWidth()<<",";
      sstream << r.Chi2();
    }

    virtual void clear() override { baseClear(); _rows.clear();}

    private:
      std::vector<Row> _rows;
      //std::map<int,std::size_t> _chanIndex;
  };

}
#endif
