#ifndef DbTables_CalCosmicEnergyCalib_hh
#define DbTables_CalCosmicEnergyCalib_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"

namespace mu2e {

  class CalCosmicEnergyCalib : public DbTable {
  public:

    class Row {
    public:
      Row(CaloSiPMId roid, float EPeak, float ErrEPeak, float Width, float ErrWidth, float chisq):_roid(roid),_EPeak(EPeak), _ErrEPeak(ErrEPeak), _Width(Width), _ErrWidth(ErrWidth), _chisq(chisq){}
      CaloSiPMId  roid() const { return _roid;}
      float EPeak() const { return _EPeak; }
      float ErrEPeak() const { return _ErrEPeak; }
      float Width() const { return _Width; }
      float ErrWidth() const { return _ErrWidth; }
      float chisq() const { return _chisq; }

    private:
      CaloSiPMId  _roid;
      float _EPeak;
      float _ErrEPeak;
      float _Width;
      float _ErrWidth;
      float _chisq;
    };

    constexpr static const char* cxname = "CalCosmicEnergyCalib";

    CalCosmicEnergyCalib():DbTable(cxname,"cal.cosmicenergycalib","roid,EPeak,ErrEPeak,Width,ErrWidth,chisq"){}

    const Row& row(std::uint16_t roid) const { return _rows.at(roid); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize()  + nrow()*sizeof(Row); };
    virtual std::size_t nrowFix() const override { return CaloConst::_nChannel; };
    const std::string orderBy() const { return std::string("roid"); }
    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index >= CaloConst::_nChannel  || index != _rows.size()) {
        throw cet::exception("CALOCOSMICCALIB_BAD_INDEX")
        << "CalCosmicEnergyCalib::addRow found index out of order: "
        <<index<< " != " << _rows.size() <<"\n";

      }
   _rows.emplace_back(CaloSiPMId(index),std::stof(columns[1]),std::stof(columns[2]),std::stof(columns[3]),std::stof(columns[4]),std::stof(columns[5]));

  }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid() <<","; //TODO
      sstream << r.EPeak()<<",";
      sstream << r.ErrEPeak()<<",";
      sstream << r.Width()<<",";
      sstream << r.ErrWidth()<<",";
      sstream << r.chisq();
    }

    virtual void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;
    //std::map<int,std::size_t> _chanIndex;
  };
}
#endif
