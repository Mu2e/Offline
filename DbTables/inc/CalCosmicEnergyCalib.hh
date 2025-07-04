#ifndef DbTables_CalCosmicEnergyCalib_hh
#define DbTables_CalCosmicEnergyCalib_hh


#include <string>
#include <iomanip>
#include <sstream>
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"

namespace mu2e {

  class CalCosmicEnergyCalib : public DbTable {
  public:
    typedef std::shared_ptr<CalCosmicEnergyCalib> ptr_t;
    typedef std::shared_ptr<const CalCosmicEnergyCalib> cptr_t;

    class Row {
    public:
      Row(CaloSiPMId roid, float Peak, float ErrPeak, float Width, float ErrWidth, float Sigma, float ErrSigma, float chisq, int Nhits):
        _roid(roid),
        _Peak(Peak),
        _ErrPeak(ErrPeak),
        _Width(Width),
        _ErrWidth(ErrWidth),
        _Sigma(Sigma),
        _ErrSigma(ErrSigma),
        _chisq(chisq),
        _Nhits(Nhits){}
      CaloSiPMId  roid() const { return _roid;}
      float Peak() const { return _Peak; }
      float ErrPeak() const { return _ErrPeak; }
      float Width() const { return _Width; }
      float ErrWidth() const { return _ErrWidth; }
      float Sigma() const { return _Sigma; }
      float ErrSigma() const { return _ErrSigma; }
      float chisq() const { return _chisq; }
      int Nhits() const { return _Nhits; }

    private:
      CaloSiPMId _roid;
      float _Peak;
      float _ErrPeak;
      float _Width;
      float _ErrWidth;
      float _Sigma;
      float _ErrSigma;
      float _chisq;
      int _Nhits;
    };

    constexpr static const char* cxname = "CalCosmicEnergyCalib";

    CalCosmicEnergyCalib():DbTable(cxname,"cal.cosmicenergycalib","roid,peak,errpeak,width,errwidth,sigma,errsigma,chisq,nhits"){}

    const Row& row(CaloSiPMId  roid) const {
                return _rows.at(roid.id()); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize()  + nrow()*sizeof(Row); };
    virtual std::size_t nrowFix() const override { return CaloConst::_nCrystalChannel; };
    const std::string orderBy() const { return std::string("roid"); }
    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index >= CaloConst::_nChannel  || index != _rows.size()) {
        throw cet::exception("CALOCOSMICCALIB_BAD_INDEX")
        << "CalCosmicEnergyCalib::addRow found index out of order: "
        <<index<< " != " << _rows.size() <<"\n";

      }
   _rows.emplace_back(CaloSiPMId(index),
      std::stof(columns[1]),
      std::stof(columns[2]),
      std::stof(columns[3]),
      std::stof(columns[4]),
      std::stof(columns[5]),
      std::stof(columns[6]),
      std::stof(columns[7]),
      std::stoi(columns[8]));

  }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << r.roid() <<",";
      sstream << std::fixed << std::setprecision(5);
      sstream << r.Peak()<<",";
      sstream << r.ErrPeak()<<",";
      sstream << r.Width()<<",";
      sstream << r.ErrWidth()<<",";
      sstream << r.Sigma()<<",";
      sstream << r.ErrSigma()<<",";
      sstream << r.chisq()<<",";
      sstream << r.Nhits();
    }

    virtual void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;
  };
}
#endif
