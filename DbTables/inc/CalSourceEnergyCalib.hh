#ifndef DbTables_CalSourceEnergyCalib_hh
#define DbTables_CalSourceEnergyCalib_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"

namespace mu2e {

  class CalSourceEnergyCalib : public DbTable {
  public:

    class Row {
    public:
      Row(CaloSiPMId roid, float fullEPeak,float fullErrEPeak,float fullWidth,float fullErrWidth,float firstescEPeak,float firstescErrEPeak,float firstescWidth,float firstescErrWidth,float secescEPeak,float secescErrEPeak,float secescWidth,float secescErrWidth, float frFull, float frFirst, float frSecond, float chisq, int ndf): _roid(roid), _fullEPeak(fullEPeak), _fullErrEPeak(fullErrEPeak), _fullWidth(fullWidth), _fullErrWidth(fullErrWidth),_firstescEPeak(firstescEPeak), _firstescErrEPeak(firstescErrEPeak), _firstescWidth(firstescWidth), _firstescErrWidth(firstescErrWidth),_secescEPeak(secescEPeak), _secescErrEPeak(secescErrEPeak), _secescWidth(secescWidth), _secescErrWidth(secescErrWidth), _frFull(frFull), _frFirst(frFirst), _frSecond(frSecond), _chisq(chisq), _ndf(ndf){}

      CaloSiPMId  roid() const { return _roid;}
      float fullEPeak() const { return _fullEPeak; }
      float fullErrEPeak() const { return _fullErrEPeak; }
      float fullWidth() const { return _fullWidth; }
      float fullErrWidth() const { return _fullErrWidth; }

      float firstescEPeak() const { return _firstescEPeak; }
      float firstescErrEPeak() const { return _firstescErrEPeak; }
      float firstescWidth() const { return _firstescWidth; }
      float firstescErrWidth() const { return _firstescErrWidth; }

      float secescEPeak() const { return _secescEPeak; }
      float secescErrEPeak() const { return _secescErrEPeak; }
      float secescWidth() const { return _secescWidth; }
      float secescErrWidth() const { return _secescErrWidth; }

      float frFull() const { return _frFull; }
      float frFirst() const { return _frFirst; }
      float frSecond() const { return _frSecond; }

      float chisq() const { return _chisq; }
      int ndf() const { return _ndf; }

    private:
      CaloSiPMId  _roid;
      float _fullEPeak;
      float _fullErrEPeak;
      float _fullWidth;
      float _fullErrWidth;
      float _firstescEPeak;
      float _firstescErrEPeak;
      float _firstescWidth;
      float _firstescErrWidth;
      float _secescEPeak;
      float _secescErrEPeak;
      float _secescWidth;
      float _secescErrWidth;
      float _frFull;
      float _frFirst;
      float _frSecond;
      float _chisq;
      int _ndf;
    };

    constexpr static const char* cxname = "CalSourceEnergyCalib";

    CalSourceEnergyCalib():DbTable(cxname,"cal.sourceenergycalib","roid,fullepeak,fullerrepeak,fullwidth,fullerrwidth,firstescepeak,firstescerrepeak,firstescwidth,firstescerrwidth,secescepeak,secescerrepeak,secescwidth,secescerrwidth,frfull,frfirst,frsecond,chisq,ndf"){}

    const Row& row(CaloSiPMId  roid) const {
                return _rows.at(roid.id()); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*sizeof(Row); };
    virtual std::size_t nrowFix() const override { return CaloConst::_nChannelDB; };
    const std::string orderBy() const { return std::string("roid"); }
    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index >= CaloConst::_nChannelDB  || index != _rows.size()) {
      throw cet::exception("CALOSOURCECALIB_BAD_INDEX")
      << "CalSourceEnergyCalib::addRow found index out of order: "
      <<index << " != " << _rows.size() <<"\n";
      }
       _rows.emplace_back(CaloSiPMId(index),std::stof(columns[1]),std::stof(columns[2]),std::stof(columns[3]),
       std::stof(columns[4]),std::stof(columns[5]),std::stof(columns[6]),std::stof(columns[7]),
       std::stof(columns[8]),std::stof(columns[9]),std::stof(columns[10]),std::stof(columns[11]),
       std::stof(columns[12]),std::stof(columns[13]),std::stof(columns[14]),std::stof(columns[15]),std::stof(columns[16]),std::stoi(columns[17]));

    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid()<<",";
      sstream << r.fullEPeak()<<",";
      sstream << r.fullErrEPeak()<<",";
      sstream << r.fullWidth()<<",";
      sstream << r.fullErrWidth()<<",";
      sstream << r.firstescEPeak()<<",";
      sstream << r.firstescErrEPeak()<<",";
      sstream << r.firstescWidth()<<",";
      sstream << r.firstescErrWidth()<<",";
      sstream << r.secescEPeak()<<",";
      sstream << r.secescErrEPeak()<<",";
      sstream << r.secescWidth()<<",";
      sstream << r.secescErrWidth()<<",";
      sstream << r.frFull()<<",";
      sstream << r.frFirst()<<",";
      sstream << r.frSecond()<<",";
      sstream << r.chisq()<<",";
      sstream << r.ndf();
    }

    virtual void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;
  };

}
#endif
