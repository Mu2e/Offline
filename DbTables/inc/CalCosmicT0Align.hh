#ifndef DbTables_CalCosmicT0Align_hh
#define DbTables_CalCosmicT0Align_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "Offline/DbTables/inc/DbTable.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"

namespace mu2e {

  class CalCosmicT0Align : public DbTable {
  public:

    class Row {
    public:
      Row(CaloSiPMId  roid, float t0val, float t0err, float t0width,
          float chi2, int nev):
        _roid(roid), _t0val(t0val), _t0err(t0err), _t0width(t0width),
        _chi2(chi2), _nev(nev) {}
      CaloSiPMId   roid() const { return _roid;}
      float t0val() const { return _t0val; }
      float t0err() const { return _t0err; }
      float t0width() const { return _t0width; }
      float chi2() const { return _chi2; }
      float nev() const { return _nev; }

    private:
      CaloSiPMId   _roid;
      float _t0val;
      float _t0err;
      float _t0width;
      float _chi2;
      int _nev;
    };

    constexpr static const char* cxname = "CalCosmicT0Align";

    CalCosmicT0Align():DbTable(cxname,"cal.cosmict0align",
                               "roid,t0val,t0err,t0width,chi2,nev"){}

    const Row& row(CaloSiPMId  roid) const {
                return _rows.at(roid.id()); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*sizeof(Row); };
    std::size_t nrowFix() const override { return CaloConst::_nChannelDB; };
    const std::string orderBy() const override { return std::string("roid"); }
    void addRow(const std::vector<std::string>& columns) override {
      std::uint16_t index = std::stoul(columns[0]);
    // enforce order, so channels can be looked up by index
    if (index >= CaloConst::_nChannelDB  || index != _rows.size()) {
        throw cet::exception("CALCOSMICT0ALIGN_BAD_INDEX")
        << "CalCosmicT0Align::addRow found index out of order: "
        <<index<< " != " << _rows.size() <<"\n";
      }

    _rows.emplace_back(CaloSiPMId(index), std::stof(columns[1]),
                       std::stof(columns[2]), std::stof(columns[3]),
                       std::stof(columns[4]), std::stoi(columns[5]) );

    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid() <<",";
      sstream << r.t0val()<<",";
      sstream << r.t0err() <<",";
      sstream << r.t0width() <<",";
      sstream << r.chi2()<<",";
      sstream << r.nev();
    }

    void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;
  };
}
#endif
