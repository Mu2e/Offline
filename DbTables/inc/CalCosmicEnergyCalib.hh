#ifndef DbTables_CaloCosmicsEneTable_hh
#define DbTables_CaloCosmicsEneTable_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "Offline/DbTables/inc/DbTable.hh"

namespace mu2e {

  class CaloCosmicsEneTable : public DbTable {
    public:

    class Row {
      public:
      Row(int roid, double Peak, double ErrPeak, double Width, double ErrWidth, double Chi2):
      _roid(roid),_Peak(Peak),_ErrPeak(ErrPeak),_Width(Width),_ErrWidth(ErrWidth),_Chi2(Chi2) {}
      int       roid()     const { return _roid;} // Offline ID
      float     Peak()     const { return _Peak; }
      float     ErrPeak()  const { return _ErrPeak; }
      float     Width()    const { return _Width; }
      float     ErrWidth() const { return _ErrWidth; }
      float     Chi2()     const { return _Chi2; }

      private:
      int   _roid;
      float _Peak;
      float _ErrPeak;
      float _Width;
      float _ErrWidth;
      float _Chi2;
    };

    constexpr static const char* cxname = "CaloCosmicsEneTable";

    CaloCosmicsEneTable():DbTable(cxname,"calocosmicsenecalib",
    "roid,Peak,ErrPeak,Width,ErrWidth,Chi2") {}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int roid) const { return _rows.at(_chanIndex.at(roid));}
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*nrow()/2 + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) override {
      int roid = std::stoi(columns[0]);
      if(roid!=int(_rows.size())) {
        throw cet::exception("CaloCosmicsEnergyCalib_BAD_INDEX")
        << "CaloCosmicsEneTable::addRow found index out of order: "
        <<roid << " != " << _rows.back().roid()+1 <<"\n";
      }
      _rows.emplace_back(roid,
      std::stoi(columns[1]),
      std::stof(columns[2]),
      std::stoi(columns[3]),
      std::stof(columns[4]),
      std::stof(columns[5]) );
      //_chanIndex[_rows.back().roid()] = _rows.size()-1;
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
    std::map<int,std::size_t> _chanIndex;
  };

};
#endif
