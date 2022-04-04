#ifndef DbTables_CalSourceCalibTable_hh
#define DbTables_CalSourceCalibTable_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "Offline/DbTables/inc/DbTable.hh"

namespace mu2e {

  class CalSourceCalibTable : public DbTable {
  public:

    class Row {
    public:
      Row(uint16_t roid, float EPeak, float Esource, float chisq):
	_roid(roid),_EPeak(EPeak),_Esource(Esource), _chisq(chisq) {}
      uint16_t  roid() const { return _roid;} 
      float     EPeak() const { return _EPeak; }
      float     Esource() const { return _Esource; }
      float     chisq() const { return _chisq; }

    private:
      uint16_t  _roid; 
      float _EPeak;
      float _Esource;
      float _chisq;
    };

    constexpr static const char* cxname = "CalSourceCalibTable";

    CalSourceCalibTable():DbTable(cxname,"cal.sourcecalib",
			"roid,EPeak,Esource,chisq") {}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int roid) const { 
                return _rows.at(_chanIndex.at(roid)); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*nrow()/2 + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) override {
      int roid = std::stoi(columns[0]);
      // enforce a strict sequential order - optional
      if(roid!=int(_rows.size())) {
	      throw cet::exception("CALOSOURCECALIB_BAD_INDEX") 
	        << "CalSourceCalibTable::addRow found index out of order: " 
	        <<roid << " != " << _rows.back().roid()+1 <<"\n";
      }
       _rows.emplace_back(roid,
			 std::stoi(columns[1]),
			 std::stof(columns[2]),
			 std::stof(columns[3]) );
      // add this channel to the map index - optional
      _chanIndex[_rows.back().roid()] = _rows.size()-1;
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid()<<",";
      sstream << r.EPeak()<<",";
      sstream << r.Esource()<<",";
      sstream << r.chisq();
    }

    virtual void clear() override { baseClear(); _rows.clear();}

  private:
    std::vector<Row> _rows;
    std::map<int,std::size_t> _chanIndex;
  };
  
};
#endif
