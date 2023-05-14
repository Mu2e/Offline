#ifndef DbTables_CalSourceEnergyCalib_hh
#define DbTables_CalSourceEnergyCalib_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "Offline/DbTables/inc/DbTable.hh"

namespace mu2e {

  class CalSourceEnergyCalib : public DbTable {
  public:

    class Row {
    public:
      Row(uint16_t roid, float EPeak, float ErrEPeak, float Width, float ErrWidth, float chisq):_roid(roid),_EPeak(EPeak), _ErrEPeak(ErrEPeak), _Width(Width), _ErrWidth(ErrWidth), _chisq(chisq){}
      uint16_t  roid() const { return _roid;} 
      float EPeak() const { return _EPeak; }
      float ErrEPeak() const { return _ErrEPeak; }
      float Width() const { return _Width; }
      float ErrWidth() const { return _ErrWidth; }
      float chisq() const { return _chisq; }

    private:
      uint16_t  _roid; 
      float _EPeak;
      float _ErrEPeak;
      float _Width;
      float _ErrWidth;
      float _chisq;
    };

    constexpr static const char* cxname = "CalSourceEnergyCalib";

    CalSourceEnergyCalib():DbTable(cxname,"cal.sourcecalib","roid,EPeak,ErrEPeak,Width,ErrWidth,chisq"){}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int roid) const { 
                return _rows.at(roid); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const override { return _rows.size(); };
    size_t size() const override { return baseSize() + nrow()*nrow()/2 + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) override {
      int roid = std::stoi(columns[0]);
      // enforce a strict sequential order - optional
      if(roid!=int(_rows.size())) {
	      throw cet::exception("CALOSOURCECALIB_BAD_INDEX") 
	        << "CalSourceEnergyCalib::addRow found index out of order: " 
	        <<roid << " != " << _rows.back().roid()+1 <<"\n";
      }
       _rows.emplace_back(roid,std::stoi(columns[1]),std::stof(columns[2]),std::stof(columns[3]),std::stof(columns[4]),std::stof(columns[5]));
      // add this channel to the map index - optional
      //_chanIndex[_rows.back().roid()] = _rows.size()-1;
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
      Row const& r = _rows.at(irow);
      sstream << std::fixed << std::setprecision(5);
      sstream << r.roid()<<",";
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
  
};
#endif
