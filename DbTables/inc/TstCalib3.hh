#ifndef DbTables_TstCalib3_hh
#define DbTables_TstCalib3_hh


#include <string>
#include <iomanip>
#include <sstream>
#include <map>
#include "DbTables/inc/DbTable.hh"

namespace mu2e {

  class TstCalib3 : public DbTable {
  public:

    class Row {
    public:
      Row(int channel, float v0, float v1, float v2, float v3, float v4,
	  float v5, float v6, float v7, float v8, float v9):
	_channel(channel),_v0(v0),_v1(v1),_v2(v2),_v3(v3),_v4(v4),
      _v5(v5),_v6(v6),_v7(v7),_v8(v8),_v9(v9){}
      int  channel() const { return _channel;}
      float     v0() const { return _v0; }
      float     v1() const { return _v1; }
      float     v2() const { return _v2; }
      float     v3() const { return _v3; }
      float     v4() const { return _v4; }
      float     v5() const { return _v5; }
      float     v6() const { return _v6; }
      float     v7() const { return _v7; }
      float     v8() const { return _v8; }
      float     v9() const { return _v9; }
    private:
      int _channel;
      float _v0;
      float _v1;
      float _v2;
      float _v3;
      float _v4;
      float _v5;
      float _v6;
      float _v7;
      float _v8;
      float _v9;
    };


    TstCalib3():DbTable("TstCalib3","tst.calib3",
			"channel,v0,v1,v2,v3,v4,v5,v6,v7,v8,v9") {}

    const Row& rowAt(const std::size_t index) const { return _rows.at(index);}
    const Row& row(const int channel) const { 
                return _rows.at(_chanIndex.at(channel)); }
    std::vector<Row> const& rows() const {return _rows;}
    std::size_t nrow() const { return _rows.size(); };
    size_t size() const { return _csv.capacity() 
	+ nrow()*nrow()/2 + nrow()*sizeof(Row); };

    void addRow(const std::vector<std::string>& columns) {
      _rows.emplace_back(std::stoul(columns[0]),
			 std::stof(columns[1]),
			 std::stof(columns[2]),
			 std::stof(columns[3]),
			 std::stof(columns[4]),
			 std::stof(columns[5]),
			 std::stof(columns[6]),
			 std::stof(columns[7]),
			 std::stof(columns[8]),
			 std::stof(columns[9]),
			 std::stof(columns[10]) );
      _chanIndex[_rows.back().channel()] = _rows.size()-1;
    }

    void rowToCsv(std::ostringstream& sstream, std::size_t irow) const {
      Row const& r = _rows.at(irow);
      sstream << r.channel()<<",";
      sstream << std::fixed << std::setprecision(5);
      sstream << r.v0()<<",";
      sstream << r.v1()<<",";
      sstream << r.v2()<<",";
      sstream << r.v3()<<",";
      sstream << r.v4()<<",";
      sstream << r.v5()<<",";
      sstream << r.v6()<<",";
      sstream << r.v7()<<",";
      sstream << r.v8()<<",";
      sstream << r.v9();
    }

    virtual void clear() { _csv.clear(); _rows.clear(); _chanIndex.clear();}

  private:
    std::vector<Row> _rows;
    std::map<int,std::size_t> _chanIndex;
  };
  
};
#endif
