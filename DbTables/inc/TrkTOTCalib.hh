#ifndef DbTables_TrkTOTCalib_hh
#define DbTables_TrkTOTCalib_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

  class TrkTOTCalibParams : public DbTable {
    public:
      typedef std::shared_ptr<TrkTOTCalibParams> ptr_t;
      typedef std::shared_ptr<const TrkTOTCalibParams> cptr_t;

      class Row {
        public:
          Row(int totTBins, float totTBinWidth, int totEBins, float totEBinWidth)
            : _totTBins(totTBins), _totTBinWidth(totTBinWidth),
            _totEBins(totEBins), _totEBinWidth(totEBinWidth) {}
          int   totTBins()    const { return _totTBins; }
          float totTBinWidth() const { return _totTBinWidth; }
          int   totEBins()    const { return _totEBins; }
          float totEBinWidth() const { return _totEBinWidth; }
        private:
          int   _totTBins;
          float _totTBinWidth;
          int   _totEBins;
          float _totEBinWidth;
      };

      constexpr static const char* cxname = "TrkTOTCalibParams";

      TrkTOTCalibParams() :
        DbTable(cxname, "trk.totcalibparams",
            "tottbins,tottbinwidth,totebins,totebinwidth") {}
      const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
      std::vector<Row> const& rows() const { return _rows; }
      std::size_t nrow() const override { return _rows.size(); };
      virtual std::size_t nrowFix() const override { return 1; }
      size_t size() const override { return baseSize() + nrow() * sizeof(Row); };
      const std::string orderBy() const { return std::string("tottbins"); }

      void addRow(const std::vector<std::string>& columns) override {
        if (_rows.size() != 0)
          throw cet::exception("TRKTOTCALIBPARAMS_BAD_INDEX")
            << "TrkTOTCalibParams::addRow adding more than one row\n";
        _rows.emplace_back(std::stoi(columns[0]), std::stof(columns[1]),
            std::stoi(columns[2]), std::stof(columns[3]));
      }

      void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
        Row const& r = _rows.at(irow);
        sstream << std::fixed << std::setprecision(1);
        sstream << r.totTBins() << ",";
        sstream << r.totTBinWidth() << ",";
        sstream << r.totEBins() << ",";
        sstream << r.totEBinWidth();
      }

      virtual void clear() override {
        baseClear();
        _rows.clear();
      }

    private:
      std::vector<Row> _rows;
  };

  class TrkTOTCalib : public DbTable {
    public:
      typedef std::shared_ptr<TrkTOTCalib> ptr_t;
      typedef std::shared_ptr<const TrkTOTCalib> cptr_t;

      class Row {
        public:
          Row(int index, float driftTime, float driftError)
            : _index(index), _driftTime(driftTime), _driftError(driftError) {}
          int   index()      const { return _index; }
          float driftTime()  const { return _driftTime; }
          float driftError() const { return _driftError; }
        private:
          int   _index;
          float _driftTime;
          float _driftError;
      };

      constexpr static const char* cxname = "TrkTOTCalib";

      TrkTOTCalib() :
        DbTable(cxname, "trk.totcalib",
            "index,drifttime,drifterror") {}
      const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
      std::vector<Row> const& rows() const { return _rows; }
      std::size_t nrow() const override { return _rows.size(); };
      size_t size() const override { return baseSize() + nrow() * sizeof(Row); };
      const std::string orderBy() const { return std::string("index"); }

      void addRow(const std::vector<std::string>& columns) override {
        int index = std::stoi(columns[0]);
        // enforce a strict sequential order
        if (index != int(_rows.size())) {
          throw cet::exception("TRKTOTCALIB_BAD_INDEX")
            << "TrkTOTCalib::addRow found index out of order: " << index
            << " != " << _rows.size() << "\n";
        }
        _rows.emplace_back(index,std::stof(columns[1]), std::stof(columns[2]));
      }

      void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
        Row const& r = _rows.at(irow);
        sstream << std::fixed << std::setprecision(1);
        sstream << r.driftTime() << ",";
        sstream << r.driftError();
      }

      virtual void clear() override {
        baseClear();
        _rows.clear();
      }

    private:
      std::vector<Row> _rows;
  };

};  // namespace mu2e
#endif
