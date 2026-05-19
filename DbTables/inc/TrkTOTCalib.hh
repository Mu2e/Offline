#ifndef DbTables_TrkTOTCalib_hh
#define DbTables_TrkTOTCalib_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

  class TrkTOTCalib : public DbTable {
    public:
      typedef std::shared_ptr<TrkTOTCalib> ptr_t;
      typedef std::shared_ptr<const TrkTOTCalib> cptr_t;

      class Row {
        public:
          Row(int totTBins, float totTBinWidth, int totEBins, float totEBinWidth, std::vector<float> driftTimes, std::vector<float> driftErrors)
            : _totTBins(totTBins), _totTBinWidth(totTBinWidth),
            _totEBins(totEBins), _totEBinWidth(totEBinWidth),
            _driftTimes(std::move(driftTimes)), _driftErrors(std::move(driftErrors)) {}
          int   totTBins()    const { return _totTBins; }
          float totTBinWidth() const { return _totTBinWidth; }
          int   totEBins()    const { return _totEBins; }
          float totEBinWidth() const { return _totEBinWidth; }
          std::vector<float> const& driftTimes() const { return _driftTimes; }
          std::vector<float> const& driftErrors() const { return _driftErrors; }
        private:
          int   _totTBins;
          float _totTBinWidth;
          int   _totEBins;
          float _totEBinWidth;
          std::vector<float> _driftTimes;
          std::vector<float> _driftErrors;
      };

      constexpr static const char* cxname = "TrkTOTCalib";

      TrkTOTCalib() :
        DbTable(cxname, "trk.totcalib",
            "tottbins,tottbinwidth,totebins,totebinwidth,drifttimes,drifterrors") {}
      const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
      std::vector<Row> const& rows() const { return _rows; }
      std::size_t nrow() const override { return _rows.size(); };
      virtual std::size_t nrowFix() const override { return 1; }
      size_t size() const override { return baseSize() + nrow() * sizeof(Row); };
      const std::string orderBy() const { return std::string("tottbins"); }

      void addRow(const std::vector<std::string>& columns) override {
        if (_rows.size() != 0)
          throw cet::exception("TRKTOTCALIB_BAD_INDEX")
            << "TrkTOTCalib::addRow adding more than one row\n";
        int totTBins = std::stoi(columns[0]);
        int totEBins = std::stoi(columns[2]);
        std::istringstream is(columns[4]); // drift times
        std::vector<float> driftTimes;
        driftTimes.reserve(totTBins*totEBins);
        float x;
        while (is >> x) driftTimes.push_back(x);
        if (!is.eof()){
          throw cet::exception("TRKTOTCALIB_BAD_INDEX")
            << "TrkTOTCalib::addRow bad driftTime\n";
        }
        std::istringstream is2(columns[5]); // drift errors
        std::vector<float> driftErrors;
        driftErrors.reserve(totTBins*totEBins);
        while (is2 >> x) driftErrors.push_back(x);
        if (!is2.eof()){
          throw cet::exception("TRKTOTCALIB_BAD_INDEX")
            << "TrkTOTCalib::addRow bad driftError\n";
        }

        _rows.emplace_back(totTBins, std::stof(columns[1]),
            totEBins, std::stof(columns[3]), driftTimes, driftErrors);
      }

      void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
        Row const& r = _rows.at(irow);
        sstream << std::fixed << std::setprecision(1);
        sstream << r.totTBins() << ",";
        sstream << r.totTBinWidth() << ",";
        sstream << r.totEBins() << ",";
        sstream << r.totEBinWidth() << ",";
        for (size_t k=0;k<r.driftTimes().size();k++){
          if (k) sstream << " ";
          sstream << r.driftTimes()[k];
        }
        sstream << ",";
        for (size_t k=0;k<r.driftErrors().size();k++){
          if (k) sstream << " ";
          sstream << r.driftErrors()[k];
        }
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
