//
//  A base class for classes which print different types of products
//  It serves to simplify the module code a bit
//
#ifndef Print_inc_ProductPrinter_hh
#define Print_inc_ProductPrinter_hh

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "CLHEP/Matrix/SymMatrix.h"
#include <boost/io/ios_state.hpp>
#include <iomanip>
#include <ios>
#include <ostream>

namespace mu2e {

class ProductPrinter {
 public:
  typedef std::vector<art::InputTag> vectag;

  // simple config common to all products
  struct Config {
    fhicl::Atom<int> verbose{fhicl::Name("verbose"),
                             fhicl::Comment("verbose flag, 0 to 1 or 2"), 0};
    fhicl::Sequence<art::InputTag> inputTags{
        fhicl::Name("inputTags"), fhicl::Comment("tags of products to print"),
        std::vector<art::InputTag>()};
  };

  // add one parameter, ecut
  struct ConfigE : public Config {
    fhicl::Atom<float> eCut{fhicl::Name("eCut"),
                            fhicl::Comment("minimum energy cut"), -1};
  };

  ProductPrinter() : _verbose(1) {}
  ProductPrinter(const Config& conf) {
    _verbose = conf.verbose();
    _tags = conf.inputTags();
  }
  virtual ~ProductPrinter() {}

  // 0 = none, 1 (default), 2, ..
  void setVerbose(int i) { _verbose = i; }
  int verbose() const { return _verbose; }
  void setTags(const vectag& tags) { _tags = tags; }
  const vectag& tags() const { return _tags; }

  virtual void Print(art::Event const& event, std::ostream& os = std::cout) {}

  virtual void PrintSubRun(art::SubRun const& subrun,
                           std::ostream& os = std::cout) {}

  virtual void PrintEndJob(std::ostream& os = std::cout) {}

  void PrintMatrix(const CLHEP::HepSymMatrix& matrix, std::ostream& os,
                   int mode = 0) {
    // when this destructs, it restores the flag state
    boost::io::ios_flags_saver ifs(os);
    // print fixed or scientific
    os.flags(std::ios::right);
    for (int r = 1; r <= matrix.num_row(); r++) {
      os << "   ";
      for (int c = 1; c <= matrix.num_col(); c++) {
        if (mode == 0) {
          os << " " << std::setw(13) << std::setprecision(6) << matrix(r, c);
        } else {
          double value = sqrt(fabs(matrix(r, r) * matrix(c, c)));
          if (value != 0.0) value = matrix(r, c) / value;
          os << " " << std::setw(13) << std::setprecision(6) << value;
        }
      }
      os << "\n";
    }
  }

 private:
  int _verbose;
  vectag _tags;
};

}  // namespace mu2e
#endif
