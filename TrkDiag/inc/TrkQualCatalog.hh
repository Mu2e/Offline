#ifndef TrkDiag_TrkQualCatalog_hh
#define TrkDiag_TrkQualCatalog_hh

//
// TrkQualCatalog is the ProditionsEntry for TrkQualDb
//

// Mu2e includes
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  struct TrkQualEntry {
    TrkQualEntry(std::string trainName, std::string xmlFileName, bool calibrated = false) : 
      _trainName(trainName), _xmlFileName(xmlFileName), _calibrated(calibrated)
    { }
    ~TrkQualEntry() {
      delete _mvaTool;
    }

    void initializeMVA() {
      _mvaTool = new MVATools(_xmlFileName);
      _mvaTool->initMVA();

      // create the MVA mask in case we have removed variables
      const auto& labels = _mvaTool->labels();
      _mvaMask = 0;
      for (int i_var = 0; i_var < TrkQual::n_vars; ++i_var) {
	for (const auto& i_label : labels) {
	  std::string i_varName = TrkQual::varName(static_cast<TrkQual::MVA_varindex>(i_var));
	  if (i_label.find(i_varName) != std::string::npos) {
	    _mvaMask ^= (1 << i_var);
	    break;
	  }
	}
      }

      if (_calibrated) {
	_mvaTool->getCalib(_effCalib);
      }
    }
    
    std::string _trainName;
    std::string _xmlFileName;
    bool _calibrated;

    MVATools* _mvaTool;
    MVAMask _mvaMask;
    std::map<Float_t, Float_t> _effCalib;
  };
  typedef std::vector<TrkQualEntry> TrkQualEntries;

  class TrkQualCatalog : virtual public ProditionsEntity {

  public:
    typedef std::shared_ptr<TrkQualCatalog> ptr_t;
    typedef std::shared_ptr<const TrkQualCatalog> cptr_t;

    TrkQualCatalog() : _name("TrkQualCatalog") {}
    TrkQualCatalog(TrkQualEntries entries) : _name("TrkQualCatalog"), _entries(entries) {}

    // accessors
    void print(std::ostream& os) const;
    std::string const& name() const { return _name; }
    TrkQualEntries& modifiableEntries() { return _entries;} 
    TrkQualEntries const& entries() const { return _entries;} 

  private:
    // data
    std::string _name;
    TrkQualEntries _entries;
  };
}

#endif
