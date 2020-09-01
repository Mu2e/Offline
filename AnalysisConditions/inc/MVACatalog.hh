#ifndef AnalysisConditions_MVACatalog_hh
#define AnalysisConditions_MVACatalog_hh

//
// MVACatalog is the ProditionsEntry for MVAToolDb
// This is a templated class that takes the MVAStruct<DETAIL> class
// See AnalysisConditions/inc/TrkQualCatalog.hh for an example use
//

// Mu2e includes
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  template <class T>
  struct MVAEntry {
    MVAEntry(std::string trainName, std::string xmlFileName, bool calibrated = false) : 
      _trainName(trainName), _xmlFileName(xmlFileName), _calibrated(calibrated)
    { }
    ~MVAEntry() {
      delete _mvaTool;
    }

    void initializeMVA() {
      _mvaTool = new MVATools(_xmlFileName);
      _mvaTool->initMVA();

      // create the MVA mask in case we have removed variables
      const auto& labels = _mvaTool->labels();
      _mvaMask = 0;
      for (int i_var = 0; i_var < T::n_vars; ++i_var) {
	for (const auto& i_label : labels) {
	  std::string i_varName = T::varName(static_cast<typename T::MVA_varindex>(i_var));
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

    float getCutVal(float calib_value) const {
      // Get the lower and upper bound of the requestd value
      // NB std::map::lower_bound and std::map::upper_bound return the same iterator 
      // when we search for an interpolated value...
      const auto& upper = _effCalib.upper_bound(calib_value);
      const auto& lower = std::prev(upper); // ...so we have to do this to get the bracketed values

      // In C++20 we will have std::lerp for liner interpolations (a + t(b-a)), for now have to do it ourselves
      float t = ((calib_value - lower->first)/(upper->first - lower->first));
      float result = lower->second + t*(upper->second - lower->second);
      // std::cout << "Request = " << calib_value << std::endl;
      // std::cout << "Lower = " << lower->first << ", " << lower->second << std::endl;
      // std::cout << "Upper = " << upper->first << ", " << upper->second << std::endl;
      // std::cout << "Request - Lower = " << calib_value - lower->first << std::endl;
      // std::cout << "Interpolated Value = " << result << std::endl;

      return result;
    }
    
    std::string _trainName;
    std::string _xmlFileName;
    bool _calibrated;

    MVATools* _mvaTool;
    MVAMask _mvaMask;
    std::map<Float_t, Float_t> _effCalib;
  };

  template <class T>
  using MVAEntries = std::vector<MVAEntry<T> >;

  template <class T>
  class MVACatalog : virtual public ProditionsEntity {

  public:
    MVACatalog() : _name("MVACatalog") { }
    MVACatalog(MVAEntries<T> entries) : _name("MVACatalog"), _entries(entries) { }

    // accessors
    std::string const& name() const { return _name; }
    MVAEntries<T>& modifiableEntries() { return _entries;} 
    MVAEntries<T> const& entries() const { return _entries;} 

    const std::string print() const {
      std::stringstream out;
      print(out);
      return out.str();
    }

    void print(std::ostream& os) const {
      os << "Entries in " << _name << ":" << std::endl;
      for (const auto& i_entry : _entries) {
	os << "Train Name: " << i_entry._trainName << ", XML File: " << i_entry._xmlFileName << std::endl;
      }
      os << std::endl;
    }

    MVAEntry<T> const& find(std::string name) const {
      for (const auto& i_entry : entries()) {
	if (i_entry._trainName == name) {
	  return i_entry;
	}
      }
      // if we get this far, then we haven't found the MVAEntry requested
      throw cet::exception("MVACatalog") << "MVA training with name " << name << " was not found in this MVACatalog: " << std::endl << print();

      return entries().at(0);
    }

    // typedefs
    typedef std::shared_ptr<MVACatalog<T> > ptr_t;
    typedef std::shared_ptr<const MVACatalog<T> > cptr_t;
    
  private:
    // data
    std::string _name;
    MVAEntries<T> _entries;
  };
}

#endif
