#ifndef GeneralUtilities_MVAStruct_hh
#define GeneralUtilities_MVAStruct_hh
//
// Template used to instantiate an MVA Structure.  Stolen heavily from Rob Kutschke's BitMap class.  D.Brown, LBNL 1/30/2017
//
// The user must supply a detail class with the following requirements:
//
//   1) The detail class must contain an enum named MVA_varindex
//      The legal values of the enums are indices into the array of input values.
//      The values must be contiguous, start at 0, and the last entry must be named n_vars
//
//   2) The detail class must contain these static functions:
//       static std::string const& typeName();
//       static std::map<std::string,MVA_varindex> const& varNames();
//
//      The first function returns one string that holds the name of detail class, usually the same name as the .hh file.
//      The template uses it decorating some printed output.
//
//      The second function returns an std::map that implements the cross reference between each value and
//      its string representation.
//
#include <string>
#include <vector>
#include <array>
#include <map>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace mu2e {

// why is inheritence needed? to provide access to protected members?
  template < class DETAIL > class MVAStruct : public DETAIL {

  public:

    typedef typename DETAIL::MVA_varindex MVA_varindex;
    typedef typename DETAIL::map_type map_type;
//    typedef std::array<double,DETAIL::n_vars> vcon_type;
    typedef std::vector<double> vcon_type;

    enum MVA_status {unset=0,calculated,failed}; 
    // default Constructor 
    explicit MVAStruct() : _values(DETAIL::n_vars,0.0){
      reset();
    // check for self-consistency
      if(DETAIL::n_vars != DETAIL::varNames().size()){
	std::ostringstream os;
	throw std::out_of_range( os.str() );
      }
    }

    // Accept compiler generated d'tor, copy c'tor and copy assignment, equivalence, etc.
    //
    // indexed accessor
    double varValue(MVA_varindex vind) const { return _values[vind]; }
    // named accessor
    double varValue(std::string vname) const {
      MVA_varindex vind = findIndexByNameOrThrow(vname);
      return varValue(vind);
    }
    double operator [] (size_t index) const {
      return _values[index];
    }
    double& operator [] (size_t index) {
      return _values[index];
    }

    // struct accessor
    const vcon_type& values() const { return _values; }
    // output accessor
    double MVAOutput() const { return _mvaout; }
    MVA_status status() const { return _status; }
    // allow setting values.  This should really be in a dedicated
    // interface associated with running the MVA FIXME!!
    void setValue(MVA_varindex vind,double val) { _values[vind] = val; }
    void setValue(std::string vname,double val) {
      MVA_varindex vind = findIndexByNameOrThrow(vname);
      setValue(vind,val);
    }
    void setMVAValue(double value) { _mvaout = value; }

    void reset(){
      for(size_t ivar=0;ivar < DETAIL::n_vars; ++ivar)
	_values[ivar] = 0.0;
      _mvaout = -1.0;
      _status = unset;
    }

    // The number of input variables
    static size_t size(){ return DETAIL::varNames().size(); }

   static map_type const& varNames() { return DETAIL::varNames(); }

  private:
    // array of MVA input values, indexed as appropriate
    vcon_type _values;
    // flag whether the variable is used in the MVA calculation or not.  Needed?  FIXME!
//    std::array<bool,DETAIL::n_vars> _used;
    // MVA output
    double _mvaout;
    MVA_status _status; // status of MVA calclulation

    MVA_varindex findIndexByNameOrThrow( std::string const& name ) const {
      typename map_type::const_iterator j = varNames().find(name);
      if ( j == varNames().end() ){
        std::ostringstream os;
        os << DETAIL::varName() << " invalid variable name : " << name;
        throw std::out_of_range( os.str() );
      }
      return j->second;
    }

  };

  template < class DETAIL >
  inline
  std::ostream& operator<<(std::ostream& ost,
                           const MVAStruct<DETAIL> & mva ){
    ost << DETAIL::typeName() << " has the following content: ";
    typename DETAIL::map_type const& varNames = DETAIL::varNames();
    for ( typename DETAIL::map_type::const_iterator i=varNames.begin(), e=varNames.end(); i!= e; ++i){
      ost << " ,Variable " << i->first << " = " << mva.varValue(i->second);
    }
    return ost;
  }

} // end namespace mu2e

#endif /* GeneralUtilities_MVAStruct_hh */
