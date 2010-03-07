#ifndef ToyDP_RandomEngineState_hh
#define ToyDP_RandomEngineState_hh
//
// Persistent data for the state of one CLHEP::HepRandomEngine.
//
// $Id: RandomEngineState.hh,v 1.1 2010/03/07 22:01:00 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/03/07 22:01:00 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) CLHEP specifies that state will be returned as vector<unsigned long>.  
//    The size of a long is machine dependent. If unsigned long is an 8 byte 
//    variable, only the least significant 4 bytes are filled and the most 
//    significant 4 bytes are zero.  We need to store the state as
//    with a machine indpendent size, which we choose to be uint32_t.  
//    So this class has constructors, setters and getters that do the
//    required type conversions.
//
// 2) Normally I would not supply setters for this class and would require
//    all state to be supplied in the c'tor.  But I want to avoid the
//    double copy implied by push_back on a container; in that case I need
//    to allow two phase construction.
//

#include <vector>
#include <string>

#include "GeneralUtilities/inc/vectorTransform.hh"

namespace mu2e {

  class RandomEngineState {

    // The data type we persist.
    typedef uint32_t StorageType;
    
  public:
    
    RandomEngineState(){}
    RandomEngineState( const std::string&              label,
		       const std::vector<StorageType>& state,
		       const std::vector<StorageType>& seed):
      label_(label),
      state_(state),
      seed_(seed){
    }
    RandomEngineState( const std::string&                label,
		       const std::vector<unsigned long>& state,
		       const std::vector<unsigned long>& seed):
      label_(label){
      setState(state);
      setSeed(seed);
    }
    ~RandomEngineState(){}
    
    const std::string&              getLabel() const { return label_; }
    const std::vector<StorageType>& getState() const { return state_; }
    const std::vector<StorageType>& getSeed()  const { return seed_;  }
    
    void setLabel(const std::string& value)              { label_ = value; }
    void setState(const std::vector<StorageType>& value) { state_ = value; }
    void setSeed (const std::vector<StorageType>& value) { seed_  = value; }
    
    // Type converting accesors and setters.
    void getState( std::vector<unsigned long>& out) { 
      vectorTransform<uint32_t, unsigned long>( state_, out );
    }

    void getSeed( std::vector<unsigned long>& out) { 
      vectorTransform<uint32_t, unsigned long>( seed_, out );
    }

    void setState(const std::vector<unsigned long>& value){
      vectorTransform<unsigned long, uint32_t>( value, state_ );
    }
    void setSeed (const std::vector<unsigned long>& value){
      vectorTransform<unsigned long, uint32_t>( value, seed_ );
    }
    
  private:
    
    std::string label_;
    std::vector<StorageType> state_;
    std::vector<StorageType> seed_;
    
  };

} // end namespace mu2e

#endif
