//
//  checkConfigRelics.hh
//
//  Checks the configuration of a geomtry file.  Throws an exception
//  if an "old style" geomtry file that is no longer supported is used.
//
//  Author: Lisa Goodenough
//  Date: 2017/4/19
//
//

#ifndef Mu2eG4_checkConfigRelics_hh
#define Mu2eG4_checkConfigRelics_hh


namespace mu2e {
    
    class SimpleConfig;
    
    void checkConfigRelics(const SimpleConfig& config);
    
}// end namespace mu2e


#endif /* Mu2eG4_checkConfigRelics_hh */
