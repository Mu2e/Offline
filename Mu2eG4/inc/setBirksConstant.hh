#ifndef Mu2eG4_setBirksConstant_hh
#define Mu2eG4_setBirksConstant_hh
//
// Set the G4 BirksConstant as specified in the fhicl file.
//
// $Id: setBirksConstant.hh,v 1.1 2015/11/04 19:28:01 genser Exp $
// $Author: genser $
// $Date: 2015/11/04 19:28:01 $
//
//-----------------------------------------------------------------------------

namespace fhicl { class ParameterSet; }

namespace mu2e{

  void setBirksConstant(const fhicl::ParameterSet& pset);

}  // end namespace mu2e

#endif /* Mu2eG4_setBirksConstant_hh */
