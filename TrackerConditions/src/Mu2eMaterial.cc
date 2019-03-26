
#include <iostream>

#include "TrackerConditions/inc/Mu2eMaterial.hh"
#include <sstream>

using namespace std;

namespace mu2e {

  void Mu2eMaterial::print( ostream& out) const{
    out << "gasMaterial " 
	<<  *( _strawtype->gasMaterial()->name() ) << endl;
    out << "wallMaterial " 
	<<  *( _strawtype->wallMaterial()->name() ) << endl;
    //out << "wireMaterial " 
    //	<<  *( _strawtype->wireMaterial()->name() ) << endl;
    out << "offset " <<  _strawtype->offset() << endl;
    out << "tolerance " <<  _strawtype->tolerance() << endl;
    out << "maxRadiusFraction " <<  _strawtype->maxRadiusFraction() << endl;
  }
  
} // namespace mu2e
