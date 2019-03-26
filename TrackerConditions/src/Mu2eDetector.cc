
#include <iostream>

#include "TrackerConditions/inc/Mu2eDetector.hh"
#include <sstream>
#include "cetlib_except/exception.h"

using namespace std;

namespace mu2e {


  const DetStrawElem* Mu2eDetector::strawElem(StrawId const& istraw) const{
    const DetStrawElem* retval(0);
    auto ifnd = _strawmap.find(istraw);
    if(ifnd != _strawmap.end())
      retval = ifnd->second;
    else
      throw cet::exception("RECO_NO_ELEMENT")
	<<"mu2e::Mu2eDetector: no element associated to straw " 
	<< istraw << std::endl;
    return retval;
  }

  Mu2eDetector::~Mu2eDetector() {
    for(auto istraw : _strawmap) {
      delete istraw.second;
    }
  }

  void Mu2eDetector::print( ostream& out) const{
    out << "Mu2eDetector has "<<_strawmap.size() << " elements" << endl;
  }
  
} // namespace mu2e
