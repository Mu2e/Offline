//
//  $Id: 
//  $Author: 
//  $Date: 
//
//  Original author Pavel Murat
//

#include "RecoDataProducts/inc/AlgorithmID.hh"

using namespace std;

namespace mu2e {

  // Constructors
  AlgorithmID::AlgorithmID() {
    _bestID  = 0;
    _algMask = 0;
  }


  AlgorithmID::AlgorithmID(const AlgorithmID & A) {
    _bestID  = A._bestID;
    _algMask = A._algMask;
  }

  AlgorithmID::~AlgorithmID() {
  }

  // operator overloading
  AlgorithmID & AlgorithmID::operator = (const AlgorithmID & A) {
    _bestID  = A._bestID;
    _algMask = A._algMask;
    return (*this);
  }


  void AlgorithmID::clear () {
    _bestID  = 0;
    _algMask = 0;
  }




} // end namespace mu2e


