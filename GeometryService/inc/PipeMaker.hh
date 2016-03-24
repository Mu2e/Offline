#ifndef GeometryServices_PipeMaker_hh
#define GeometryServices_PipeMaker_hh
//
// Class to construct and return Pipe
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class Pipe;
  class SimpleConfig;

  class PipeMaker {
  public:

    static std::unique_ptr<Pipe>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* GeometryServices_PipeMaker_hh */
