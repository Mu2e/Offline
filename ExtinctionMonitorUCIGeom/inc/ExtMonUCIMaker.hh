// $Id: ExtMonUCIMaker.hh,v 1.4 2012/02/16 20:25:46 youzy Exp $
// $Author: youzy $
// $Date: 2012/02/16 20:25:46 $

#ifndef EXTMONUCIMAKER_HH
#define EXTMONUCIMAKER_HH

#include <memory>

namespace mu2e { class SimpleConfig; }
namespace mu2e { namespace ExtMonUCI { class ExtMon; } }

namespace mu2e {

  namespace ExtMonUCI {

    //forward declarations
    class ExtMon;

    class ExtMonMaker {

    public:
      ExtMonMaker(const SimpleConfig& config);

      // interface to GeometryService
      std::auto_ptr<ExtMon> _det;
      std::auto_ptr<ExtMon> getDetectorPtr() { return _det; }

    private:
      void MakeCols();
      void MakeMags();
      void MakeTofs();
      void MakeShds();

    };
  }
}

#endif/*EXTMONUCIMAKER_HH*/
