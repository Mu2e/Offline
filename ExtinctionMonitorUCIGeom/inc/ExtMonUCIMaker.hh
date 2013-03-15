// $Id: ExtMonUCIMaker.hh,v 1.5 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $

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
      std::unique_ptr<ExtMon> _det;
      std::unique_ptr<ExtMon> getDetectorPtr() { return std::move(_det); }

    private:
      void MakeCols();
      void MakeMags();
      void MakeTofs();
      void MakeShds();

    };
  }
}

#endif/*EXTMONUCIMAKER_HH*/
