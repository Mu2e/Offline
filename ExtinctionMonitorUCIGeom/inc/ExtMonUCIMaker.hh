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
      void MakeTofs();

    };
  }
}

#endif/*EXTMONUCIMAKER_HH*/
