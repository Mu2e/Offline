#ifndef GeometryService_PTMonMaker_hh
#define GeometryService_PTMonMaker_hh
//
// construct and return a PTargetMon
//
// original author Helenka Casler
//

namespace mu2e {

  class SimpleConfig;
  // TODO: class PTMon (or something)

  class PTMonMaker {

  public:
    PTMonMaker( SimpleConfig const& config );

  private:
    // Extract info from the config file.
    void parseConfig( const SimpleConfig& config );
    // TODO all other methods should also return void
    // TODO internal numbers, etc to hold config/geom data
  };

} // namespace mu2e

#endif