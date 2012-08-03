// Magnet parameters, used by both filter and detector magnets.
//
// Andrei Gaponenko, 2011

#ifndef EXTMONFNALMAGNET_HH
#define EXTMONFNALMAGNET_HH

#include <vector>

namespace mu2e {

  class ExtMonFNALMagnetMaker;

  class ExtMonFNALMagnet {
    std::vector<double> _outerHalfSize;
    double _apertureWidth;
    double _apertureHeight;
    double _fieldStrength;


  public:

    ExtMonFNALMagnet();
    // An initialized instance of this class should be obtained via ExtMonFNALMagnetMaker
    friend class ExtMonFNALMagnetMaker;

    const std::vector<double> &outerHalfSize() const { return _outerHalfSize; }
    double apertureWidth() const { return _apertureWidth; }
    double apertureHeight() const { return _apertureHeight; }
    double fieldStrength() const { return _fieldStrength; }

    // derived:
    double trackBendRadius(double momentum) const;
    double trackBendHalfAngle(double momentum) const;
  };

}// namespace mu2e

#endif/*EXTMONFNALMAGNET_HH*/
