#ifndef GEOM_DET_H
#define GEOM_DET_H

//
// A base class for detector components.
//
// $Id: Detector.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//


#include <string>

namespace mu2e
{
  class Detector
  {
  public:
    virtual ~Detector();
    virtual void update();
    virtual std::string name() const;
  };
}

#endif
