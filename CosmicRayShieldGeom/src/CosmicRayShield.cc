//
// Class to represent CosmicRayShield
//
// $Id: CosmicRayShield.cc,v 1.4 2013/09/13 06:42:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/09/13 06:42:44 $
//
// Original author KLG
//

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

namespace mu2e 
{
    void CosmicRayShield::getMinMaxPoints(const std::string &sectorName, std::vector<double> &minPoint, std::vector<double> &maxPoint) const
    {
      bool first=true;
      minPoint.resize(3);
      maxPoint.resize(3);
      size_t nSectors=_scintillatorShields.size();
      for(size_t s=0; s<nSectors; s++)
      {
        const CRSScintillatorShield &sector = _scintillatorShields[s];
        if(sector.getName().find(sectorName,4)==4)
        {
          size_t nModules = sector.nModules();
          for(size_t m=0; m<nModules; m++)
          {
            const CRSScintillatorModule &module = sector.getModule(m);
            for(size_t l=0; l<2; l++)  //the two aluminum sheets are the outermost layers of each CRV module
            {
              const CRSAluminumSheet &aluminumSheet = module.getAluminumSheet(l);
              const CLHEP::Hep3Vector &position = aluminumSheet.getPosition();
              const std::vector<double> &halfLengths = aluminumSheet.getHalfLengths();
              if(first)
              {
                for(int i=0; i<3; i++)
                {
                  minPoint[i]=position[i]-halfLengths[i];
                  maxPoint[i]=position[i]+halfLengths[i];
                }
                first=false;
              }
              else
              {
                for(int i=0; i<3; i++)
                {
                  double minPointTmp=position[i]-halfLengths[i];
                  double maxPointTmp=position[i]+halfLengths[i];
                  if(minPoint[i]>minPointTmp) minPoint[i]=minPointTmp;
                  if(maxPoint[i]<maxPointTmp) maxPoint[i]=maxPointTmp;
                }
              }
            }
          }
        }
      }
    }

    std::vector<double> CosmicRayShield::getSectorHalfLengths(const std::string &sectorName) const
    {
      std::vector<double> minPoint, maxPoint;
      CosmicRayShield::getMinMaxPoints(sectorName, minPoint, maxPoint);

      std::vector<double> halfLengths;
      halfLengths.resize(3);
      for(int i=0; i<3; i++)
      {
        halfLengths[i]=(maxPoint[i]-minPoint[i])/2.0;
      }
      return halfLengths;
    }

    CLHEP::Hep3Vector CosmicRayShield::getSectorPosition(const std::string &sectorName) const
    {
      std::vector<double> minPoint, maxPoint;
      CosmicRayShield::getMinMaxPoints(sectorName, minPoint, maxPoint);

      CLHEP::Hep3Vector position;
      for(int i=0; i<3; i++)
      {
        position[i]=(maxPoint[i]+minPoint[i])/2.0;
      }
      return position;
    }
}


