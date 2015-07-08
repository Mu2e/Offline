//
// Utility to compute the cluster moments
//  
// Original author B. Echenard
//

// Mu2e includes
#include "BTrk/BaBar/Constants.hh"
#include "CaloCluster/inc/CaloClusterMoments.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"


// C++ includes
#include <iostream>
#include <list>

using CLHEP::Hep3Vector;


namespace mu2e {


       void CaloClusterMoments::calculate(cogtype mode) 
       {
	   double sxi(0),sxi2(0),syi(0),syi2(0),sxyi(0),szi(0),szi2(0),swi(0);

       	   auto const& main     = _caloCluster.caloCrystalHitsPtrVector();
	   //           double clusterEnergy = _caloCluster.energyDep();
	   	   
	   
	   //calculate first and second moments in one pass
	   for (auto it = main.begin(); it !=main.end(); ++it)
	   {        
		int    crId((*it)->id());
		double energy((*it)->energyDep());

		if (_cal.crystal(crId).sectionId() != _iSection) continue;

		double xCrystal = _cal.crystalOrigin(crId).x();
		double yCrystal = _cal.crystalOrigin(crId).y();
		double zCrystal = _cal.crystalOrigin(crId).z();

		double weight = energy;
		if (mode == LinearMod)   weight = 1.01*energy - 6.25;
		if (mode == Logarithm)   weight = log(energy);
		if (mode == Sqrt)        weight = sqrt(energy);

		sxi  += xCrystal*weight;
		sxi2 += xCrystal*xCrystal*weight;
		syi  += yCrystal*weight;
		syi2 += yCrystal*yCrystal*weight;
		sxyi += xCrystal*yCrystal*weight;
		szi  += zCrystal*weight;
		szi2 += zCrystal*yCrystal*weight;
		swi  += weight;
	  }
          
	  
	  CLHEP::Hep3Vector cogMu2eFrame(sxi/swi,syi/swi,szi/swi);
	  _cog = _cal.toSectionFrameFF(_iSection,cogMu2eFrame);
	  
	  _secondMoment = sxi2 -sxi*sxi/swi + syi2 -syi*syi/swi; 
         
	 
	 double beta  = (sxyi - sxi*syi/swi)/(sxi2-sxi*sxi/swi);
         _angle = atan(beta);
         if (_angle < 0) _angle += Constants::pi;
	 
       }

}


       /*
       void CaloClusterMoments::calculate(cogtype mode) 
       {
	   double sumWeights(0),secondMoment(0);
	   CLHEP::Hep3Vector cogMu2eFrame(0,0,0);


       	   auto const& main     = _caloCluster.caloCrystalHitsPtrVector();
           double clusterEnergy = _caloCluster.energyDep();
	   	   

	   //we need to calculate the cog first
	   for (auto it = main.begin(); it !=main.end(); ++it)
	   {        
		int    crId((*it)->id());
		double energy((*it)->energyDep());

		if (_cal.crystal(crId).sectionId() != _iSection) continue;
		if (energy/clusterEnergy < 0.05)                 continue;


		double xCrystal = _cal.crystalOrigin(crId).x();
		double yCrystal = _cal.crystalOrigin(crId).y();
		double zCrystal = _cal.crystalOrigin(crId).z();

		double weight = -5.26 + 1.24*energy;
		//if (mode == Logarithm)   weight = -1.73+log(energy);
		//if (mode == Sqrt) weight = sqrt(energy);

		cogMu2eFrame[0] += xCrystal*weight;
		cogMu2eFrame[1] += yCrystal*weight;
		cogMu2eFrame[2] += zCrystal*weight;
		sumWeights += weight;
	  }

	  cogMu2eFrame[0] /= sumWeights;
	  cogMu2eFrame[1] /= sumWeights;
	  cogMu2eFrame[2] /= sumWeights;

	  _cog = _cal.toSectionFrameFF(_iSection,cogMu2eFrame);
	  
	  

          //calculate second moment and major axis of cluster
          for (auto it = main.begin(); it !=main.end(); ++it)
	  {        
	      int crId((*it)->id());
	      double energy = (*it)->energyDep();
	      CLHEP::Hep3Vector crystalPos = _cal.crystalOrigin(crId);

	      if (_cal.crystal(crId).sectionId() != _iSection) continue;
	      if (energy/clusterEnergy < 0.05) continue;
	      
	      double dx = crystalPos.x() - _cog.x();
	      double dy = crystalPos.y() - _cog.y();
	      secondMoment += energy*(dx*dx + dy*dy);	      
	  }
          _secondMoment = secondMoment/sumWeights;
         
       }

      */ 
 


