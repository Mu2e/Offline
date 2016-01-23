//
// Utility to study the MC content of a calo cluster. Browser over SimParticles of each crystal, and keep only distinct
// entries, updating the total energy, time and position
//
//
// Original author B. Echenard
//

// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "CaloCluster/inc/CaloContentMC.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

#include "CLHEP/Vector/ThreeVector.h"

// C++ includes
#include <map>
#include<iostream>



namespace mu2e {


       CaloContentMC::CaloContentMC(CaloHitMCNavigator const& navi, CaloCluster const& cluster) :
        _navi(&navi), _cluster(&cluster), _simUtilMap(), _simBase(), _simPart(), _edep(), _time(), _mom(), _pos()
       {
           fillCluster();
       }




       void CaloContentMC::fillCluster()
       {

            std::set<SimParticlePtr > simBaseSet;


            for (int i=0;i<_cluster->size();++i)
            {

                //run only over the first readout, since the other are duplicating the StepPointMc in the crystal
                CaloCrystalHit const& hit      = *(_cluster->caloCrystalHitsPtrVector().at(i));
                CaloHit const& caloHit         = *(hit.readouts().at(0));
                CaloHitSimPartMC const& hitSim = _navi->sim(caloHit);

                fillMaps(hitSim,simBaseSet);
            }


            for (auto const& it: _simUtilMap)
            {
               _simPart.push_back(it.first);
               _edep.push_back(it.second.edepTot());
               _time.push_back(it.second.time());
               _mom.push_back(it.second.momentum());
               _pos.push_back(it.second.position());
            }


            for (std::set<SimParticlePtr >::const_iterator it=simBaseSet.begin();it!=simBaseSet.end();++it)
               _simBase.push_back(*it);


       }


       void CaloContentMC::fillMaps(CaloHitSimPartMC const& hitSim, std::set<SimParticlePtr >& simBaseSet)
       {

            std::vector<SimParticlePtr >    const& sim  = hitSim.simParticles();
            std::vector<double>             const& edep = hitSim.eDep();
            std::vector<double>             const& time = hitSim.time();
            std::vector<double>             const& mom  = hitSim.momentum();
            std::vector<CLHEP::Hep3Vector>  const& pos  = hitSim.position();

            for (unsigned int i=0;i<sim.size();++i)
            {
                SimUtilMap::iterator mfind = _simUtilMap.find(sim[i]);
                if (mfind != _simUtilMap.end()) mfind->second.update(edep[i],time[i],mom[i],pos[i]);
                else _simUtilMap.insert(std::pair<SimParticlePtr, CaloClusterMCUtil>(sim[i],CaloClusterMCUtil(edep[i],time[i],mom[i],pos[i])) );

                SimParticlePtr mother = sim[i];
                while ( mother->hasParent() ) mother = mother->parent();
                simBaseSet.insert(mother);
            }
       }





       bool CaloContentMC::hasConversion()
       {
             for (auto const& sim:_simBase)
             {
                 if (sim->genParticle() && sim->genParticle()->generatorId()==GenId::conversionGun) return true;
             }
             return false;
       }


}
