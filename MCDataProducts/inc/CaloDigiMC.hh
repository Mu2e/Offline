#ifndef MCDataProducts_CaloDigiMC_hh
#define MCDataProducts_CaloDigiMC_hh

// Original author G. Pezzullo

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/CaloShowerStepMC.hh"

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes

namespace mu2e {

  class CaloDigiMC{

  public:

    CaloDigiMC():
      _nParticles(0),
      _totalEDep (0.),
      _meanTime  (0.){}
    
    CaloDigiMC(int                   NParticles, 
	       std::vector<double>   EDep,
	       std::vector<double>   Time,
	       std::vector<int>      PdgId, 
	       std::vector<art::Ptr<SimParticle> >  SimCol):
      _nParticles(NParticles)
    {
      _totalEDep  = 0.;
      _meanTime   = 0.;
      
      for (int i=0; i< NParticles; ++i){
	_totalEDep    += EDep.  at(i);
	_meanTime     += Time.  at(i) * EDep.  at(i);
	_eDep.        push_back(EDep.  at(i));
	_time.        push_back(Time.  at(i));
	_pdgId.       push_back(PdgId. at(i));
	_simParticle. push_back(SimCol.at(i));
      }
      _meanTime  /= _totalEDep;

    }

    CaloDigiMC(const CaloDigiMC &CaloMC) {
      int     size = CaloMC.nParticles();
      
      _nParticles = size;
      _totalEDep  = CaloMC.totalEDep();
      _meanTime   = CaloMC.meanTime();

      for (int i=0; i<size; ++i){
	_eDep.        push_back(CaloMC.eDep(i));
	_time.        push_back(CaloMC.time(i));
	_pdgId.       push_back(CaloMC.pdgId(i));
	_simParticle. push_back(CaloMC.simParticle(i));
      }
      
    }


    // Accessors
    int                      nParticles  ()          const {return _nParticles; }
    double                   totalEDep   ()          const {return _totalEDep;  }
    double                   meanTime    ()          const {return _meanTime;   }
    double                   timeFirst   ()          const {return _timeFirst;   }

    double                   eDep        (int Index) const {return _eDep.at(Index); }
    double                   time        (int Index) const {return _time.at(Index); }
    int                      pdgId       (int Index) const {return _pdgId.at(Index);}

    art::Ptr<SimParticle>    simParticle (int Index)const {return _simParticle.at(Index); }

    void                     addCaloShower(const CaloShowerStepMC* CaloShower);
    void                     init();

  private:
    int                           _nParticles;

    double                        _totalEDep;
    double                        _meanTime;
    double                        _timeFirst;

    std::vector<double>           _eDep;
    std::vector<double>           _time;
    std::vector<double>           _pdgId;

    std::vector<art::Ptr<SimParticle> >       _simParticle;

  };

} // namespace mu2e

#endif /* MCDataProducts_CaloDigiMC_hh */
