// interface with Andrei's code
// volume indexing at the constructor level

//
//
// An EDProducer Module that creates a compressed representation of Calorimeter StepPointMCs
//
// Original author Bertrand Echenard
//
// Note: Replace std::map<art:Ptr<SimParticle, XXXXX> > by std::unordered_map<art::Ptr<SimParticle>, xxx> when the art::Ptr hashing 
//       function is implemented in a future release (recommended to not implement your own, it will clash with art)
//
//
// This modules compresses calorimeter StepPointMCs into CaloShower objects, and flags "interesting" SimParticles. 
// We call a SimParticle entering the calorimeter an Ancestore SimParticle. This ancestor will generate a shower of SimParticles 
// and StepPointMcs in the crystal. 
// Note: if a SimParticvle enters the calorimeter, generates a secondary SimParticle that hit another section of the calorimeter 
// (e.g. a leakage photon in the3 first disk hits the second disk), then the SimParticle hitting the second section is trated 
// as an ancestor SimParticle
//
// Basic idea: We recollect all StepPointMCS attached to a SimParticle ancestor. If the SimParticle is compressible, all StepPointsMC 
// are collapsed into CaloShowerStepMC objects. If not, we create a CaloShowerStepMC object for each SimParticle created by the ancestor
// SimParticle (basically, compress the StepPointMC for each SimParticle). We do the same thing for readout hits but without
// compressing them. At the end of the modules, all StepPointMCs can be dropped, as well as a large fraction of SimParticles 
// with almost no loss of information.

// The compressibility is determined by looking at the interaction codes of the StepPointMCs. 
// Particles are compressed in small intervals of time and crystal longitudinal slices.
//
// See doc-db XXXX for more details
//


// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <unordered_set>
#include <utility>
#include <queue>


// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "HitMakers/inc/CaloShowerStepMCUtil.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepMCCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"

#include "Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"
#include "TH2F.h"





namespace mu2e {



    // Anonymous namespace to hold some helper classes.
    namespace {

	 class caloCompressUtil {

	     public:

		 caloCompressUtil() : _steps(), _sims(), _procs() {}

		 std::vector<StepPointMC const* > const&     steps()        const  {return _steps;}
		 std::vector<art::Ptr<SimParticle> > const&  sims()         const  {return _sims;}
		 std::unordered_set<int> const&              processCodes() const  {return _procs;}

		 void fill(StepPointMC const* step, std::vector<art::Ptr<SimParticle> > sims)
		 {
	             _steps.push_back(step);
		     _procs.insert(step->endProcessCode());	       
		     for (auto const& sim: sims) _sims.push_back(sim); 
		 }
		 
	       
	      private:
	      
	         std::vector<StepPointMC const* >      _steps;
	         std::vector<art::Ptr<SimParticle> >   _sims;
	         std::unordered_set<int>               _procs;	       	       
	 }; 

    } 














  class MakeCaloCompressedHits : public art::EDProducer {

       public:

	 explicit MakeCaloCompressedHits(fhicl::ParameterSet const& pset) :
	    _numZSlices(              pset.get<int>("numZSlices") ),
	    _deltaTime(               pset.get<double>("deltaTime") ),
	    _calorimeterStepPoints(   pset.get<std::string>("calorimeterStepPoints") ),
	    _calorimeterROStepPoints( pset.get<std::string>("calorimeterROStepPoints") ),
	    _physVolInfoInput(        pset.get<std::string>("physVolInfoInput") ),
	    _caloMaterial(            pset.get<std::vector<std::string> >("caloMaterial") ),
	    _diagLevel(               pset.get<int>("diagLevel",0) ),
	    _messageCategory("CaloCompressHits"),
	    _vols(),
	    _procCodes(),
	    _zSliceSize(0)
	 {
	     produces<CaloShowerStepMCCollection>("calorimeter");
	     produces<CaloShowerStepMCCollection>("calorimeterRO");
	     produces<SimParticlePtrCollection>();
	 }

	 virtual      ~MakeCaloCompressedHits() {}
	 virtual void beginJob();
        	 void produce( art::Event& e);
        	 void beginSubRun(art::SubRun& sr) override;



       private:


	 typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;
	 typedef art::Ptr<SimParticle> SimPtr;



	 int                                            _numZSlices;
	 double                                         _deltaTime;
	 std::string                                    _stepPointsCrystal;
	 std::string                                    _calorimeterStepPoints;
	 std::string                                    _calorimeterROStepPoints;
	 art::InputTag                                  _physVolInfoInput;
	 std::vector<std::string>                       _caloMaterial;
	 int                                            _diagLevel;  
	 const std::string                              _messageCategory;
	 const PhysicalVolumeInfoMultiCollection*       _vols;
	 std::map<int,std::set<int> >                   _procCodes;
	 double                                         _zSliceSize;
         std::unordered_set<PhysicalVolumeInfo const*>  _mapPhysVol;
         TH2F *_hStartPos;
         TH2F *_hStopPos;
         TH2F *_hStopPos2;
         TH2F *_hStartPos2;
	 TH1F *_hEtot;



	 void makeCompressedHits(HandleVector const&, HandleVector const&, CaloShowerStepMCCollection&,  
                        	 CaloShowerStepMCCollection&, SimParticlePtrCollection&);

	 void collectStepBySimAncestor(Calorimeter const&, PhysicalVolumeMultiHelper const& , 
                                       HandleVector const&, std::map<SimPtr,caloCompressUtil>&);
	 void collectStepBySim(HandleVector const&, std::map<SimPtr,std::vector<StepPointMC const*> >&);
	 bool isInsideCalorimeter(PhysicalVolumeMultiHelper const&, art::Ptr<SimParticle> const&);

	 bool isCompressible(int simPdgId, std::unordered_set<int> const&, SimParticlePtrCollection const&);
	 void compressSteps(Calorimeter const&, CaloShowerStepMCCollection &, bool isCrystal,
                            int volId, SimPtr const&, std::vector<StepPointMC const*>&);

    };






    //--------------------------------------------------------------------
    void MakeCaloCompressedHits::beginJob()
    {            
	// procCodes are the process codes for the StepPointMC.
	// see MCDataProducts/inc/ProcessCode.hh for code numbering scheme
	//
	// --- These are hardcoded to make sure changes are intended and carefully considered ---
	//      
	_procCodes[11].insert(   {2,16,17,21,23,29,40,49,58} );    // electron
	_procCodes[-11].insert(  {2,16,17,21,23,29,40,49,58} );    // positron
	_procCodes[22].insert(   {2,12,16,17,21,23,29,40,49,58} ); // photon
	_procCodes[2112].insert( {2,16,17,21,23,29,40,49,58,74} ); // neutron
	_procCodes[2212].insert( {16,17,21,23,29,40,45,49,58} );   // proton      

        art::ServiceHandle<art::TFileService> tfs;
        _hStartPos = tfs->make<TH2F>("hStartPos", "Sim start position",  1000,  5000, 15000, 200, 0, 1000);
        _hStopPos  = tfs->make<TH2F>("hStopPos",  "Sim stop position",   1000,  5000, 15000, 200, 0, 1000);
        _hStopPos2 = tfs->make<TH2F>("hStopPos2", "Sim stop position",   1000,  5000, 15000, 200, 0, 1000);
        _hStartPos2 = tfs->make<TH2F>("hStartPos2", "Sim stop position", 1000,  5000, 15000, 200, 0, 1000);
        _hEtot      = tfs->make<TH1F>("hEtot", "Sim stop position",       150,     0, 150);
   
    }

    void MakeCaloCompressedHits::beginSubRun(art::SubRun& sr) 
    {
	art::Handle<PhysicalVolumeInfoMultiCollection> volh;
	sr.getByLabel(_physVolInfoInput, volh);
	_vols = &*volh;
	
	_mapPhysVol.clear();
	for (auto const& vol : *volh)
	{	
	   for (auto const& mv : vol.second) 
	   {
	       std::string material = mv.second.materialName();
	       for (std::string const& caloMaterial : _caloMaterial) 
	           if (material == caloMaterial) _mapPhysVol.insert(&mv.second);
	   }  		
	}
	
    }




    //--------------------------------------------------------------------
    void  MakeCaloCompressedHits::produce(art::Event& event)
    {

	if ( _diagLevel > 0 ) std::cout << "MakeCaloCompressedHits: produce() begin" << std::endl;


	// Check that calorimeter geometry description exists
	art::ServiceHandle<GeometryService> geom;    
	if( !(geom->hasElement<Calorimeter>()) ) return;

	// A container to hold the output hits.
	std::unique_ptr<CaloShowerStepMCCollection> caloShowerStepMCs(new CaloShowerStepMCCollection);
	std::unique_ptr<CaloShowerStepMCCollection> caloROShowerStepMCs(new CaloShowerStepMCCollection);
        std::unique_ptr<SimParticlePtrCollection>   simsToKeep(new SimParticlePtrCollection());

	// These selectors will select data products with the given instance name, and ignore all other fields of the product ID.
	art::ProductInstanceNameSelector getCrystalSteps(_calorimeterStepPoints);
	art::ProductInstanceNameSelector getReadoutSteps(_calorimeterROStepPoints);

	// Get the StepPointMCs from the event.
	HandleVector crystalStepsHandles, readoutStepsHandles;
	event.getMany( getCrystalSteps, crystalStepsHandles);
	event.getMany( getReadoutSteps, readoutStepsHandles);


	makeCompressedHits(crystalStepsHandles,readoutStepsHandles,*caloShowerStepMCs,*caloROShowerStepMCs,*simsToKeep);

	// Add the output hit collection to the event
	event.put(std::move(caloShowerStepMCs),"calorimeter");
	event.put(std::move(caloROShowerStepMCs),"calorimeterRO");
	event.put(std::move(simsToKeep));

	if ( _diagLevel > 0 ) std::cout << "MakeCaloCompressedHits: produce() end" << std::endl;

    } 



    void MakeCaloCompressedHits::makeCompressedHits(HandleVector const& crystalStepsHandles, 
                                                    HandleVector const& readoutStepsHandles, 
                                                    CaloShowerStepMCCollection& caloShowerStepMCs, 
						    CaloShowerStepMCCollection& caloROShowerStepMCs,
						    SimParticlePtrCollection& simsToKeep)
    {

	 PhysicalVolumeMultiHelper vi(*_vols);

	 Calorimeter const & cal = *(GeomHandle<Calorimeter>());
	 _zSliceSize             = (2.0*cal.caloGeomInfo().crystalHalfLength()+0.01)/float(_numZSlices);





	 // Collect the StepPointMC's produced by each SimParticle Ancestor
	 //-----------------------------------------------------------------

	 std::map<SimPtr,caloCompressUtil> crystalAncestorsMap; 
	 collectStepBySimAncestor(cal,vi,crystalStepsHandles,crystalAncestorsMap);



	 //Loop over ancestor simParticles, check if they are compressible, and produce the corresponding caloShowerStepMC
	 //---------------------------------------------------------------------------------------------------------------
	 int nCompress(0),nCompressAll(0);

	 for (auto const& iter : crystalAncestorsMap )
	 {
              SimPtr const&      sim  = iter.first;
              caloCompressUtil const& info = iter.second;

              bool doCompress = isCompressible(sim->pdgId(),info.processCodes(),info.sims());

	      std::map<int,std::vector<StepPointMC const*> > crystalMap;	    
              for (StepPointMC const* step : info.steps()) crystalMap[step->volumeId()].push_back(step);

	      for (auto const& iterCrystal : crystalMap)
	      {
		   int crid = iterCrystal.first;
		   std::vector<StepPointMC const*> steps = iterCrystal.second;

		   if ( doCompress )
		   {
			simsToKeep.push_back(sim);
			compressSteps(cal, caloShowerStepMCs, true, crid, sim, steps);

			if (_diagLevel > 0)
			{
			    _hStartPos->Fill(sim->startPosition().z(),sqrt(pow(sim->startPosition().x()+3904,2)+pow(sim->startPosition().y(),2) ));
                	    _hStopPos->Fill(sim->endPosition().z(),sqrt(pow(sim->endPosition().x()+3904,2)+pow(sim->endPosition().y(),2) ));
			    for (auto const& simD: info.sims()) 
				_hStartPos2->Fill(simD->startPosition().z(),sqrt(pow(simD->startPosition().x()+3904,2)+pow(simD->startPosition().y(),2)) );
			}	      
		   }
		   
		   else 
		   {
			std::map<SimPtr, std::vector<StepPointMC const*> > newSimStepMap;
			for (StepPointMC const* step : steps) newSimStepMap[step->simParticle()].push_back(step);
			for (auto& iter : newSimStepMap)
			{ 
			   compressSteps(cal, caloShowerStepMCs, true, crid, iter.first, iter.second);
			   simsToKeep.push_back(iter.first);
			}   
		   }
	      }   
              ++nCompressAll;
	      if (doCompress) ++nCompress;
	 }



         // Do the same for the readouts, but there is no need to compress
         //---------------------------------------------------------------

	 std::map<SimPtr,std::vector<StepPointMC const*> > simStepROMap;
	 collectStepBySim(readoutStepsHandles, simStepROMap);

	 for (auto const& iter : simStepROMap )
	 {
              SimPtr   const& sim   = iter.first;
              std::vector<StepPointMC const*> const& steps = iter.second;

	      std::map<int,std::vector<StepPointMC const*> > crystalMap;	    
              for (StepPointMC const* step : steps) crystalMap[step->volumeId()].push_back(step);

	      for (auto const& iterCrystal : crystalMap)
	      {
		   int roid = iterCrystal.first;
		   std::vector<StepPointMC const*> steps = iterCrystal.second;              		   
		   compressSteps(cal, caloROShowerStepMCs, false, roid, sim, steps);
	      }   
	 }
         
	 
	 
	 //Final statistics
	 if (_diagLevel > 0) std::cout<<"MakeCaloCompressedHits compressed "<<nCompress<<" / "<<nCompressAll<<" incoming SimParticles"<<std::endl;
         if (_diagLevel > 0) std::cout<<"MakeCaloCompressedHits keeping "<<simsToKeep.size()<<" CaloShowerSteps"<<std::endl;

    } 




















  //------------------------------------------------------------------------------------------------------------------
    void MakeCaloCompressedHits::collectStepBySimAncestor(Calorimeter const& cal, PhysicalVolumeMultiHelper const& vi,
                                                          HandleVector const& stepsHandles,
                                                          std::map<SimPtr,caloCompressUtil>& ancestorsMap)
    {       

	// double crystalLength = 2.0*cal.caloGeomInfo().crystalHalfLength();
	 std::map<SimPtr,SimPtr> simToAncestorMap;

	 for ( HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end(); i != e; ++i )
	 {     
	      art::Handle<StepPointMCCollection> const& handle(*i);
	      StepPointMCCollection const& steps(*handle);

	      for ( auto const& step : steps )
	      {
		    SimPtr sim = step.simParticle();

		    SimParticlePtrCollection inspectedSims;	      
		    while (sim->hasParent() && isInsideCalorimeter(vi,sim) )  
		    //while (sim->hasParent() && cal.isInsideCalorimeter(sim->startPosition()) )  
		    {
			//simparticle starting in one section and ending in another one see note above
			if (!cal.isContainedSection(sim->startPosition(),sim->endPosition()) ) break; 
			
			auto const alreadyInspected = simToAncestorMap.find(sim);
			if (alreadyInspected != simToAncestorMap.end()) {sim = alreadyInspected->second; break;}

			inspectedSims.push_back(sim);
			sim = sim->parent();
		    }
		    
                    for (SimPtr const& inspectedSim : inspectedSims) simToAncestorMap[inspectedSim] = sim;		    
		    ancestorsMap[sim].fill(&step,inspectedSims);
	      }
	 }
	 
         if (_diagLevel > 0) std::cout<<"MakeCaloCompressedHits starting with "<<simToAncestorMap.size()+ancestorsMap.size()<<" simParticles"<<std::endl;
    }


  //--------------------------------------------------------------------------------------------------------------
    bool MakeCaloCompressedHits::isInsideCalorimeter(PhysicalVolumeMultiHelper const& vi, art::Ptr<SimParticle> const& thisSimPtr)
    {	
        return 	_mapPhysVol.find(&vi.startVolume(*thisSimPtr)) != _mapPhysVol.end();
	
	//const std::string material = vi.startVolume(*thisSimPtr).materialName();
	//for (std::string const& caloMaterial : _caloMaterial) if (material== caloMaterial) return true;
	//return false;
    }



  //-------------------------------------------------------------------------------------------------------------
    void MakeCaloCompressedHits::collectStepBySim(HandleVector const& stepsHandles, 
                                                  std::map<SimPtr,std::vector<StepPointMC const*> >& simStepMap)
    {       
	 for ( HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end(); i != e; ++i )
	 {     
	      art::Handle<StepPointMCCollection> const& handle(*i);
	      StepPointMCCollection const& steps(*handle);

	      for ( auto const& step : steps ) simStepMap[step.simParticle()].push_back(&step);
	 }
    }


  //-------------------------------------------------------------------------------------------------------------------------------
    void MakeCaloCompressedHits::compressSteps(Calorimeter const& cal, CaloShowerStepMCCollection &caloShowerStepMCs, 
                                               bool isCrystal, int volId, SimPtr const& sim, std::vector<StepPointMC const*>& steps)
    {

	  std::sort( steps.begin(), steps.end(), [](const StepPointMC* a, const StepPointMC* b) {return a->time() < b->time();} );      

	  CaloShowerStepMCUtil buffer(_numZSlices,CaloShowerStepMCUtil::weight_type::energy );	
	  for (StepPointMC const* step : steps)
	  {
	       CLHEP::Hep3Vector pos  = (isCrystal) ? cal.toCrystalFrame(volId,step->position()) : cal.toCrystalFrame(cal.crystalByRO(volId),step->position());
	       int               idx  = (isCrystal) ? int(pos.z()/_zSliceSize) : 0;

	       if (buffer.entries(idx) == 0) buffer.setTinit(idx,step->time()); 

	       if (step->time()-buffer.t0(idx) > _deltaTime)
	       { 
		  caloShowerStepMCs.push_back(CaloShowerStepMC(volId, sim, buffer.entries(idx), buffer.time(idx), buffer.eDep(idx), buffer.pos(idx), buffer.covPos(idx)));
		  if (_diagLevel > 1) {std::cout<<"MakeCaloCompressedHits  inserted     "; buffer.printBucket(idx);}
		  buffer.reset(idx);   
		  buffer.setTinit(idx,step->time());
	       }

	       buffer.add( idx, step->eDep(), step->time(), pos);
	 }

	 //do not forget to flush final buffer :-)
	 for (int i=0;i<buffer.nBuckets();++i)
	 {
	      if (buffer.entries(i) == 0) continue;
	      caloShowerStepMCs.push_back(CaloShowerStepMC(volId, sim,  buffer.entries(i), buffer.time(i), buffer.eDep(i), buffer.pos(i), buffer.covPos(i)));
	      if (_diagLevel > 1) {std::cout<<"MakeCaloCompressedHits  inserted     ";  buffer.printBucket(i);}
	 }	 	 	 
    }


  //----------------------------------------------------------------------------------------------------------------------------
    bool MakeCaloCompressedHits::isCompressible(int simPdgId, std::unordered_set<int> const& processCodes, SimParticlePtrCollection const& sims)
    {      

	if (simPdgId > 1000000000) return true;                            //ions are always compressed
	if (_procCodes.find(simPdgId) == _procCodes.end()) return false; 

	std::set<int> const& proc = _procCodes[simPdgId];      
	for (int code : processCodes) if ( proc.find(code) == proc.end() ) return false;

	if (simPdgId==2212 || simPdgId==2112)
	{	 
	  for ( auto const& sim : sims) 
	    if (sim->pdgId()==22 || sim->pdgId()==11 || sim->pdgId()==-11) {if (sim->startMomentum().mag() > 1) return false;}
	}	 

	return true;    
    }

	

 

}



using mu2e::MakeCaloCompressedHits;
DEFINE_ART_MODULE(MakeCaloCompressedHits);


