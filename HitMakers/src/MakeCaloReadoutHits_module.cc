//
// An EDProducer Module that reads StepPointMC objects and turns them into
// CaloHit, CaloHitMCTruth, CaloHitMCPtr, CrystalOnlyHit, CrystalHitMCPtr
// objects.
//
// Original author Ivan Logashenko; major changes by Krzysztof Genser and Rob Kutschke.
//
// Notes
// 1) The CrystalOnlyHitCollection is a form of MC truth on a per crystal basis.
//    It represents an idea per crystal response if no readouts were hit directly.
//
// 2) We still have two times the crystal StepPointMC hits saved per hit in a crystal. Note sure we can 
//    keep only a single copy before simulating the whole readout chain, so keep them for now. 


// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <unordered_map>
#include <vector>
#include <utility>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/sort_functors.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "HitMakers/inc/CaloCrystalMCUtil.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "SeedService/inc/SeedService.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandGaussQ.h"


using namespace std;

namespace mu2e {

  // Anonymous namespace to hold some helper classes.
  namespace {


       // A helper class to hold temporary information.
       class ROHit {
	   
	   public:

	     // Is this StepPointMC from a hit in the crystal or a hit in the readout plane?
	     enum step_type {crystal, readout};

	     art::Ptr<StepPointMC> _step;    // Ptr back to the StepPointMC.
	     double    _edep;                // copy of energy from the StepPointMC.
	     double    _edep_corr;           // Energy corrected for absorption within the crystal.
	     step_type _type;                // Is this a hit in the crystal or the readout?
	     double    _time;                // copy of the time from the StepPointMC.

	     ROHit(art::Ptr<StepPointMC> const& step, double edep, double edep1, step_type charged, double time):
               _step(step), _edep(edep), _edep_corr(edep1),
               _type(charged), _time(time) { }

	     // This operator is overloaded in order to time-sort the hits
	     bool operator <(const ROHit& b) const { return (_time < b._time); }

       }; 


       // A helper class to add Ptr's to the appropriate PtrStepPointMCVector collection.
       class PtrAdder{

	   public:
	      PtrAdder( PtrStepPointMCVector& crystals, 
                	PtrStepPointMCVector& readouts ) : 
			_crystals(crystals),_readouts(readouts){}

	      // Choose the appropriate collection and add the Ptr for this hit.
	      void operator()( ROHit const& hit )
	      {
        	if ( hit._type == ROHit::crystal) _crystals.push_back( hit._step );
        	else                              _readouts.push_back( hit._step );
	      }


	   private:

	     PtrStepPointMCVector& _crystals;
	     PtrStepPointMCVector& _readouts;

       }; 


  } // end anonymous namespace





  //--------------------------------------------------------------------
  //
  //
  class MakeCaloReadoutHits : public art::EDProducer {
  private:

 
    typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;
    typedef art::Ptr<StepPointMC>   StepPtr;
    typedef std::vector<StepPtr >   StepPtrs;
    typedef art::Ptr<SimParticle>   SimPtr;
    typedef std::vector<SimPtr >    SimPtrs;
    typedef std::map<int,StepPtrs > HitMap;

	 
    int   _diagLevel;  
    int   _maxFullPrint;
    bool  _fillDetailedHit;
    double                _timeGap;
    double                _blindTime;
    bool  _caloLRUcorrection;
    bool  _caloNonLinCorrection;
    CLHEP::RandGaussQ _randGauss;

    std::string _stepPoints;
    std::string _rostepPoints;
    std::string _g4ModuleLabel;   // Name of the module that made these hits.

    SimParticleTimeOffset _toff;  // time offset smearing
    double _mbtime;               // period of 1 microbunch
    double _mbbuffer;             // buffer on that for ghost hits (wrapping)

    const std::string _messageCategory;
 
  public:

       // First vector is list of crystal steps, associated with particular readout element.
       // Second vector is list of readout steps, associated with particular readout element.

    explicit MakeCaloReadoutHits(fhicl::ParameterSet const& pset) :

      // Parameters
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _fillDetailedHit(pset.get<bool>("fillDetailedHit",0)),
      _timeGap               (pset.get<double>("timeGap"               , 2.)),
      _blindTime             (pset.get<double>("blindTime"             ,300.)),
      _caloLRUcorrection(pset.get<bool>("caloLRUcorrection",0)),
      _caloNonLinCorrection(pset.get<bool>("caloNonLinCorrection",0)),
      _randGauss( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
      _stepPoints(pset.get<string>("calorimeterStepPoints","calorimeter")),
      _rostepPoints(pset.get<string>("calorimeterROStepPoints","calorimeterRO")),
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
      _toff(pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
      _mbbuffer(pset.get<double>("TimeFoldingBuffer",100.0)), // nsec
      _messageCategory("CaloReadoutHitsMakerNew"){
      
      // Tell the framework what we make.
      produces<CaloHitCollection>();
      produces<CaloHitMCTruthCollection>();
      produces<CaloCrystalOnlyHitCollection>();
      produces<PtrStepPointMCVectorCollection>("CaloHitMCCrystalPtr");
      produces<PtrStepPointMCVectorCollection>("CaloHitMCReadoutPtr");
      produces<CaloHitSimPartMCCollection>();
    }

    virtual ~MakeCaloReadoutHits() { }
    virtual void beginJob();
    void produce( art::Event& e);

  private:

    void makeCalorimeterHits (HandleVector const& crystalStepsHandles,
			      HandleVector const& readoutStepsHandles,
			      CaloHitCollection &,
			      CaloHitMCTruthCollection&,
			      CaloCrystalOnlyHitCollection&,
			      PtrStepPointMCVectorCollection&,
			      PtrStepPointMCVectorCollection&,
			      CaloHitSimPartMCCollection& );


    void nonLinearityCorrection(StepPointMC const& h, 
				double& energy, 
				int cryId, 
				ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations,
				GlobalConstantsHandle<ParticleDataTable>& pdt );

    void longitudinalResponseUniformityCorrection(double posZ, 
						  double cryhalflength, 
						  double& energy, 
						  int crid,
						  ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);

    void fillMapById(HitMap& map,HandleVector const& crystalStepsHandles);
    void printStepPointDebug(Calorimeter const & cal,StepPointMC const& h, int crid);


    // Print information about the data products found by the selector functions. 
    void printDataProductInfo( HandleVector const& crystalStepsHandles,
			       HandleVector const& readoutStepsHandles );
  };

  //--------------------------------------------------------------------

  void MakeCaloReadoutHits::beginJob(){}


  void  MakeCaloReadoutHits::produce(art::Event& event) {

    if ( _diagLevel > 0 ) cout << "MakeCaloReadoutHits: produce() begin" << endl;

    static int ncalls(0);
    ++ncalls;

    // Check that calorimeter geometry description exists
    art::ServiceHandle<GeometryService> geom;    
    if( !(geom->hasElement<Calorimeter>()) ) return;
   

    //update condition cache
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime = accPar->deBuncherPeriod;
    _toff.updateMap(event);
   
    // A container to hold the output hits.
    unique_ptr<CaloHitCollection>               caloHits         (new CaloHitCollection);
    unique_ptr<CaloHitMCTruthCollection>        caloMCHits       (new CaloHitMCTruthCollection);
    unique_ptr<CaloCrystalOnlyHitCollection>    caloCrystalMCHits(new CaloCrystalOnlyHitCollection);
    unique_ptr<PtrStepPointMCVectorCollection>  caloMCptrHits    (new PtrStepPointMCVectorCollection);
    unique_ptr<PtrStepPointMCVectorCollection>  caloMCroptrHits  (new PtrStepPointMCVectorCollection);
    unique_ptr<CaloHitSimPartMCCollection>      caloMCSimParts   (new CaloHitSimPartMCCollection);


    // These selectors will select data products with the given instance name, and ignore
    // all other fields of the product ID.
    art::ProductInstanceNameSelector getCrystalSteps(_stepPoints);
    art::ProductInstanceNameSelector getReadoutSteps(_rostepPoints);

    // Get the StepPointMCs from the event.
    HandleVector crystalStepsHandles, readoutStepsHandles;
    event.getMany( getCrystalSteps, crystalStepsHandles);
    event.getMany( getReadoutSteps, readoutStepsHandles);


    static bool firstEvent(true);
    if ( firstEvent ) {
      printDataProductInfo( crystalStepsHandles, readoutStepsHandles);
      firstEvent = false;
      std::cout<<"Timegap is set to "<<_timeGap<<std::endl;
    }
    

    makeCalorimeterHits(crystalStepsHandles, readoutStepsHandles,
                        *caloHits, *caloMCHits, *caloCrystalMCHits,
                        *caloMCptrHits, *caloMCroptrHits, *caloMCSimParts);


    if ( ncalls < _maxFullPrint && _diagLevel > 0 ) {
      cout << "MakeCaloReadoutHits: Total number of calorimeter hits = "
           << caloHits->size()
           << endl;
      cout << "MakeCaloReadoutHits: Total number of crystal MC hits = "
           << caloCrystalMCHits->size()
           << endl;
    }

    
    
    // Add the output hit collection to the event
    event.put(std::move(caloHits));
    event.put(std::move(caloMCHits));
    event.put(std::move(caloCrystalMCHits));
    event.put(std::move(caloMCptrHits),"CaloHitMCCrystalPtr");
    event.put(std::move(caloMCroptrHits),"CaloHitMCReadoutPtr");
    event.put(std::move(caloMCSimParts));

    if ( _diagLevel > 0 ) cout << "MakeCaloReadoutHits: produce() end" << endl;

  } 










  // The simulation is expected to follow the data reconstruction, i.e. read the n readouts for each crystal 
  // (n=2 as of Oct 2013), then combine these n readouts to form a crystal hit.
  // The simulation collects the hits in a crystal, not a readout, so we create first readout hits (this module)


  void MakeCaloReadoutHits::makeCalorimeterHits (HandleVector const& crystalStepsHandles,
                                                 HandleVector const& readoutStepsHandles,
                                                 CaloHitCollection &caloHits,
                                                 CaloHitMCTruthCollection& caloHitsMCTruth,
                                                 CaloCrystalOnlyHitCollection& caloCrystalHitsMCTruth,
                                                 PtrStepPointMCVectorCollection& caloHitsMCCrystalPtr,
                                                 PtrStepPointMCVectorCollection& caloHitsMCReadoutPtr,
						 CaloHitSimPartMCCollection& caloMCSimParts ){
    
    
    
    GlobalConstantsHandle<ParticleDataTable> pdt;
    // Handle to the conditions service
    ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
    
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());
    //    double timeGap    = cal.caloGeomInfo().timeGap();
    double addEdep    = cal.caloGeomInfo().electronEdep();
    double cryhalflength = cal.caloGeomInfo().crystalHalfLength();
    
    
    
    // Fill map: crystal id -> StepMCPtr and RO id -> StepMCPtr
    HitMap hitmapCrystal,hitmapRO;        
    fillMapById( hitmapCrystal, crystalStepsHandles);
    fillMapById( hitmapRO, readoutStepsHandles);

    CaloCrystalMCUtil readoutUtil;
    std::unordered_map<const StepPointMC*, double> tmap;
    
    // First step, we loop over the StepPoints of each crystal and create a Hit for each stepPoint.
    // Faster than looping over the readouts and associating the same crystal hits several times
    for (auto const& crIter : hitmapCrystal) {      	
      
      vector<ROHit> cr_hits;
      int crid = crIter.first;
      StepPtrs const& isteps = crIter.second;
      
      
	
      // Loop over steps inside the crystal for a given crystal id
      for (StepPtrs::const_iterator i=isteps.begin(), e=isteps.end(); i != e; ++i ) {
	
	StepPointMC const& h = **i;

	if ( h.eDep()<=0.0 ) continue;
	double edep_corr = h.eDep();

	// Calculate correction for edep
	CLHEP::Hep3Vector const& posInMu2e = h.position();
        double posZ = cal.toCrystalFrame(crid,posInMu2e).z();
	
	//should probably go to another class when the code for these corrections grows larger
	if (_caloNonLinCorrection && h.simParticle().isNonnull()) {
	  nonLinearityCorrection(h, edep_corr, crid, calorimeterCalibrations, pdt);
	}

	if (_caloLRUcorrection) {
	  longitudinalResponseUniformityCorrection(posZ, cryhalflength, edep_corr,crid, calorimeterCalibrations);
	}
        
	// time folding and Adding ghost hits to properly treat boundary conditions with folding, see docdb-3425
	double hitTimeUnfolded = _toff.timeWithOffsetsApplied(h);
	double hitTime         = fmod(hitTimeUnfolded,_mbtime);
        
        const StepPointMC* hptr = &h;
	tmap[hptr]=hitTime;
	
	// 2015-09-30 P.Murat: given that hits are 'point-like' things, need only the first if

	if (hitTime > _blindTime) {
	  cr_hits.push_back(ROHit(*i,h.eDep(),edep_corr,ROHit::crystal,hitTime));
	}

// 	if (hitTime < _mbbuffer) {
// 	  if (hitTime+_mbtime > _blindTime) {
// 	    cr_hits.push_back(ROHit(*i,h.eDep(),edep_corr,ROHit::crystal,hitTime + _mbtime));
// 	  }
// 	}

// for _mbtime=1695, _mbbuffer=100, and _blindTime=300 the next two if's together are always false 
// no need to comment anything out
	if (hitTime > (_mbtime-_mbbuffer)) {
	  if (hitTime-_mbtime > _blindTime) {
	    cr_hits.push_back(ROHit(*i,h.eDep(),edep_corr,ROHit::crystal,hitTime - _mbtime));
	  }
	}
      }
//-----------------------------------------------------------------------------
// Second, loop over all RO for a given crystal id (Roid = Roidbase + j, j=0,nROPerCrystal-1)
// for each readout, assign hits in the crystal + hits in the readout
//-----------------------------------------------------------------------------
      int ROidBase = cal.ROBaseByCrystal(crid);
      int ROidEnd  = ROidBase+cal.caloGeomInfo().nROPerCrystal();
	
      for (int roid=ROidBase;roid<ROidEnd;++roid) {
	if (cr_hits.size() == 0) continue;

	// copy hits from the crystal
	vector<ROHit> ro_hits(cr_hits);

	//find the entry in RO map and add the RO hits
	HitMap::const_iterator irIter = hitmapRO.find(roid);

	if (irIter != hitmapRO.end() )    {
	  StepPtrs const& irosteps = irIter->second;
	  for( StepPtrs::const_iterator i=irosteps.begin(); i!=irosteps.end(); ++i ) {
	    StepPointMC const& h = **i;
//-----------------------------------------------------------------------------
// There is no cut on energy deposition here - may be, we need to add one?
// time folding and Adding ghost hits to properly treat boundary conditions with folding, see docdb-3425
//-----------------------------------------------------------------------------
	    double hitTimeUnfolded = _toff.timeWithOffsetsApplied(h);
	    double hitTime         = fmod(hitTimeUnfolded,_mbtime);

	    if (hitTime > _blindTime) {
	      ro_hits.push_back(ROHit(*i,0.,0.,ROHit::readout,hitTime));
	    }

// 	    if (hitTime < _mbbuffer) {
// 	      if (hitTime+_mbtime > _blindTime) {
// 		ro_hits.push_back(ROHit(*i,0.,0.,ROHit::readout,hitTime + _mbtime));
// 	      }
// 	    }

// for _mbtime=1695, _mbbuffer=100, and _blindTime=300 the next two if's together are always false 
// no need to comment anything out
	    if (hitTime > (_mbtime-_mbbuffer)) {
	      if (hitTime-_mbtime > _blindTime) {
		ro_hits.push_back(ROHit(*i,0.,0.,ROHit::readout,hitTime - _mbtime));
	      }
	    }
	  }
	}
//-----------------------------------------------------------------------------
// Sort hits by time
//-----------------------------------------------------------------------------
	sort(ro_hits.begin(), ro_hits.end());

	// A buffer to collect output.
	PtrStepPointMCVector mcptr_crystal;
	PtrStepPointMCVector mcptr_readout;
	     
	// A tool to add the art::Ptr to the correct output collection.
	PtrAdder addPtr( mcptr_crystal, mcptr_readout );
	
        

	// Loop over sorted hits and form complete ro/calorimeter hits
	// We will need to digitize the hits here, simulate APD response
	
	double h_time    = ro_hits[0]._time;
	double h_edep    = ro_hits[0]._edep;
	double h_edepc   = ro_hits[0]._edep_corr;
	ROHit::step_type h_type = ro_hits[0]._type;
	addPtr( ro_hits[0] );
	
	for( size_t i=1; i<ro_hits.size(); ++i )  {
	  if( (ro_hits[i]._time - h_time) > _timeGap )   {
	    // Save current hit , remove ghosts
	    if (h_time>0 && h_time < _mbtime)  {
	      if(_diagLevel > 1) std::cout<<"MakeCaloReadoutHits:: insert hit roid="<<roid<<" with time="<< h_time
					  <<"and energy = "<<h_edepc+h_type*addEdep<<std::endl;		      		    
	      
	      CaloHitSimPartMC  caloHitSimPartMC;
	      if (_fillDetailedHit) readoutUtil.fillSimMother(cal,mcptr_crystal,caloHitSimPartMC,tmap);
	      
	      caloHits.push_back(       CaloHit(       roid,h_time,h_edepc+h_type*addEdep));
	      caloHitsMCTruth.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_type));
	      caloHitsMCCrystalPtr.push_back(mcptr_crystal);
	      caloHitsMCReadoutPtr.push_back(mcptr_readout);
	      caloMCSimParts.push_back(caloHitSimPartMC);	
	    }   
	    
	    // ...and create new hit
	    mcptr_crystal.clear();
	    mcptr_readout.clear();
	    addPtr(ro_hits[i]);
	    
	    h_time    = ro_hits[i]._time;
	    h_edep    = ro_hits[i]._edep;
	    h_edepc   = ro_hits[i]._edep_corr;
	    h_type    = ro_hits[i]._type;
	  } 
	  else {
	    // Append data to hit
	    h_edep  += ro_hits[i]._edep;
	    h_edepc += ro_hits[i]._edep_corr;
	    if( ro_hits[i]._type != ROHit::crystal ) h_type = ROHit::readout; // this does not count the charge...
	    addPtr(ro_hits[i]);
	  }
	}

	if(_diagLevel > 1) std::cout<<"MakeCaloReadoutHits:: insert hit roid="<<roid<<" with time="<< h_time
				    <<"and energy = "<<h_edepc+h_type*addEdep<<std::endl;
	
	CaloHitSimPartMC  caloHitSimPartMC;
	if (_fillDetailedHit) readoutUtil.fillSimMother(cal,mcptr_crystal,caloHitSimPartMC,tmap);
	
	
	caloHits.push_back(       CaloHit(roid,h_time,h_edepc+h_type*addEdep));
	caloHitsMCTruth.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_type));
	caloHitsMCCrystalPtr.push_back(mcptr_crystal);
	caloHitsMCReadoutPtr.push_back(mcptr_readout);
	caloMCSimParts.push_back(caloHitSimPartMC);	   
      }
    }
//-----------------------------------------------------------------------------
// Finally, collect the unsmeared crystal hits for each crystal
//-----------------------------------------------------------------------------
    for (auto const& cr: hitmapCrystal){
      CaloCrystalOnlyHitCollection cr_hits;
      int cid = cr.first;
      StepPtrs const& steps = cr.second;

      for( size_t i=0; i<steps.size(); i++ )	{
	StepPointMC const& h2 = *(steps[i]);
	cr_hits.push_back(CaloCrystalOnlyHit(cid,h2.time(),h2.eDep()));
      }

      sort(cr_hits.begin(), cr_hits.end(), lessByTime<CaloCrystalOnlyHitCollection::value_type>());

      // now form final hits if they are close enough in time
      CaloCrystalOnlyHitCollection::value_type cHitMCTruth  = cr_hits[0];
      for ( size_t i=1; i<cr_hits.size(); ++i )	{
	if ( (cr_hits[i].time()-cr_hits[i-1].time()) > _timeGap ) 	   {
	  // Save current hit and create new onw, reject ghost as well	    
	  if (cr_hits[i-1].time()>0) caloCrystalHitsMCTruth.push_back(cHitMCTruth);
	  cHitMCTruth  = cr_hits[i];
	} 
	else {
	  // Add energy to the old hit (keep the "earlier" time)
	  cHitMCTruth.setEnergyDep(cHitMCTruth.energyDep()+cr_hits[i].energyDep());
	}
      }
      caloCrystalHitsMCTruth.push_back(cHitMCTruth);
    }
  } // end makeCalorimeterHits


//-----------------------------------------------------------------------------
  void MakeCaloReadoutHits::fillMapById( HitMap& hitmap,HandleVector const& crystalStepsHandles) {
    for ( HandleVector::const_iterator i=crystalStepsHandles.begin(), e=crystalStepsHandles.end(); i != e; ++i ){
      
      art::Handle<StepPointMCCollection> const& handle(*i);
      StepPointMCCollection const& steps(*handle);
      
      StepPointMCCollection::const_iterator j0=steps.begin();
      for ( StepPointMCCollection::const_iterator j=j0, je=steps.end(); j != je; ++j ){
	StepPointMC const& step(*j);
	hitmap[step.volumeId()].push_back( StepPtr(handle,j-j0) );
      }
    }
  } 


//-----------------------------------------------------------------------------
  void MakeCaloReadoutHits::printDataProductInfo( HandleVector const& crystalStepsHandles,
						  HandleVector const& readoutStepsHandles )
  {
    mf::LogInfo log(_messageCategory);
    log << "MakeCaloReadoutHit::produce will use StepPointMCs from: \n";
    for ( HandleVector::const_iterator i=crystalStepsHandles.begin(), e=crystalStepsHandles.end();
	  i != e; ++i ){
      art::Provenance const& prov(*(i->provenance()));
      log  << "   " << prov.branchName() << "\n";
    }
    log << "\nMakeCaloReadoutHit::produce will use StepPointMCs from: \n";
    for ( HandleVector::const_iterator i=readoutStepsHandles.begin(), e=readoutStepsHandles.end();
	  i != e; ++i ){
      art::Provenance const& prov(*(i->provenance()));
      log  << "   " << prov.branchName() << "\n";
    }
  } 
  
  
//-----------------------------------------------------------------------------
  void MakeCaloReadoutHits::nonLinearityCorrection(StepPointMC const& h, 
						   double& energy, int cryId, 
                                                   ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations, 
						   GlobalConstantsHandle<ParticleDataTable>& pdt  )
  { 
    double edep_save(energy);
    double MeV2keV = 1000.0;
      
    int particleCode = std::abs(h.simParticle()->pdgId());
    if (particleCode!= 11 && particleCode != 22 ) return;
    
    const HepPDT::ParticleData& data = pdt->particle(particleCode).ref();
    double mass = data.mass().value();
    double kinetic_energy = std::sqrt(h.momentum().mag2() + mass*mass) - mass;
    if (kinetic_energy > 1.0) return;
    
    //the formula used returns positive values if tmpEnergy>(approximately) 3keV, see Mu2e doc 1748-v1
    kinetic_energy = MeV2keV*kinetic_energy;
    double trackKine = kinetic_energy;
    kinetic_energy = std::log10(kinetic_energy);
    kinetic_energy *= _randGauss.fire(calorimeterCalibrations->LINpar2(cryId), calorimeterCalibrations->LINpar2Err(cryId) );
    kinetic_energy += _randGauss.fire(calorimeterCalibrations->LINpar1(cryId), calorimeterCalibrations->LINpar1Err(cryId) );
    kinetic_energy = std::log10(kinetic_energy);
    kinetic_energy *= _randGauss.fire(calorimeterCalibrations->LINpar0(cryId), calorimeterCalibrations->LINpar0Err(cryId) );
    kinetic_energy /= _randGauss.fire(calorimeterCalibrations->LINpar3(cryId), calorimeterCalibrations->LINpar3Err(cryId) );
    kinetic_energy /= std::log10(trackKine);
    
    if (kinetic_energy>0) energy *= kinetic_energy;		
    
    if (_diagLevel > 2) {
      std::cout<<"************************** BEFORE / AFTER NON-LINEARITY EFFECT-> edep_corr = "
	       << edep_save<<"  /  "<<energy<<std::endl
	       <<", energyKin = "
	       << trackKine 
	       << ", mass = "<< mass
	       << ", momentum.mag2() = "<< h.momentum().mag2()
	       <<std::endl;
    }
  }
  
//-----------------------------------------------------------------------------
  void MakeCaloReadoutHits::longitudinalResponseUniformityCorrection(double posZ, 
								     double cryhalflength, 
								     double& energy, 
								     int crid, 
								     ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations )
  {
    double edep_save(energy);
    
    posZ *= -_randGauss.fire(calorimeterCalibrations->LRUpar0(crid), calorimeterCalibrations->LRUpar0Err(crid) );
    posZ += 1.0;
    energy *= posZ;
    
    if (_diagLevel > 2) { 
      std::cout<<"***************BEFORE /  AFTER LRU EFFECT-> edep_corr = "
	       << edep_save<<"  /  "<<energy
	       << std::endl;	  
    }
  }


//-----------------------------------------------------------------------------
  void MakeCaloReadoutHits::printStepPointDebug(Calorimeter const & cal,StepPointMC const& h, int crid)
  {
     CLHEP::Hep3Vector testpos = h.position();
     std::cout<<"Reco "<<crid<<"   "<<h.volumeId()<<std::endl;
     std::cout<<"reco position Mu2e    "<<testpos<<std::endl;
     std::cout<<"reco position disk    "<<cal.toSectionFrame(cal.crystal(crid).sectionId(),testpos)<<std::endl;
     std::cout<<"reco position diskFF  "<<cal.toSectionFrameFF(cal.crystal(crid).sectionId(),testpos)<<std::endl;
     std::cout<<"reco position local   "<<cal.toCrystalFrame(crid,testpos)<<std::endl;
     std::cout<<"reco position disk    "<<cal.fromSectionFrame(cal.crystal(crid).sectionId(),cal.toSectionFrame(cal.crystal(crid).sectionId(),testpos))<<std::endl;
     std::cout<<"reco position diskFF  "<<cal.fromSectionFrameFF(cal.crystal(crid).sectionId(),cal.toSectionFrameFF(cal.crystal(crid).sectionId(),testpos))<<std::endl;
     std::cout<<"reco position local   "<<cal.fromCrystalFrame(crid,cal.toCrystalFrame(crid,testpos))<<std::endl;

     std::cout<<"reco Crystal orig     "<<cal.crystal(crid).position()<<std::endl;
     std::cout<<"reco Crystal orig sec "<<cal.crystal(crid).localPosition()<<std::endl;
     std::cout<<"Is inside I           "<<cal.isInsideCalorimeter(testpos)<<std::endl;
     std::cout<<"Is inside II          "<<cal.isInsideCalorimeter(testpos+CLHEP::Hep3Vector(0,-1,0))<<std::endl;
     std::cout<<"Is inside III         "<<cal.isInsideCalorimeter(testpos+CLHEP::Hep3Vector(0,0,-1))<<std::endl;
  }

} // end namespace mu2e

using mu2e::MakeCaloReadoutHits;
DEFINE_ART_MODULE(MakeCaloReadoutHits);



