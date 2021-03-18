// ---------------------------------------------------
//  Filter for selecting Calo Cosmic MIP-like events
//  From CaloCosmicCalib by clusters to cells level
//  V1 version: S.Miscetti 28/January/2021
//  from original version of G. Pezzullo
// -------------------------------------------------
// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"
#include "fhiclcpp/ParameterSet.h"
//#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
//#include "GeometryService/inc/DetectorSystem.hh"
//
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

using namespace CLHEP;
// c++
#include <string>
#include <vector>
#include <iostream>
#include <memory>

using namespace std;

namespace mu2e
{
  class CaloCosmicCalib : public art::EDFilter
  {
  public:
    struct Config
    {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<std::string>  cltag            { Name("CaloClusterModuleLabel"), Comment("Calo cluster collection name ") }; 
      fhicl::Atom<int>          minncrystalhits  { Name("MinNCrystalHits"),        Comment("Minimum Number of crystal hits") }; 
      fhicl::Atom<double>       minenergy        { Name("MinEnergy"),              Comment("Minimum Energy for the clusters") }; 
      fhicl::Atom<double>       maxenergy        { Name("MaxEnergy"),              Comment("Maximum Energy for the clusters") }; 
      fhicl::Atom<std::string>  trigPath         { Name("triggerPath"),            Comment("Name of the trigger path") }; 
      fhicl::Atom<int>          debug            { Name("debugLevel"),             Comment("Debug Level") }; 
      fhicl::Atom<int>          mincelinout      { Name("MinNumCelinout"),         Comment("Minimum Number of crystals for IN-OUT algo") }; 
      fhicl::Atom<int>          minceldiag       { Name("MinNumCeldiagver"),       Comment("Minimum Number of crystals for DIAG-VER algo") }; 
      fhicl::Atom<float>        mincelcut        { Name("MinCelEneCut"),           Comment("Minimum Energy for the crystals") }; 
      fhicl::Atom<float>        outradius        { Name("OutRadius"),              Comment("Value of the Outer radius considered") }; 
      fhicl::Atom<float>        innradius        { Name("InnerRadius"),            Comment("Value of the Inner radius considered") };
      fhicl::Atom<float>        diagminval       { Name("diagminval"),             Comment("Min slope val for DIAG algo") }; 
      fhicl::Atom<float>        diagmaxval       { Name("diagmaxval"),             Comment("Max slope value for DIAG algo") };
      fhicl::Atom<float>        dyminval         { Name("dyminval"),               Comment("Min value of Y-depths for OUT-OUT algo   ") };     
    };

    explicit CaloCosmicCalib(const art::EDFilter::Table<Config>& config);
    virtual bool filter(art::Event& event) override;
    virtual bool beginRun(art::Run&   run   );
    virtual bool endRun( art::Run& run ) override;

  private:
    art::InputTag _cltag;
    int           _minncrystalhits;
    double        _minenergy, _maxenergy;
    std::string   _trigPath;
    int           _debug;
    int           _mincelinout;
    int           _minceldiag;
    float         _mincelcut;
    float         _outradius;
    float         _innradius;
    float         _crystalsize;
    float         _diagminval;
    float         _diagmaxval;
    float         _dyminval;
    // counters
    unsigned _nevt, _npass;
    const Calorimeter* _calogeom;
  };

  CaloCosmicCalib::CaloCosmicCalib(const art::EDFilter::Table<Config>& config) :
     EDFilter{config},
    _cltag          (config().cltag()),
    _minncrystalhits(config().minncrystalhits()),      
    _minenergy      (config().minenergy()),
    _maxenergy      (config().maxenergy()),
    _trigPath       (config().trigPath()),
    _debug          (config().debug()),
    _mincelinout    (config().mincelinout()),
    _minceldiag     (config().minceldiag()),
    _mincelcut      (config().mincelcut()),
    _outradius      (config().outradius()),
    _innradius      (config().innradius()),
    _diagminval     (config().diagminval()),
    _diagmaxval     (config().diagmaxval()),
    _dyminval       (config().dyminval()),     
    _nevt(0), _npass(0)
  {
    produces<TriggerInfo>();
  }

  bool CaloCosmicCalib::beginRun(art::Run & run){
    GeomHandle<Calorimeter> ch;    
    _calogeom    = ch.get();
    _crystalsize = _calogeom->caloInfo().getDouble("crystalXYLength") +  2.*_calogeom->caloInfo().getDouble("wrapperThickness");

    return true;
  }

  bool CaloCosmicCalib::filter(art::Event& evt){
    // create output
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool retval(false); // preset to fail
    // find the collection
    auto clH = evt.getValidHandle<CaloClusterCollection>(_cltag);
    const CaloClusterCollection* clcol = clH.product();
    size_t trig_ind(0);
    
    const CaloHit* hitcalo(0);
    const CaloHitPtrVector* caloClusterHits(0);
    
    // loop over the collection: if any pass the selection, pass this event
    
    int numcluster=0;
    
    for(auto icl = clcol->begin();icl != clcol->end(); ++icl) {
      auto const& cl = *icl;
      
      // get energy and cluster size
      float energy     = cl.energyDep();
      int   clsize     = cl.size();
      
      if(_debug > 2){
        std::cout << moduleDescription().moduleLabel() << " nhits = " << cl.size() << " energy = " << energy << std::endl;
      }
      if( (energy >= _minenergy) && 
	  (energy <= _maxenergy) && 
	  (clsize >= _minncrystalhits) ) {
	// ----------------------------------------------
	// Loop over the hits of the selected cluster
	// -----------------------------------------------
	caloClusterHits = &cl.caloHitsPtrVector();

	int selflag(0), NumMipCel(0);
	float Eouter(0), Einner(0);
	float DXmax, DYmax;
	float RMax(-999), RMin(-999), XMin(999.), XMax(-999.), YoutMin(999.), YoutMax(-999.);
	
	float xval,yval,rval;
	for( int icry=0; icry< clsize; icry++){
	  hitcalo = &(*caloClusterHits->at(icry));	  
	  float ene   = hitcalo->energyDep();
	  float time  = hitcalo->time();	  
	  int _cryid = hitcalo->crystalID();	  
	  int DiskId = _calogeom->crystal(_cryid).diskID();	
	  
	  CLHEP::Hep3Vector crystalpos = _calogeom->geomUtil().mu2eToDiskFF(DiskId,_calogeom->crystal(_cryid).position());
	  
	  if( ene> _mincelcut){         // Use cells with energy above a MIP-like threshold
	    xval = crystalpos.x();
	    yval = crystalpos.y();
	    rval = sqrt(pow(xval,2)+ pow(yval,2));
	    // ----  find the cell with the highest and lowest Y ---------
	    if( yval >  YoutMax ){ // Ymax
	      YoutMax = yval;
	      XMax = xval;
	      RMax = rval; } 	   
	    if( yval < YoutMin ){  // Ymin
	      YoutMin = yval;
	      XMin = xval;
	      RMin = rval; }
	    // calculate Outer/Inner Ring Energies
	    
	    if( rval > _outradius ) Eouter = Eouter+ene; // Oouter ring 
	    if( rval < _innradius ) Einner = Einner+ene; // Inner  ring 
	    NumMipCel++; // Increment number of MIP cells
	    
	  }
	  
	  if(_debug>2) std::cout << moduleDescription().moduleLabel() << " Nc " << numcluster << " cel " << icry << " T-E " << time  << " " << ene << " icry " <<   _cryid  <<  " " << DiskId << " X-Y " << xval << " " << yval << std::endl;
	  
	  
	} // end loop over hits/cluster
	//====================================
	// Now fill up the selection flag  ...
	//====================================
	if( Eouter>_mincelcut && Einner>_mincelcut && NumMipCel>_mincelinout ) selflag=selflag+1;   // Out-IN or IN-Out
	
	if( NumMipCel>_minceldiag ){ // 	  
	  DXmax = abs(XMax-XMin);
	  DYmax = abs(YoutMax-YoutMin);
	  float DiagVal=0.;
	  if( DXmax> 0 ) DiagVal = DYmax/DXmax;	  
	  if( DXmax  < _crystalsize  ) selflag = selflag+100;  // Vertical selection, cell size+1mm
	  if( DiagVal > _diagminval && DiagVal < _diagmaxval )     // vertical slope boundaries (1.5, 2.5)
	    selflag = selflag+200;                     // Diagonal selection
	  if( DYmax > _dyminval && RMax > _outradius && RMin > _outradius ) // at least 4 crystals vertical depth
	    selflag=selflag+10;                        // Outer-Outer
	}
	
	if( selflag > 0){
	  if(_debug>1) std::cout << moduleDescription().moduleLabel() << " SELFLAG:  " << selflag << " EO-EI-NM " << Eouter << " " << Einner << " " << NumMipCel << std::endl;
	}
	    	  
	//---------------------------------------------------
	if(selflag>0){
	  retval = true;
	  ++_npass;
	  // Fill the trigger info object
	  if (trig_ind == 0){
	    triginfo->_triggerBits.merge(TriggerFlag::caloCalib);
	    triginfo->_triggerPath = _trigPath;
	  }
	  //
	  // associate to the caloClusters which trigger
	  //
	  size_t index = std::distance(clcol->begin(),icl);
	  triginfo->_caloClusters.push_back(art::Ptr<CaloCluster>(clH,index));   
	  ++trig_ind;
	  
	  if(_debug > 1){
	    std::cout << moduleDescription().moduleLabel() << " passed event " << evt.id() << std::endl;
	  }	  
	}  // End Loop over selflag by cell-cuts
      }    // End Loop over Cosmic-cluster "selected" by averag
      numcluster++;
    }    // End Loop over cluster objects
    
    evt.put(std::move(triginfo));
    return retval;
  }
  
  bool CaloCosmicCalib::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      std::cout << moduleDescription().moduleLabel() << " passed " <<  _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << std::endl;
    }
    return true;
  }
}
using mu2e::CaloCosmicCalib;
DEFINE_ART_MODULE(CaloCosmicCalib);
