#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// Root includes
#include "TDirectory.h"
#include "TTree.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

// C++ includes
#include <vector>

//using namespace std;

namespace mu2e {

  struct EcalTrigDigi{
    int cryId;
    float x;
    float y;
    float phi;
    float tpeak;
    float amp;
  };

  struct ecalPeak{
    int disk;
    float x;
    float y;
    float t;
    float ene;
    float E10;
    float E11;
    float E20;
  };

  struct strawHit{
    float x;
    float y;
    float z;
    float t;
    float dadz;
    float dfdz;
    float dt;
  };
  
  class FilterEcalStrawHitTrigger : public art::EDFilter {
  public:
    explicit FilterEcalStrawHitTrigger(fhicl::ParameterSet const& pset);
    virtual ~FilterEcalStrawHitTrigger() { }

    bool filter( art::Event& event);

    //virtual bool beginRun(art::Run &run) override;

    virtual void beginJob() override;
    virtual void endJob() override;

  private:

    int _diagLevel;

    std::string _MVAMethodLabel;
    std::string _g4ModuleLabel;
    std::string _virtualDetectorLabel;
    std::string _caloDigiModuleLabel;
    std::string _ecalweightsfile;
    std::string _mixedweightsfile;
    art::InputTag _shTag;   
    art::InputTag _shpTag;

    SimParticleTimeOffset _toff;  // time offset smearing

    float                      _ADC2MeV;
    float                      _PEAK2E;
    float                      _DIGITHRESHOLD;
    float                      _T0MIN;
    float                      _PEAKmin;
    float                      _DTmax;
    float                      _wcry;
    float                      _MVArpivot;
    float                      _ecalMVAhighcut0;
    float                      _ecalMVApivotcut0;
    float                      _ecalMVAlowcut0;
    float                      _ecalMVAhighcut1;
    float                      _ecalMVApivotcut1;
    float                      _ecalMVAlowcut1;
    float                      _mixedMVAhighcut0;
    float                      _mixedMVApivotcut0;
    float                      _mixedMVAlowcut0;
    float                      _mixedMVAhighcut1;
    float                      _mixedMVApivotcut1;
    float                      _mixedMVAlowcut1;
    int                        _downscale500_factor;
    int                        _step;
    int                        _nProcessed;

    std::string _MVAmethod;
    static const int nECALDISKs=2;
    float _ecalMVAcutA[nECALDISKs];
    float _ecalMVAcutB[nECALDISKs];
    float _ecalMVAlowcut[nECALDISKs];

    float _mixedMVAcutA[nECALDISKs];
    float _mixedMVAcutB[nECALDISKs];
    float _mixedMVAlowcut[nECALDISKs];

    float R1SQ,R2SQ;
    float dadzbinwidth;
    float dfdzbinwidth;

    TMVA::Reader *ecalreader;
    TMVA::Reader *mixedreader;

    float _rpeak,_tpeak;
    float _Epeak,_E10,_E11,_E20;
    float _fndphi0,_fndphi1;
    float _meanz,_sigmaz;
    float _meanr,_sigmar;
    float _sigmadphi;
    float _fevt,_fdiskpeak; // needed by TMVA reader

    // Virtual hits
    float ptvd,rvd,xcvd,ycvd,pvd,sthvd,cthvd,dphidzvd,betazvd;

    // Calo Digi unpacking
    std::vector < int > _DigiIds;
    std::vector < double > _DigiT0s;
    std::vector < std::vector < int > > _DigiWaves;
  
    // Calo Digi processing
    std::vector <int> _cryIDs;
    std::vector <double> _TPEAKs;
    std::vector <double> _EPEAKs;
    int digi_index;
    int cryID;
    float tnow,tmin;

    // Energy ordered calo digi map
    std::multimap<float,EcalTrigDigi> _DigiMapE[nECALDISKs];
    EcalTrigDigi _EcalTrigDigi;
    int cry_index;
    float amp;
    int disk;

    // Good shower peaks finding
    std::vector <ecalPeak> _peaklist;
    float xpeak,ypeak;
    std::multimap<float,EcalTrigDigi> _Ring1MapE;
    std::multimap<float,EcalTrigDigi> _Ring2MapE;
    float x,y;
    float dist;
    float _ecalMVA;
    float ecalMVAcut;
    ecalPeak peak;
  
    // Straw hits unpacking
    std::vector < float > _StrawTimes;
    std::vector < float > _StrawEnes;
    std::vector < float > _StrawXs;
    std::vector < float > _StrawYs;
    std::vector < float > _StrawZs;

    // Straw hit analysis
    float zpeak;
    float DTPAR0,DTPAR1;
    float DAPAR0,DAPAR1;
    // circle centers close to peak
    float BC,h,x_C,y_C;
    float x10,y10,x11,y11;
    // angle at circumference
    float alphapeak;
    static const int _nDADZBINs= 30;
    int _ndadz[_nDADZBINs];
    // Straw Hits time and position selection
    std::vector <strawHit> _shitlist;
    strawHit shit;
    int straw_index;
    float strawtime;
    float dt;
    int circleok;
    float xstraw,ystraw,zstraw,tstraw;
    // circle radius first estimate
    float alphastraw;
    float dalpha;
    int n2pi;
    float dz;
    float dalphadz;
    int dadzbin;  
    int ndadzmax,dadzbinmax;
    float dalphadzbest;
    float rfirst;
    // circle center
    static const int _nXBINs= 28;
    static const int _nYBINs= 28;
    int _nxy[_nXBINs][_nYBINs];
    float xo,yo;
    float k;
    float xci,yci,rci;
    int xrbin,yrbin;
    int nxymax;
    float xmax,ymax;
    // Straw hits on circle
    std::vector <strawHit> _shitcircle;
    static const int _nDFDZBINs= 60;
    int _ndfdz[_nDFDZBINs];
    float phipeak;
    float r;
    float phistraw;
    float dphi;
    float dphidz;
    int dfdzbin;
    int ndfdzmax,dfdzbinmax;
    float dphidzbest;
    // find circle radius and beta_z
    float rbest;
    float betaz_est;
    // Straw Hit refined time selection
    std::vector <strawHit> _shittime;
    // Find average and sigma for z,r,phi
    int   ndphicpeak;
    int   ndphicmed;
    float meanz2; 
    float meanr2;
    float meandphi; 
    float meandphi2; 
    float rstraw;
    float dphipeak;
    int   ndphigood;
    // Calculate mixed MVA
    float _mixedMVA;
    float mixedMVAcut;
  };

  FilterEcalStrawHitTrigger::FilterEcalStrawHitTrigger(fhicl::ParameterSet const& pset):
    _diagLevel(pset.get<int>("diagLevel",0)),

    _MVAMethodLabel(pset.get<std::string>("MVAMethod","BDT")), 
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
    _virtualDetectorLabel(pset.get<std::string>("virtualDetectorName")),
    _caloDigiModuleLabel      (pset.get<std::string>("caloDigiModuleLabel","CaloDigiFromShower")), 
    _ecalweightsfile               (pset.get<std::string>("ecalweightsfile")),
    _mixedweightsfile               (pset.get<std::string>("mixedweightsfile")),
    _shTag(pset.get<art::InputTag>("StrawHitCollectionTag","makeSH")),
    _shpTag(pset.get<art::InputTag>("StrawHitPositionCollectionTag","MakeStrawHitPositions")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
    _ADC2MeV                   (pset.get<float>("ADC2MeV",0.0076)),
    _PEAK2E                    (pset.get<float>("PEAK2E",6.04)),
    _DIGITHRESHOLD             (pset.get<float>("DigiThreshold",1.)),
    _T0MIN                     (pset.get<float>("T0min",500.)),
    _PEAKmin                   (pset.get<float>("PEAKmin",20.)),
    _DTmax                     (pset.get<float>("DTmax",6.)),
    _wcry                      (pset.get<float>("wcry",34.3)),
    _MVArpivot                 (pset.get<float>("MVArpivot",445.)),
    _ecalMVAhighcut0                (pset.get<float>("ecalMVAhighcut0",0.5)),
    _ecalMVApivotcut0               (pset.get<float>("ecalMVApivotcut0",0.2)),
    _ecalMVAlowcut0                 (pset.get<float>("ecalMVAlowcut0",0.2)),
    _ecalMVAhighcut1                (pset.get<float>("ecalMVAhighcut1",0.5)),
    _ecalMVApivotcut1               (pset.get<float>("ecalMVApivotcut1",0.2)),
    _ecalMVAlowcut1                 (pset.get<float>("ecalMVAlowcut1",0.2)),
    _mixedMVAhighcut0                (pset.get<float>("mixedMVAhighcut0",0.5)),
    _mixedMVApivotcut0               (pset.get<float>("mixedMVApivotcut0",0.2)),
    _mixedMVAlowcut0                 (pset.get<float>("mixedMVAlowcut0",0.2)),
    _mixedMVAhighcut1                (pset.get<float>("mixedMVAhighcut1",0.5)),
    _mixedMVApivotcut1               (pset.get<float>("mixedMVApivotcut1",0.2)),
    _mixedMVAlowcut1                 (pset.get<float>("mixedMVAlowcut1",0.2)),
    _downscale500_factor        (pset.get<float>("downscale500factor",0)),
    _step                       (pset.get<float>("step",10)),
    _nProcessed(0)  
  { 
    _MVAmethod= _MVAMethodLabel + " method"; 
    //
    _ecalMVAcutB[0]=(_ecalMVAhighcut0-_ecalMVApivotcut0)/(395.-_MVArpivot); // high cut is at r=395 mm
    _ecalMVAcutA[0]=_ecalMVApivotcut0-_ecalMVAcutB[0]*_MVArpivot;
    _ecalMVAlowcut[0]= _ecalMVAlowcut0;
    _ecalMVAcutB[1]=(_ecalMVAhighcut1-_ecalMVApivotcut1)/(395.-_MVArpivot); // high cut is at r=395 mm
    _ecalMVAcutA[1]=_ecalMVApivotcut1-_ecalMVAcutB[1]*_MVArpivot;
    _ecalMVAlowcut[1]= _ecalMVAlowcut1;
    //
    _mixedMVAcutB[0]=(_mixedMVAhighcut0-_mixedMVApivotcut0)/(395.-_MVArpivot); // high cut is at r=395 mm
    _mixedMVAcutA[0]=_mixedMVApivotcut0-_mixedMVAcutB[0]*_MVArpivot;
    _mixedMVAlowcut[0]= _mixedMVAlowcut0;
    _mixedMVAcutB[1]=(_mixedMVAhighcut1-_mixedMVApivotcut1)/(395.-_MVArpivot); // high cut is at r=395 mm
    _mixedMVAcutA[1]=_mixedMVApivotcut1-_mixedMVAcutB[1]*_MVArpivot;
    _mixedMVAlowcut[1]= _mixedMVAlowcut1;

    R1SQ=99225.; // 315.*315.
    R2SQ=156025.; //395.*395.
    dadzbinwidth=0.003/_nDADZBINs;
    dfdzbinwidth=0.003/_nDFDZBINs;

  }


  void FilterEcalStrawHitTrigger::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("diag");
    
    ConfigFileLookupPolicy configFile;
    // this loads the TMVA library
    TMVA::Tools::Instance();
    //
    ecalreader = new TMVA::Reader( "!Color:!Silent" );  
    ecalreader->AddVariable("rpeak",&_rpeak);
    ecalreader->AddVariable("tpeak",&_tpeak);
    ecalreader->AddVariable("Epeak",&_Epeak);
    ecalreader->AddVariable("E10",&_E10);
    ecalreader->AddVariable("E11",&_E11);
    ecalreader->AddVariable("E20",&_E20);
    ecalreader->AddSpectator("fdiskpeak",&_fdiskpeak);
    ecalreader->AddSpectator("fevt",&_fevt);  
    ecalreader->BookMVA(_MVAmethod ,configFile(_ecalweightsfile));
    //
    mixedreader = new TMVA::Reader( "!Color:!Silent" );  
    mixedreader->AddVariable("rpeak",&_rpeak);
    mixedreader->AddVariable("tpeak",&_tpeak);
    mixedreader->AddVariable("Epeak",&_Epeak);
    mixedreader->AddVariable("E10",&_E10);
    mixedreader->AddVariable("E11",&_E11);
    mixedreader->AddVariable("E20",&_E20);
    mixedreader->AddVariable("fndphi0",&_fndphi0);
    mixedreader->AddVariable("fndphi1",&_fndphi1);
    mixedreader->AddVariable("meanz",&_meanz);
    mixedreader->AddVariable("sigmaz",&_sigmaz);
    mixedreader->AddVariable("meanr",&_meanr);
    mixedreader->AddVariable("sigmar",&_sigmar);
    mixedreader->AddVariable("sigmadphi",&_sigmadphi);
    mixedreader->AddSpectator("fdiskpeak",&_fdiskpeak);
    mixedreader->AddSpectator("fevt",&_fevt);  
    mixedreader->BookMVA(_MVAmethod ,configFile(_mixedweightsfile));
  }

  bool FilterEcalStrawHitTrigger::filter(art::Event& event) {

    if (_step==0) return false;

    _fevt=(float)event.id().event();

    ++_nProcessed;
    if (_nProcessed%10==0 && _diagLevel > 0) std::cout<<"Processing event from FilterEcalStrawHitTrigger =  "<<_nProcessed <<std::endl;

    //Handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return false;
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    // Microbunch time
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    double _mbtime = accPar->deBuncherPeriod;
    _toff.updateMap(event);
    //--------------------------  Virtual hits --------------------------------  
    if (_diagLevel > 1){
      art::Handle<StepPointMCCollection> vdhits;
      event.getByLabel(_g4ModuleLabel,_virtualDetectorLabel,vdhits);
      
      if (vdhits.isValid()){
	for (auto iter=vdhits->begin(), ie=vdhits->end(); iter!=ie; ++iter){
	  // the virtual hit has the same structure of the Step Point
	  const StepPointMC& hit = *iter;
	  // HITS in CALORIMETER VIRTUAL DETECTORS               
	  if (hit.volumeId()<VirtualDetectorId::EMC_Disk_0_SurfIn || hit.volumeId()>VirtualDetectorId::EMC_Disk_1_EdgeOut) continue;
	  
	  double hitTimeUnfolded = _toff.timeWithOffsetsApplied(hit);
	  double hitTime         = fmod(hitTimeUnfolded,_mbtime);
	  
	  CLHEP::Hep3Vector VDPos = cal.geomUtil().mu2eToTracker(hit.position());
	  
	  
	  if (hit.simParticle()->pdgId()==11 &&  hit.momentum().mag()>90.){
	    std::cout << "*** EVENT " << event.id().event() << " CE VIRTUAL HIT FOUND: x=" << VDPos.x() << " y=" << VDPos.y() << " z=" << VDPos.z() << " t=" << hitTime;
	      
	    ptvd=sqrt(hit.momentum().x()*hit.momentum().x()+hit.momentum().y()*hit.momentum().y());
	    rvd= ptvd/0.3;
	    xcvd=VDPos.x()-rvd*hit.momentum().y()/ptvd;
	    ycvd=VDPos.y()+rvd*hit.momentum().x()/ptvd;
	    pvd=sqrt(hit.momentum().x()*hit.momentum().x()+hit.momentum().y()*hit.momentum().y()+hit.momentum().z()*hit.momentum().z());
	    sthvd=ptvd/pvd;
	    cthvd=sqrt(1-sthvd*sthvd);
	    dphidzvd=sthvd/cthvd/rvd;
	    betazvd=cthvd;
	    std::cout<< " r=" << rvd << " xc=" << xcvd << " yc=" << ycvd << " dphidz=" << dphidzvd << " betaz=" << betazvd << std::endl; 
	
	    break;
	  }
	}
      }   
    }
    //--------------------------  Unpack Calo Digis  --------------------------------    
    //art::Handle<CaloDigiCollection> caloDigisHandle;
    //event.getByLabel(_caloDigiModuleLabel, caloDigisHandle);
    auto caloDigiFlag = event.getValidHandle<mu2e::CaloDigiCollection>(_caloDigiModuleLabel);
    if ( caloDigiFlag.product() == 0){
      //if (! caloDigisHandle.isValid() ){
      std::cout << "FILTERECALMVATRIGGER: CaloDigiHandle NOT VALID" << std::endl;
      return false;
    }
    const CaloDigiCollection& caloDigis = *caloDigiFlag.product();
    //const CaloDigiCollection& caloDigis(*caloDigisHandle);
    _DigiIds.clear();
    _DigiT0s.clear();
    _DigiWaves.clear();
    for (const auto& caloDigi : caloDigis){
      _DigiIds.push_back(caloDigi.roId());
      _DigiT0s.push_back(caloDigi.t0());
      _DigiWaves.push_back(caloDigi.waveform());
    }


    //******************************************************************************
    // Process Calo Digi
 
    CLHEP::Hep3Vector  crystalPos(0);
    _cryIDs.clear();
    _TPEAKs.clear();
    _EPEAKs.clear();
    digi_index=0;
    for (std::vector<std::vector<int> >::iterator it=_DigiWaves.begin(); it!=_DigiWaves.end(); ++it){
      if (_DigiIds.at(digi_index)%2==1){
	digi_index++;
	continue; // consider only even sensors
      }
      cryID=_DigiIds.at(digi_index)/2;
      tmin=(float) _DigiT0s.at(digi_index);
      for (unsigned int i=2;i<it->size()-2;++i){
	if ( it->at(i)*_ADC2MeV*_PEAK2E > _DIGITHRESHOLD){
	  tnow=tmin+(i+0.5)*5.; // 5 ns tick 
	  if ( tnow > _T0MIN - _DTmax) {
	    if ( it->at(i) > it->at(i-1) &&
		 it->at(i) > it->at(i-2) &&
		 it->at(i) > it->at(i+1) &&
		 it->at(i) > it->at(i+2)){
	      _cryIDs.push_back(cryID); 
	      _TPEAKs.push_back((double)tnow);
	      _EPEAKs.push_back((double)it->at(i)*_ADC2MeV*_PEAK2E);
	    }
	  }
	}
      }
      digi_index++;
    }
  
    if (_step==1) return false;
    // Iterate on reconstructed digi
    for (int idisk=0;idisk<nECALDISKs;idisk++){
      _DigiMapE[idisk].clear();
    }
    cry_index=0;
    for (std::vector<int>::iterator it=_cryIDs.begin(); it!=_cryIDs.end(); ++it){
      _tpeak= _TPEAKs.at(cry_index);
      if (_tpeak<_T0MIN){
	cry_index++;
	continue;
      }
      cryID=*it;
      amp= (float) _EPEAKs.at(cry_index);
      //
      _EcalTrigDigi.cryId= cryID;
      crystalPos = cal.crystal(cryID).localPosition();
      _EcalTrigDigi.x=crystalPos.x();
      _EcalTrigDigi.y=crystalPos.y();
      _EcalTrigDigi.phi = atan2(_EcalTrigDigi.y,_EcalTrigDigi.x);
      _EcalTrigDigi.tpeak = _tpeak; 
      _EcalTrigDigi.amp = amp; 
      disk = cal.crystal(cryID).diskId();
      _DigiMapE[disk].insert(std::make_pair(amp,_EcalTrigDigi));
      cry_index++;
    }
    if (_step==2) return false;


    //******************************************************************************
    // Find Good shower peaks
    _peaklist.clear();
    // Loop on disks
    for (int idisk=0;idisk<nECALDISKs;idisk++){
      // Loop on digi map for this disk (amp decreasing ordered)
      for (std::multimap<float,EcalTrigDigi>::reverse_iterator rit=_DigiMapE[idisk].rbegin(); rit!=_DigiMapE[idisk].rend(); ++rit){
	if ( rit->second.amp < _PEAKmin) break; // peak amplitude threshold
	xpeak=rit->second.x;
	ypeak=rit->second.y;
	_rpeak=sqrt(xpeak*xpeak+ypeak*ypeak);
	_tpeak=rit->second.tpeak;
	_Epeak=rit->first;
	// Search crystals in close rings
	_Ring1MapE.clear();
	_Ring2MapE.clear();
	for (std::multimap<float,EcalTrigDigi>::reverse_iterator rit2=_DigiMapE[idisk].rbegin(); rit2!=_DigiMapE[idisk].rend(); ++rit2){
	  //
	  x=rit2->second.x;
	  y=rit2->second.y;
	  //
	  dist=sqrt((x-xpeak)*(x-xpeak)+(y-ypeak)*(y-ypeak));
	  if (dist>2.5*_wcry) continue; // too far
	  if (dist<0.5*_wcry) continue; // same as peak
	  //
	  if (abs(rit2->second.tpeak-_tpeak)>_DTmax) continue; // too early or too late
	  if (dist<1.5*_wcry) _Ring1MapE.insert(std::make_pair(rit2->second.amp,rit2->second));
	  else _Ring2MapE.insert(std::make_pair(rit2->second.amp,rit2->second));
	}
	if (_Ring1MapE.size()>0){
	  std::multimap<float,EcalTrigDigi>::reverse_iterator ring1it=_Ring1MapE.rbegin();
	  _E10=ring1it->first;
	  if (_Ring1MapE.size()>1){
	    std::advance(ring1it,1);
	    _E11=ring1it->first;
	  }
	  else _E11=0.;
	}
	else{
	  _E10=0.;
	  _E11=0.;
	}
	if (_Ring2MapE.size()>0){
	  _E20=_Ring2MapE.rbegin()->first; 
	}
	else _E20=0.;
	_ecalMVA= ecalreader->EvaluateMVA(_MVAmethod);
	if (_rpeak>_MVArpivot){
	  ecalMVAcut=_ecalMVAlowcut[idisk];
	}
	else{
	  ecalMVAcut=_ecalMVAcutA[idisk]+_ecalMVAcutB[idisk]*_rpeak;
	}
	if (_diagLevel > 0){
	  std::cout<<"EVENT " << event.id().event() << " DISK " << idisk << " Epeak=" << _Epeak << " tpeak=" << _tpeak << " E10=" << _E10 << " E11=" << _E11 << " E20=" << _E20 << " ecalMVA=" << _ecalMVA << " rpeak=" << _rpeak ;
	  if (_rpeak>_MVArpivot){
	    std::cout << " cut at " << _ecalMVAlowcut[idisk] << std::endl;
	  }
	  else{
	    ecalMVAcut=_ecalMVAcutA[idisk]+_ecalMVAcutB[idisk]*_rpeak;
	    std::cout << " cut at " << ecalMVAcut <<  std::endl;
	  }
	}
	if (_rpeak>_MVArpivot){
	  if (_ecalMVA<_ecalMVAlowcut[idisk]) continue;
	}
	else{
	  ecalMVAcut=_ecalMVAcutA[idisk]+_ecalMVAcutB[idisk]*_rpeak;
	  if (_ecalMVA<ecalMVAcut) continue;
	}
	// good peak
	peak.disk=idisk;
	peak.x=xpeak;
	peak.y=ypeak;
	peak.t=_tpeak;
	peak.ene=_Epeak;
	peak.E10=_E10;
	peak.E11=_E11;
	peak.E20=_E20;
	_peaklist.push_back(peak);
      }
    }
    if (_peaklist.size()<1) return false;


    // ------------------------- Unpack Straw Hits -------------------------------
    // (from ReadStrawHit_module.cc

    /*
      art::Handle<StrawHitCollection> shHandle;
      event.getByLabel(_shTag,shHandle);
      art::Handle<mu2e::StrawHitPositionCollection> shpHandle;
      event.getByLabel(_shpTag,shpHandle);
      if (! (shHandle.isValid() && shpHandle.isValid())){
      if (! shHandle.isValid()) std::cout << "FILTERECALMVATRIGGER: StrawHit Handle NOT VALID" << std::endl;
      if (! shpHandle.isValid()) std::cout << "FILTERECALMVATRIGGER: StrawHitPos Handle NOT VALID" << std::endl;
      return false;
      }
      const StrawHitCollection& strawHits(*shHandle);
      const StrawHitPositionCollection& strawHitPositions(*shpHandle);
    */
    auto strawHitFlag = event.getValidHandle<mu2e::StrawHitCollection>(_shTag);
    auto strawHitPosFlag = event.getValidHandle<mu2e::StrawHitPositionCollection>(_shpTag);
    if (strawHitFlag.product() == 0 || strawHitPosFlag.product() == 0 ){
      if (strawHitFlag.product() == 0) std::cout << "FILTERECALMVATRIGGER: StrawHit Handle NOT VALID" << std::endl;
      if (strawHitPosFlag.product() == 0) std::cout << "FILTERECALMVATRIGGER: StrawHitPos Handle NOT VALID" << std::endl;
      return false;
    }
    const StrawHitCollection& strawHits = *strawHitFlag.product();
    const StrawHitPositionCollection& strawHitPositions = *strawHitPosFlag.product();
    _StrawTimes.clear();
    _StrawEnes.clear();
    _StrawXs.clear();
    _StrawYs.clear();
    _StrawZs.clear();
    if ( strawHits.size()==strawHitPositions.size()){
      for (const auto& strawHit : strawHits){
	_StrawTimes.push_back(strawHit.time());
	_StrawEnes.push_back(strawHit.energyDep());
      }
      for (const auto& strawHitPosition : strawHitPositions){
	_StrawXs.push_back(strawHitPosition.pos().x());
	_StrawYs.push_back(strawHitPosition.pos().y());
	_StrawZs.push_back(strawHitPosition.pos().z());
      }
    }
      
    //-------------------------- STRAW HIT ANALYSIS --------------------------------
  
    // Loop on good shower peaks
    for (std::vector<ecalPeak>::iterator it=_peaklist.begin(); it!=_peaklist.end(); ++it){
      peak = *it;
      _fdiskpeak=(float) peak.disk;
      if (peak.disk==0){
	zpeak=1645;
	DTPAR0=6.08;
	DTPAR1=-0.00370;
	DAPAR0=-8.01;
	DAPAR1=0.0048;
      }
      else{
	zpeak=2345.;
	DTPAR0=9.08;
	DTPAR1=-0.00430;
	DAPAR0=-11.4;
	DAPAR1=0.0048;
      }
      //
      // Find circle centers close to peak 
      xpeak=peak.x;
      ypeak=peak.y;
      _rpeak=sqrt(xpeak*xpeak+ypeak*ypeak);
      BC= (R2SQ-R1SQ+_rpeak*_rpeak)/2./_rpeak;
      h= sqrt(R2SQ-BC*BC);
      x_C=xpeak/_rpeak*BC;
      y_C=ypeak/_rpeak*BC;
      x10=x_C-h*ypeak/_rpeak;
      y10=y_C+h*xpeak/_rpeak;
      x11=x_C+h*ypeak/_rpeak;
      y11=y_C-h*xpeak/_rpeak;
      // angle at circumference
      alphapeak=atan2(ypeak,xpeak);
      for (int ibin=0;ibin<_nDADZBINs;ibin++){
	_ndadz[ibin]=0;
      }
      // peak time
      _tpeak=peak.t;
      // recover peak variables needed for MVA calulation
      _Epeak=peak.ene;
      _E10=peak.E10;
      _E11=peak.E11;
      _E20=peak.E20;
      // Straw Hits time and position selection
      _shitlist.clear();
      straw_index=0;
      for (std::vector<float>::iterator it=_StrawTimes.begin(); it!=_StrawTimes.end(); ++it){
	strawtime=_StrawTimes.at(straw_index);
	dt = _tpeak - strawtime;
	zstraw=_StrawZs.at(straw_index);
	if (abs(dt-DTPAR0-DTPAR1*zstraw)<25.){
	  // POSITION SELECTION
	  circleok=0;
	  xstraw=_StrawXs.at(straw_index);
	  ystraw=_StrawYs.at(straw_index);
	  if (sqrt((xstraw-x10)*(xstraw-x10)+(ystraw-y10)*(ystraw-y10))<315.) circleok++;
	  if (sqrt((xstraw-x11)*(xstraw-x11)+(ystraw-y11)*(ystraw-y11))<315.) circleok++;
	  if (circleok==1){
	    shit.x=xstraw;
	    shit.y=ystraw;
	    shit.z=zstraw;
	    shit.t=strawtime;
	    // ANGULAR VELOCITY SELECTION (ANGLES AT THE CIRCUMFERENCE)
	    alphastraw=atan2(ystraw,xstraw);
	    dalpha=2.*(alphastraw-alphapeak);
	    if (dalpha>0) dalpha-=6.2832;
	    n2pi=(int)((dalpha-DAPAR0-DAPAR1*zstraw)/6.2832+1.5);
	    dalpha-=6.2832*(float)(n2pi-1);
	      //
	    dz=zstraw-zpeak;
	    dalphadz=dalpha/dz;
	    dadzbin=(int)((dalphadz-0.003)/0.0001);
	    if (dadzbin>=0 && dadzbin<_nDADZBINs) _ndadz[dadzbin]++;
	    shit.dadz=dalphadz;
	    shit.dfdz=0.;
	    shit.dt=100.;
	    _shitlist.push_back(shit);
	  }
	}
	straw_index++;
      }
      if (_diagLevel > 1) std::cout<< "hits in space-time=" << _shitlist.size() << std::endl;
      // Radius First Estimate
      ndadzmax=0;
      dadzbinmax=-1;
      for (int ibin=0;ibin<_nDADZBINs;ibin++){
	if (_ndadz[ibin]>ndadzmax){
	  dadzbinmax=ibin;
	  ndadzmax=_ndadz[ibin];
	}
      }
      if (_diagLevel > 1) std::cout<< "ndadzmax=" << ndadzmax << std::endl;
      if (ndadzmax>0){
	dalphadzbest=0.003+(dadzbinmax+0.5)*dadzbinwidth;
	rfirst= sqrt(103.3*103.3/0.3/0.3-1./dalphadzbest/dalphadzbest);
      }
      else continue;
      if (_diagLevel > 1) std::cout<< "dalphadzbest=" <<  dalphadzbest << " rfirst=" << rfirst << std::endl;
      // Find Center
      for (int ibin=0;ibin<_nXBINs;ibin++){
	for (int jbin=0;jbin<_nYBINs;jbin++){
	  _nxy[ibin][jbin]=0;
	}
      }
      for (std::vector<strawHit>::iterator it=_shitlist.begin(); it!=_shitlist.end(); ++it){
	shit=*it;
	dalphadz=shit.dadz;
	// preselect the straw hits to determine the best center
	if (abs(dalphadz-dalphadzbest)<0.00015){
	  xstraw=shit.x;
	  ystraw=shit.y;
	  // FIND CENTER CANDIDATES
	  xo=(xpeak+xstraw)/2.;
	  yo=(ypeak+ystraw)/2.;
	  dist=sqrt((xstraw-xpeak)*(xstraw-xpeak)+(ystraw-ypeak)*(ystraw-ypeak));
	  k=rfirst*rfirst-dist*dist/4.; 
	  if (k>0.0001){
	    k=sqrt(k);
	    if (xo*(ystraw-ypeak)-yo*(xstraw-xpeak)>0){
	      xci=xo-k*(ystraw-ypeak)/dist;
	      yci=yo+k*(xstraw-xpeak)/dist;
	    }
	    else{
	      xci=xo+k*(ystraw-ypeak)/dist;
	      yci=yo-k*(xstraw-xpeak)/dist;
	    }
	    rci=sqrt(xci*xci+yci*yci);
	    if (abs(rci-rfirst)<90.){ // 1 cm TOLERANCE
	      if (_diagLevel > 2)  std::cout << "xci=" << xci << " yci=" << yci << std::endl;
	      xrbin=(int)((xci+700.)/50.);
	      if (xrbin>=0 && xrbin<_nXBINs){
		yrbin=(int)((yci+700.)/50.);
		if (yrbin>=0 && yrbin<_nYBINs){
		  _nxy[xrbin][yrbin]++;
		}
	      }
	    }
	  }
	}
      }
      //  
      nxymax=0;
      xrbin=-1;
      yrbin=-1;
      for (int ibin=0;ibin<_nXBINs;ibin++){
	for (int jbin=0;jbin<_nYBINs;jbin++){
	  if (_nxy[ibin][jbin]>nxymax){
	    nxymax=_nxy[ibin][jbin];
	    xrbin=ibin;
	    yrbin=jbin;
	  }
	}
      }
      if (nxymax>0){
	xmax=-675.+xrbin*50.;
	ymax=-675.+yrbin*50.;
      }
      else continue;
      if (_diagLevel > 1) std::cout << "xrbin=" << xrbin << " yrbin=" << yrbin << " xmax=" << xmax << " ymax=" << ymax << std::endl; 
	   
      // Add Straw Hits To Circle List (with dphi/dz)
      _shitcircle.clear();
      for (int ibin=0;ibin<_nDFDZBINs;ibin++){
	_ndfdz[ibin]=0;
      }
      phipeak=atan2(ypeak-ymax,xpeak-xmax);
      if (_diagLevel > 1) std::cout << "phipeak=" << phipeak << std::endl;
      for (std::vector<strawHit>::iterator it=_shitlist.begin(); it!=_shitlist.end(); ++it){
	shit=*it;
	xstraw=shit.x;
	ystraw=shit.y;
	r=sqrt((xstraw-xmax)*(xstraw-xmax)+(ystraw-ymax)*(ystraw-ymax));
	// DELTA R SELECTION
	if (abs(r-rfirst)<50.){ // circle selection
	  //
	  // ANGULAR VELOCITY SELECTION (ANGLES AT THE CENTER)
	  zstraw=shit.z;
	  phistraw=atan2(ystraw-ymax,xstraw-xmax);
	  dphi=phistraw-phipeak;
	  if (dphi>0) dphi-=6.2832;
	  n2pi=(int)((dphi-DAPAR0-dalphadzbest*zstraw)/6.2832+1.5);
	  dphi-=6.2832*(float)(n2pi-1);
	  dz=zstraw-zpeak;
	  dphidz=dphi/dz;
	  dfdzbin=(int)((dphidz-0.003)/dfdzbinwidth);
	  if (dfdzbin>=0 && dfdzbin<_nDFDZBINs) _ndfdz[dfdzbin]++;
	  shit.dfdz=dphidz;
	  shit.dt=100.; // dummy value
	  _shitcircle.push_back(shit);
	}
      }    
      if (_diagLevel > 1) std::cout << "hits in circle:" <<  _shitcircle.size() << std::endl;
      if (_shitcircle.size()<1) continue;
      // Find best dphi/dz
      ndfdzmax=0;
      dfdzbinmax=-1;
      for (int ibin=0;ibin<_nDFDZBINs;ibin++){
	if (_ndfdz[ibin]>ndfdzmax){
	  dfdzbinmax=ibin;
	  ndfdzmax=_ndfdz[ibin];
	}
      }
      if (_diagLevel > 1) std::cout << "ndfdzmax=" <<  ndfdzmax << std::endl;
      if (ndfdzmax>0){
	dphidzbest=0.003025+dfdzbinmax*0.00005;
      }
      else continue;
      if (_diagLevel > 1) std::cout << "dphidzbest=" <<  dphidzbest << std::endl;
      //
      // Estimate beta_Z
      rbest= sqrt(103.3*103.3/0.3/0.3-1./dphidzbest/dphidzbest);
      betaz_est=0.3/(dphidzbest)/103.3;
      if (peak.disk==0) betaz_est+=0.02;
      else betaz_est+=0.015;
      if (_diagLevel > 1) std::cout << "rbest=" <<  rbest << " betaz_est=" << betaz_est <<  std::endl;
      // Refined time selection
      _shittime.clear();
      for (std::vector<strawHit>::iterator it=_shitcircle.begin(); it!=_shitcircle.end(); ++it){
	shit=*it;
	zstraw=shit.z;
	dz=zpeak-zstraw;
	tstraw=shit.t;
	dt= tstraw-(_tpeak-dz/betaz_est/300.);
	if (dt>-10. && dt<20.){
	  shit.dt=dt;
	  _shittime.push_back(shit);
	}
      }
      if (_diagLevel > 1) std::cout << "hits in time:" <<  _shittime.size() << std::endl;
      // Find average and sigma for z,r,phi
      ndphicpeak=0;
      ndphicmed=0;
      _meanz=0.;
      meanz2=0.;
      _meanr=0.;
      meanr2=0.;
      meandphi=0.;
      meandphi2=0.;
      for (std::vector<strawHit>::iterator it=_shittime.begin(); it!=_shittime.end(); ++it){
	shit=*it;
	dphidz=shit.dfdz;
	xstraw=shit.x;
	ystraw=shit.y;
	rstraw=sqrt(xstraw*xstraw+ystraw*ystraw);
	zstraw=shit.z;
	dphipeak=(dphidz-dphidzbest)*(zstraw-zpeak);
	if (dphipeak<-3.1416) dphipeak+=3.1416;
	if (dphipeak>3.1416) dphipeak-=3.1416;
	dphi=dphidz*(zstraw-zpeak);
	n2pi=(int)(dphi/6.2832+7.5);
	dphi-=6.2832*(float)(n2pi-7);
	if (abs(dphipeak)< 0.08){
	  ndphicpeak++;
	}
	if (abs(dphipeak)> 0.08 && abs(dphipeak)<0.2){
	  ndphicmed++;
	}
	if (abs(dphipeak)<0.2){
	  _meanz+=zstraw;
	  meanz2+=zstraw*zstraw;
	  _meanr+=rstraw;
	  meanr2+=rstraw*rstraw;
	  meandphi+=dphi;
	  meandphi2+=dphi*dphi;
	}
      }
      ndphigood=ndphicpeak+ndphicmed;
      if (_diagLevel > 1) std::cout << "ndphigood=" <<  ndphigood << std::endl;
      
      if (ndphigood>0){
	_meanz/=ndphigood;
	_meanr/=ndphigood;
	meandphi/=ndphigood;
      }
      else continue;
      if (ndphigood>1 && 
	  meanz2>_meanz*_meanz*ndphigood &&
	  meanr2>_meanr*_meanr*ndphigood &&
	  meandphi2>meandphi*meandphi*ndphigood ){
	_sigmaz=sqrt((meanz2-_meanz*_meanz*ndphigood)/(ndphigood-1));
	_sigmar=sqrt((meanr2-_meanr*_meanr*ndphigood)/(ndphigood-1));
	_sigmadphi=sqrt((meandphi2-meandphi*meandphi*ndphigood)/(ndphigood-1));
      }
      else continue;
      _fndphi0=(float) ndphicpeak;
      _fndphi1=(float) ndphicmed;
      //
      if (_step>3){
	_mixedMVA= mixedreader->EvaluateMVA(_MVAmethod);
	if (_diagLevel > 0){
	  std::cout<<"EVENT " << event.id().event() << " DISK " << peak.disk << " Epeak=" << _Epeak << " tpeak=" << _tpeak << " E10=" << _E10 << " E11=" << _E11 << " E20=" << _E20 << " meanz=" << _meanz << " sigmaz=" << _sigmaz << " meanr=" << _meanr << " sigmar=" << _sigmar << " sigmadphi=" << _sigmadphi << " ndphicpeak=" << ndphicpeak << " ndphicmed=" << ndphicmed << " mixedMVA=" << _mixedMVA << " rpeak=" << _rpeak;
	  if (_rpeak>_MVArpivot){
	    std::cout << " cut at " << _mixedMVAlowcut[peak.disk] << std::endl;
	  }
	  else{
	    mixedMVAcut=_mixedMVAcutA[peak.disk]+_mixedMVAcutB[peak.disk]*_rpeak;
	    std::cout << " cut at " << mixedMVAcut <<  std::endl;
	  }
	}
	if (_rpeak>_MVArpivot){
	  if (_mixedMVA>_mixedMVAlowcut[peak.disk]) return true;
	}
	else{
	  mixedMVAcut=_mixedMVAcutA[peak.disk]+_mixedMVAcutB[peak.disk]*_rpeak;
	  if (_mixedMVA>mixedMVAcut) return true;
	}
      }
    }
    return false;
  }
  void FilterEcalStrawHitTrigger::endJob(){
    std::cout << "FilterEcalStrawHitTrigger filter end job:" << _nProcessed << " events processed" << std::endl;
  }
  
}

using mu2e::FilterEcalStrawHitTrigger;
DEFINE_ART_MODULE(FilterEcalStrawHitTrigger);
