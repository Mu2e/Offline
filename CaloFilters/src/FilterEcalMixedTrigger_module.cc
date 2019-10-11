#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"

//#include "MCDataProducts/inc/SimParticleCollection.hh"
//#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"

//#include "MCDataProducts/inc/GenParticleCollection.hh"
//#include "DataProducts/inc/VirtualDetectorId.hh"

#include "RecoDataProducts/inc/CaloTrigSeedCollection.hh"

#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "RecoDataProducts/inc/TriggerInfo.hh"


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
  
  class FilterEcalMixedTrigger : public art::EDFilter {
  public:
    explicit FilterEcalMixedTrigger(fhicl::ParameterSet const& pset);
    virtual ~FilterEcalMixedTrigger() { }

    bool filter( art::Event& event);

    //virtual bool beginRun(art::Run &run) override;

    virtual void beginJob() override;
    virtual void endJob() override;

  private:

    int _diagLevel;
    std::string    _trigPath;

    std::string _MVAMethodLabel;
    std::string _caloTrigSeedModuleLabel;
    std::string _ecalweightsfile;
    std::string _mixedweightsfile;
    art::InputTag _shTag;   
    
    float _TOFF; // time offset to align fast clustering with tracker

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

    int nch; // number of combo hits
    const ComboHit*     comboHit;
    const ComboHitCollection* _chcol;


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
    int disk;
  
    // Straw hits unpacking
    std::vector < float > _StrawTimes;
    std::vector < float > _StrawEnes;
    std::vector < float > _StrawXs;
    std::vector < float > _StrawYs;
    std::vector < float > _StrawZs;

    // Peak variables
    // Good shower peaks finding
    std::vector <ecalPeak> _peaklist;
    float xpeak,ypeak;
    float _ecalMVA;
    float ecalMVAcut;
    ecalPeak peak;

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
    float dist;
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

  FilterEcalMixedTrigger::FilterEcalMixedTrigger(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    _diagLevel(pset.get<int>("diagLevel",0)),
    _trigPath(pset.get<std::string>("triggerPath")),
    _MVAMethodLabel(pset.get<std::string>("MVAMethod","BDT")), 
    _caloTrigSeedModuleLabel(pset.get<std::string>("caloTrigSeedModuleLabel")), 
    _ecalweightsfile               (pset.get<std::string>("ecalweightsfile")),
    _mixedweightsfile               (pset.get<std::string>("mixedweightsfile")),
    _shTag(pset.get<art::InputTag>("StrawHitCollectionTag","makePH")),
    _TOFF (pset.get<float>("TimeOFFSET",22.5)),
    _MVArpivot                 (pset.get<float>("MVArpivot",445.)),
    _ecalMVAhighcut0                (pset.get<float>("ecalMVAhighcut0",-0.3)),
    _ecalMVApivotcut0               (pset.get<float>("ecalMVApivotcut0",0.3)),
    _ecalMVAlowcut0                 (pset.get<float>("ecalMVAlowcut0",-0.3)),
    _ecalMVAhighcut1                (pset.get<float>("ecalMVAhighcut1",-0.3)),
    _ecalMVApivotcut1               (pset.get<float>("ecalMVApivotcut1",-0.3)),
    _ecalMVAlowcut1                 (pset.get<float>("ecalMVAlowcut1",-0.3)),
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
    produces<TriggerInfo>();

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


  void FilterEcalMixedTrigger::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("diag");
    
    ConfigFileLookupPolicy configFile;
    // this loads the TMVA library
    TMVA::Tools::Instance();
    //
    ecalreader = new TMVA::Reader( "!Color:Silent" );  
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
    mixedreader = new TMVA::Reader( "!Color:Silent" );  
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

  bool FilterEcalMixedTrigger::filter(art::Event& event) {
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);

    if (_step==0) return false;

    _fevt=(float)event.id().event();

    ++_nProcessed;
    if (_nProcessed%10==0 && _diagLevel > 0) std::cout<<"Processing event from FilterEcalMixedTrigger =  "<<_nProcessed <<std::endl;


    //Handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return false;
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    // ------------------------- Unpack Calo Trig Seeds -------------------------------
    art::Handle<CaloTrigSeedCollection> caloTrigSeedsHandle;
    event.getByLabel(_caloTrigSeedModuleLabel, caloTrigSeedsHandle);
    CaloTrigSeedCollection const& caloTrigSeeds(*caloTrigSeedsHandle);
    if (_step==1) return false;
    // Find Good shower peaks
    _peaklist.clear();

    for (CaloTrigSeedCollection::const_iterator seedIt = caloTrigSeeds.begin(); seedIt != caloTrigSeeds.end(); ++seedIt){
      disk= cal.crystal((int)seedIt->crystalid()).diskId();
      _fdiskpeak   = (float) disk;
      _Epeak   = seedIt->epeak();
      _tpeak   = seedIt->tpeak()+_TOFF;
      _rpeak   = seedIt->rpeak();
      _E10   = seedIt->ring1max();
      _E11   = seedIt->ring1max2();
      _E20   = seedIt->ring2max();
      //
      _ecalMVA= ecalreader->EvaluateMVA(_MVAmethod);
      if (_rpeak>_MVArpivot){
	ecalMVAcut=_ecalMVAlowcut[disk];
      }
      else{
	ecalMVAcut=_ecalMVAcutA[disk]+_ecalMVAcutB[disk]*_rpeak;
      }
      if (_diagLevel > 0){
	std::cout<<"EVENT " << event.id().event() << " DISK " << disk << " Epeak=" << _Epeak << " tpeak=" << _tpeak << " E10=" << _E10 << " E11=" << _E11 << " E20=" << _E20 << " ecalMVA=" << _ecalMVA << " rpeak=" << _rpeak ;
	if (_rpeak>_MVArpivot){
	  std::cout << " cut at " << _ecalMVAlowcut[disk] << std::endl;
	}
	else{
	    ecalMVAcut=_ecalMVAcutA[disk]+_ecalMVAcutB[disk]*_rpeak;
	    std::cout << " cut at " << ecalMVAcut <<  std::endl;
	}
      }
      if (_rpeak>_MVArpivot){
	if (_ecalMVA<_ecalMVAlowcut[disk]) continue;
      }
      else{
	ecalMVAcut=_ecalMVAcutA[disk]+_ecalMVAcutB[disk]*_rpeak;
	if (_ecalMVA<ecalMVAcut) continue;
      }
      // good peak
      peak.disk=disk;
      xpeak= cal.crystal((int)seedIt->crystalid()).localPosition().x();
      ypeak= cal.crystal((int)seedIt->crystalid()).localPosition().y();
      peak.x=xpeak;
      peak.y=ypeak;
      peak.t=_tpeak;
      peak.ene=_Epeak;
      peak.E10=_E10;
      peak.E11=_E11;
      peak.E20=_E20;
      _peaklist.push_back(peak);
    }
  
    if (_peaklist.size()<1) return false;
    
    if (_step==2) return false;


    // ------------------------- Unpack Straw Hits -------------------------------
    // (from ReadStrawHit_module.cc

    auto chcolH = event.getValidHandle<mu2e::ComboHitCollection>(_shTag);
    if (chcolH.product() != 0){
      _chcol= chcolH.product();
      _StrawTimes.clear();
      _StrawEnes.clear();
      _StrawXs.clear();
      _StrawYs.clear();
      _StrawZs.clear();
      nch=_chcol->size();
      for(int istr=0; istr<nch;++istr) {
	comboHit=&_chcol->at(istr);
	_StrawTimes.push_back(comboHit->time());
	_StrawEnes.push_back(comboHit->energyDep());
	_StrawXs.push_back(comboHit->pos().x());
	_StrawYs.push_back(comboHit->pos().y());
	_StrawZs.push_back(comboHit->pos().z());
      }
    }
    else{
      std::cout << "myntuple: ComboHitCollection Handle NOT VALID" << std::endl;
      return false;
    }
    
    if (_step==3) return false;
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
      
      xpeak= peak.x;
      ypeak= peak.y;
      //
      // Find circle centers close to peak 
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
      _Epeak = peak.ene;
      _E10   = peak.E10;
      _E11   = peak.E11;
      _E20   = peak.E20;
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
      if (_step>4){
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
	  if (_mixedMVA>_mixedMVAlowcut[peak.disk]) {
	    //FIX ME!!!!
	    triginfo->_triggerBits.merge(TriggerFlag::caloTrigSeed);
	    triginfo->_triggerBits.merge(TriggerFlag::hitCluster);	    
	    triginfo->_triggerPath = _trigPath;
	    event.put(std::move(triginfo));
    	    return true;
	  }
	}
	else{
	  mixedMVAcut=_mixedMVAcutA[peak.disk]+_mixedMVAcutB[peak.disk]*_rpeak;
	  if (_mixedMVA>mixedMVAcut) {
	    //FIX ME!!!!
	    triginfo->_triggerBits.merge(TriggerFlag::caloTrigSeed);
	    triginfo->_triggerBits.merge(TriggerFlag::hitCluster);	    
	    triginfo->_triggerPath = _trigPath;
	    event.put(std::move(triginfo));	    
	    return true;
	  }
	}
      }
    }
    return false;
  }
  void FilterEcalMixedTrigger::endJob(){
    std::cout << "FilterEcalMixedTrigger filter end job:" << _nProcessed << " events processed" << std::endl;
  }
  
}

using mu2e::FilterEcalMixedTrigger;
DEFINE_ART_MODULE(FilterEcalMixedTrigger);
