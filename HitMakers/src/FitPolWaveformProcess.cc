///////////////////////////////////////////////////////////////////////////////
// class to process calo-digitized waveform
// 
///////////////////////////////////////////////////////////////////////////////
#include "HitMakers/inc/FitPolWaveformProcess.hh"
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Run.h"


#include "ConditionsService/inc/ConditionsHandle.hh"

#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
// #include "TROOT.h"
// #include "TFolder.h"
// #include "TDirectory.h"


namespace mu2e {

//-----------------------------------------------------------------------------
  FitPolWaveformProcess::FitPolWaveformProcess(fhicl::ParameterSet const& PSet) :

    WaveformProcess(PSet),
    _debugLevel       (PSet.get<int>   ("debugLevel")),
    _acquisitionLenght(PSet.get<int>   ("acquisitionLenght")),
    _digiSampling     (PSet.get<int>   ("digiSampling")),
    _wave_point_error (PSet.get<double>("pulseErrorSample")),
    _psdThreshold     (PSet.get<double>("psdThreshold")),
    _timeFraction     (PSet.get<double>("timeFraction")),
    _fitThresholdMin  (PSet.get<double>("fitThresholdMin")),
    _fitThresholdMax  (PSet.get<double>("fitThresholdMax"))
  {
    TFile*f = TFile::Open("/mu2e/data/users/gianipez/test-CsI-2015-10-06.root");
    _refPulse = (TH1F*) f->Get("CsIPulse");
    f->Close();

    _fWave    = new TF1("fWave" , "pol04" , 0., 2000.);
    _flogn    = new TF1("flogn" , logn , 0., 2000., 4);
    _fdlogn   = new TF1("fdlogn", dlogn, 0., 2000., 6);


    _wave     = new TH1F("wave", "wave", int(_acquisitionLenght/_digiSampling), 0, _acquisitionLenght);

    _counter        = 0;
    _debugHistIndex = 0;
  }


//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  FitPolWaveformProcess::~FitPolWaveformProcess() {}

  void        FitPolWaveformProcess::psdFromChi2(){
    _psd  = 0;

    double   tpeak     = _refPulse->GetBinCenter(_refPulse->GetMaximumBin());
    double   tshift    =  _wave->GetBinCenter(_wave->GetMaximumBin()) - tpeak;

    double   wfMaximum = _wave->GetMaximum();

    double   pcChi2 = 0; 
    double   errChi2(10.); 

    //    double   wnorm = _wave->Integral();
    double   shapeAmp, wfAmp;
    double   binWidth = _refPulse->GetBinWidth(1);
    int      shapeBin;

    int      nSamples = _wave->GetNbinsX();     

    for(int ibin = 0; ibin<nSamples; ibin++){
      wfAmp    = _wave ->GetBinContent(ibin+1);
      wfAmp    /= wfMaximum;

      shapeBin = (_wave ->GetBinCenter(ibin+1) - tshift)/binWidth;
      if(shapeBin < 0) {
	shapeAmp = 0.;
      } else {
	shapeAmp = _refPulse->GetBinContent(shapeBin+1);
      }
     
      if(wfAmp <= 0.)                     continue;

      pcChi2 += pow( (wfAmp - shapeAmp), 2);
    }

    //normalize the chi2
    pcChi2 /= errChi2;
    
    //if (pcChi2 > 0.5) _wfWithPileUp = 1;
    _psd = pcChi2;
  }

  void        FitPolWaveformProcess::psdFromAmplitude(){
    double    amplitude = _wave->GetMaximum();
    double    charge    = _wave->Integral();
    
    _psd  =  amplitude/(charge*_ADCToMeV);
  }
  

//--------------------------------------------------------------------------------//
  void        FitPolWaveformProcess::findMaxima(){
  
    double   maxCont   [20] = {0};
    double   min   [20];
    
    for (int i=0; i<20; ++i){
      min[i] = 1e10;
    } 
  
    int      nMin(0);
    
    _nMax = 0;
    
    int    nBins = _wave->GetNbinsX();
  
    double   gradient, tmpGradient(1), content, reference(0);
  
    for (int i=0; i<nBins; ++i){
      content = _wave->GetBinContent(i+1);
      if (content <= 0 )          continue;
    
      gradient = content - reference;

      if (gradient > 0){
	if (gradient*tmpGradient > 0){
	  if (content > maxCont[_nMax]) {
	    maxCont    [_nMax]= content;
	    _max [_nMax]= i;
	  }
	  reference = maxCont[_nMax];
	}else{
	  reference = min[nMin];
	  ++_nMax;
	  if (content > maxCont[_nMax]) {
	    maxCont    [_nMax]= content;
	    _max [_nMax]= i;
	  }
	}
      }else{
	if (gradient*tmpGradient > 0){
	  if (content < min[nMin]){
	    min   [nMin] = content;
	  }
	  reference = min[nMin];
	}else{
	  reference = maxCont[_nMax];
	  ++nMin;
	  if (content < min[nMin]){
	    min   [nMin] = content;
	  }
	}
      }

      tmpGradient = gradient;

      //      printf("_nMax = %i nMin = %i grad = %4.3f min[1] = %3.2f max[1] = %3.2f\n", _nMax, nMin, gradient, min[1], max[1]);
    }
    
    _nMax +=1;
  }




//--------------------------------------------------------------------------------//
  void        FitPolWaveformProcess::separeteWaveform(){
    //search for the local maxima
    findMaxima();

    //set function parameters
    double   peak0 = _wave->GetBinCenter(_max[0]+1);
    double   amp0  = _wave->GetBinContent(_max[0]+1);
    double   peak1 = _wave->GetBinCenter(_max[1]+1);
    double   amp1  = _wave->GetBinContent(_max[1]+1);

    _fdlogn->SetParameters(-0.5, 22, peak0, peak1, amp0*_digiSampling, amp1*_digiSampling);
    _fdlogn->SetParLimits (0, -1, -0.2);
    _fdlogn->SetParLimits (1, 10, 100);
    _fdlogn->SetParLimits (2, peak0-10, peak0+10);
    _fdlogn->SetParLimits (3, peak1-10, peak1+10);
    _fdlogn->SetParLimits (4, amp0, amp0*1e3);
    _fdlogn->SetParLimits (5, amp1, amp1*1e3);


    //set range
    double   min    = _wave->GetBinCenter(_max[0]+1) - 50;//FindFirstBinAbove(10);
    double   max    = _wave->GetBinCenter(_max[1]+1) + 50;//_wave->FindLastBinAbove(10);

    _fdlogn->SetRange(min, max);
    
    //perform the fit
    _wave->Fit("fdlogn","RQ");
    _wave->Fit("fdlogn","RQ");
    _wave->Fit("fdlogn","RQ");

    //now calculate the reconstructed energies
    double   eta   = _fdlogn->GetParameter(0);
    double   sigma = _fdlogn->GetParameter(1);

    peak0 = _fdlogn->GetParameter(2);
    peak1 = _fdlogn->GetParameter(3);

    amp0  = _fdlogn->GetParameter(4);
    amp1  = _fdlogn->GetParameter(5);
    
    _flogn->SetParameters(eta, sigma, peak0, amp0);
    _eDep[0] =  _flogn->Integral(min, max)*_ADCToMeV;
    
    _flogn->SetParameters(eta, sigma, peak1, amp1);
    _eDep[1] = _flogn->Integral(min, max)*_ADCToMeV;
    //------------------------------------------------------------//
    

    //now extract the times
    calculateTime  (_wave->GetBinCenter(_max[0]+1), _time[0], _timeChi2[0]);
    
    _flogn->SetParameters(eta, sigma, peak1, amp1);

    _timeChi2[1] =  _fdlogn->GetChisquare()/_fdlogn->GetNDF();

    _flogn->SetParameters(eta, sigma, peak1, amp1);
    //now extract the time with the digital constant fraction
    //    _timeFraction =  0.05;
    _time    [1] = _flogn->GetX(_flogn->GetMaximum()*_timeFraction, 0, _flogn->GetMaximumX());
   
  }

//--------------------------------------------------------------------------------//

  void        FitPolWaveformProcess::calculateEnergy(double &Edep){
    Edep = _wave->Integral("width")*_ADCToMeV; 
  }

  void        FitPolWaveformProcess::calculateTime  (double WaveMax, double &Time, double &Chi2){
   
    //initilize function parameters
    double histPeak = WaveMax;//_wave->GetBinCenter(_wave->GetMaximumBin());
    //    _flogn->SetParameter(0, -1.);

    // _flogn->SetParLimits(0,  -10.,  0.);
    // //_flogn->FixParameter(0,  -0.21);
    // _flogn->SetParLimits(1,   1., 1e2);
    // //    _flogn->FixParameter(1,  11.7);
    // _flogn->SetParLimits(2, histPeak - 20.,  histPeak + 20.);
    // _flogn->SetParLimits(3,   0.,   _wave->Integral("width")*1000.);

    // _flogn->SetParameter(0, -0.21);
    // _flogn->SetParameter(1,  12.);
    // _flogn->SetParameter(2, _wave->GetBinCenter(_wave->GetMaximumBin()));
    // _flogn->SetParameter(3, _wave->Integral("width"));

    //define the fit range
    double           xmin(0),xmax(2000.);
    double           ymax = _wave->GetMaximum();
    //double binmin(0), binpsd(0);

 
    //double offset   = -1.;
    //    _timeFraction = 0.02;//0.05;
  
    //double           xmaxFrac(0.85);
    //    double           xmaxFrac(0.95);
    for (int ibin=_wave->GetMaximumBin();ibin>0;ibin--){
      //      if (_wave->GetBinContent(ibin)<xmaxFrac*ymax){
      if (_wave->GetBinContent(ibin+1)<=_fitThresholdMax*ymax){
	xmax = _wave->GetBinCenter(ibin+1);
	break;
      }
    }
    
    //search the minimum x of the fitting range
    //    double    fracmin(2.5e-04);//1e-02);//(5e-03);
    double    content;
    for (int ibin=_wave->GetMaximumBin();ibin>0;ibin--){
      content  = _wave->GetBinContent(ibin+1);
      //      if ( content <= fracmin*ymax ){
      if ( content <= _fitThresholdMin ){
	xmin = _wave->GetBinCenter(ibin+1);//(ibin);
	break;
      }
    }
    
    // if (_debugLevel > 0){
    //   int    contentMin =  _wave->GetBinContent(xmin/_wave->GetBinWidth(1)+1);
    //   printf("[FitPolWaveformProcess::calculateTime] time-fit range: xmin = %5.3f pulseHightMin = %i xmax = %.3f\n", 
    // 	     xmin, contentMin, xmax);
    // }

    //perform the fit
    _fWave->SetRange(xmin,xmax);
    _wave->Fit("fWave","RQ");
    _wave->Fit("fWave","RQ");
  
    double    nDof = _fWave->GetNDF();
    Chi2 = _fWave->GetChisquare()/nDof;

    
    //now extract the time with the digital constant fraction
    Time = _fWave->GetX(_wave->GetMaximum()*_timeFraction, histPeak-50., histPeak);//0, _wave->GetBinCenter(_wave->GetMaximumBin()));
    
    
  }
  
//--------------------------------------------------------------------------------//
  double    FitPolWaveformProcess::logn(double *x, double *par){
    double Epeak, sigma, eta, norm;
    double Aterm;
    double logterms0,s0;
    double logn,logterm;
    double expterm;
    double pigreco=3.14159265;
      
    eta = par[0];
    sigma = par[1];
    Epeak = par[2];
    norm = par[3];
      
    logterms0 = eta*2.35/2+sqrt(1+pow((eta*2.35/2),2));
    s0 = (2/2.35)*log(logterms0);
    
    Aterm = eta/(sqrt(2*pigreco)*sigma*s0);

    logterm = 1-(eta/sigma)*(x[0]-Epeak);  

    if(logterm<0){
      logterm = 0.0001;
    }
    expterm = log(logterm)/s0;
    expterm = -0.5*pow(expterm,2);

    logn = norm*Aterm *exp(expterm);      
    return logn;
  }

//--------------------------------------------------------------------------------//
  double    FitPolWaveformProcess::dlogn(double *x, double *par){
    double Epeak0, Epeak1, sigma, eta, norm0, norm1;
    double Aterm;
    double logterms0,s0;
    double logn,logterm0, logterm1;
    double expterm0, expterm1;
    double pigreco=3.14159265;
      
    eta    = par[0];
    sigma  = par[1];
    Epeak0 = par[2];
    Epeak1 = par[3];
    norm0  = par[4];
    norm1  = par[5];
      
    logterms0 = eta*2.35/2+sqrt(1+pow((eta*2.35/2),2));
    s0 = (2/2.35)*log(logterms0);
    
    Aterm = eta/(sqrt(2*pigreco)*sigma*s0);

    logterm0 = 1-(eta/sigma)*(x[0]-Epeak0);  
    logterm1 = 1-(eta/sigma)*(x[0]-Epeak1);  
      
    if(logterm0<0){
      logterm0 = 0.0001;
    }
    if(logterm1<0){
      logterm1 = 0.0001;
    }
    expterm0 = log(logterm0)/s0;
    expterm0 = -0.5*pow(expterm0,2);

    expterm1 = log(logterm1)/s0;
    expterm1 = -0.5*pow(expterm1,2);

    logn = Aterm *( norm0*exp(expterm0) + norm1*exp(expterm1));      

    return logn;

  }


//------------------------------------------------------------------------------------------//
  void        FitPolWaveformProcess::processWaveform(double   ADCToMeV, 
						  CaloDigi CaloHit , 
						  RecoCaloDigiCollection &RecoCaloHits) {
    _ADCToMeV = ADCToMeV;
    _nMax     = 1;
    //fill the histogram

    //time of the first digitized timestamp
    double         t0       = CaloHit.t0();

    //lenght of the waveform
    int            wfSize   = CaloHit.waveform().size();
    
    double         content, time, timeBin;
    double         binWidth = _wave->GetBinWidth(1);

    _wave->Reset();

    for (int i=0; i<wfSize; ++i){
      timeBin = i*_digiSampling;
      time    = t0 + timeBin;
      content = CaloHit.waveform().at(i);

      _wave->Fill(time, content);
      _wave->SetBinError(time/binWidth + 1 , _wave_point_error);
    }

    //--------------------------------------------------------------------------------//


    //understand how many signals are within the pulse
    psdFromAmplitude();
    //    psdFromChi2();



    if (_psd < _psdThreshold){//found pile up
      separeteWaveform();
      
      for(int i=0; i<_nMax; ++i){
	
	RecoCaloHits.push_back( RecoCaloDigi(CaloDigi(CaloHit),
					     _eDep    [i],
					     _wave->GetBinContent(_max[i]+1), 
					     _time    [i], 
					     _timeChi2[i],
					     _psd        )); 
      }
      
    }else{                   //single pulse
      
      double   eDep(0), waveMax, time, chi2;
      waveMax = _wave->GetMaximum();
      
      calculateEnergy(eDep);

      // if (_debugLevel > 0){
      // 	printf("[FitPolWaveformProcess::processWaveform] pulse charge = %5.3f\n", eDep);
      //      }
      calculateTime  (_wave->GetBinCenter(_wave->GetMaximumBin()+1), time, chi2);

      RecoCaloHits.push_back( RecoCaloDigi(CaloDigi(CaloHit),
					   eDep,
					   waveMax, 
					   time, 
					   chi2,
					   _psd    )); 
    }

  }
//------------------------------------------------------------------------------------------//

  void        FitPolWaveformProcess::book(){
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("CaloDigiDiag");
    
    _hist._hEdep         = tfdir.make<TH1F>("hEdep","Hit energy deposition",200,0.,500);
    _hist._hAmplitude    = tfdir.make<TH1F>("hAmplitude","Waveform amplitude", 2500,0.,2500);

    _hist._hTime         = tfdir.make<TH1F>("hTime","Hit time ", 12000, 0., 2000);

    _hist._hTimeVsE      = tfdir.make<TH2F>("hTimeVsE",
					   "Hit time resolution vs reco charge",
					   500 ,   0., 500,
					   1000, -50.,  50);
    _hist._hTimeVsAmp    = tfdir.make<TH2F>("hTimeVsAmplitude",
					      "Hit amplitude resolution vs reco charge",
					   1250 ,   0., 2500,
					   1000, -20.,  20);

    _hist._hAmpVsE       = tfdir.make<TH2F>("hAmpVsE","Amplitude versus charge; Amplitue [mV]; Energy [MeV]",
					   2500,0.,2500,
					   500 ,   0., 500);
    _hist._hPSDVsAmp     = tfdir.make<TH2F>("hPSDVsAmp","Pile-up disriminator versus amplitude",
					    500 ,   0.,  100,
					    2500,   0., 2500);
    

    _hist._hPSD          = tfdir.make<TH1F>("hPSD","PSD ", 600, 0., 300);
    _hist._hChi2Time     = tfdir.make<TH1F>("hChi2Time","#Chi^2 from time reco alg:#Chi2^2/ndof"   , 100, 0., 10);
    _hist._hNDofTime     = tfdir.make<TH1F>("hNDofTime","nDof from time reconstruction fit; nDof [#]", 200, 0., 10);
    _hist._hDtBinTime    = tfdir.make<TH1F>("hDtBinTime","#Delta t time reco; t_{reco} - t_{bin-edge} [ns]", 400, 0.,_digiSampling);
    
    //histograms for the time reconstruction fit results
    _hist._hFitEta       = tfdir.make<TH1F>("hFitEta", "Fit - #eta distribution; #eta", 100, -5, 5);
    _hist._hFitPeak      = tfdir.make<TH1F>("hFitPeak","Fit - peak distribution; peak [ns]",4000, 0.,  2000);
    _hist._hFitSigma     = tfdir.make<TH1F>("hFitSigma","Fit - #sigma distribution; #sigma [ns]", 100, 0, 100);
    _hist._hFitNorm      = tfdir.make<TH1F>("hFitNorm","Fit - Norm distribution; Norm [a.u.]", 5000, 0, 10000);



    _hist._hDt           = tfdir.make<TH1F>("hDt","#Deltat distribution; #Delta t = t_{reco} t_{MC} [ns]", 4000, -100., 100);

    for (int i=0; i<20; ++i){
      _hist._debugWf[i] = tfdir.make<TH1F>(Form("hDWf%i", i), Form("waveform %i",i),  _acquisitionLenght/_digiSampling, 0, _acquisitionLenght);
    }
    //create the  TTree      
    _tree  = tfdir.make<TTree>("Calo", "Calo");
       
    _tree->Branch("counter",      &_counter ,    "counter/I");

    _tree->Branch("time",         &_timeWf  ,        "time/F");
    _tree->Branch("mcTime",       &_mcTime ,      "mcTime/F");
    _tree->Branch("mcEDep",       &_mcEDep ,      "mcEDep/F");

    _tree->Branch("Chi2Time",     &_Chi2Time ,    "Chi2Time/F");
    _tree->Branch("psdWf   ",     &_psdWf      ,    "psdWf/F");

    _tree->Branch("nDof",         &_nDof ,        "nDof/F");
    _tree->Branch("fitEta",       &_fitEta ,      "fitEta/F");
    _tree->Branch("fitNorm",      &_fitNorm ,     "fitNorm/F");
    _tree->Branch("fitSigma",     &_fitSigma ,    "fitSigma/F");
    _tree->Branch("fitPeak",      &_fitPeak ,     "fitPeak/F");
    _tree->Branch("fitDlogn",     &_fitDlogn ,    "fitDlogn/F");

    _tree->Branch("nParticles",   &_nParticles,   "nParticles/I");
    _tree->Branch("pdgId"     ,   &_pdgId     ,   "pdgId[nParticles]/I");

    _tree->Branch("amp",          &_amp ,         "amp/F");
    _tree->Branch("nPeaks",       &_nPeaks ,      "nPeaks/I");
    _tree->Branch("charge",       &_charge ,      "chargexs/F");
  }

  //--------------------------------------------------------------------------------//


  void      FitPolWaveformProcess::fillDiag(CaloDigiMC*DigiMC, RecoCaloDigiCollection* RecoCaloHits){
    RecoCaloDigi *recoHit;
    int          size = RecoCaloHits->size();

    double       eDep, content(0), recoTimeBest(0), timeMCBest(0), eDepMax(0);
    
    for (int i=0; i<_nMax; ++i){
      recoHit = &RecoCaloHits->at(size - 1 -i);

      _counter     = _hitCounter;
      _psdWf       = _psd;
      _charge      = recoHit->edep();
      
      double     recoTime = recoHit->time();
      double     nDof     = _fWave->GetNDF();
      double     eta      = _fWave->GetParameter(0);
      double     sigma    = _fWave->GetParameter(1);
      double     peak     = _fWave->GetParameter(2);
      double     norm     = _fWave->GetParameter(3);

      _timeWf      = recoTime;
      _Chi2Time    = recoHit->tChi2();
      _nDof        = nDof;
      _amp         = recoHit->amplitude();
      _nPeaks      = _nMax;
      _fitEta      = eta   ;
      _fitSigma    = sigma ;
      _fitPeak     = peak  ;
      _fitNorm     = norm  ;
      if (_psd < _psdThreshold){//pile up
	_fitDlogn = _fdlogn->GetChisquare()/_fdlogn->GetNDF();
      }else{
	_fitDlogn = -9999;
      }

      double        timeMC     = DigiMC->timeFirst();//meanTime();
      int           nParticles(0);// = DigiMC->nParticles();

      _mcTime     = timeMC;
      _mcEDep     = DigiMC->totalEDep();


      for (int j=0; j<DigiMC->nParticles(); ++j){
	if (DigiMC->eDep(j) > 1.){
	_pdgId[j] = DigiMC->pdgId(j);
	++nParticles;
	}

      }
      _nParticles = nParticles;
      
      _tree->Fill();

      ++_hitCounter;


      //now fill the histograms
      _hist._hPSD       ->Fill(_psd);
      _hist._hPSDVsAmp  ->Fill(_psd, recoHit->amplitude());

      _hist._hChi2Time  ->Fill(recoHit->tChi2());
      _hist._hNDofTime  ->Fill(nDof  );
      _hist._hFitEta    ->Fill(eta   );
      _hist._hFitSigma  ->Fill(sigma );
      _hist._hFitPeak   ->Fill(peak  );
      _hist._hFitNorm   ->Fill(norm  );

      _hist._hEdep      ->Fill(recoHit->edep());
      _hist._hAmplitude ->Fill(recoHit->amplitude());

      _hist._hTime      ->Fill(recoTime);
      _hist._hTimeVsAmp ->Fill(recoTime, recoHit->amplitude());
      

      eDep = recoHit->edep();

      if (eDep > eDepMax){
	eDepMax      = eDep;
	recoTimeBest = recoTime;
	timeMCBest   = timeMC;
      }


      if (_debugHistIndex < 20){
	double      dt = (recoTime - timeMC);
	
	//	if (recoTime > 600.){
	if (dt > -5.3 && recoHit->edep() > 30.){
	  for (int i=0; i<_wave->GetNbinsX(); ++i){
	    content   = _wave->GetBinContent(i+1);
	    _hist._debugWf[_debugHistIndex] ->SetBinContent(i+1, content);
	    _hist._debugWf[_debugHistIndex] ->SetBinError  (i+1, _wave_point_error);
	  }
	  ++_debugHistIndex;
	}
	
	// if (eDep > 15.){
	//   for(int i=0; i< _wave->GetNbinsX(); ++i) {
	//       content   = _wave->GetBinContent(i+1);
	//       _hist._debugWf[_debugHistIndex] ->SetBinContent(i+1, content);
	//       _hist._debugWf[_debugHistIndex] ->SetBinError  (i+1, _wave_point_error);
	//   }
		 
	//   ++_debugHistIndex;
	// }
      }//end filling pulses
    
    }//end loop on the RecoCaloDigi
    
    
    // if (_debugLevel > 0){
    //   double      dt = (recoTimeBest - timeMCBest);
    //   if ( (dt < 23.4) && (eDepMax > 30)){
    // 	printf("[FitPolWaveformProcess::fillDiag] MCTime = %5.3f MCEDep = %5.3f dt = %5.3f\n", 
    // 	       DigiMC->meanTime(), DigiMC->totalEDep(), (recoTimeBest - timeMCBest));
    //   }
    // }
 
    
    _hist._hDtBinTime ->Fill( (recoTimeBest/_digiSampling - int(recoTimeBest/_digiSampling))*_digiSampling);
    _hist._hDt        ->Fill(recoTimeBest - timeMCBest);
    _hist._hTimeVsE   ->Fill(eDepMax, recoTimeBest- timeMCBest);


  }
  //------------------------------------------------------------------------------------------//

}
