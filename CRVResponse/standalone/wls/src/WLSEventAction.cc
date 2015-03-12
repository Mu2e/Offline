#include "WLSEventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "WLSSteppingAction.hh"
#include "WLSDetectorConstruction.hh"

#include "Randomize.hh"
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TFile.h>

#include <TStyle.h>
#include <TText.h>
#include <TGraph.h>
#include <TMarker.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TPaveStats.h>

#include "MakeCrvSiPMResponses.hh"
#include "MakeCrvWaveforms.hh"
#include "MakeCrvRecoPulses.hh"

#include <stdexcept>

#include "CLHEP/Random/Randomize.h"

WLSEventAction* WLSEventAction::_fgInstance = NULL;

WLSEventAction::WLSEventAction(int mode, int id) : _mode(mode), _storeConstants(false)
{
  _fgInstance = this;
  _fileLookupTable = NULL;

  if(_mode==-1)
  {
    std::ostringstream filename;
    filename<<"CRVLookupTable_"<<id<<".root";
    _fileLookupTable = new TFile(filename.str().c_str(),"RECREATE");
    if(_fileLookupTable==NULL) throw std::logic_error("CRVLookupTable.root could not be created.");

    WLSDetectorConstruction *detector = WLSDetectorConstruction::Instance();
    std::vector<double> xbinsVector = detector->GetXBins();
    std::vector<double> ybinsVector = detector->GetYBins();
    std::vector<double> zbinsVector = detector->GetZBins();
    double *xbins = new double[xbinsVector.size()];
    double *ybins = new double[ybinsVector.size()];
    double *zbins = new double[zbinsVector.size()];
    for(unsigned int v=0; v<xbinsVector.size(); v++) xbins[v]=xbinsVector[v];
    for(unsigned int v=0; v<ybinsVector.size(); v++) ybins[v]=ybinsVector[v];
    for(unsigned int v=0; v<zbinsVector.size(); v++) zbins[v]=zbinsVector[v];

    double *tbins = new double[2001];
    double *ebins = new double[101];
    for(unsigned int v=0; v<2001; v++) tbins[v]=0.1*v;
    for(unsigned int v=0; v<101; v++) ebins[v]=v;

    for(int table=0; table<4; table++)  //0...scintillation, 1...cerenkov, 
    {                                   //2...cerenkovFiber0, 3...cerenkovFiber1
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        std::string tableTag[4]={"scintillation","cerenkov","cerenkovFiber0","cerenkovFiber1"};
        std::stringstream s1, s2, s3;
        s1<<tableTag[table]<<"_SurvivalProb_SiPM_"<<SiPM;
        s2<<tableTag[table]<<"_TimeDifference_SiPM_"<<SiPM;
        s3<<tableTag[table]<<"_FiberEmissions_SiPM_"<<SiPM;

        _histSurvivalProb[table][SiPM] = new TH3D(s1.str().c_str(),s1.str().c_str(),
                        xbinsVector.size()-1, xbins,
                        ybinsVector.size()-1, ybins,
                        zbinsVector.size()-1, zbins);

        _histTimeDifference[table][SiPM] = new TH3D(s2.str().c_str(),s2.str().c_str(),
                        ybinsVector.size()-1, ybins,
                        zbinsVector.size()-1, zbins,
                        2000, tbins);

        _histFiberEmissions[table][SiPM] = new TH3D(s3.str().c_str(),s3.str().c_str(),
                        ybinsVector.size()-1, ybins,
                        zbinsVector.size()-1, zbins,
                        100, ebins);
      }
    }
    delete[] xbins;
    delete[] ybins;
    delete[] zbins;
    delete[] tbins;
    delete[] ebins;
  }

  if(_mode==0)
  {
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      std::stringstream s0, s1, title;
      s0<<"Photons_Mode_"<<0<<"__SiPM_"<<SiPM;
      s1<<"Photons_Mode_"<<1<<"__SiPM_"<<SiPM;
      title<<"Fiber: "<<SiPM/2<<",  Side: "<<SiPM%2;
      _histP[0][SiPM] = new TH1D(s0.str().c_str(),title.str().c_str(),700,0,700);
      _histP[1][SiPM] = new TH1D(s1.str().c_str(),title.str().c_str(),700,0,700);
      _histP[0][SiPM]->GetXaxis()->SetTitle("Photons");
      _histP[0][SiPM]->SetLineColor(1);
      _histP[1][SiPM]->GetXaxis()->SetTitle("Photons");
      _histP[1][SiPM]->SetLineColor(2);
    }

    for(int SiPM=0; SiPM<4; SiPM++)
    {
      std::stringstream s0, s1, title;
      s0<<"ArrivalTimes_Mode_"<<0<<"__SiPM_"<<SiPM;
      s1<<"ArrivalTimes_Mode_"<<1<<"__SiPM_"<<SiPM;
      title<<"Fiber: "<<SiPM/2<<",  Side: "<<SiPM%2;
      _histT[0][SiPM] = new TH1D(s0.str().c_str(),title.str().c_str(),250,0,250);
      _histT[1][SiPM] = new TH1D(s1.str().c_str(),title.str().c_str(),250,0,250);
      _histT[0][SiPM]->GetXaxis()->SetTitle("t [ns]");
      _histT[0][SiPM]->SetLineColor(1);
      _histT[1][SiPM]->GetXaxis()->SetTitle("t [ns]");
      _histT[1][SiPM]->SetLineColor(2);
    }

    _photonsVsIntegral = new TH2D("PhotonsVsIntegral","PhotonsVsIntegral", 200,0,2.0, 100,0,100);
    _photonsVsIntegral->SetYTitle("Photons");
    _photonsVsIntegral->SetXTitle("Integral");
    _photonsVsPulseHeight = new TH2D("PhotonsVsPulseHeight","PhotonVsPulseHeight", 200,0,0.5, 100,0,100);
    _photonsVsPulseHeight->SetYTitle("Photons");
    _photonsVsPulseHeight->SetXTitle("PulseHeight [mV]");
    _PEsVsIntegral = new TH2D("PEsVsIntegral","PEsVsIntegral", 200,0,2.0, 100,0,100);
    _PEsVsIntegral->SetYTitle("PEs");
    _PEsVsIntegral->SetXTitle("Integral");
    _PEsVsPulseHeight = new TH2D("PEsVsPulseHeight","PEsVsPulseHeight", 200,0,0.5, 100,0,100);
    _PEsVsPulseHeight->SetYTitle("PEs");
    _PEsVsPulseHeight->SetXTitle("PulseHeight [mV]");
  }
}

WLSEventAction::~WLSEventAction()
{
  if(_mode==-1)
  {
    for(int table=0; table<4; table++)
    {
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        delete _histSurvivalProb[table][SiPM];
        delete _histTimeDifference[table][SiPM];
        delete _histFiberEmissions[table][SiPM];
      }
    }
    _fileLookupTable->Close();
    delete _fileLookupTable;
  }

  if(_mode==0)
  {
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      delete _histP[0][SiPM];
      delete _histP[1][SiPM];
      delete _histT[0][SiPM];
      delete _histT[1][SiPM];
    }
  
    delete _photonsVsIntegral;
    delete _photonsVsPulseHeight;
  }

}

G4ThreeVector WLSEventAction::GetHistBinCenter(int binx, int biny, int binz) 
{
  double x=_histSurvivalProb[0][0]->GetXaxis()->GetBinCenter(binx);
  double y=_histSurvivalProb[0][0]->GetYaxis()->GetBinCenter(biny);
  double z=_histSurvivalProb[0][0]->GetZaxis()->GetBinCenter(binz);
  G4ThreeVector point(x,y,z);
  return point;
}

double WLSEventAction::GetHistBinWidthX(int binx) {return _histSurvivalProb[0][0]->GetXaxis()->GetBinWidth(binx);} 
double WLSEventAction::GetHistBinWidthY(int biny) {return _histSurvivalProb[0][0]->GetYaxis()->GetBinWidth(biny);} 
double WLSEventAction::GetHistBinWidthZ(int binz) {return _histSurvivalProb[0][0]->GetZaxis()->GetBinWidth(binz);} 

void WLSEventAction::BeginOfEventAction(const G4Event* evt)
{
  std::cout<<"Event # "<<evt->GetEventID()<<std::endl;

  WLSSteppingAction::Instance()->Reset();
}

void WLSEventAction::EndOfEventAction(const G4Event* evt)
{
  if(_mode==-1)
  {
    int table = evt->GetEventID();

    G4Material *material = G4Material::GetMaterial("PMMA",true); //fiber
    G4MaterialPropertiesTable* materialPropertiesTable = material->GetMaterialPropertiesTable();
    double refractiveIndexFiber = (*materialPropertiesTable->GetProperty("RINDEX"))[0];
    double speedOfLightFiber = CLHEP::c_light/refractiveIndexFiber;

    WLSDetectorConstruction *detector = WLSDetectorConstruction::Instance();

    for(int SiPM=0; SiPM<4; SiPM++)
    {
      double p=(double)WLSSteppingAction::Instance()->GetArrivalTimes(0,SiPM).size()/(double)_generatedPhotons;
      _histSurvivalProb[table][SiPM]->Fill(_start.x(),_start.y(),_start.z(), p);

      double zSiPM=(SiPM%2==0?-detector->GetScintillatorHalfLength():detector->GetScintillatorHalfLength());
      double straightLineTravelTime=abs(_start.z()-zSiPM)/speedOfLightFiber;
      const std::vector<double> &arrivalTimes = WLSSteppingAction::Instance()->GetArrivalTimes(0,SiPM);
      for(size_t i=0; i<arrivalTimes.size(); i++)
      {
        double t=arrivalTimes[i]-straightLineTravelTime;  //the decay times have been set to 0 for mode==-1
        if(t<0) t=0;
        _histTimeDifference[table][SiPM]->Fill(_start.y(),_start.z(),t);
      }

      const std::vector<int> &fiberEmissions = WLSSteppingAction::Instance()->GetFiberEmissions(SiPM);
      for(size_t i=0; i<fiberEmissions.size(); i++)
      {
        int n=fiberEmissions[i];
        _histFiberEmissions[table][SiPM]->Fill(_start.y(),_start.z(),n);
      }

      _histSurvivalProb[table][SiPM]->Write();
      _histTimeDifference[table][SiPM]->Write();
      _histFiberEmissions[table][SiPM]->Write();
    }

  //write some constants to the file
    if(table==0 && _storeConstants)
    {
      G4Material *material = G4Material::GetMaterial("Polystyrene",true); //scintillator
      G4MaterialPropertiesTable *materialPropertiesTable = material->GetMaterialPropertiesTable();
      G4MaterialPropertyVector *rindex = materialPropertiesTable->GetProperty("RINDEX");
      double cerenkovRindex = (*rindex)[0];
      double cerenkovEnergyMin = rindex->GetMinLowEdgeEnergy();
      double cerenkovEnergyMax = rindex->GetMaxLowEdgeEnergy();
      double cerenkovEnergyInterval = cerenkovEnergyMax - cerenkovEnergyMin;
      double ratioFastSlow = materialPropertiesTable->GetConstProperty("YIELDRATIO");
 
      double scintillatorDensity = material->GetDensity();
      double scintillatorBirksConstant = material->GetIonisation()->GetBirksConstant();

      TNtuple *ntuple = new TNtuple("CRVbarConstants","CRVbarConstants","halfThickness:halfWidth:halfLength:speedOfLightFiber:cerenkovRindex:cerenkovEinterval:ratioFastSlow:scintillatorDensity:scintillatorBirksConstant:fiberSeparation:holeRadius:fiberRadius");

      ntuple->Fill(detector->GetScintillatorHalfThickness(),
                   detector->GetScintillatorHalfWidth(), 
                   detector->GetScintillatorHalfLength(),
                   speedOfLightFiber,
                   cerenkovRindex,
                   cerenkovEnergyInterval,
                   ratioFastSlow,
                   scintillatorDensity,
                   scintillatorBirksConstant,
                   detector->GetFiberSeparation(),
                   detector->GetHoleRadius(),
                   detector->GetClad2Radius());
      ntuple->Write();
      delete ntuple;
    }
  }

  if(_mode==0)
  {
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      _histP[0][SiPM]->Fill(WLSSteppingAction::Instance()->GetArrivalTimes(0,SiPM).size());
      _histP[1][SiPM]->Fill(WLSSteppingAction::Instance()->GetArrivalTimes(1,SiPM).size());

      const std::vector<double> &arrivalTimes0 = WLSSteppingAction::Instance()->GetArrivalTimes(0,SiPM);
      const std::vector<double> &arrivalTimes1 = WLSSteppingAction::Instance()->GetArrivalTimes(1,SiPM);
      for(size_t i=0; i<arrivalTimes0.size(); i++) _histT[0][SiPM]->Fill(arrivalTimes0[i]);
      for(size_t i=0; i<arrivalTimes1.size(); i++) _histT[1][SiPM]->Fill(arrivalTimes1[i]);
    }

    for(int SiPM=0; SiPM<4; SiPM++) std::cout<<_histP[0][SiPM]->GetMean()<<"/"<<_histP[1][SiPM]->GetMean()<<"  ";
    std::cout<<std::endl;
    for(int SiPM=0; SiPM<4; SiPM++) std::cout<<_histT[0][SiPM]->GetMean()<<"/"<<_histT[1][SiPM]->GetMean()<<"  ";
    std::cout<<std::endl;

    Draw(evt);
  }
}

void WLSEventAction::Draw(const G4Event* evt) const
{
  MakeCrvSiPMResponses::ProbabilitiesStruct probabilities;
  probabilities._constGeigerProbCoef = 2;
  probabilities._constGeigerProbVoltScale = 3;
  probabilities._constTrapType0Prob = 0.14;   //trap_prob*trap_type0_prob=0.2*0.7
  probabilities._constTrapType1Prob = 0.06;   //trap_prob*trap_type1_pron=0.2*0.3
  probabilities._constTrapType0Lifetime = 5;
  probabilities._constTrapType1Lifetime = 50;
  probabilities._constThermalProb = 6.25e-7; //1MHz at SiPM --> 1e-3/(#pixel*t[ns])  //exp(-E_th/T)=1.6e-6
  probabilities._constPhotonProduction = 0.1; //0.4;

  CLHEP::HepJamesRandom engine(1);
  CLHEP::RandFlat randFlat(engine);
  CLHEP::RandPoissonQ randPoissonQ(engine);
  MakeCrvSiPMResponses sim(randFlat,randPoissonQ);
  sim.SetSiPMConstants(1600, 2.5, 0, 1695, 0.08, probabilities);

  MakeCrvWaveforms makeCrvWaveform, makeCrvWaveform2;
  double binWidth = 12.5; //ns
  double binWidth2 = 1.0; //ns
  gStyle->SetOptStat(0);
  makeCrvWaveform.LoadSinglePEWaveform("singlePEWaveform.txt", 1.0, 200);
  makeCrvWaveform2.LoadSinglePEWaveform("singlePEWaveform.txt", 1.0, 200);

  double startTime=-G4UniformRand()*binWidth;
  std::vector<double> siPMtimes[4], siPMcharges[4];
  std::vector<double> waveform[4], waveform2[4];
  for(int SiPM=0; SiPM<4; SiPM++)
  {
    const std::vector<double> &photonTimes = WLSSteppingAction::Instance()->GetArrivalTimes(1,SiPM);

    std::vector<SiPMresponse> SiPMresponseVector;
    sim.Simulate(photonTimes, SiPMresponseVector);
    for(unsigned int i=0; i<SiPMresponseVector.size(); i++)
    {
      siPMtimes[SiPM].push_back(SiPMresponseVector[i]._time);
      siPMcharges[SiPM].push_back(SiPMresponseVector[i]._charge);
    }

    makeCrvWaveform.MakeWaveform(siPMtimes[SiPM], siPMcharges[SiPM], waveform[SiPM], startTime, binWidth);
    makeCrvWaveform2.MakeWaveform(siPMtimes[SiPM], siPMcharges[SiPM], waveform2[SiPM], startTime, binWidth2);
  }

  std::ostringstream s1;
  s1<<"waveform_"<<evt->GetEventID();

  gStyle->SetOptStat(0);
  TCanvas c(s1.str().c_str(),s1.str().c_str(),1000,1000);
  c.Divide(2,2);
  TGraph *graph[4]={NULL};
  TGraph *graph2[4]={NULL};
  TH1D *hist[4], *histSiPMResponse[4];
  std::vector<TGraph*> graphVector;
  std::vector<TMarker*> markerVector;
  std::vector<TGaxis*> axisVector;
  std::vector<TLine*> lineVector;
  std::vector<TText*> textVector;

  for(int SiPM=0; SiPM<4; SiPM++)
  {
    c.cd(SiPM+1);

//Photon Arrival Times
    std::ostringstream s2, s3;
    s2<<"Photons_"<<evt->GetEventID()<<"__"<<SiPM;
    s3<<"Fiber: "<<SiPM/2<<",  Side: "<<SiPM%2;
    hist[SiPM]=new TH1D(s2.str().c_str(),s3.str().c_str(),100,0,200);
    const std::vector<double> &photonTimes = WLSSteppingAction::Instance()->GetArrivalTimes(1,SiPM);
    for(unsigned int i=0; i<photonTimes.size(); i++)
    {
      hist[SiPM]->Fill(photonTimes[i]);
    }

    hist[SiPM]->SetLineColor(kBlue);
    hist[SiPM]->GetXaxis()->SetTitle("t [ns]");
    hist[SiPM]->GetYaxis()->SetTitle("Photons");
    hist[SiPM]->GetYaxis()->SetTitleOffset(0.5);
    hist[SiPM]->GetYaxis()->SetAxisColor(kBlue);
    hist[SiPM]->GetYaxis()->SetTitleColor(kBlue);
    hist[SiPM]->GetYaxis()->SetLabelColor(kBlue);
    hist[SiPM]->Draw();

//SiPM response
    double scaleSiPMResponse = 0.25;
    histSiPMResponse[SiPM]=new TH1D((s2.str()+"SiPMResponse").c_str(),"",100,0,200);
    for(unsigned int j=0; j<siPMtimes[SiPM].size(); j++)
    {
      histSiPMResponse[SiPM]->Fill(siPMtimes[SiPM][j], siPMcharges[SiPM][j]*scaleSiPMResponse);
    }
    histSiPMResponse[SiPM]->SetLineColor(kOrange);
    histSiPMResponse[SiPM]->Draw("same");

//waveforms with 1 ns bin width
    unsigned int n2 = waveform2[SiPM].size();
    if(n2==0) continue;
    double *t2 = new double[n2];
    double *v2 = new double[n2];
    double histMax = hist[SiPM]->GetMaximum();
    double waveformMax = *std::max_element(waveform2[SiPM].begin(),waveform2[SiPM].end());
    double scale = histMax/waveformMax;
    for(unsigned int j=0; j<n2; j++)
    {
      t2[j]=startTime+j*binWidth2;
      v2[j]=waveform2[SiPM][j];
      v2[j]*=scale;
    }
    graph2[SiPM]=new TGraph();
    graph2[SiPM]->SetTitle("");
    graph2[SiPM]->SetLineWidth(1);
    graph2[SiPM]->SetLineColor(kRed);
    graph2[SiPM]->DrawGraph(n2,t2,v2,"same");

    delete[] t2;
    delete[] v2;

//waveforms with 12.5 ns bin width
    unsigned int n = waveform[SiPM].size();
    if(n==0) continue;
    double *t = new double[n];
    double *v = new double[n];
    for(unsigned int j=0; j<n; j++)
    {
      t[j]=startTime+j*binWidth;
      v[j]=waveform[SiPM][j];
      v[j]*=scale;
    }
    graph[SiPM]=new TGraph(n,t,v);
    graph[SiPM]->SetTitle("");
    graph[SiPM]->SetMarkerStyle(24);
    graph[SiPM]->SetMarkerSize(1.5);
    graph[SiPM]->SetMarkerColor(kRed);
    graph[SiPM]->Draw("sameP");

    delete[] t;
    delete[] v;

//waveforms with 12.5 ns bin width for points above the threshold
//used for the integral
    for(unsigned int j=0; j<n; j++)
    {
      double tI=startTime+j*binWidth;
      double vI=waveform[SiPM][j];
      if(tI>200) break;
      if(vI<=0.015) continue;
      TMarker *marker = new TMarker(tI, vI*scale, markerVector.size());
      markerVector.push_back(marker);
      marker->SetMarkerStyle(20);
      marker->SetMarkerSize(1.5);
      marker->SetMarkerColor(kRed);
      marker->Draw("same");
    }

//Landau fit
    MakeCrvRecoPulses makeRecoPulses(0.015,0.2, 3.2,33.6);
    makeRecoPulses.SetWaveform(waveform[SiPM], startTime, binWidth);
    unsigned int nPulse = makeRecoPulses.GetNPulses();
    for(unsigned int pulse=0; pulse<nPulse; pulse++)
    {
      double tL1=makeRecoPulses.GetT1(pulse);
      double tL2=makeRecoPulses.GetT2(pulse);
      int nL=(tL2-tL1)/1.0 + 1;
      double *tL = new double[nL];
      double *vL = new double[nL];
      for(int iL=0; iL<nL; iL++)
      {
        double p0 = makeRecoPulses.GetLandauParam0(pulse);
        double p1 = makeRecoPulses.GetLandauParam1(pulse);
        double p2 = makeRecoPulses.GetLandauParam2(pulse);
        tL[iL] = tL1 + iL*1.0;
        vL[iL] = p0*TMath::Landau(tL[iL], p1, p2);
        vL[iL]*=scale;
        if(isnan(vL[iL])) nL=0;
      }
      if(nL>0)
      {
        TGraph *graphL=new TGraph();
        graphVector.push_back(graphL);
        graphL->SetTitle("");
        graphL->SetLineWidth(2);
        graphL->SetLineColor(kGreen);
        graphL->DrawGraph(nL,tL,vL,"same");
      }

      double leadingEdge=makeRecoPulses.GetLeadingEdge(pulse);
      if(!isnan(leadingEdge) && leadingEdge<200.0)
      {
        TMarker *marker = new TMarker(leadingEdge,
                                      0.2*makeRecoPulses.GetPulseHeight(pulse)*scale,
                                      markerVector.size());
        markerVector.push_back(marker);
        marker->SetMarkerStyle(21);
        marker->SetMarkerSize(1.5);
        marker->SetMarkerColor(kGreen);
        marker->Draw("same");
      }
      delete[] tL;
      delete[] vL;
    }

    if(makeRecoPulses.GetNPulses()==1)
    {
      int photons = photonTimes.size();
      std::cout<<"Photons/integral: "<<photons/makeRecoPulses.GetIntegral(0)<<"         ";
      std::cout<<"Photons/maxBin: "<<photons/makeRecoPulses.GetPulseHeight(0)<<std::endl;
      _photonsVsIntegral->Fill(makeRecoPulses.GetIntegral(0),photons);
      _photonsVsPulseHeight->Fill(makeRecoPulses.GetPulseHeight(0),photons);
      double PEs=0;
      for(unsigned int j=0; j<siPMtimes[SiPM].size(); j++) PEs+=siPMcharges[SiPM][j];
      std::cout<<"PEs/integral: "<<PEs/makeRecoPulses.GetIntegral(0)<<"         ";
      std::cout<<"PEs/maxBin: "<<PEs/makeRecoPulses.GetPulseHeight(0)<<std::endl;
      _PEsVsIntegral->Fill(makeRecoPulses.GetIntegral(0),PEs);
      _PEsVsPulseHeight->Fill(makeRecoPulses.GetPulseHeight(0),PEs);
    }

    TGaxis *axis = new TGaxis(170.0,0,170.0,histMax,0,histMax/scale,10,"+L");
    axisVector.push_back(axis);
    axis->SetTitle("voltage [V]");
    axis->SetTitleOffset(-0.5);
    axis->SetTitleColor(kRed);
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->Draw("same");

    TGaxis *axisSiPMResponse = new TGaxis(200.0,0,200.0,histMax,0,histMax/scaleSiPMResponse,10,"+L");
    axisVector.push_back(axisSiPMResponse);
    axisSiPMResponse->SetTitle("SiPM output [PE]");
    axisSiPMResponse->SetTitleOffset(-0.5);
    axisSiPMResponse->SetTitleColor(kOrange);
    axisSiPMResponse->SetLineColor(kOrange);
    axisSiPMResponse->SetLabelColor(kOrange);
    axisSiPMResponse->Draw("same");

    TLine *landauLine = new TLine(100, histMax*0.9, 130, histMax*0.9);
    lineVector.push_back(landauLine);
    landauLine->SetLineWidth(2);
    landauLine->SetLineColor(kGreen);
    landauLine->Draw("same");
    TMarker *landauMarker = new TMarker(110, histMax*0.9, markerVector.size());
    markerVector.push_back(landauMarker);
    landauMarker->SetMarkerStyle(21);
    landauMarker->SetMarkerSize(1.5);
    landauMarker->SetMarkerColor(kGreen);
    landauMarker->Draw("same");
    TText *landauText1 = new TText(100, histMax*0.82, "Landau fit with");
    TText *landauText2 = new TText(100, histMax*0.76, "leading edge");
    textVector.push_back(landauText1);
    textVector.push_back(landauText2);
    landauText1->SetTextSize(0.03);
    landauText2->SetTextSize(0.03);
    landauText1->SetTextColor(kGreen);
    landauText2->SetTextColor(kGreen);
    landauText1->Draw("same");
    landauText2->Draw("same");
  }
  if(evt->GetEventID()<20) c.SaveAs((s1.str()+".C").c_str());
  _photonsVsIntegral->SaveAs("PhotonsVsIntegral.C");
  _photonsVsPulseHeight->SaveAs("PhotonsVsPulseHeight.C");
  _PEsVsIntegral->SaveAs("PEsVsIntegral.C");
  _PEsVsPulseHeight->SaveAs("PEsVsPulseHeight.C");


  gStyle->SetOptStat(1111);
  TCanvas c1("Photons","Photons",1000,1000);
  c1.Divide(2,2);
  for(int SiPM=0; SiPM<4; SiPM++)
  {
    c1.cd(SiPM+1);
    gPad->SetLogy();
    _histP[0][SiPM]->Draw();
    gPad->Update();
    TPaveStats *stats0 = (TPaveStats*)_histP[0][SiPM]->FindObject("stats");
    stats0->SetTextColor(1);
    stats0->SetLineColor(1);
    double X1 = stats0->GetX1NDC();
    double Y1 = stats0->GetY1NDC();
    double X2 = stats0->GetX2NDC();
    double Y2 = stats0->GetY2NDC();
    _histP[1][SiPM]->Draw();
    gPad->Update();
    TPaveStats *stats1 = (TPaveStats*)_histP[1][SiPM]->FindObject("stats");
    stats1->SetTextColor(2);
    stats1->SetLineColor(2);
    stats1->SetX1NDC(X1);
    stats1->SetY1NDC(Y1-(Y2-Y1));
    stats1->SetX2NDC(X2);
    stats1->SetY2NDC(Y1);
    _histP[0][SiPM]->Draw();
    _histP[1][SiPM]->Draw("same");
  }      
  c1.SaveAs("Photons.C");

  TCanvas c2("ArrivalTimes","ArrivalTimes",1000,1000);
  c2.Divide(2,2);
  for(int SiPM=0; SiPM<4; SiPM++)
  {
    c2.cd(SiPM+1);
    gPad->SetLogy();
    _histT[0][SiPM]->Draw();
    gPad->Update();
    TPaveStats *stats0 = (TPaveStats*)_histT[0][SiPM]->FindObject("stats");
    stats0->SetTextColor(1);
    stats0->SetLineColor(1);
    double X1 = stats0->GetX1NDC();
    double Y1 = stats0->GetY1NDC();
    double X2 = stats0->GetX2NDC();
    double Y2 = stats0->GetY2NDC();
    _histT[1][SiPM]->Draw();
    gPad->Update();
    TPaveStats *stats1 = (TPaveStats*)_histT[1][SiPM]->FindObject("stats");
    stats1->SetTextColor(2);
    stats1->SetLineColor(2);
    stats1->SetX1NDC(X1);
    stats1->SetY1NDC(Y1-(Y2-Y1));
    stats1->SetX2NDC(X2);
    stats1->SetY2NDC(Y1);
    _histT[0][SiPM]->Draw();
    _histT[1][SiPM]->Draw("same");
  }      
  c2.SaveAs("ArrivalTimes.C");

  for(int SiPM=0; SiPM<4; SiPM++)
  {
    delete hist[SiPM];
    delete histSiPMResponse[SiPM];
    if(graph[SiPM]) delete graph[SiPM];
    if(graph2[SiPM]) delete graph2[SiPM];
  }
  for(size_t i=0; i<graphVector.size(); i++) delete graphVector[i];
  for(size_t i=0; i<markerVector.size(); i++) delete markerVector[i];
  for(size_t i=0; i<axisVector.size(); i++) delete axisVector[i];
  for(size_t i=0; i<textVector.size(); i++) delete textVector[i];
  for(size_t i=0; i<lineVector.size(); i++) delete lineVector[i];
}

