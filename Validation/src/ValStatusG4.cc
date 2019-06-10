
#include "Validation/inc/ValStatusG4.hh"


int mu2e::ValStatusG4::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hstat = tfs.make<TH1D>( "stat", "Status", 50, -0.5, 50.5);
  _hnTrk = tfs.make<TH1D>( "Ntrk", "N gtrk", 101, -0.5, 100.5);
  _hnTrk2 = tfs.make<TH1D>( "Ntrk2", "N gtrk", 100, 0.0, 1000.0);
  _hnTrk3 = tfs.make<TH1D>( "Ntrk3", "log10(N gtrk)", 100, 0.0, 7.0);
  _hover = tfs.make<TH1D>( "over", "SimPart overflow", 2, -0.5, 1.5);
  _hkill = tfs.make<TH1D>( "kill", "N trk killed", 11, -0.5, 10.5);
  _hkillfp = tfs.make<TH1D>( "killfp", "N trk killed", 11, -0.5, 10.5);
  _hCPU1 = tfs.make<TH1D>( "CPU1", "CPU [info]", 100, 0.0, 1.0);
  _hCPU2 = tfs.make<TH1D>( "CPU2", "CPU [info]", 100, 0.0, 100.0);
  _hCPU3 = tfs.make<TH1D>( "CPU3", "log10(CPU) [info]", 100, -3.0, 4.0);
  _hWall1 = tfs.make<TH1D>( "Wall1", "Wall [info]", 100, 0.0, 1.0);
  _hWall2 = tfs.make<TH1D>( "Wall2", "Wall [info]", 100, 0.0, 100.0);
  _hWall3 = tfs.make<TH1D>( "Wall3", "log10(Wall) [info]", 100, -3.0, 4.0);

  return 0;
}

int mu2e::ValStatusG4::fill(const mu2e::StatusG4 & obj,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(1.0);

  double logx;
  _hstat->Fill(obj.status());
  _hnTrk->Fill(obj.nG4Tracks());
  _hnTrk2->Fill(obj.nG4Tracks());
  logx = (obj.nG4Tracks()>0 ? log10(obj.nG4Tracks()) : -1.0);
  _hnTrk3->Fill(logx);
  _hover->Fill(obj.overflowSimParticles());
  _hkill->Fill(obj.nKilledStepLimit());
  _hkillfp->Fill(obj.nKilledByFieldPropagator());
  _hCPU1->Fill(obj.cpuTime());
  _hCPU2->Fill(obj.cpuTime());
  logx = (obj.cpuTime()>0 ? log10(obj.cpuTime()) : -20.0);
  _hCPU3->Fill(logx);
  _hWall1->Fill(obj.realTime());
  _hWall2->Fill(obj.realTime());
  logx = (obj.realTime()>0 ? log10(obj.realTime()) : -20.0);
  _hWall3->Fill(logx);

  return 0;
}
