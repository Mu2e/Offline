#include "TrkDiag/test/FillChain.C+"
#include "TrkDiag/test/TrkRecoEff.C"
#include "TrkDiag/test/TrkRecoTrig.C"
void TRD() {
  TChain* mtce(0);
  FillChain(mtce,"cemix.txt","TrkRecoDiag/trdiag");
  TChain* mtbk(0);
  FillChain(mtbk,"bkgonly.txt","TrkRecoDiag/trdiag");
  TrkRecoEff trdce(mtce,1e5);
  trdce.Loop();
  trdce.drawHistos();
  TrkRecoTrig trdbk(mtbk);
  trdbk.Loop();
  trdbk.drawHistos();
}
