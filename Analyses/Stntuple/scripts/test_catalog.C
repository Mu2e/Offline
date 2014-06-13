///////////////////////////////////////////////////////////////////////////////
//
// root [1] .x Stntuple/scripts/test_catalog.C("stntuple/dev_240:etau08-taumet",151900,152700)
// root [2] chain->GetListOfFiles()->Print()
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0020.PR4779.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0021.PR4792.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0022.PR4840.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0023.PR5518.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0060.PR7480.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0061.PR7492.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0062.PR7501.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0063.PR7591.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0064.PR7665.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0065.PR7670.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0066.PR7675.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0067.PR7681.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0069.PR7718.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0124.PS6918.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0125.PS6933.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0126.PS7084.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0127.PS7086.0.s.01
// OBJ: TChainElement      STNTUPLE        root://fcdfdata030.fnal.gov//cdf/scratch/ewk/etau08-taumet/etau08-taumet.0128.PS7087.0.s.01
//
///////////////////////////////////////////////////////////////////////////////
TStnCatalog*         catalog = NULL;
TChain*              chain   = NULL;

int test_catalog(const char* Dataset, Int_t Run1=0, Int_t Run2=1000000) {
  catalog = new TStnCatalog();
  chain   = new TChain("STNTUPLE");
  catalog->InitChain(chain,Dataset,Run1,Run2);
}
