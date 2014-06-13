int copy() {

  TChain* chain = new TChain("STNTUPLE");

  // chain->Add("/cdf/scratch/murat/etme01-141218-143297-s232.stn");

  chain->Add("/cdf/scratch/murat/etme01-139339-141216-s232.stn");

  //  chain->Add("/cdf/scratch/murat/etme01-138425-139339-s232.stn");

  TStnAna* x = new TStnAna(chain);

  TStnOutputModule* om;

  om = new TStnOutputModule("./results/junk.root");
  om->SetMaxFileSize(1000);
  // om->DropDataBlock("PhotonBlock");
  // om->DropDataBlock("MuonBlock");

  x->SetOutputModule(om);
  x->Run(500);
}
