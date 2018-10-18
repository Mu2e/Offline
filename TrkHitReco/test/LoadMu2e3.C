{
  TF1 fun("fun","x*x",0);
  gSystem->AddIncludePath("-I/home/online1/sethhirsh/localMu2e/Offline");
  gSystem->AddDynamicPath("/home/online1/sethhirsh/localMu2e/Offline/lib/");
  gSystem->Load("libmu2e_TrkChargeReco");

}
