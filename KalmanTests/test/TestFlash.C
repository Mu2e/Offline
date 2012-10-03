void FlashVD(TTree* ntvd,const char* vdname) {
	
  TH2F* xy = new TH2F("xy","Y vs X position at VD",100,-4200,-3600,100,-300,300);
  TH1F* pdg = new TH1F("pdg","PDG code",200,-250,2250);
  TH1F* emt = new TH1F("emt","Time at VD",200,0,2000);
  TH1F* ept = new TH1F("ept","Time at VD",200,0,2000);
  TH1F* pt = new TH1F("pt","Time at VD",200,0,2000);
  TH1F* mumt = new TH1F("mumt","Time at VD",200,0,2000);
  TH1F* mupt = new TH1F("mupt","Time at VD",200,0,2000);
  TH1F* pimt = new TH1F("pimt","Time at VD",200,0,2000);
  TH1F* pht = new TH1F("pht","Time at VD",200,0,2000);
  TH1F* nt = new TH1F("nt","Time at VD",200,0,2000);
  emt->SetLineColor(kRed);
  ept->SetLineColor(kBlack);
  pt->SetLineColor(kOrange);
  mumt->SetLineColor(kBlue);
  mupt->SetLineColor(kYellow);
  pimt->SetLineColor(kGreen);
  pht->SetLineColor(kMagenta);
  nt->SetLineColor(kCyan);
  emt->SetStats(0);

  TH1F* emke = new TH1F("emke","Kinetic Energy at VD",100,0,100);
  TH1F* epke = new TH1F("epke","Kinetic Energy at VD",100,0,100);
  TH1F* pke = new TH1F("pke","Kinetic Energy at VD",100,0,100);
  TH1F* mumke = new TH1F("mumke","Kinetic Energy at VD",100,0,100);
  TH1F* mupke = new TH1F("mupke","Kinetic Energy at VD",100,0,100);
  TH1F* pimke = new TH1F("pimke","Kinetic Energy at VD",100,0,100);
  TH1F* phke = new TH1F("phke","Kinetic Energy at VD",100,0,100);
  TH1F* nke = new TH1F("nke","Kinetic Energy at VD",100,0,100);

  emke->SetLineColor(kRed);
  epke->SetLineColor(kBlack);
  pke->SetLineColor(kOrange);
  mumke->SetLineColor(kBlue);
  mupke->SetLineColor(kYellow);
  pimke->SetLineColor(kGreen);
  phke->SetLineColor(kMagenta);
  nke->SetLineColor(kCyan);
  emke->SetStats(0);

  TH1F* emct = new TH1F("emct","Cos(#theta) at VD",100,-1.1,1.1);
  TH1F* epct = new TH1F("epct","Cos(#theta) at VD",100,-1.1,1.1);
  TH1F* pct = new TH1F("pct","Cos(#theta) at VD",100,-1.1,1.1);
  TH1F* mumct = new TH1F("mumct","Cos(#theta) at VD",100,-1.1,1.1);
  TH1F* mupct = new TH1F("mupct","Cos(#theta) at VD",100,-1.1,1.1);
  TH1F* pimct = new TH1F("pimct","Cos(#theta) at VD",100,-1.1,1.1);
  TH1F* phct = new TH1F("phct","Cos(#theta) at VD",100,-1.1,1.1);
  TH1F* nct = new TH1F("nct","Cos(#theta) at VD",100,-1.1,1.1);

  emct->SetLineColor(kRed);
  epct->SetLineColor(kBlack);
  pct->SetLineColor(kOrange);
  mumct->SetLineColor(kBlue);
  mupct->SetLineColor(kYellow);
  pimct->SetLineColor(kGreen);
  phct->SetLineColor(kMagenta);
  nct->SetLineColor(kCyan);
  emct->SetStats(0);

  ntvd->Project("xy","y:x");
  ntvd->Project("pdg","pdg");

  ntvd->Project("emt","time","pdg==11");
  ntvd->Project("ept","time","pdg==-11");
  ntvd->Project("pt","time","pdg==2212");
  ntvd->Project("mumt","time","pdg==13");
  ntvd->Project("mupt","time","pdg==-13");
  ntvd->Project("pimt","time","pdg==-211");
  ntvd->Project("pht","time","pdg==22");
  ntvd->Project("nt","time","pdg==2112");
  emt->SetMinimum(1);

  ntvd->Project("emke","ke","pdg==11");
  ntvd->Project("epke","ke","pdg==-11");
  ntvd->Project("pke","ke","pdg==2212");
  ntvd->Project("mumke","ke","pdg==13");
  ntvd->Project("mupke","ke","pdg==-13");
  ntvd->Project("pimke","ke","pdg==-211");
  ntvd->Project("phke","ke","pdg==22");
  ntvd->Project("nke","ke","pdg==2112");
  emke->SetMinimum(1);

  ntvd->Project("emct","pz/sqrt(px^2+py^2+pz^2)","pdg==11");
  ntvd->Project("epct","pz/sqrt(px^2+py^2+pz^2)","pdg==-11");
  ntvd->Project("pct","pz/sqrt(px^2+py^2+pz^2)","pdg==2212");
  ntvd->Project("mumct","pz/sqrt(px^2+py^2+pz^2)","pdg==13");
  ntvd->Project("mupct","pz/sqrt(px^2+py^2+pz^2)","pdg==-13");
  ntvd->Project("pimct","pz/sqrt(px^2+py^2+pz^2)","pdg==-211");
  ntvd->Project("phct","pz/sqrt(px^2+py^2+pz^2)","pdg==22");
  ntvd->Project("nct","pz/sqrt(px^2+py^2+pz^2)","pdg==2112");

  TCanvas* flashcan = new TCanvas("flashcan","flashcan",1200,800);
  flashcan->Divide(2,2);
  flashcan->cd(1);
  xy->Draw("box");
  
	flashcan->cd(2);
  emct->Draw();
  epct->Draw("same");
  mumct->Draw("same");
	mupct->Draw("same");
  pct->Draw("same");
  pimct->Draw("same");
  phct->Draw("same");
  nct->Draw("same");
  
	flashcan->cd(3);
  gPad->SetLogy();
  emt->Draw();
  ept->Draw("same");
  mumt->Draw("same");
  mupt->Draw("same");
  pt->Draw("same");
  pimt->Draw("same");
  pht->Draw("same");
  nt->Draw("same");
  TLegend* leg = new TLegend(0.6,0.4,0.8,0.9);
  leg->AddEntry(emke,"e^{-}","l");
  leg->AddEntry(epke,"e^{+}","l");
  leg->AddEntry(mumke,"#mu^{-}","l");
  leg->AddEntry(mupke,"#mu^{+}","l");
  leg->AddEntry(pke,"P^{+}","l");
  leg->AddEntry(pimke,"#pi^{-}","l");
  leg->AddEntry(phke,"#gamma","l");
  leg->AddEntry(nke,"neutron","l");
  leg->Draw();

  flashcan->cd(4);
  gPad->SetLogy();
  emke->Draw();
  epke->Draw("same");
  mumke->Draw("same");
  mupke->Draw("same");
  pke->Draw("same");
  pimke->Draw("same");
  phke->Draw("same");
  nke->Draw("same");
}

void TestHits(TTree* sh) {
  TH1F* emhtime = new TH1F("emhtime","Hit time;time(ns)",200,0,2000);
  TH1F* ephtime = new TH1F("ephtime","Hit time;time(ns)",200,0,2000);
  TH1F* phtime = new TH1F("phtime","Hit time;time(ns)",200,0,2000);
  TH1F* mumhtime = new TH1F("mumhtime","Hit time;time(ns)",200,0,2000);
  TH1F* muphtime = new TH1F("muphtime","Hit time;time(ns)",200,0,2000);
  TH1F* pimhtime = new TH1F("pimhtime","Hit time;time(ns)",200,0,2000);
  TH1F* phhtime = new TH1F("phhtime","Hit time;time(ns)",200,0,2000);
  TH1F* nhtime = new TH1F("nhtime","Hit time;time(ns)",200,0,2000);
  TH1F* ahtime = new TH1F("ahtime","Hit time;time(ns)",200,0,2000);

  emhtime->SetLineColor(kRed);
  ephtime->SetLineColor(kBlack);
  phtime->SetLineColor(kOrange);
  mumhtime->SetLineColor(kBlue);
  muphtime->SetLineColor(kYellow);
  pimhtime->SetLineColor(kGreen);
  phhtime->SetLineColor(kMagenta);
  nhtime->SetLineColor(kCyan);
  ahtime->SetLineColor(kSpring);
  emhtime->SetStats(0);

  sh->Project("emhtime","hittime","mcppdg==11");
  sh->Project("ephtime","hittime","mcppdg==-11");
  sh->Project("phtime","hittime","mcppdg==2212");
  sh->Project("mumhtime","hittime","mcppdg==13");
  sh->Project("muphtime","hittime","mcppdg==-13");
  sh->Project("pimhtime","hittime","mcppdg==-211");
  sh->Project("phhtime","hittime","mcppdg==22");
  sh->Project("nhtime","hittime","mcppdg==2112");
  sh->Project("ahtime","hittime");
  ahtime->SetMinimum(1);

	TH1F* emstrawz = new TH1F("emstrawz","Hit z;z(mm)",200,-1600,1600);
  TH1F* epstrawz = new TH1F("epstrawz","Hit z;z(mm)",200,-1600,1600);
  TH1F* pstrawz = new TH1F("pstrawz","Hit z;z(mm)",200,-1600,1600);
  TH1F* mumstrawz = new TH1F("mumstrawz","Hit z;z(mm)",200,-1600,1600);
  TH1F* mupstrawz = new TH1F("mupstrawz","Hit z;z(mm)",200,-1600,1600);
  TH1F* pimstrawz = new TH1F("pimstrawz","Hit z;z(mm)",200,-1600,1600);
  TH1F* phstrawz = new TH1F("phstrawz","Hit z;z(mm)",200,-1600,1600);
  TH1F* nstrawz = new TH1F("nstrawz","Hit z;z(mm)",200,-1600,1600);
  TH1F* astrawz = new TH1F("astrawz","Hit z;z(mm)",200,-1600,1600);

  emstrawz->SetLineColor(kRed);
  epstrawz->SetLineColor(kBlack);
  pstrawz->SetLineColor(kOrange);
  mumstrawz->SetLineColor(kBlue);
  mupstrawz->SetLineColor(kYellow);
  pimstrawz->SetLineColor(kGreen);
  phstrawz->SetLineColor(kMagenta);
  nstrawz->SetLineColor(kCyan);
  astrawz->SetLineColor(kSpring);
  astrawz->SetStats(0);

  sh->Project("emstrawz","strawz","mcppdg==11");
  sh->Project("epstrawz","strawz","mcppdg==-11");
  sh->Project("pstrawz","strawz","mcppdg==2212");
  sh->Project("mumstrawz","strawz","mcppdg==13");
  sh->Project("mupstrawz","strawz","mcppdg==-13");
  sh->Project("pimstrawz","strawz","mcppdg==-211");
  sh->Project("phstrawz","strawz","mcppdg==22");
  sh->Project("nstrawz","strawz","mcppdg==2112");
  sh->Project("astrawz","strawz");

	TH1F* emstrawr = new TH1F("emstrawr","Hit #rho;#rho(mm)",200,350,700);
  TH1F* epstrawr = new TH1F("epstrawr","Hit #rho;#rho(mm)",200,350,700);
  TH1F* pstrawr = new TH1F("pstrawr","Hit #rho;#rho(mm)",200,350,700);
  TH1F* mumstrawr = new TH1F("mumstrawr","Hit #rho;#rho(mm)",200,350,700);
  TH1F* mupstrawr = new TH1F("mupstrawr","Hit #rho;#rho(mm)",200,350,700);
  TH1F* pimstrawr = new TH1F("pimstrawr","Hit #rho;#rho(mm)",200,350,700);
  TH1F* phstrawr = new TH1F("phstrawr","Hit #rho;#rho(mm)",200,350,700);
  TH1F* nstrawr = new TH1F("nstrawr","Hit #rho;#rho(mm)",200,350,700);
  TH1F* astrawr = new TH1F("astrawr","Hit #rho;#rho(mm)",200,350,700);

  emstrawr->SetLineColor(kRed);
  epstrawr->SetLineColor(kBlack);
  pstrawr->SetLineColor(kOrange);
  mumstrawr->SetLineColor(kBlue);
  mupstrawr->SetLineColor(kYellow);
  pimstrawr->SetLineColor(kGreen);
  phstrawr->SetLineColor(kMagenta);
  nstrawr->SetLineColor(kCyan);
  astrawr->SetLineColor(kSpring);
  astrawr->SetStats(0);

  sh->Project("emstrawr","strawrho","mcppdg==11");
  sh->Project("epstrawr","strawrho","mcppdg==-11");
  sh->Project("pstrawr","strawrho","mcppdg==2212");
  sh->Project("mumstrawr","strawrho","mcppdg==13");
  sh->Project("mupstrawr","strawrho","mcppdg==-13");
  sh->Project("pimstrawr","strawrho","mcppdg==-211");
  sh->Project("phstrawr","strawrho","mcppdg==22");
  sh->Project("nstrawr","strawrho","mcppdg==2112");
  sh->Project("astrawr","strawrho");

	TCut late("hittime>700");

	TH1F* emlatestrawz = new TH1F("emlatestrawz","Hit z, t>700 ns;z(mm)",50,-1600,1600);
  TH1F* eplatestrawz = new TH1F("eplatestrawz","Hit z, t>700 ns;z(mm)",50,-1600,1600);
  TH1F* platestrawz = new TH1F("platestrawz","Hit z, t>700 ns;z(mm)",50,-1600,1600);
  TH1F* mumlatestrawz = new TH1F("mumlatestrawz","Hit z, t>700 ns;z(mm)",50,-1600,1600);
  TH1F* muplatestrawz = new TH1F("muplatestrawz","Hit z, t>700 ns;z(mm)",50,-1600,1600);
  TH1F* pimlatestrawz = new TH1F("pimlatestrawz","Hit z, t>700 ns;z(mm)",50,-1600,1600);
  TH1F* phlatestrawz = new TH1F("phlatestrawz","Hit z, t>700 ns;z(mm)",50,-1600,1600);
  TH1F* nlatestrawz = new TH1F("nlatestrawz","Hit z, t>700 ns;z(mm)",50,-1600,1600);
  TH1F* alatestrawz = new TH1F("alatestrawz","Hit z, t>700 ns;z(mm)",50,-1600,1600);

  emlatestrawz->SetLineColor(kRed);
  eplatestrawz->SetLineColor(kBlack);
  platestrawz->SetLineColor(kOrange);
  mumlatestrawz->SetLineColor(kBlue);
  muplatestrawz->SetLineColor(kYellow);
  pimlatestrawz->SetLineColor(kGreen);
  phlatestrawz->SetLineColor(kMagenta);
  nlatestrawz->SetLineColor(kCyan);
  alatestrawz->SetLineColor(kSpring);
  alatestrawz->SetStats(0);

  sh->Project("emlatestrawz","strawz",late+"mcppdg==11");
  sh->Project("eplatestrawz","strawz",late+"mcppdg==-11");
  sh->Project("platestrawz","strawz",late+"mcppdg==2212");
  sh->Project("mumlatestrawz","strawz",late+"mcppdg==13");
  sh->Project("muplatestrawz","strawz",late+"mcppdg==-13");
  sh->Project("pimlatestrawz","strawz",late+"mcppdg==-211");
  sh->Project("phlatestrawz","strawz",late+"mcppdg==22");
  sh->Project("nlatestrawz","strawz",late+"mcppdg==2112");
  sh->Project("alatestrawz","strawz",late);

	TH1F* emlatestrawr = new TH1F("emlatestrawr","Hit #rho, t>700 ns;#rho(mm)",50,350,700);
  TH1F* eplatestrawr = new TH1F("eplatestrawr","Hit #rho, t>700 ns;#rho(mm)",50,350,700);
  TH1F* platestrawr = new TH1F("platestrawr","Hit #rho, t>700 ns;#rho(mm)",50,350,700);
  TH1F* mumlatestrawr = new TH1F("mumlatestrawr","Hit #rho, t>700 ns;#rho(mm)",50,350,700);
  TH1F* muplatestrawr = new TH1F("muplatestrawr","Hit #rho, t>700 ns;#rho(mm)",50,350,700);
  TH1F* pimlatestrawr = new TH1F("pimlatestrawr","Hit #rho, t>700 ns;#rho(mm)",50,350,700);
  TH1F* phlatestrawr = new TH1F("phlatestrawr","Hit #rho, t>700 ns;#rho(mm)",50,350,700);
  TH1F* nlatestrawr = new TH1F("nlatestrawr","Hit #rho, t>700 ns;#rho(mm)",50,350,700);
  TH1F* alatestrawr = new TH1F("alatestrawr","Hit #rho, t>700 ns;#rho(mm)",50,350,700);

  emlatestrawr->SetLineColor(kRed);
  eplatestrawr->SetLineColor(kBlack);
  platestrawr->SetLineColor(kOrange);
  mumlatestrawr->SetLineColor(kBlue);
  muplatestrawr->SetLineColor(kYellow);
  pimlatestrawr->SetLineColor(kGreen);
  phlatestrawr->SetLineColor(kMagenta);
  nlatestrawr->SetLineColor(kCyan);
  alatestrawr->SetLineColor(kSpring);
  alatestrawr->SetStats(0);

  sh->Project("emlatestrawr","strawrho",late+"mcppdg==11");
  sh->Project("eplatestrawr","strawrho",late+"mcppdg==-11");
  sh->Project("platestrawr","strawrho",late+"mcppdg==2212");
  sh->Project("mumlatestrawr","strawrho",late+"mcppdg==13");
  sh->Project("muplatestrawr","strawrho",late+"mcppdg==-13");
  sh->Project("pimlatestrawr","strawrho",late+"mcppdg==-211");
  sh->Project("phlatestrawr","strawrho",late+"mcppdg==22");
  sh->Project("nlatestrawr","strawrho",late+"mcppdg==2112");
  sh->Project("alatestrawr","strawrho",late);

	TH2F* pdg = new TH2F("pdg","Direct PDG Id vs Parent PDG Id",51,-25.5,25.5,51,-25.5,25.5);
	pdg->SetStats(0);
	sh->Project("pdg","mcpdg:mcppdg");


  TCanvas* hcan = new TCanvas("hcan","hcan",1200,800);
	hcan->Clear();
	hcan->Divide(2,3);
	
	hcan->cd(1);
	gPad->SetLogy();
  ahtime->Draw();
  emhtime->Draw("same");
  ephtime->Draw("same");
  mumhtime->Draw("same");
  muphtime->Draw("same");
  phtime->Draw("same");
  pimhtime->Draw("same");
  phhtime->Draw("same");
  nhtime->Draw("same");

	hcan->cd(2);
	pdg->Draw("box");

	hcan->cd(3);
  astrawz->Draw();
  emstrawz->Draw("same");
  epstrawz->Draw("same");
  mumstrawz->Draw("same");
  mupstrawz->Draw("same");
  pstrawz->Draw("same");
  pimstrawz->Draw("same");
  phstrawz->Draw("same");
  nstrawz->Draw("same");

	hcan->cd(4);
  alatestrawz->Draw();
  emlatestrawz->Draw("same");
  phlatestrawz->Draw("same");
  eplatestrawz->Draw("same");
  mumlatestrawz->Draw("same");
  muplatestrawz->Draw("same");
  pimlatestrawz->Draw("same");
  platestrawz->Draw("same");
  nlatestrawz->Draw("same");

	hcan->cd(5);
  astrawr->Draw();
  emstrawr->Draw("same");
  epstrawr->Draw("same");
  mumstrawr->Draw("same");
  mupstrawr->Draw("same");
  pstrawr->Draw("same");
  pimstrawr->Draw("same");
  phstrawr->Draw("same");
  nstrawr->Draw("same");

	TLegend* leg = new TLegend(0.6,0.35,0.8,0.9);
  leg->AddEntry(ahtime,"All parents","l");
  leg->AddEntry(emhtime,"parent e^{-}","l");
  leg->AddEntry(ephtime,"parent e^{+}","l");
  leg->AddEntry(mumhtime,"parent #mu^{-}","l");
  leg->AddEntry(muphtime,"parent #mu^{+}","l");
  leg->AddEntry(phtime,"parent P^{+}","l");
  leg->AddEntry(pimhtime,"parent #pi^{-}","l");
  leg->AddEntry(phhtime,"parent #gamma","l");
  leg->AddEntry(nhtime,"parent n","l");
  leg->Draw();
	
	hcan->cd(6);
  alatestrawr->Draw();
  emlatestrawr->Draw("same");
  phlatestrawr->Draw("same");
  eplatestrawr->Draw("same");
  mumlatestrawr->Draw("same");
  muplatestrawr->Draw("same");
  pimlatestrawr->Draw("same");
  platestrawr->Draw("same");
  nlatestrawr->Draw("same");

}

