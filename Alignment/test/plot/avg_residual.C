

void avg_residual(std::string filename)
{
    TFile f(filename.c_str());
    TH1F* residusum = (TH1F*)f.Get("AlignTrackCollector/plane_residusum");
    TH1F* plane_trackcount = (TH1F*)f.Get("AlignTrackCollector/plane_trackcount");

    residusum->SetDirectory(0);
    plane_trackcount->SetDirectory(0);

    f.Close();


    TCanvas* t = new TCanvas();
    gStyle->SetOptStat(0);

    residusum->Divide(plane_trackcount);

    residusum->SetTitle("Average hit residual per plane");
    residusum->GetYaxis()->SetTitle("Average hit residual");
    residusum->GetXaxis()->SetTitle("Plane ID [0,35]");

    residusum->GetYaxis()->SetRangeUser(-0.5,0.5);
    
    residusum->Draw();

    TLine* l = new TLine(0,0,36,0);
    l->Draw();
    t->Draw();
    
}
