// Generates overlapping spectra from the .root files generated with MakeTree.fcl
// Original author: Pawel Plesniak
#include <boost/lexical_cast.hpp>

void collectData(std::string filename, int pdgId, std::vector<float> &energy_full, std::vector<float> &energy_reduced, double reduced_spectrum_max_E){
    // Extracts the energy data from the TTrees generated with MakeTree.fcl. The parameters are as follows
    // filename - as per function PlotOverlapSpectra()
    // pdgId - as per function PlotOverlapSpectra()
    // energy_full - pointer to the vector that contains all the particles used in the spectrum without an upper limit
    // energy_reduced- pointer to the vector that contains all the particles used in the spectrum with an upper limit defined by reduced_spectrum_max_E
    // reduced_spectrum_max_E - as per function PlotOverlapSpectra()
    TFile *input = new TFile(filename.c_str());
    TTree *tree = (TTree*)input->Get("Tree/ttree");
    float data_E; tree->SetBranchAddress("E", &data_E);
    int data_pdgID; tree->SetBranchAddress("pdgId", &data_pdgID);
    int entries = tree->GetEntries();
    for (int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);
        if (pdgId == data_pdgID)
        {
            energy_full.push_back(data_E);
            if ((reduced_spectrum_max_E != 0) && (data_E < reduced_spectrum_max_E))
            energy_reduced.push_back(data_E);
        }
    }
    double max_energy = *std::max_element(energy_full.begin(), energy_full.end());
    if (max_energy < reduced_spectrum_max_E)
        std::cout << "WARNING! The maximum particle energy in this spectrum (" << max_energy << ") is smaller than the requested cut (" << reduced_spectrum_max_E << ")." << std::endl;
    return;
}

void generateOverlapPlot(std::vector<float> &electron_energy, std::vector<float> &muon_energy, double bin_width, std::string particle_name, std::string plot_title_extension, std::string output_file_name_extension, const Double_t scale_factor = 1){
    // Generates and saves the spectrum plots. The parameters are as follows
    // electron_energy - vector of electron energies extracted from the TTree to be plotted
    // muon_energy - vector of muon energies extracted from the TTree to be plotted
    // bin_width - as per function PlotOverlapSpectra()
    // particle_name - as per function PlotOverlapSpectra()
    // plot_title_extension - used for the plot_title, which becomes "<particle_name><plot_title_extension>"
    // output_file_name_extension - used for the output filename, which becomes "<particle_name><output_file_name_extension>.jpg"
    // scale_factor - as per function PlotOverlapSpectra()
    std::string plot_title       = particle_name + plot_title_extension;
    std::string output_file_name = particle_name + output_file_name_extension;
    double max_energy_electron   = *std::max_element(electron_energy.begin(), electron_energy.end());
    double max_muon_electron     = *std::max_element(muon_energy.begin(),     muon_energy.end());
    double max_energy            = std::max(max_energy_electron, max_muon_electron);
    int n_bins = max_energy/bin_width;
    double plot_text_size = 0.06;
    plot_title = plot_title + ";Energy [MeV];Count";

    auto hs = new THStack("hs", plot_title.c_str());
    gStyle->SetTitleFontSize(plot_text_size);
    gStyle->SetTitleAlign(33);
    gStyle->SetTitleX(0.99);

    auto h1 = new TH1F("EleBeamCat", "EleBeamCat", n_bins, 0, max_energy);
    for (auto e1 : electron_energy)
        h1->Fill(e1);
    h1->Scale(scale_factor);

    auto h2 = new TH1F("MuBeamCat", "MuBeamCat", n_bins, 0, max_energy);
    for (auto e3 : muon_energy){
        h2->Fill(e3);
    }
    h1->SetFillColor(kRed);
    h1->SetLineWidth(4);
    h1->SetLineColor(kRed);
    h2->SetFillColor(kBlue);
    h2->SetLineWidth(4);
    h2->SetLineColor(kBlue);


    // auto c = new TCanvas("c", "c", 800, 500);
    auto c = new TCanvas("c", "c", 1500, 900);
    auto h = new TH1F("Combined", "Combined", n_bins, 0, max_energy);
    h->SetStats(0);
    hs->SetHistogram(h);
    hs->GetHistogram()->GetXaxis()->SetLabelSize(plot_text_size);
    hs->GetHistogram()->GetXaxis()->SetTitleSize(plot_text_size);
    hs->GetHistogram()->GetYaxis()->SetLabelSize(plot_text_size);
    hs->GetHistogram()->GetYaxis()->SetTitleSize(plot_text_size);

    c->SetBottomMargin(0.14);
    c->SetLeftMargin(0.175);
    c->SetRightMargin(0.02);
    hs->Add(h1);
    hs->Add(h2);
    hs->Draw("HIST");

    TLegend *legend = new TLegend(0.6, 0.5, 0.98, 0.9, "Input dataset");
    legend->SetHeader("Dataset", "C");
    legend->AddEntry("EleBeamCat", "EleBeamCat", "l");
    legend->AddEntry("MuBeamCat", "MuBeamCat", "l");
    legend->Draw();

    c->SaveAs(output_file_name.c_str());

    h1->Delete();
    h2->Delete();
    hs->Delete();
    c->Close();
    return;
}

void PlotOverlapSpectra(){
    // Script to generate spectra using the results of MakeTree.fcl. This will generate a full and reduced spectrum for each file named as <particle_name><spectrum_type>.jpg such that <particle_name> is the name of the particle being plotted and <spectrum_type> will be either "Full" or "Red" (reduced to the limit specified by reduced_spectrum_max_E). The parameters are as follows
    // reduced_spectrum_max_E - upper energy threshold of the spectrum to generate
    // bin_width - the bin width in the generated histograms
    // scale_factor - scales the electron spectrum by this constant. Used for scaling to an equivalent number of POTs
    int reduced_spectrum_max_E = 2.0;
    double bin_width = 0.05;
    double scale_factor = 42;
    std::vector<int> pdgIDs{11, -11, 22, 2112};
    std::vector<std::string> pdgIDNames{"Electron", "Positron", "Photon", "Neutron"};

    std::vector<float> electron_energy_full, electron_energy_reduced;
    std::vector<float> muon_energy_full, muon_energy_reduced;

    int pdgID = pdgIDs[0];
    std::string pdgIDName = pdgIDNames[0];
    for (int i = 0; i < pdgIDs.size(); i++)
    {
        int pdgID = pdgIDs[i];
        std::string pdgIDName = pdgIDNames[i];
        collectData("VD101Ele.root", pdgID, electron_energy_full, electron_energy_reduced, reduced_spectrum_max_E);
        collectData("VD101Mu.root",  pdgID, muon_energy_full,     muon_energy_reduced,     reduced_spectrum_max_E);
        collectData("VD101Mu3.root", pdgID, muon_energy_full,     muon_energy_reduced,     reduced_spectrum_max_E);
        if (scale_factor == 1)
        {
            generateOverlapPlot(electron_energy_full,    muon_energy_full,    bin_width, pdgIDName, " full spectrum stacked",    "Full.jpg", scale_factor);
            generateOverlapPlot(electron_energy_reduced, muon_energy_reduced, bin_width, pdgIDName, " reduced spectrum stacked", "Red.jpg" , scale_factor);
        }
        else
        {
            generateOverlapPlot(electron_energy_full,    muon_energy_full,    bin_width, pdgIDName, " full spectrum stacked and scaled",    "Full.jpg", scale_factor);
            generateOverlapPlot(electron_energy_reduced, muon_energy_reduced, bin_width, pdgIDName, " reduced spectrum stacked and scaled", "Red.jpg" , scale_factor);
        }
        electron_energy_full.clear();
        electron_energy_reduced.clear();
        muon_energy_full.clear();
        muon_energy_reduced.clear();
    }
    return;
}
