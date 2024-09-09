// Generates spectra from the .root files generated with MakeTree.fcl
// Original author: Pawel Plesniak
#include <boost/lexical_cast.hpp>

void collectData(std::string filename, int pdgId, std::vector<float> &energy_full, std::vector<float> &energy_reduced, double reduced_spectrum_max_E){
    // Extracts the energy data from the TTrees generated with MakeTree.fcl. The parameters are as follows
    // filename - as per function PlotSpectra()
    // pdgId - as per function PlotSpectra()
    // energy_full - pointer to the vector that contains all the particles used in the spectrum without an upper limit
    // energy_reduced- pointer to the vector that contains all the particles used in the spectrum with an upper limit defined by reduced_spectrum_max_E
    // reduced_spectrum_max_E - as per function PlotSpectra()
    TFile *input = new TFile(filename.c_str());
    TTree *tree = (TTree*)input->Get("Tree/ttree");
    float data_E;   tree->SetBranchAddress("E", &data_E);
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
    input->Close();
    return;
}

void generatePlot(std::vector<float> &energy_data, double bin_width, double scale_factor, std::string particle_name, std::string plot_title_extension, std::string output_file_name_extension){
    // Generates and saves the spectrum plots. The parameters are as follows
    // energy_data - vector of energies extracted from the TTree to be plotted
    // bin_width - bin width for the histogram spectra
    // scale_factor - as per function plotParticle()
    // particle_name - as per function plotParticle()
    // plot_title_extension - used for the plot_title, which becomes "<particle_name><plot_title_extension>"
    // output_file_name_extension - used for the output filename, which becomes "<particle_name><output_file_name_extension>.jpg"
    std::string plot_title       = particle_name + plot_title_extension;
    std::string output_file_name = particle_name + output_file_name_extension;
    double max_energy = *std::max_element(energy_data.begin(), energy_data.end());
    int n_bins = max_energy/bin_width;
    double plot_text_size = 0.06;

    auto c = new TCanvas("c1", "c1", 800, 500);
    auto h = new TH1F(particle_name.c_str(), plot_title.c_str(), n_bins, 0, max_energy);
    gStyle->SetTitleFontSize(plot_text_size);
    gStyle->SetTitleAlign(33);
    gStyle->SetTitleX(0.99);

    c->SetBottomMargin(0.14);
    c->SetLeftMargin(0.175);
    c->SetRightMargin(0.02);

    h->GetXaxis()->SetTitle("Energy [MeV]");
    h->GetXaxis()->SetLabelSize(plot_text_size);
    h->GetXaxis()->SetTitleSize(plot_text_size);
    h->GetYaxis()->SetTitle("Count");
    h->GetYaxis()->SetLabelSize(plot_text_size);
    h->GetYaxis()->SetTitleSize(plot_text_size);

    for (auto e : energy_data)
        h->Fill(e);
    h->Scale(scale_factor);
    h->Draw("HIST");
    gPad->Update();

    TPaveStats *st = (TPaveStats*)h->FindObject("stats");
    st->SetX1NDC(0.7);
    st->SetX2NDC(0.98);
    st->SetY1NDC(0.6);
    st->SetY2NDC(0.9);

    h->Draw("HIST");
    c->SaveAs(output_file_name.c_str());

    c->Close();
    h->Delete();
    return;
}

void plotParticle(std::vector<std::string> filenames, int pdgId, std::string particle_name, double reduced_spectrum_max_E = 2.0, double bin_width = 0.05, const Double_t scale_factor = 1){
    // Generates both the full and reduced spectra for a single particle type. The parameters are as follows
    // filenames - vector of strings containing the path and name of the ROOT TTrees generated with MakeTree.fcl
    // pdgId - pdgID of the particle being plotted
    // particle_name - name of the particle being plotted. This will be used in the plot title and the generated spectrum filename
    // reduced_spectrum_max_E - as per function PlotSpectra()
    // bin_width - as per function PlotSpectra()
    // scale_factor - as per function PlotSpectra()
    std::vector<float> energy_full, energy_reduced;
    for (auto file : filenames)
    {
        std::cout << "Reading file " << file << std::endl;
        collectData(file, pdgId, energy_full, energy_reduced, reduced_spectrum_max_E);
    }
    generatePlot(energy_full,    bin_width, scale_factor, particle_name, " full spectrum",    "Full.jpg");
    generatePlot(energy_reduced, bin_width, scale_factor, particle_name, " reduced spectrum", "Red.jpg");
    return;
}

void PlotSpectra(){
    // Script to generate spectra using the results of MakeTree.fcl. This will generate a full and reduced spectrum for each file named as <particle_name><spectrum_type>.jpg such that <particle_name> is the name of the particle being plotted and <spectrum_type> will be either "Full" or "Red" (reduced to the limit specified by reduced_spectrum_max_E). The parameters are as follows
    // reduced_spectrum_max_E - upper energy threshold of the spectrum to generate
    // bin_width - the bin width in the generated histograms
    // scale_factor - scales the spectrum by this constant. Used for scaling to an equivalent number of POTs
    double reduced_spectrum_max_E = 2.0;
    double bin_width = 0.05;
    double scale_factor = 1;
    std::vector<int> pdgIDs{11, -11, 22, 2112};
    std::vector<std::string> pdgIDNames{"Electron", "Positron", "Photon", "Neutron"};
    std::vector<std::string> filenames{"VD101Mu.root", "VD101Mu3.root"};

    for (int i = 0; i < pdgIDNames.size(); i++)
        plotParticle(filenames, pdgIDs[i], pdgIDNames[i], reduced_spectrum_max_E, bin_width, scale_factor);
    return;
}
