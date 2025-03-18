// Generates a plot of the digitized waveforms
// See plotDigitizedWaveform.sh for usage examples
// Original author - Pawel Plesniak

#include <TError.h>
#include <iostream>
#include <limits>
#include <math.h>


void customErrorHandler(int level, Bool_t abort, const char* location, const char* message) {
    /*
        Description
            Define a custom error handler that won't print the stack trace but will print an error message and exit.
    */
    std::cerr << message << std::endl;
    if (level > kInfo)
        exit(1);
};

void collectData(const std::string fileName, const std::string treeName, std::vector<int16_t> &ADCs, std::vector<uint32_t> &times, double tMin, double tMax) {
    /*
        Description
            Collects all the required data from virtual detector TTrees

        Arguments
            fileName - as documented in function "plotWaveform"
            treeName - as documented in function "plotWaveform"
            ADCs - as documented in function "plotWaveform"
            times - as documented in function "plotWaveform"
            tMin - as documented in function "plotWaveform"
            tMax - as documented in function "plotWaveform"

        Variables
            file - ROOT TFile interface
            branches - ROOT TFile interface to TTree branches
            branchNames - vector of branch names
            dataADC - ADC from file
            dataTime - time from file
            entries - number of entries in the TTree
    */
    std::cout << "Processing file " << fileName << std::endl;

    // Get the branch
    TFile *file = new TFile(fileName.c_str());
    if (!file || file->IsZombie()) {
        Fatal("collectData", "Failed to open the file.");
    };
    TTree *tree = (TTree*)file->Get(treeName.c_str());
    if (!tree)
        Fatal("collectData", "Requested tree does not exist in the file.");

    // Get the list of branches to check if they exist
    TObjArray *branches = tree->GetListOfBranches();
    std::vector<std::string> branchNames;
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch *branch = dynamic_cast<TBranch*>(branches->At(i));
        if (branch)
            branchNames.push_back(branch->GetName());
    };

    // If the branches exist, assign them to the appropriate variables
    int16_t dataADC;
    uint32_t dataTime;
    if (std::find(branchNames.begin(), branchNames.end(), "ADC") != branchNames.end())
        tree->SetBranchAddress("ADC", &dataADC);
    else
        Fatal("collectData", "Requested branch 'ADC' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "time") != branchNames.end())
        tree->SetBranchAddress("time", &dataTime);
    else
        Fatal("collectData", "Requested branch 'time' does not exist in the file.");

    // Get the number of entries
    int entries = tree->GetEntries();

    // Set the limit variables
    const double tADC = 3.125;
    const uint32_t tMinADC = static_cast<uint32_t>(tMin/tADC), tMaxADC = static_cast<uint32_t>(tMax/tADC);

    // Collect the data
    for (int i = 0; i < entries; i++) {
        // Get the corresponding entry
        tree->GetEntry(i);

        // Collect the data
        if ((tMinADC == 0 || dataTime > tMinADC) && (tMaxADC == 0 || dataTime < tMaxADC)){
            ADCs.push_back(dataADC);
            times.push_back(dataTime);
        };
    };
    // Check if data has been collected
    if (ADCs.size() == 0)
        Fatal("collectData", "No data was collected from this file");

    // Clear
    file->Close();
    delete file;
    std::cout << "Finished processing file " << fileName << std::endl;
    return;
};

void plotDigitizedWaveforms(const std::string fileName, const std::string treeName, double tMin = 0.0, double tMax = 0.0) {
    /*
        Description
            Generates plots of the digitized waveforms generated with HPGeWaveformsFromStepPointMCs

        Arguments
            fileName - name of ROOT file generated with HPGeWaveformsFromStepPointMCs
            treeName - name of tree containing all the data in the file
            tMin - if non-zero, applies a minimum time cut to the data used for plotting
            tMin - if non-zero, applies a maximum time cut to the data used for plotting

        Variables
            ADCs - vector of ADCs from the file
            times - vector of times from the file
            plotADCs - vector of ADCs for plotting cast to the correct data type
            plotTimes - vector of times cast from ADC clock ticks to ns
            tADC - ADC clock tick in ns
            ADC - iterator variable
            c - canvas used to save the digitized waveform plot
            g - graph of digitized waveform
            cFileName - file name of generated plot
    */

    // Update global parameters
    SetErrorHandler(customErrorHandler);
    gROOT->SetBatch(kTRUE);

    // Construct the variables used to collect data from the file
    std::vector<int16_t> ADCs;
    std::vector<uint32_t> times;

    // Collect the data
    collectData(fileName, treeName, ADCs, times, tMin, tMax);

    // Set up the variables used to store the data used to generate the plot
    std::vector<double> plotADCs, plotTimes;

    // Convert the collected data to the correct types for the plots
    const double tADC = 3.125;
    for (auto ADC : ADCs)
        plotADCs.push_back(static_cast<double>(ADC));
    ADCs.clear();
    ADCs.shrink_to_fit();
    for (auto time : times)
        plotTimes.push_back(time * tADC);
    times.clear();
    times.shrink_to_fit();

    // Sanity checks
    if (plotADCs.empty())
        Fatal("plotConcatenated", "ADCs vector is empty");
    if (plotTimes.empty())
        Fatal("plotConcatenated", "Time vector is empty");

    // Construct the canvas used to save the plot
    TCanvas* c = new TCanvas("c", "c", 1500, 1000);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.08);

    // Construct the graph
    TGraph *g = new TGraph(plotADCs.size(), plotTimes.data(), plotADCs.data());
    g->Draw("APL");
    g->SetTitle("Concatenated waveform;Time [ns];ADC [arb. unit]");

    // Save the file
    std::string cFileName = "plotDigitizedWaveform.png";
    c->SaveAs(cFileName.c_str());

    // Cleanup
    g->Delete();
    c->Close();

    return;
};
