// Generates plots of the analysis progression from HPGeWaveformsFromStepPointMCs
// See plotDigitization.sh for usage examples
// Original author: Pawel Plesniak

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

void collectData(const std::string fileName, const std::string treeName, std::vector<double> &chargeCollected, std::vector<double> &chargeDecayed, std::vector<int16_t> &ADCs, std::vector<unsigned int> &eventIds) {
    /*
        Description
            Collects all the required data from virtual detector TTrees

        Arguments
            fileName - as documented in function "plotWaveform"
            treeName - as documented in function "plotWaveform"
            charge - as documented in function "plotWaveform"
            chargeCollected - as documented in function "plotWaveform"
            chargeDecayed - as documented in function "plotWaveform"
            adcs - as documented in function "plotWaveform"
            eventIDs - as documented in function "plotWaveform"

        Variables
            file - ROOT TFile interface
            branches - ROOT TFile interface to TTree branches
            branchNames - vector of branch names
            dataChargeCollected - charge collected from file
            dataChargeDecayed - charge decayed from file
            dataADC - ADC from file
            dataEventId - event ID from file
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
    double dataChargeCollected, dataChargeDecayed;
    int16_t dataADC;
    unsigned int dataEventId;
    if (std::find(branchNames.begin(), branchNames.end(), "chargeCollected") != branchNames.end())
        tree->SetBranchAddress("chargeCollected", &dataChargeCollected);
    else
        Fatal("collectData", "Requested branch 'chargeCollected' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "chargeDecayed") != branchNames.end())
        tree->SetBranchAddress("chargeDecayed", &dataChargeDecayed);
    else
        Fatal("collectData", "Requested branch 'chargeDecayed' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "ADC") != branchNames.end())
        tree->SetBranchAddress("ADC", &dataADC);
    else
        Fatal("collectData", "Requested branch 'ADC' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "eventId") != branchNames.end())
        tree->SetBranchAddress("eventId", &dataEventId);
    else
        Fatal("collectData", "Requested branch 'eventId' does not exist in the file.");

    // Get the number of entries
    int entries = tree->GetEntries();

    // Collect the data
    for (int i = 0; i < entries; i++) {
        // Get the corresponding entry
        tree->GetEntry(i);

        // Collect the data
        chargeCollected.push_back(dataChargeCollected);
        chargeDecayed.push_back(dataChargeDecayed);
        ADCs.push_back(dataADC);
        eventIds.push_back(dataEventId);
    };

    // Check if data has been collected
    if (eventIds.size() == 0)
        Fatal("collectData", "No data was collected from this file");

    // Clear
    file->Close();
    delete file;
    std::cout << "Finished processing file " << fileName << std::endl;
    return;
};

void plot(std::vector<double> &chargeCollected, std::vector<double> &chargeDecayed, std::vector<int16_t> &ADCs, unsigned int eventId) {
    /*
        Description
            Generates the plots

        Arguments
            chargeCollected - vector of charge collected used to generate the plot
            chargeDecayed - vector of charges decayed used to generate the plot
            ADCs - vector of ADCs used to generate the plott
            eventId - eventID used to generate the plot

        Variables
            N - number of data points in the plot
            chargeCollectedA - array of collected charges
            chargeDecayedA - array of decayed charges
            ADCsA - array of ADCs
            i - iterator variable
            time - vector of times to use as x axis in plots
            timeA - array of time
            cChargeCollected - canvas for charge collected plot
            gChargeCollected - plot of charge collectted
            cChargeCollectedFileName - charge collected plot file name
            cChargeDecayed - canvas for charge decayed plot
            gChargeDecayed - plot of charge decayed
            cChargeDecayedFileName - charge decayed plot file name
            cADCs - canvas for ADC plot
            gADCs - plot of ADCs
            cCADCsFileName - ADC plot file name
    */
    // Sanity checks
    const int N = ADCs.size();
    if (chargeCollected.size() != N)
        Fatal("plot", "Incorrect number of chargeCollected entries to generate this plot");
    if (chargeDecayed.size() != N)
        Fatal("plot", "Incorrect number of chargeDecayed entries to generate this plot");

    // Convert the data to the correct type for plots
    double* chargeCollectedA = chargeCollected.data();
    double* chargeDecayedA = chargeDecayed.data();
    double ADCsA[N];
    for (int i = 0; i < N; i++)
        ADCsA[i] = static_cast<double>(ADCs[i]);

    // Construct the time vector for plotting
    std::vector<double> time;
    const double tADC = 3.125;
    for (int i = 0; i < N; i++) time.push_back(i * tADC);
    double* timeA = time.data();

    // Generate the plot for the collected charge
    if (!chargeCollected.empty()) {
        TCanvas* cChargeCollected = new TCanvas("cChargeCollected", "chargeCollected", 800, 600);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.02);
        TGraph *gChargeCollected = new TGraph(N, timeA, chargeCollectedA);
        gChargeCollected->SetTitle("Charge collected;Time [ns];Charge [e]");
        gChargeCollected->Draw("APL");
        std::string cChargeCollectedFileName = "ChargeCollectedEvent" + std::to_string(eventId) + ".png";
        cChargeCollected->SaveAs(cChargeCollectedFileName.c_str());
        gChargeCollected->Delete();
        cChargeCollected->Close();
    };

    // Generate the plot for the decayed charge
    if (!chargeDecayed.empty()) {
        TCanvas* cChargeDecayed = new TCanvas("cChargeDecayed", "chargeDecayed", 800, 600);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.02);
        TGraph *gChargeDecayed = new TGraph(N, timeA, chargeDecayedA);
        gChargeDecayed->SetTitle("Charge deacyed;Time [ns];Charge [e]");
        gChargeDecayed->Draw("APL");
        std::string cChargeDecayedFileName = "ChargeDecayedEvent" + std::to_string(eventId) + ".png";
        cChargeDecayed->SaveAs(cChargeDecayedFileName.c_str());
        gChargeDecayed->Delete();
        cChargeDecayed->Close();
    };

    // Generate the plot for the digitized waveform
    if (!ADCs.empty()) {
        TCanvas* cADCs = new TCanvas("cADCs", "ADCs", 800, 600);
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.02);
        TGraph *gADCs = new TGraph(N, timeA, ADCsA);
        gADCs->SetTitle("Digitized waveform;Time [ns];ADC [arb. unit]");
        gADCs->Draw("APL");
        std::string cADCsFileName = "ADCsEvent" + std::to_string(eventId) + ".png";
        cADCs->SaveAs(cADCsFileName.c_str());
        gADCs->Delete();
        cADCs->Close();
    };

    return;
};

std::vector<unsigned int> makeVectorUnique(std::vector<unsigned int>& v) {
    /*
        Description
            Generate a unique list of event IDs

        Argument
            v - vector to get unique entries from

        Variables
            vSet - set of unique entries from vector v
            vUnique - contains same elements as vSet but as a vector
    */
    // Make the entries unique
    std::set<unsigned int> vSet(v.begin(), v.end());

    // Convert back into a vector
    std::vector<unsigned int> vUnique(vSet.begin(), vSet.end());

    return vUnique;
};

void plotDigitization(const std::string fileName, const std::string treeName, const unsigned int eventID = 0){
    /*
        Description
            Plots the results from HPGeWaveformsFromStepPointMCs

        Arguments
            fileName - root file generated with HPGeWaveformsFromStepPointMCs
            treeName - name of ttree in fileName
            eventID - eventID to plot. If left as 0, plots all the events, otherwise pltos the selected event ID
    */
    // Update global parameters
    SetErrorHandler(customErrorHandler);
    gROOT->SetBatch(kTRUE);

    // Construct the variables to hold the data from the file
    std::vector<double> chargeCollected, chargeDecayed;
    std::vector<int16_t> ADCs;
    std::vector<unsigned int> eventIds;

    // Read in the data from the file
    collectData(fileName, treeName, chargeCollected, chargeDecayed, ADCs, eventIds);
    int N = ADCs.size();

    // Construct the variables to hold tthe selected data for plotting
    std::vector<double> plotChargeCollected, plotChargeDecayed;
    std::vector<int16_t> plotADCs;
    std::vector<unsigned int> plotEventIds = makeVectorUnique(eventIds);

    // Select the relevant event ID if suittable
    if (eventID != 0) {
        if (std::find(plotEventIds.begin(), plotEventIds.end(), eventID) == plotEventIds.end())
            std::cout << "Available Event IDs: ";
            for (unsigned int eventId : plotEventIds)
                std::cout << eventId << ", ";
            std::cout << std::endl;
            Fatal("plotDigitization", "Requested event ID not in the list, exiting.");
    };

    // Select the relevant datta and generate the plots
    for (unsigned int plotEventId : plotEventIds) {
        std::cout << "Generating plot for event with ID " << plotEventId << std::endl;

        // Select the relevant data
        for (int i = 0; i < N; i++) {
            if (eventIds[i] == plotEventId) {
                plotChargeCollected.push_back(chargeCollected[i]);
                plotChargeDecayed.push_back(chargeDecayed[i]);
                plotADCs.push_back(ADCs[i]);
            };
        };

        // Generate the plots
        plot(plotChargeCollected, plotChargeDecayed, plotADCs, plotEventId);

        // Prepare the storage variables for the next event ID
        plotChargeCollected.clear();
        plotChargeDecayed.clear();
        plotADCs.clear();
    };

    return;
};
