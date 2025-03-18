// Plots the MWD analysis chain
// see plotMWDAnalysis.sh for usage example
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

void collectMWDData(const std::string fileName, const std::string treeName, std::vector<int16_t> &ADCs, std::vector<double> &deconvoluted, std::vector<double> &differentiated, std::vector<double> &averaged, std::vector<uint32_t> &times, std::vector<uint> &eventIds) {
    /*
        Description
            Collects all the required data from virtual detector TTrees

        Arguments
            fileName - as documented in function "plotMWDResults"
            treeName - as documented in function "plotMWDResults"
            ADCs - as documented in function "plotMWDResults"
            deconvoluted - as documented in function "plotMWDResults"
            differentiated - as documented in function "plotMWDResults"
            averaged - as documented in function "plotMWDResults"
            times - as documented in function "plotMWDResults"
            eventIds - as documented in function "plotMWDResults"

        Variables
            file - ROOT TFile interface
            branches - ROOT TFile interface to TTree branches
            branchNames - vector of branch names
            dataADC - ADC from file
            dataDeconvoluted - deconvoluted data from file
            dataDifferentiated - differentiated data from file
            dataAveraged - averaged data from file
            dataTime - time from file
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
    int16_t dataADC;
    double dataDeconvoluted, dataDifferentiated, dataAveraged;
    uint32_t dataTime;
    uint dataEventId;
    if (std::find(branchNames.begin(), branchNames.end(), "ADC") != branchNames.end())
        tree->SetBranchAddress("ADC", &dataADC);
    else
        Fatal("collectData", "Requested branch 'ADC' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "deconvoluted") != branchNames.end())
        tree->SetBranchAddress("deconvoluted", &dataDeconvoluted);
    else
        Fatal("collectData", "Requested branch 'deconvoluted' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "differentiated") != branchNames.end())
        tree->SetBranchAddress("differentiated", &dataDifferentiated);
    else
        Fatal("collectData", "Requested branch 'differentiated' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "averaged") != branchNames.end())
        tree->SetBranchAddress("averaged", &dataAveraged);
    else
        Fatal("collectData", "Requested branch 'averaged' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "time") != branchNames.end())
        tree->SetBranchAddress("time", &dataTime);
    else
        Fatal("collectData", "Requested branch 'time' does not exist in the file.");
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
        ADCs.push_back(dataADC);
        deconvoluted.push_back(dataDeconvoluted);
        differentiated.push_back(dataDifferentiated);
        averaged.push_back(dataAveraged);
        times.push_back(dataTime);
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

uint32_t normalizeTimeToFirstMicrospill(uint32_t time, bool subtract = false) {
    /*
        Description
            Converts the time from continuous time to modulated within the first microspill, returning it as a number of ADC clock ticks
            Note - this currnetly only  works for 320MHz sampling

        Arguments
            time - input time to convert from absolute timing to microspill timing, with the start of the microspill starting at t=0
            subtract - subtract one from the timing. TODO - find out why this is needed and resolve the issue

        Variables
            returnTime - time to return
            fineTimeRecoQuot - number of microspills that fit in a group of 5. This is used to get the reconstructed time as close to the real time as possible, minimizing the drift
    */
    uint32_t returnTime = time % 2712;
    int fineTimeRecoQuot = returnTime / 542;
    returnTime = returnTime % 542;
    if (subtract) {
        returnTime--;
    };
    if (fineTimeRecoQuot == 0 || fineTimeRecoQuot == 3)
        returnTime++;
    return returnTime;
};

void plot(std::vector<int16_t> &ADCs, std::vector<double> &deconvoluted, std::vector<double> &differentiated, std::vector<double> &averaged, std::vector<uint32_t> &times, const unsigned int eventId, const double threshold, const bool plotMicrospill) {
    /*
        Description
            Generates the plots

        Arguments
            ADCs - vector of ADC values to plot
            deconvoluted - vector of deconvoluted data to plot
            differentiated - vector of differentiated data to plot
            averaged - vector of averaged data to plot
            times - vector of times in ADC clock ticks to plot
            eventId - event ID being plotted
            threshold - threshold for averaged gradients used in peak finding
            plotMicrospill - if true, sets the start time of the event to zero and applies a modulo to the time setting it to the first microspill

        Variables
            tADC - ADC clock tick
            plotTimes - time to plot in [ns]
            plotADCs - ADCs to plot cast to double
            c - canvas used for plotting
            gADCs - ADCs graph
            gDeconvoluted - deconvoluted data graph
            gDifferentiated - differentiated data graph
            gAveraged - averaged data graph
            thresholdLine - averaged gradient threshold line for peak finding
            legend - plot legend
    */

    // Convert the time from ADC clock ticks to time in [ns]
    const double tADC = 3.125;
    std::vector<double> plotTimes;
    for (auto i : times)
        plotTimes.push_back((plotMicrospill ? normalizeTimeToFirstMicrospill(i) : i) * tADC);

    // Conver the ADCs to doubles for plotting
    std::vector<double> plotADCs;
    for (auto i : ADCs)
        plotADCs.push_back(static_cast<double>(i));

    // Generate the canvas and format it
    TCanvas* c = new TCanvas("c", "c", 1500, 1000);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.02);

    // Generate the ADCs plot and format the general plot
    TGraph* gADCs = new TGraph(ADCs.size(), plotTimes.data(), plotADCs.data());
    gADCs->Draw("APL");
    gADCs->SetLineColor(kRed);
    gADCs->SetLineWidth(2);
    gADCs->GetXaxis()->SetLimits(*std::min_element(plotTimes.begin(), plotTimes.end()), *std::max_element(plotTimes.begin(), plotTimes.end()));
    gADCs->SetTitle("Moving Window Deconvolution algorithm;Time [ns];ADC [arb. unit]");

    // Generate the deconvoluted plot
    TGraph* gDeconvoluted = new TGraph(ADCs.size(), plotTimes.data(), deconvoluted.data());
    gDeconvoluted->Draw("PL SAME");
    gDeconvoluted->SetLineWidth(2);

    // Generate the differentiated plot
    TGraph* gDifferentiated = new TGraph(ADCs.size(), plotTimes.data(), differentiated.data());z
    gDifferentiated->Draw("PL SAME");
    gDifferentiated->SetLineColor(kBlue);
    gDifferentiated->SetLineWidth(2);

    // Generate the averaged gradient plot
    TGraph* gAveraged = new TGraph(ADCs.size(), plotTimes.data(), averaged.data());
    gAveraged->Draw("PL SAME");
    gAveraged->SetLineColor(kGreen+2);
    gAveraged->SetLineWidth(2);

    // Generate the threshold boundary line
    TLine *thresholdLine = new TLine(*std::min_element(plotTimes.begin(), plotTimes.end()), threshold, *std::max_element(plotTimes.begin(), plotTimes.end()), threshold);
    thresholdLine->SetLineColor(kGreen + 2);
    thresholdLine->SetLineWidth(2);
    thresholdLine->Draw();

    // Generate the legend
    TLegend *legend = new TLegend(0.2, 0.2, 0.5, 0.4);
    legend->SetBorderSize(1);
    legend->SetFillColor(0);
    legend->SetTextSize(0.04);
    legend->AddEntry(gADCs, "ADCs", "l")->SetTextColor(kRed);
    legend->AddEntry(gDeconvoluted, "Deconvoluted", "l");
    legend->AddEntry(gDifferentiated, "Differentiated", "l")->SetTextColor(kBlue);
    legend->AddEntry(gAveraged, "Averaged", "l")->SetTextColor(kGreen + 2);
    legend->Draw();
    c->Update();

    // Save the plot
    std::string cFileName = "MWDAnalysis.Event" + std::to_string(eventId) + ".png";
    c->SaveAs(cFileName.c_str());

    // Cleanup
    gADCs->Delete();
    gDeconvoluted->Delete();
    gDifferentiated->Delete();
    gAveraged->Delete();
    thresholdLine->Delete();
    legend->Delete();
    c->Close();

    // Done
    return;
};

std::vector<unsigned int> makeUniqueEventIds(std::vector<unsigned int> &eventIds) {
    /*
        Description
            Generates a unique intersection vector of the event IDs retrieved from the input files

        Argument
            eventIds - as documented in function "plotZSAnalysis"

        Variables
            eventIdsSet - set of event IDs
            uniqueInputEventIds - vector of unique event IDs
    */
    std::cout << "\nMaking the event IDs unique\n" << std::endl;

    // Get a unqiue set of input event IDs
    std::set<unsigned int> eventIdsSet(eventIds.begin(), eventIds.end());

    // Convert the set to a vector
    std::vector<unsigned int> uniqueEventIds(eventIdsSet.begin(), eventIdsSet.end());

    // Done
    return uniqueEventIds;
};

void plotMWDAnalysis(const std::string fileName, const std::string treeName, const unsigned int eventID = 0, const double threshold = -100, const bool plotMicrospill = false) {
    /*
        Description
            Plots the MWD analysis result chain

        Arguments
            fileName - name of ROOT file generated with MWDTree
            treeName - name of ttree containing the data in MWDfileNameFileName
            eventID - event ID to plot. If 0, plots all the event IDs
            threshold - averaged gradient threshold, used in plot
            plotMicrospill - if true, resets all times to the first microspill

        Variables
            ADCs - vector of ADCs in the file
            deconvoluted - vector of deconvolution data in the file
            differentiated - vector of differentiation data in the file
            averaged - vector of averaged data in the file
            times - vector of times in the file
            eventIds - vector of event IDs in the file
            plotADCs - vector of ADCs used for plotting
            plotDeconvoluted - vector of deconvolution data used for plotting
            plotDifferentiated - vector of differentiation data used for plotting
            plotAveraged - vector of averaged data used for plotting
            plotTimes - vector of times used for plotting
            plotEventIds - vector of event IDs to be plotted
            nWaveformADCs - number of entries in the file
            plotEventId - event ID of plot being generated
            i - iterator variable
    */

    // Update global parameters
    SetErrorHandler(customErrorHandler);
    gROOT->SetBatch(kTRUE);

    // Construct the variables used to collect the data from the input file
    std::vector<int16_t> ADCs;
    std::vector<double> deconvoluted, differentiated, averaged;
    std::vector<uint32_t> times;
    std::vector<uint> eventIds;

    // Collect the data
    collectMWDData(fileName, treeName, ADCs, deconvoluted, differentiated, averaged, times, eventIds);

    // Construct the variables that will hold the filtered data for plotting
    std::vector<int16_t> plotADCs;
    std::vector<double> plotDeconvoluted, plotDifferentiated, plotAveraged;
    std::vector<uint32_t> plotTimes;
    std::vector<uint> plotEventIds = makeUniqueEventIds(eventIds);

    // Update the event IDs that will be plotted if selected to
    if (eventID != 0) {
        // Validate that the chosen event ID is found in the data
        if (std::find(plotEventIds.begin(), plotEventIds.end(), eventID) == plotEventIds.end())
            Fatal("plotMWDAnalysis", "Requested event ID is not found in the data, exiting");

        // Update the event ID to plot
        plotEventIds.clear();
        plotEventId.emplace_back(eventID);
    };

    // Select the data by event ID, plot it
    int nWaveformADCs = ADCs.size();
    for (unsigned int plotEventId : plotEventIds) {
        std::cout << "Generating plot for event with ID " << plotEventId << std::endl;

        // Select the data to plot
        for (int i = 0; i < nWaveformADCs; i++) {
            if (eventIds[i] == plotEventId) {
                plotADCs.push_back(ADCs[i]);
                plotDeconvoluted.push_back(deconvoluted[i]);
                plotDifferentiated.push_back(differentiated[i]);
                plotAveraged.push_back(averaged[i]);
                plotTimes.push_back(times[i]);
            };
        };

        // Plot the selected data
        plot(plotADCs, plotDeconvoluted, plotDifferentiated, plotAveraged, plotTimes, plotEventId, threshold, plotMicrospill);

        // Prepare the variables for the next event ID
        plotADCs.clear();
        plotDeconvoluted.clear();
        plotDifferentiated.clear();
        plotAveraged.clear();
        plotTimes.clear();
    };

    // Done
    return;
};
