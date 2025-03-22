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

void collectMWDData(const std::string fileName, const std::string treeName, std::vector<int16_t> &ADCs, std::vector<double> &deconvoluted, std::vector<double> &differentiated, std::vector<double> &averaged, std::vector<double> &times, std::vector<uint> &eventIds, std::vector<uint> &waveformIds) {
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
            waveformIds - as documented in function "plotMWDResults"

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
            dataWaveformId - waveform ID from file
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
    uint dataEventId, dataWaveformId;
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
    if (std::find(branchNames.begin(), branchNames.end(), "waveformID") != branchNames.end())
        tree->SetBranchAddress("waveformID", &dataWaveformId);
    else
        Fatal("collectData", "Requested branch 'waveformID' does not exist in the file.");

    // Get the number of entries
    int entries = tree->GetEntries();

    // Construct the variable casting the data type to the correct variable
    const double tADC = 3.125;

    // Collect the data
    for (int i = 0; i < entries; i++) {
        // Get the corresponding entry
        tree->GetEntry(i);

        // Collect the data
        ADCs.push_back(dataADC);
        deconvoluted.push_back(dataDeconvoluted);
        differentiated.push_back(dataDifferentiated);
        averaged.push_back(dataAveraged);
        times.push_back(dataTime * tADC);
        eventIds.push_back(dataEventId);
        waveformIds.push_back(dataWaveformId);
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

void plot(std::vector<int16_t> &ADCs, std::vector<double> &deconvoluted, std::vector<double> &differentiated, std::vector<double> &averaged, std::vector<double> &times, const unsigned int eventId, const unsigned int waveformId, const double threshold, const bool highResolution) {
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
            waveformId - waveform ID being plotted
            threshold - threshold for averaged gradients used in peak finding
            highResolution - controls the resolution of the generated plots

        Variables
            plotADCs - ADCs to plot cast to double
            px - number of x pixels for canvases
            py - number of y pixels for canvases
            canvasLeftMargin - padding for plots on cavases from left border
            canvasRightMargin - padding for plots on canvases from right border
            plotFileNamePrefix - prefix for all plot file names
            plotFileNameSuffix - suffix for all plot file names
            plotTitleSuffix - suffix for all plot titles, defining the x and y axis labels
            cCombined - combined canvas
            legend - legend for combined plot
            gADCs - ADCs graph
            plotTitle - title of ploted including x and y axis labels
            gDeconvoluted - deconvoluted data graph
            gDifferentiated - differentiated data graph
            gAveraged - averaged data graph
            thresholdLine - averaged gradient threshold line for peak finding. If the plot is not in range, it stays as a nullptr
            legend - plot legend
            cADC - canvas for ADC only plot
            cDeconv - canvas for deconvoluted only plot
            cDiff - canvas for differentiated only plot
            cAvg - canvas for averaged only plot
            canvases - vector of canvases to use with plots
            graphs - vector of all plots
            lowResAxisTextSize - updated text size for low resolution plots
            canvas - iterator for canvases
            graph - iterator for graphs
            fileName - name of plot file
    */

    // Convert the ADCs to doubles for plotting
    std::vector<double> plotADCs;
    for (auto i : ADCs)
        plotADCs.push_back(static_cast<double>(i));

    // Generate the canvas properties and apply them
    const int px = highResolution ? 1500 : 800;
    const int py = highResolution ? 1000 : 500;
    const double canvasLeftMargin = 0.14, canvasRightMargin = 0.05;
    gStyle->SetPadLeftMargin(canvasLeftMargin);
    gStyle->SetPadRightMargin(canvasRightMargin);

    // Define plot IO variables
    const std::string plotFileNamePrefix = "MWD.Analysis.";
    const std::string plotFileNameSuffix = ".Event" + std::to_string(eventId) + ".Waveform" + std::to_string(waveformId) + "." + (highResolution ? "high" : "low") + ".png";
    const std::string plotTitleSuffix = ";Time [ns];ADC [arb. unit]";

    // Generate the canvas for the combined plot
    TCanvas* cCombined = new TCanvas("cCombined", "cCombined", px, py);

    // Generate the legend
    TLegend *legend = new TLegend(0.6, 0.5, 0.85, 0.75);
    legend->SetBorderSize(1);
    legend->SetFillColor(0);
    legend->SetTextSize(0.04);

    // Generate the ADCs plot and format it
    TGraph* gADCs = new TGraph(ADCs.size(), times.data(), plotADCs.data());
    gADCs->Draw("APL");
    gADCs->SetLineColor(kRed);
    gADCs->SetLineWidth(2);
    gADCs->GetXaxis()->SetLimits(*std::min_element(times.begin(), times.end()), *std::max_element(times.begin(), times.end()));
    std::string plotTitle = "Moving Window Deconvolution algorithm" + plotTitleSuffix;
    gADCs->SetTitle(plotTitle.c_str());
    legend->AddEntry(gADCs, "ADCs", "l")->SetTextColor(kRed);

    // Generate the deconvoluted plot
    TGraph* gDeconvoluted = new TGraph(ADCs.size(), times.data(), deconvoluted.data());
    gDeconvoluted->Draw("PL SAME");
    gDeconvoluted->SetLineWidth(2);
    legend->AddEntry(gDeconvoluted, "Deconvoluted", "l");

    // Generate the differentiated plot
    TGraph* gDifferentiated = new TGraph(ADCs.size(), times.data(), differentiated.data());
    gDifferentiated->Draw("PL SAME");
    gDifferentiated->SetLineColor(kBlue);
    gDifferentiated->SetLineWidth(2);
    legend->AddEntry(gDifferentiated, "Differentiated", "l")->SetTextColor(kBlue);

    // Generate the averaged gradient plot
    TGraph* gAveraged = new TGraph(ADCs.size(), times.data(), averaged.data());
    gAveraged->Draw("PL SAME");
    gAveraged->SetLineColor(kGreen+2);
    gAveraged->SetLineWidth(2);
    legend->AddEntry(gAveraged, "Averaged", "l")->SetTextColor(kGreen + 2);

    // Generate the threshold boundary line
    TLine *thresholdLine = nullptr;
    if (threshold < *std::max_element(ADCs.begin(), ADCs.end())) {
        thresholdLine = new TLine(gADCs->GetXaxis()->GetXmin(), threshold, gADCs->GetXaxis()->GetXmax(), threshold);
        thresholdLine->SetLineColor(kGreen + 2);
        thresholdLine->SetLineWidth(2);
        thresholdLine->SetLineStyle(2); // Make the line dashed
        thresholdLine->Draw();
        legend->AddEntry(thresholdLine, "Threshold", "l")->SetTextColor(kGreen + 2);
    };

    // Draw the legend, update the canvas
    legend->Draw();
    cCombined->Update();


    // Generate the canvases for the remining plots
    TCanvas* cADC = new TCanvas("cADC", "cADC", px, py);
    TCanvas* cDeconv = new TCanvas("cDeconv", "cDeconv", px, py);
    TCanvas* cDiff = new TCanvas("cDiff", "cDiff", px, py);
    TCanvas* cAvg = new TCanvas("cAvg", "cAvg", px, py);

    // Format the canvas text sizes for low resolution plots
    const double lowResAxisTextSize = 0.04;
    std::vector<TCanvas*> canvases = {cCombined, cADC, cDeconv, cDiff, cAvg};
    std::vector<TGraph*> graphs = {gADCs, gDeconvoluted, gDifferentiated, gAveraged};
    gStyle->SetTitleFontSize(lowResAxisTextSize);
    if (!highResolution) {
        for (TGraph* graph : graphs) {
            graph->GetXaxis()->SetLabelSize(lowResAxisTextSize);
            graph->GetXaxis()->SetTitleSize(lowResAxisTextSize);
            graph->GetYaxis()->SetLabelSize(lowResAxisTextSize);
            graph->GetYaxis()->SetTitleSize(lowResAxisTextSize);
        };
        for (TCanvas* canvas : canvases) {
            canvas->SetBottomMargin(0.14);
            canvas->SetLeftMargin(0.175);
            canvas->Update();
        };
    };


    // Save the combined plot
    cCombined->cd();
    std::string fileName = plotFileNamePrefix + "comb" + plotFileNameSuffix;
    cCombined->SaveAs(fileName.c_str());

    // Generate the ADC plot
    cADC->cd();
    plotTitle = "MWD: Input ADCs" + plotTitleSuffix;
    gADCs->SetTitle(plotTitle.c_str());
    gADCs->Draw("APL");
    fileName = plotFileNamePrefix + "ADC" + plotFileNameSuffix;
    cADC->SaveAs(fileName.c_str());

    // Generate the deconvoluted plot
    cDeconv->cd();
    plotTitle = "MWD: Deconvoluted ADCs" + plotTitleSuffix;
    gDeconvoluted->SetTitle(plotTitle.c_str());
    gDeconvoluted->Draw("APL");
    fileName = plotFileNamePrefix + "deconvoluted" + plotFileNameSuffix;
    cDeconv->SaveAs(fileName.c_str());

    // Generate the deconvoluted plot
    cDiff->cd();
    plotTitle = "MWD: Differentiated deconvoluted ADCs" + plotTitleSuffix;
    gDifferentiated->SetTitle(plotTitle.c_str());
    gDifferentiated->Draw("APL");
    fileName = plotFileNamePrefix + "differentiated" + plotFileNameSuffix;
    cDiff->SaveAs(fileName.c_str());

    // Generate the averaged plot
    cAvg->cd();
    plotTitle = "MWD: Averaged differentiated deconvoluted ADCs" + plotTitleSuffix;
    gAveraged->SetTitle(plotTitle.c_str());
    gAveraged->Draw("APL");
    if (thresholdLine != nullptr)
        thresholdLine->Draw();
    fileName = plotFileNamePrefix + "averaged" + plotFileNameSuffix;
    cAvg->SaveAs(fileName.c_str());

    // Cleanup
    gADCs->Delete();
    gDeconvoluted->Delete();
    gDifferentiated->Delete();
    gAveraged->Delete();
    legend->Delete();
    cCombined->Close();
    cADC->Close();
    cDeconv->Close();
    cDiff->Close();
    cAvg->Close();

    // Done
    return;
};

std::vector<unsigned int> makeUnique(std::vector<unsigned int> &v) {
    /*
        Description
            Generates a unique vector from a vector of repeated numbers

        Argument
            v - vector to generate unique entries for

        Variables
            vSet - set of entries in v
            vUnique - vector of unique entries from vector v
    */

    // Get a unqiue set of input event IDs
    std::set<unsigned int> vSet(v.begin(), v.end());

    // Convert the set to a vector
    std::vector<unsigned int> vUnique(vSet.begin(), vSet.end());

    // Done
    return vUnique;
};

void plotMWDAnalysis(const std::string fileName, const std::string treeName, const unsigned int eventID = 0, const unsigned int waveformID = 0, const double threshold = -100) {
    // TODO - add tMin and tMax
    /*
        Description
            Plots the MWD analysis result chain, generating files as
                MWD.Analysis.<type>.Event<EventID>.Waveform<WaveformID>.png
            such that
                <type> is one of "deconv", "diff", "avg", "comb" for deconvoluted data, differentiated data, averaged data, and all data
                <EventID> is the art event ID
                <waveformID> is the counter of the waveforms generated for a specific event

        Arguments
            fileName - name of ROOT file generated with MWDTree
            treeName - name of ttree containing the data in MWDfileNameFileName
            eventID - event ID to plot. If 0, plots all the event IDs
            waveformID - waveform ID to plot. If 0, plots all the waveform IDs
            threshold - averaged gradient threshold, used in plot

        Variables
            ADCs - vector of ADCs in the file
            deconvoluted - vector of deconvolution data in the file
            differentiated - vector of differentiation data in the file
            averaged - vector of averaged data in the file
            times - vector of times in the file
            eventIds - vector of event IDs in the file
            waveformIds - vector of waveform IDs in the file
            plotADCs - vector of ADCs used for plotting
            plotDeconvoluted - vector of deconvolution data used for plotting
            plotDifferentiated - vector of differentiation data used for plotting
            plotAveraged - vector of averaged data used for plotting
            plotTimes - vector of times used for plotting
            plotEventIds - vector of event IDs to be plotted
            plotWaveformIds - vector of waveform IDs to be plotted
            nWaveformADCs - number of entries in the file
            plotEventId - event ID of plot being generated
            plotWaveformId - event ID of plot being generated
            i - iterator variable
    */

    // Update global parameters
    SetErrorHandler(customErrorHandler);
    gROOT->SetBatch(kTRUE);

    // Construct the variables used to collect the data from the input file
    std::vector<int16_t> ADCs;
    std::vector<double> deconvoluted, differentiated, averaged, times;
    std::vector<uint> eventIds, waveformIds;

    // Collect the data
    collectMWDData(fileName, treeName, ADCs, deconvoluted, differentiated, averaged, times, eventIds, waveformIds);

    // Construct the variables that will hold the filtered data for plotting
    std::vector<int16_t> plotADCs;
    std::vector<double> plotDeconvoluted, plotDifferentiated, plotAveraged, plotTimes;
    std::vector<uint> plotEventIds = makeUnique(eventIds), plotWaveformIds = makeUnique(waveformIds);

    // Update the event IDs and waveform IDs that will be plotted if selected to
    if (eventID != 0) {
        // Validate that the chosen event ID is found in the data
        if (std::find(plotEventIds.begin(), plotEventIds.end(), eventID) == plotEventIds.end())
            Fatal("plotMWDAnalysis", "Requested event ID is not found in the data, exiting");

        // Update the event ID to plot
        plotEventIds.clear();
        plotEventIds.emplace_back(eventID);
    };
    if (waveformID != 0) {
        // Validate that the chosen event ID is found in the data
        if (std::find(plotWaveformIds.begin(), plotWaveformIds.end(), waveformID) == plotWaveformIds.end())
            Fatal("plotMWDAnalysis", "Requested event ID is not found in the data, exiting");

        // Update the event ID to plot
        plotWaveformIds.clear();
        plotWaveformIds.emplace_back(waveformID);
    };

    // Construct iterator for different types of plot being generated
    std::vector<bool> boolValues = {true, false};

    // Select the data by event ID, plot it
    int nWaveformADCs = ADCs.size();
    for (unsigned int plotEventId : plotEventIds) {
        for (unsigned int plotWaveformId : plotWaveformIds) {
            // Select the data to plot
            for (int i = 0; i < nWaveformADCs; i++) {
                if (eventIds[i] == plotEventId && waveformIds[i] == plotWaveformId) {
                    plotADCs.push_back(ADCs[i]);
                    plotDeconvoluted.push_back(deconvoluted[i]);
                    plotDifferentiated.push_back(differentiated[i]);
                    plotAveraged.push_back(averaged[i]);
                    plotTimes.push_back(times[i]);
                };
            };

            // Plot the selected data
            std::cout << "Generating plot for event with ID " << plotEventId << " and waveform ID " << plotWaveformId << std::endl;
            for (bool highResolution : boolValues)
                plot(plotADCs, plotDeconvoluted, plotDifferentiated, plotAveraged, plotTimes, plotEventId, plotWaveformId, threshold, highResolution);

            // Prepare the variables for the next event ID
            plotADCs.clear();
            plotDeconvoluted.clear();
            plotDifferentiated.clear();
            plotAveraged.clear();
            plotTimes.clear();
        };
    };

    // Generate a single combined plot at the end if relevant
    if (eventID == 0 && waveformID == 0) {
        for (bool highResolution : boolValues)
            plot(ADCs, deconvoluted, differentiated, averaged, times, eventID, waveformID, threshold, highResolution);
    };

    // Done
    return;
};
