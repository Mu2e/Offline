// Generates a plot of previous STM energy fluxes
// See plotEnergyFlux.sh for usage examples
// Original author: Pawel Plesniak

void customErrorHandler(int level, Bool_t abort, const char* location, const char* message) {
    /*
        Description
            Define a custom error handler that won't print the stack trace but will print an error message and exit.
    */
    std::cerr << message << std::endl;
    if (level > kInfo)
        exit(1);
};

void getHPGeFluxes(std::vector<double> &flux, std::vector<double> &fluxUncertainty) {
    /*
        Description
            Stores all the HPGe flux data. Variables that store the fluxes are named as
                flux_HPGe_DOC-DB-ID_LocationInSource - defines the energy flux from the indicated source
                uncert_flux_HPGe_DOC-DB-ID_LocationInSource - defines the energy flux uncertainty from the indicated source

        Arguments
            flux - vector of the flux data
            fluxUncertainty - vector of the associated flux uncertainty
    */
    const double flux_HPGe_DocDB6453v4_Tab6     = 0.0463;
    const double flux_HPGe_DocDB25188v1_Fig6    = 0.3869, uncert_flux_HPGe_DocDB25188v1_Fig6 = 0.0511;
    const double flux_HPGe_DocDB34596v7_Page12  = 0.1254;
    const double flux_HPGe_DocDB36803v6_TabI    = 0.12;
    const double flux_HPGe_DocDB36575v1_Tab5    = 0.0978;
    const double flux_HPGe_DocDB37535v1_Page5   = 0.1168;
    const double flux_HPGe_DocDB42969v1_Tab3    = 0.0894, uncert_flux_HPGe_DocDB42969v1_Tab3 = 0.0013;
    const double flux_HPGe_DocDB49109_Tab3      = 0.3181, uncert_flux_HPGe_DocDB49109_Tab3 = 0.0016;
    const double flux_HPGe_DocDB45229v1_Page9   = 0.0946;
    const double flux_HPGe_DocDB49769_Tab1      = 0.0179, uncert_flux_HPGe_DocDB49769_Tab1 = 0.001;

    std::vector<double> fluxData = {
        flux_HPGe_DocDB6453v4_Tab6,
        flux_HPGe_DocDB25188v1_Fig6,
        flux_HPGe_DocDB34596v7_Page12,
        flux_HPGe_DocDB36803v6_TabI,
        flux_HPGe_DocDB36575v1_Tab5,
        flux_HPGe_DocDB37535v1_Page5,
        flux_HPGe_DocDB42969v1_Tab3,
        flux_HPGe_DocDB49109_Tab3,
        flux_HPGe_DocDB45229v1_Page9,
        flux_HPGe_DocDB49769_Tab1};
    std::vector<double> fluxUncertaintyData = {
        0.0,
        uncert_flux_HPGe_DocDB25188v1_Fig6,
        0.0,
        0.0,
        0.0,
        0.0,
        uncert_flux_HPGe_DocDB42969v1_Tab3,
        uncert_flux_HPGe_DocDB49109_Tab3,
        0.0,
        uncert_flux_HPGe_DocDB49769_Tab1};

    flux.insert(flux.end(), fluxData.begin(), fluxData.end());
    fluxUncertainty.insert(fluxUncertainty.end(), fluxUncertaintyData.begin(), fluxUncertaintyData.end());
    return;
};

void getLaBrFluxes(std::vector<double> &flux, std::vector<double> &fluxUncertainty) {
    /*
        Description
            NOTE - the data has not been collated from the sources
            Stores all the LaBr flux data. Variables that store the fluxes are named as
                flux_LaBr_DOC-DB-ID_LocationInSource - defines the energy flux from the indicated source
                uncert_flux_LaBr_DOC-DB-ID_LocationInSource - defines the energy flux uncertainty from the indicated source

        Arguments
            flux - vector of the flux data
            fluxUncertainty - vector of the associated flux uncertainty
    */
    Fatal("getLaBrFluxes", "Required data has not been collated from the appropriate sources, exiting.");
};

void plotEnergyFlux(const std::string detector, const bool highResolution, std::vector<double> referenceTableSources) {
    /*
        Description
            Plot the energy fluxes from previous STM studies

        Arguments
            detector - either "HPGe" or "LaBr"
            highResolution - determines whether the plot is generated in high resolution or in low resolution
            referenceTableSources - a vector used to index what entries in a linked publication correspond to the corresponding entry in the plot

        // Variables
            allowedDetectors - list of all the allowed STM detectors
            flux - vector of the fluxes
            fluxUncertainty - vector of the uncertainties corresponding to the entries in `flux`
            fluxMin - vector of the minimum bounds of the fluxes incorporating the associated uncertainties
            fluxMax - vector of the maximum bounds of the fluxes incorporating the associated uncertainties
            xMin - minimum x value to use for the plot
            xMax - maximum x value to use for the plot
            yMin - minimum y value to use for the plot
            yMax - maximum y value to use for the plot
            eyV - vector of uncertainties associated with the entries (all zeros)
            opcacityErrors - opacity of all the flux errors
            opacityLimit - opacity of the flux limit
            lineWidth - width of the saturation and average lines
            px - number of x pixels for the TCanvas
            py - number of y pixels for the TCanvas
            c1 - the TCanvas
            x - array of all the flux data
            y - array of all the reference table sources
            ex - array of all the flux data uncertainty
            ey - array of all the reference table source uncertainty (all zeros)
            graph - pointer to the TGraph for plotting the fluxes
            errorX - defines the x coordinates of the uncertainty shading area, going counter-clockwise from the bottom left
            errorY - as per `errorX`, but for the Y coordinates
            errorBounds - vector of TGraphs for each flux with a non-zero uncertainty
            xAvg - average flux value
            exAvg - average uncertainty value
            avgErrorBound - TGraph for the averaged flux uncertainty
            avgLine - TLine for the averaged flux
            avgLabel - TText for the averaged flux
            saturationBound - TGraph for the saturation region
            saturationLine - TLine for the saturation boundary
            saturationLabel - TText for the saturation region
            outputFileName - name of the output file
    */

    // Update global parameters
    SetErrorHandler(customErrorHandler);

    // Validate the data
    std::vector<std::string> allowedDetectors = {"HPGe", "LaBr"};
    if (std::find(allowedDetectors.begin(), allowedDetectors.end(), detector) == allowedDetectors.end())
        Fatal("plotFlux", "detector is not one of the allowed types, it must be either 'HPGe' or 'LaBr'");

    // Collect the relevant data
    std::vector<double> flux, fluxUncertainty;
    if (detector == "HPGe")
        getHPGeFluxes(flux, fluxUncertainty);
    else
        getLaBrFluxes(flux, fluxUncertainty);
    if (flux.size() != fluxUncertainty.size())
        Fatal("plotEnergyFlux", "Collected number of data elements is not equal, validate the correct number of entries have been collected");
    const int N = flux.size();

    // Collate the plotting data and determine the plotting range in x and y
    std::vector<double> fluxMin(N), fluxMax(N);
    std::transform(flux.begin(), flux.end(), fluxUncertainty.begin(), fluxMin.begin(), std::minus<double>());
    std::transform(flux.begin(), flux.end(), fluxUncertainty.begin(), fluxMax.begin(), std::plus<double>());
    double saturationLimit = (detector == "HPGe") ? 0.37 : 4; // based on docDB 36803v6
    const double xMin = *std::min_element(fluxMin.begin(), fluxMin.end()) * 0.85; // scaling to include buffer space in plot
    const double xMax = std::max(*std::max_element(fluxMax.begin(), fluxMax.end()), saturationLimit) * 1.15; // scaling to include buffer space in plot
    const double yMax = *std::max_element(referenceTableSources.begin(), referenceTableSources.end()) + 1; // + 1 to include buffer space in plot
    const double yMin = *std::min_element(referenceTableSources.begin(), referenceTableSources.end()) - 1; // - 1 to include buffer space in plot
    std::vector<double> eyV(referenceTableSources.size(), 0.0);

    // Set up the plotting variables
    const double opacityErrors = 0.15, opacityLimit = 0.4;
    const int lineWidth = 3;
    const int px = (highResolution ? 1500 : 750), py = (highResolution ? 1000 : 500);

    // Set up the TCanvas
    TCanvas *c1 = new TCanvas("Fluxes", "Fluxes", px, py);
    gStyle->SetOptStat(0); // Remove the stat box

    // Convert the vectors to arrays for ROOT
    double* x  = flux.data();
    double* y  = referenceTableSources.data();
    double* ex = fluxUncertainty.data();
    double* ey = eyV.data();

    // Set up the TGraph
    TGraphErrors *graph = new TGraphErrors(N, x, y, ex, ey);
    graph->SetTitle((std::string("Energy fluxes at ") + detector + " detector;Energy flux [TeV s^{-1}];Row number in Table 5.1").c_str());
    graph->GetYaxis()->SetRangeUser(yMin, yMax);
    graph->SetMarkerStyle(5);
    graph->SetMarkerSize(highResolution ? 3 : 1);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    // Set up the TGraphs to shade the areas of uncertainty
    std::vector<double> errorX, errorY = {yMin, yMin, yMax, yMax};
    std::vector<TGraph*> errorBounds;
    for (int i = 0; i < N; i++) {
        if (fluxUncertainty[i] > std::numeric_limits<double>::epsilon()) {
            errorX = {fluxMin[i], fluxMax[i], fluxMax[i], fluxMin[i]};
            errorBounds.emplace_back(new TGraph(4, errorX.data(), errorY.data()));
        };
    };

    // Set up the average flux measurement from literature, including the associated TGraph, TLine, and TText
    double xAvg = std::accumulate(flux.begin(), flux.end(), 0.0) / N;
    double exAvg = std::accumulate(fluxUncertainty.begin(), fluxUncertainty.end(), 0.0) / N;
    errorX = {xAvg - exAvg, xAvg + exAvg, xAvg + exAvg, xAvg - exAvg};
    TGraph *avgErrorBound = new TGraph(4, errorX.data(), errorY.data());
    avgErrorBound->SetFillColorAlpha(kBlue, opacityErrors);
    TLine *avgLine = new TLine(xAvg, yMin, xAvg, yMax);
    avgLine->SetLineColor(kBlue);
    avgLine->SetLineWidth(lineWidth);
    TText *avgLabel = new TText((xAvg + exAvg) * 1.3, (yMax - yMin)/2, "Average");
    avgLabel->SetTextColor(kBlue);
    avgLabel->SetTextAlign(22);

    // Set up the detector flux limit TGraph, TLine, and TText
    errorX = {saturationLimit, xMax, xMax, saturationLimit};
    TGraph *saturationBound = new TGraph(4, errorX.data(), errorY.data());
    saturationBound->SetFillColorAlpha(kRed, opacityLimit);
    TLine *saturationLine = new TLine(saturationLimit, yMin, saturationLimit, yMax);
    saturationLine->SetLineColor(kRed);
    saturationLine->SetLineWidth(lineWidth);
    TText *saturationLabel = new TText((xMax + saturationLimit - 0.02)/2, (yMax - yMin)/2, "Saturated");
    saturationLabel->SetTextColor(kRed);
    saturationLabel->SetTextAlign(22);

    // Gather all the TGraphs, TLines and TTexts
    graph->Draw("AP");
    for (TGraph* errorBound : errorBounds) {
        errorBound->SetFillColorAlpha(kBlue, opacityErrors);
        errorBound->SetLineColor(0);
        errorBound->Draw("F SAME");
    };
    saturationBound->Draw("F SAME");
    saturationLine->Draw();
    saturationLabel->Draw();
    avgErrorBound->Draw("F SAME");
    avgLine->Draw();
    avgLabel->Draw();

    // Save the plot
    std::string outputFileName = detector + "EnergyFlux.png";
    c1->SaveAs(outputFileName.c_str());
    return;
};
