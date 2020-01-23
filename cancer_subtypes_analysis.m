function cancer_subtypes_analysis(filename_km, outputDir, numericData, patientsNames, geneNames,...
    txtDataKM, patientsNamesKM, mainGeneIdx, colCancerSubtypes, use_cbioportal_not_lab_std, geneName)
fileOutput = fullfile(outputDir, 'cancer_subtypes_historgram_70_30.xlsx');
getHistogramSubtypes(filename_km, fileOutput, numericData, patientsNames, geneNames, txtDataKM,...
    patientsNamesKM, 0.3, 0.7, mainGeneIdx, colCancerSubtypes, use_cbioportal_not_lab_std, geneName);

fileOutput = fullfile(outputDir, 'cancer_subtypes_historgram_80_20.xlsx');
getHistogramSubtypes(filename_km, fileOutput, numericData, patientsNames, geneNames, txtDataKM,...
    patientsNamesKM, 0.2, 0.8, mainGeneIdx, colCancerSubtypes, use_cbioportal_not_lab_std, geneName);
end
    
    %% histogram function for breast cancer subtypes:

function getHistogramSubtypes(fileInput, fileOutput, numericData, patientsNames, geneNames,...
    txtDataKM, patientsNamesKM ,topPerc, lowPerc, mainGeneIdx, colCancerSubtypes, ...
    use_cbioportal_not_lab_std, geneName1)
% Getting all the different breast cancer categories
if use_cbioportal_not_lab_std
    subtypeData = txtDataKM(6:end, colCancerSubtypes);
else
    subtypeData = txtDataKM(4:end, colCancerSubtypes);
end
subtypeCategorical = categorical(subtypeData);
xLabels = unique(subtypeData);
xLabels = xLabels(xLabels ~= "");

% For the excel file
tableTitle = {'geneName', ''};
tableTitle = [tableTitle, xLabels'];
tmp = {fileInput};
tmp{end + length(xLabels) + 1} = '';
tableTitle = [tmp; tableTitle];
table = {};

% Splitting main gene into high expression vs. low expression
[~, ~, patientsNamesTopGene1, patientsNamesLowGene1]...
        = split_gene_data_by_percentage(numericData, patientsNames, topPerc, lowPerc, mainGeneIdx);

% Going over all the genes

    for geneIdx = 1:length(geneNames)
    geneName2 = geneNames{geneIdx};
    % Splitting the patients into: Patients with high expression of current
    % gene vs. patients with low expression.
        [~, ~, patientsNamesTopGene2, patientsNamesLowGene2]...
        = split_gene_data_by_percentage(numericData, patientsNames, topPerc, lowPerc, geneIdx);
    [~, idxTop] = intersect(patientsNamesKM, patientsNamesTopGene2);
    [~, idxLow] = intersect(patientsNamesKM, patientsNamesLowGene2);
    % Calculating the percentages of the no. of patients in each cancer subtype category. Comparing high/low
    % expression in current gene.
    percentageHigh = getCategoriesPercentages(subtypeCategorical, idxTop);
    percentageLow = getCategoriesPercentages(subtypeCategorical, idxLow);
    histDataOfGene = [percentageHigh, percentageLow];
    % Plotting the histogram
    plotTitle = sprintf('Histogram of cancer subtypes %d%% vs. %d%%- %s',...
        lowPerc*100, topPerc*100, geneName2);
    figure('Name', plotTitle, 'visible', 'off');
    bar(histDataOfGene);
    clear title xlabel ylabel; 
    xticklabels(xLabels);
    xtickangle(315);
    legend('High','Low',...
        'Location', 'northeastoutside');
    xlabel('Cancer subtypes');
    ylabel('Percentage of patients in each high/low category');
    title(plotTitle);
        for i = 1:length(xLabels)
        text(i-0.3,percentageHigh(i),num2str(percentageHigh(i),2),'vert','bottom','horiz','left');
        text(i+0.3,percentageLow(i),num2str(percentageLow(i),2),'vert','bottom','horiz','right');
        end
    if geneIdx ~= mainGeneIdx
    % Calculating the percentages of the no. of patients in each cancer subtype category. 
    % Comparing high-high, high-low, low-high, low-low expressions of main gene and current gene, respectively.
    patientsNamesHH = intersect(patientsNamesTopGene1, patientsNamesTopGene2);
    patientsNamesHL = intersect(patientsNamesTopGene1, patientsNamesLowGene2);
    patientsNamesLH = intersect(patientsNamesLowGene1, patientsNamesTopGene2);
    patientsNamesLL = intersect(patientsNamesLowGene1, patientsNamesLowGene2);
    [~, idxHH] = intersect(patientsNamesKM, patientsNamesHH);
    [~, idxHL] = intersect(patientsNamesKM, patientsNamesHL);
    [~, idxLH] = intersect(patientsNamesKM, patientsNamesLH);
    [~, idxLL] = intersect(patientsNamesKM, patientsNamesLL);
    percentageHH = getCategoriesPercentages(subtypeCategorical, idxHH);
    percentageHL = getCategoriesPercentages(subtypeCategorical, idxHL);
    percentageLH = getCategoriesPercentages(subtypeCategorical, idxLH);
    percentageLL = getCategoriesPercentages(subtypeCategorical, idxLL);
    histDataOfGene1AndGene2 = [percentageHH, percentageHL, percentageLH, percentageLL];
    % Plotting the histogram 
    plotTitle = sprintf('Histogram of cancer subtypes %d%% vs. %d%%- %s (G2) with %s (G1)',...
        lowPerc*100, topPerc*100, geneName2, geneName1);
    figure('Name', plotTitle, 'visible', 'off');
    bar(histDataOfGene1AndGene2);
    clear title xlabel ylabel; 
    xticklabels(xLabels);
    xtickangle(315);
    legend('High G1-High G2','High G1-Low G2','Low G1-High G2','Low G1-Low G2',...
        'Location', 'northeastoutside');
    xlabel('Cancer subtypes');
    ylabel('Percentage of patients in each category');
    title(plotTitle);
    for i = 1:length(xLabels)
        text(i-0.4,percentageHH(i),num2str(percentageHH(i),2),'vert','bottom','horiz','left', 'FontSize', 7);
        text(i-0.2,percentageHL(i),num2str(percentageHL(i),2),'vert','bottom','horiz','left', 'FontSize', 7);
        text(i+0.2,percentageLH(i),num2str(percentageLH(i),2),'vert','bottom','horiz','right', 'FontSize', 7);
        text(i+0.4,percentageLL(i),num2str(percentageLL(i),2),'vert','bottom','horiz','right', 'FontSize', 7);
    end

    % Exporting an excel file with all the percentages data and p-values
    pvaluesGene = pvalueForEachSubtype(patientsNames, subtypeCategorical, {percentageHigh,percentageLow});
    pvaluesGene1AndGene2 =  pvalueForEachSubtype(patientsNames, subtypeCategorical,...
        {percentageHH,percentageHL, percentageLH, percentageLL});
    table = [table; addGeneDataInTable(geneName2, percentageHigh, percentageLow', pvaluesGene,...
    percentageHH, percentageHL, percentageLH, percentageLL, pvaluesGene1AndGene2')];
        
    else % Add main gene data into the table. 
        % See 'addGeneDataInTable' - very similar, but specific for MET.
        mainGeneTable1 = {geneName1,'High';'','Low';'', 'p-value High vs Low'};
        mainGeneTable2 = doubleArrayToCellArray(percentageHigh);
        mainGeneTable2 = [mainGeneTable2; doubleArrayToCellArray(percentageLow)];
        pvaluesGene = pvalueForEachSubtype(patientsNames, subtypeCategorical, {percentageHigh,percentageLow});
        mainGeneTable2 = [mainGeneTable2; doubleArrayToCellArray(pvaluesGene)];
        mainGeneTable = [mainGeneTable1, mainGeneTable2];
        table = [table; mainGeneTable];
    end
    end
    table = [tableTitle;table];
    xlswrite(fileOutput, table);
end


function percentagesArray = getCategoriesPercentages(subtypeCategorical, idxOfPatients)
numOfPatients = countcats(subtypeCategorical(idxOfPatients)); % ignores undefined.
sumOfPatients = sum(numOfPatients);
percentagesArray = (numOfPatients/sumOfPatients)*100;
end

function table = addGeneDataInTable(geneName2, percentagesHigh, percentagesLow, pvaluesGene,...
    percentagesHH, percentagesHL, percentagesLH, percentagesLL, pvaluesGene1AndGene2)
table1 = {geneName2,'High';'','Low';'', 'p-value High vs Low'; '', 'High-High';...
    '', 'High-Low'; '', 'Low-High'; '', 'Low-Low'; '', 'p-value HH, HL, LH, LL'};
table2 = {};
table2 = [table2; doubleArrayToCellArray(percentagesHigh)];
table2 = [table2; doubleArrayToCellArray(percentagesLow)];
table2 = [table2; doubleArrayToCellArray(pvaluesGene)];
table2 = [table2; doubleArrayToCellArray(percentagesHH)];
table2 = [table2; doubleArrayToCellArray(percentagesHL)];
table2 = [table2; doubleArrayToCellArray(percentagesLH)];
table2 = [table2; doubleArrayToCellArray(percentagesLL)];
table2 = [table2; doubleArrayToCellArray(pvaluesGene1AndGene2)];
table = [table1,table2];
end

function cellArray = doubleArrayToCellArray(columnDoubleArray)
rowDoubleArray = columnDoubleArray';
tmp = sprintf('%f*', rowDoubleArray);
tmp(end) = [];
cellArray = regexp(tmp, '*', 'split');
end

% percentagesGroups is a cell array that includes percentages column array.
% For example: percentaesGroups = {percentageHigh, percentageLow}.
function pvalue = pvalueForEachSubtype(patientsNames, subtypeCategorical, percentagesGroups)
totNumOfPatients = length(patientsNames);
numPatientsInSubtype = countcats(subtypeCategorical);
expectedValues = numPatientsInSubtype * 100 / totNumOfPatients; % column array
chi2 = 0;
for groupIdx = 1:size(percentagesGroups, 2)
    chi2 = chi2 + (percentagesGroups{groupIdx} - expectedValues).^2./expectedValues;
end
% chi2cdf - Chi-square cumulative distribution function
pvalue = chi2cdf(chi2, size(percentagesGroups,2) - 1, 'upper'); % Column vector.
end

