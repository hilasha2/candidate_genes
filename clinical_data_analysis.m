function clinical_data_analysis(filename_km, outputDir, txtDataKM,...
    numericDataKM, patientsNamesKM, patientsNames, extraClinicalDataCellMatrix,...
    use_cbioportal_not_lab_std)
% ---------- Clinical data histograms
% If extraClinicalDataCellMatrix == [], then no calculations are done.
% extraClinicalData - [{'hello'}, {1}; {'bye'}, {2}] or [{'hello', 1}; {'bye', 2}]
% if ~isempty(extraClinicalDataCellMatrix) clinical_data_analysis end.
for i = 1:size(extraClinicalDataCellMatrix,1)
    clinicalDataName = extraClinicalDataCellMatrix{i, 1};
    clinicalDataIdx = extraClinicalDataCellMatrix{i, 2};
    str = sprintf('clinical_data_%s.xlsx', clinicalDataName);
    fileOutput = fullfile(outputDir, str);
    [patientsNamesUpdated, clinicalData] = readData(txtDataKM,...
    numericDataKM, patientsNamesKM ,use_cbioportal_not_lab_std, clinicalDataIdx);
    getHistogramClinicalData(filename_km, fileOutput, patientsNames,...
    patientsNamesUpdated, topIdx, lowIdx, clinicalData, clinicalDataName);
end
end

function [patientsNamesKM, clinicalData] = readData(txtDataKM,...
    numericDataKM, patientsNamesKM ,use_cbioportal_not_lab_std, clinicalDataIdx)
if use_cbioportal_not_lab_std
    dataType = txtDataKM(3, clinicalDataIdx);
    if isequal(dataType, 'NUMBER')
        clinicalData = numericDataKM(3:end, clinicalDataIdx - 1);
        
    elseif isequal(dataType, 'STRING')
        clinicalData = txtDataKM(6:end, clinicalDataIdx);
    end
else
    dataType = txtDataKM(1, clinicalDataIdx);
    if isequal(dataType, 'NUMBER')
        clinicalData = numericDataKM(3:end, clinicalDataIdx - 1);
    elseif isequal(dataType, 'STRING')
        clinicalData = txtDataKM(4:end, clinicalDataIdx);
    end
end
nanIndices = isnan(clinicalDataIdx); 
clinicalData(nanIndices) = [];
patientsNamesKM(nanIndices) = [];
end

function getHistogramClinicalData(filename_km, fileOutput, patientsNames,...
    patientsNamesUpdated, topIdx, lowIdx, clinicalData, clinicalDataName)

[~,idxTop] = intersect(patientsNamesKM, patientsNamesTop);
[~,idxLow] = intersect(patientsNamesKM, patientsNamesLow);
clinicalDataCategorical = categorical(clinicalData);
clinicalDataCategorical = clinicalDataCategorical(~isundefined(clinicalDataCategorical));
categories = unique(clinicalDataCategorical);
percentageHigh = getCategoriesPercentages(clinicalDataCategorical, idxTop);
percentageLow = getCategoriesPercentages(clinicalDataCategorical, idxLow);
histData = [percentageHigh, percentageLow];
% Plotting the histogram
plotTitle = sprintf('Histogram of %s clinical data, separated into high-low %s', clinicalDataName, geneName);
figure('Name', plotTitle);
bar(histData);
clear title xlabel ylabel; 
xticklabels(cellstr(categories));
legend('High','Low','Location', 'northeastoutside');
xLabel = sprintf('%s categories', clinicalDataName);
xlabel(xLabel);
ylabel('Percentage of patients in each high/low category');
title(plotTitle);
for i = 1:length(categories)
    text(i-0.3,percentageHigh(i),num2str(percentageHigh(i),2),'vert','bottom','horiz','left');
    text(i+0.3,percentageLow(i),num2str(percentageLow(i),2),'vert','bottom','horiz','right');
end

categoricalTop = clinicalDataCategorical(idxTop);
categoricalLow = clinicalDataCategorical(idxLow);
totPatientsCategorical = [categoricalTop; categoricalLow];
% Writing to an excel file the p-values
table = {fileInput, ''; 'Categories','p-value'};
% TODO - pvalues are calculated wrongly.
% TODO!!!!!!!!!! Write a specific pvalue calculating function
pvalues = pvalueForEachSubtype(totPatientsCategorical, totPatientsCategorical, {percentageHigh,percentageLow});
tmp = [cellstr(categories), num2cell(pvalues)];
table = [table; tmp];
xlswrite(fileOutput, table);
end


