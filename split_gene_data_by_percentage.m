function [topData, lowData, patientsNamesTop, patientsNamesLow]...
    = split_gene_data_by_percentage(numericData, patientsNames, topPerc, lowPerc, geneIdx)
[~, topIndices] = maxk(numericData(geneIdx,:), ceil(size(numericData, 2)*topPerc));
patientsNamesTop = patientsNames(1, topIndices);
topData = numericData(:, topIndices);
[~, lowIndices] = mink(numericData(geneIdx,:), floor(size(numericData, 2)*lowPerc));
patientsNamesLow = patientsNames(1, lowIndices); 
lowData = numericData(:, lowIndices); 
end