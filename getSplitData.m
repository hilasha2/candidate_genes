function  [structNumericDataH, structNumericDataL, ...
    structPatientsH, structPatientsL] = getSplitData(numericData, patientsNames, geneIdx)
% -------- splitting it into different expression percetages of main gene's expressions.
structNumericDataH = struct;
structNumericDataL = struct;
structPatientsH = struct;
structPatientsL = struct;
% Splitting the table into the top 30% main gene values and 70% low main gene values.
[structNumericDataH.h30, structNumericDataL.l70, structPatientsH.h30, structPatientsL.l70]...
    = split_gene_data_by_percentage(numericData, patientsNames, 0.3, 0.7, geneIdx);

% Splitting the table into the top 20% main gene values and 80% low main gene values.
[structNumericDataH.h20, structNumericDataL.l80, structPatientsH.h20, structPatientsL.l80]...
    = split_gene_data_by_percentage(numericData, patientsNames, 0.2, 0.8, geneIdx);

% Splitting the table into the top 10% main gene values and 90% low main gene values.
[structNumericDataH.h10, structNumericDataL.l90, structPatientsH.h10, structPatientsL.l90]...
    = split_gene_data_by_percentage(numericData, patientsNames, 0.1, 0.9, geneIdx);

% Splitting the table into the top 80% main gene values and 20% low main gene values.
[structNumericDataH.h80, structNumericDataL.l20, structPatientsH.h80, structPatientsL.l20]...
    = split_gene_data_by_percentage(numericData, patientsNames, 0.8, 0.2, geneIdx);

% Splitting the table into the top 70% main gene values and 30% low main gene values.
[structNumericDataH.h70, structNumericDataL.l30, structPatientsH.h70, structPatientsL.l30]...
    = split_gene_data_by_percentage(numericData, patientsNames, 0.7, 0.3, geneIdx);

end