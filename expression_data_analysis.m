function expression_data_analysis(filename, outputDir, numericData,...
    geneNames, patientsNames, geneIdx, gene_name, structNumericDataH, ...
    structNumericDataL, structPatientsH, printAllFigures)

% -------- Clustergram & Zscore

if eq(printAllFigures, 1)
    getClustergram(numericData, geneNames, patientsNames,...
        'normalized gene expression of all patients', 'clustergram on all patients');
    [titleCluster, titleZScore] = getClustergramTitles(30, 'high', gene_name);
    getClustergram(structNumericDataH.h30, geneNames, structPatientsH.h30,...
        titleZScore, titleCluster); 
    
    [titleCluster, titleZScore] = getClustergramTitles(20, 'high', gene_name);
    getClustergram(structNumericDataH.h20, geneNames, structPatientsH.h20,...
        titleZScore, titleCluster);
    
    [titleCluster, titleZScore] = getClustergramTitles(10, 'high', gene_name);
    getClustergram(structNumericDataH.h10, geneNames, structPatientsH.h10,...
        titleZScore, titleCluster);
    
    [titleCluster, titleZScore] = getClustergramTitles(20, 'low', gene_name);
    getClustergram(structNumericDataH.h80, geneNames, structPatientsH.h80,...
        titleZScore, titleCluster); 
    
    [titleCluster, titleZScore] = getClustergramTitles(30, 'low', gene_name);
    getClustergram(structNumericDataH.h70, geneNames, structPatientsH.h70,...
        titleZScore, titleCluster); 
end

% ---------- Bar and Box plots

% ------- Box plots

if eq(printAllFigures, 1)
    boxTitle = getBoxPlotTitle(30, 'high', gene_name);
    getBoxPlot(structNumericDataH.h30, structNumericDataL.l70, ...
        geneNames, boxTitle);

    boxTitle = getBoxPlotTitle(30, 'low', gene_name);
    getBoxPlot(structNumericDataH.h70, structNumericDataL.l30, ....
        geneNames, boxTitle);

    boxTitle = getBoxPlotTitle(20, 'high', gene_name);
    getBoxPlot(structNumericDataH.h20, structNumericDataL.l80, ....
        geneNames, boxTitle);

    boxTitle = getBoxPlotTitle(20, 'low', gene_name);
    getBoxPlot(structNumericDataH.h80, structNumericDataL.l20, ...
        geneNames, boxTitle);

    boxTitle = getBoxPlotTitle(10, 'high', gene_name);
    getBoxPlot(structNumericDataH.h10, structNumericDataL.l90, ....
        geneNames, boxTitle);
end

% ------- Barplots

[barTitle, barFileName] = getBarTitleAndFileName(30, 'high', gene_name, outputDir);
getBarPlot(structNumericDataH.h30, structNumericDataL.l70, geneNames, barTitle, ...
    filename, gene_name, barFileName, 'expression', printAllFigures, 30);

[barTitle, barFileName] = getBarTitleAndFileName(30, 'low', gene_name, outputDir);
getBarPlot(structNumericDataH.h70, structNumericDataL.l30, geneNames, barTitle, ...
    filename, gene_name, barFileName, 'expression', printAllFigures, 70);

[barTitle, barFileName] = getBarTitleAndFileName(20, 'high', gene_name, outputDir);
getBarPlot(structNumericDataH.h20, structNumericDataL.l80, geneNames, barTitle, ...
    filename, gene_name, barFileName, 'expression', printAllFigures, 20);

[barTitle, barFileName] = getBarTitleAndFileName(20, 'low', gene_name, outputDir);
getBarPlot(structNumericDataH.h80, structNumericDataL.l20, geneNames, barTitle,...
    filename, gene_name, barFileName, 'expression', printAllFigures, 80);

[barTitle, barFileName] = getBarTitleAndFileName(10, 'high', gene_name, outputDir);
getBarPlot(structNumericDataH.h10, structNumericDataL.l90, geneNames, barTitle, ...
    filename, gene_name, barFileName, 'expression', printAllFigures, 10);

end


%% Supporting functions

function CGobj = getClustergram(numericData, geneNames, patientsNames, titleZS, titleCG) 
figure('Name', titleZS, 'visible', 'off');
imagesc(numericData); % numericData is supposed to be already normalized 
colormap(jet);
colorbar
title(titleZS);
xlabel('patients');
set(gca,'YTickLabel', geneNames,'YTick', 1:length(geneNames));
% Rows are genes, and columns are samples.
CGobj = clustergram(numericData,'ColumnLabels',patientsNames,'RowLabels',geneNames,... 
    'ImputeFun', @knnimpute, 'Standardize', 'row', 'ColumnLabelsRotate', 30);
addTitle(CGobj, titleCG);
end


function getBoxPlot(topData, lowData, geneNames, titlePlot)
numLbls = length(geneNames);
figure('Name', titlePlot, 'visible', 'off');
boxplot(topData', geneNames, 'positions',1:1:numLbls, 'width', 0.18);
hold on
boxplot(lowData', geneNames, 'positions',1.3:1:(numLbls+0.3), 'width', 0.18, 'colors','g');
hold off
clear title xlabel ylabel; 
xticks(1:numLbls);
xticklabels(geneNames);
xtickangle(90);
xlabel('Gene Names');
ylabel('Expression of Genes');
title(titlePlot);
end


function [titleCluster, titleZScore] = getClustergramTitles(numPercentage, strLowOrHigh, geneName)
titleCluster = sprintf('clustergram on %u%% patients with %s expression of %s', ...
    numPercentage, strLowOrHigh, geneName);
titleZScore = sprintf('normalized gene expression of %u%% patients with %s %s expression', ...
    numPercentage, strLowOrHigh, geneName);
end



function boxTitle = getBoxPlotTitle(numPercentage, strLowOrHigh, geneName)
boxTitle = sprintf('Boxplot of average gene expressions of %u%% patients with %s %s expression vs rest of patients',...
    numPercentage, strLowOrHigh, geneName);
end


function [barTitle, barFileName] = getBarTitleAndFileName(numPercentage, strLowOrHigh, geneName, outputDir)
barTitle = sprintf('Barplot of average gene expressions of %u%% patients with %s %s expression vs rest of patients',...
    numPercentage, strLowOrHigh, geneName);
if strcmp(strLowOrHigh, 'high')
    fileName = sprintf('barplots_pvalues_R%u_H%u.xlsx', ...
        100-numPercentage, numPercentage);
   
elseif strcmp(strLowOrHigh, 'low')
    fileName = sprintf('barplots_pvalues_L%u_R%u.xlsx', ...
        numPercentage, 100-numPercentage);
end
 barFileName = fullfile(outputDir,fileName);
end