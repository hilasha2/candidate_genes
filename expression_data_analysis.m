function expression_data_analysis(filename, outputDir, numericData,...
    geneNames, patientsNames, geneIdx, gene_name, structNumericDataH, ...
    structNumericDataL, structPatientsH)

% -------- Clustergram & Zscore

getClustergram(numericData, geneNames, patientsNames,...
    'zscored', 'clustergram on all patients');
title = sprintf('clustergram on all 30%% top %s patients', gene_name);
getClustergram(structNumericDataH.h30, geneNames, structPatientsH.h30,...
    'zscored top 30%', title); 
title = sprintf('clustergram on all 20%% top %s patients', gene_name);
getClustergram(structNumericDataH.h20, geneNames, structPatientsH.h20,...
    'zscored top 20%', title);
title = sprintf('clustergram on all 10%% top %s patients', gene_name);
getClustergram(structNumericDataH.h10, geneNames, structPatientsH.h10,...
    'zscored top 10%', title);
title = sprintf('clustergram on all 80%% top %s patients', gene_name);
getClustergram(structNumericDataH.h80, geneNames, structPatientsH.h80,...
    'zscored top 80%', title); 
title = sprintf('clustergram on all 70%% top %s patients', gene_name);
getClustergram(structNumericDataH.h70, geneNames, structPatientsH.h70,...
    'zscored top 70%', title); 

% ---------- TSNE

% getTsne(numericData, numGroups,...
%     'TSNE','TSNE with KMeans');
% numGroupsTop30 = 4;
% getTsne(top30Data, numGroupsTop30,...
%     'TSNE top 30%','TSNE top 30% with KMeans');


% ---------- Bar and Box plots

fileOutput = fullfile(outputDir,'barplots_pvalues_L70_H30.xlsx');
title = sprintf('Mean gene expressions of patients with high 30 %s vs low 70 %s', ...
    gene_name, gene_name);
getBarPlot(structNumericDataH.h30, structNumericDataL.l70, geneNames, title, filename, 30, geneIdx, gene_name, fileOutput);
title = sprintf('Box plot of gene expressions of patients with high 30 %s vs low 70 %s',...
    gene_name, gene_name);
getBoxPlot(structNumericDataH.h30, structNumericDataL.l70, geneNames, title);

fileOutput = fullfile(outputDir,'barplots_pvalues_L80_H20.xlsx');
title = sprintf('Mean gene expressions of patients with high 20 %s vs low 80 %s', ...
    gene_name, gene_name);
getBarPlot(structNumericDataH.h20, structNumericDataL.l80, geneNames, title, filename, 20, geneIdx, gene_name, fileOutput);
title = sprintf('Box plot of gene expressions of patients with high 20 %s vs low 80 %s',...
    gene_name, gene_name);
getBoxPlot(structNumericDataH.h20, structNumericDataL.l80, geneNames, title);

fileOutput = fullfile(outputDir,'barplots_pvalues_L90_H10.xlsx');
title = sprintf('Mean gene expressions of patients with high 10 %s vs low 90 %s', ...
    gene_name, gene_name);
getBarPlot(structNumericDataH.h10, structNumericDataL.l90, geneNames, title, filename, 10, geneIdx, gene_name, fileOutput);
title = sprintf('Box plot of gene expressions of patients with high 10 %s vs low 90 %s',...
    gene_name, gene_name);
getBoxPlot(structNumericDataH.h10, structNumericDataL.l90, geneNames,...
    title);

fileOutput = fullfile(outputDir,'barplots_pvalues_L20_H80.xlsx');
title = sprintf('Mean gene expressions of patients with high 80 %s vs low 20 %s', ...
    gene_name, gene_name);
getBarPlot(structNumericDataH.h80, structNumericDataL.l20, geneNames, title, filename, 80, geneIdx, gene_name, fileOutput);
title = sprintf('Box plot of gene expressions of patients with high 80 %s vs low 20 %s',...
    gene_name, gene_name);
getBoxPlot(structNumericDataH.h80, structNumericDataL.l20, geneNames, title);

fileOutput = fullfile(outputDir,'barplots_pvalues_L30_H70.xlsx');
title = sprintf('Mean gene expressions of patients with high 70 %s vs low 30 %s', ...
    gene_name, gene_name);
getBarPlot(structNumericDataH.h70, structNumericDataL.l30, geneNames, title, filename, 70, geneIdx, gene_name, fileOutput);
title = sprintf('Box plot of gene expressions of patients with high 70 %s vs low 30 %s',...
    gene_name, gene_name);
getBoxPlot(structNumericDataH.h70, structNumericDataL.l30, geneNames, title);


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


% TODO - improve this function. 
% function getTsne(numericData, numGroups, titleTsne, titleKmeans)
% zscored = zscore(numericData'); % zscored - n * p 
% wcoeff = pca(zscored); % wcoeff - p * p, columns - coefficients for one principal component.
% tsne_z = tsne(wcoeff,'Distance','euclidean');
% cls = clusterdata(tsne_z(:,1:2),'maxclust',numGroups,'distance','euclidean');
% figure('Name', titleTsne, 'visible', 'off');
% gscatter(tsne_z(:,1),tsne_z(:,2),cls)
% title(titleTsne)
% idx = kmeans(tsne_z(:,1:2),numGroups); % kmeans clustering partition the observations, idx - n * 1.
% figure('Name', titleKmeans, 'visible', 'off');
% gscatter(tsne_z(:,1),tsne_z(:,2),idx)
% title(titleKmeans)
% end


function tbl = getBarPlot(topData, lowData, geneNames, titlePlot, inputFile, ...
    topPerc, geneIdx, gene_name, outputFile)
topMean = mean(topData, 2);
lowMean = mean(lowData, 2);
barData = cat(2, topMean, lowMean);
numBars = length(geneNames);
[~, pvalues] = ttest2(topData', lowData');
figure('Name', titlePlot, 'visible', 'off');
hBar = bar(1:numBars, barData);
bar1Coor = bsxfun(@plus, hBar(1).XData, [hBar(1).XOffset]'); 
bar2Coor = bsxfun(@plus, hBar(2).XData, [hBar(2).XOffset]');
hold on
stdTop = std(topData, 1, 2)./sqrt(length(topMean));
stdLow = std(lowData, 1, 2)./sqrt(length(lowMean));
errorbar(bar1Coor, topMean, stdTop, '.');
errorbar(bar2Coor, lowMean, stdLow, '.');
hold off
clear title xlabel ylabel; 
xticks(1:numBars);
xticklabels(geneNames);
xtickangle(90);
hGene = sprintf('High %s', gene_name);
lGene = sprintf('Low %s', gene_name);
legend(hGene, lGene,'Location', 'northeastoutside');
xlabel('Gene Names');
ylabel('Expression of Genes');
title(titlePlot);
geneName = sprintf('%s', gene_name);
% Constructing a table for excel summarising the results
% First 2 rows:
tbl = {inputFile, '', '', '', '', '';'data split according to', geneName, '', '', '', ''};
% 3rd & 4th rows
numHighPatientsStr = sprintf('no. of high-%s patients', gene_name);
numLowPatientsStr = sprintf('no. of low-%s patients', gene_name);
tbl = [tbl; {'top percentage', numHighPatientsStr, 'low percentage', numLowPatientsStr, '', ''}];
tbl = [tbl; {topPerc, size(topData, 2), 100-topPerc, size(lowData, 2), '', ''}];
% 5th row
tbl = [tbl; {'Gene Names', 'mean high expression', 'mean low expression', 'Ratio high/low', 'P-value', 'Significant'}];
% bool = (pvalues < 0.05 & pvalues <= pvalues(geneIdx));
bool = (pvalues < 0.05);
tmp = [geneNames, num2cell(topMean), num2cell(lowMean), num2cell(topMean./lowMean), num2cell(pvalues'), num2cell(bool')];
tmp = sortrows(tmp, 5); % sort by p-values
tbl = [tbl; tmp];
xlswrite(outputFile, tbl);
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