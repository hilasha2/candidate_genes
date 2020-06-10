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

fileOutput = fullfile(outputDir,'barplots_pvalues_R70_H30.xlsx');
title = sprintf('Mean gene expressions of 30 percent patients with high %s vs rest', ...
    gene_name);
getBarPlot(structNumericDataH.h30, structNumericDataL.l70, geneNames, title, ...
    filename, gene_name, fileOutput, 'expression', 30);
title = sprintf('Box plot of gene expressions of 30 percent patients with high %s vs rest',...
    gene_name);
getBoxPlot(structNumericDataH.h30, structNumericDataL.l70, geneNames, title);

fileOutput = fullfile(outputDir,'barplots_pvalues_R80_H20.xlsx');
title = sprintf('Mean gene expressions of 20 percent patients with high %s vs rest', ...
    gene_name);
getBarPlot(structNumericDataH.h20, structNumericDataL.l80, geneNames, title, ...
    filename, gene_name, fileOutput, 'expression', 20);
title = sprintf('Box plot of gene expressions of 20 percent patients with high %s vs rest',...
    gene_name);
getBoxPlot(structNumericDataH.h20, structNumericDataL.l80, geneNames, title);

fileOutput = fullfile(outputDir,'barplots_pvalues_R90_H10.xlsx');
title = sprintf('Mean gene expressions of 10 percent patients with high %s vs rest', ...
    gene_name);
getBarPlot(structNumericDataH.h10, structNumericDataL.l90, geneNames, title, ...
    filename, gene_name, fileOutput, 'expression', 10);
title = sprintf('Box plot of gene expressions of 10 percent patients with high %s vs rest',...
    gene_name);
getBoxPlot(structNumericDataH.h10, structNumericDataL.l90, geneNames,...
    title);

fileOutput = fullfile(outputDir,'barplots_pvalues_L20_R80.xlsx');
title = sprintf('Mean gene expressions of 20 percent patients with low %s vs rest', ...
    gene_name);
getBarPlot(structNumericDataH.h80, structNumericDataL.l20, geneNames, title,...
    filename, gene_name, fileOutput, 'expression', 80);
title = sprintf('Box plot of gene expressions of 20 percent patients with low %s vs rest',...
    gene_name);
getBoxPlot(structNumericDataH.h80, structNumericDataL.l20, geneNames, title);

fileOutput = fullfile(outputDir,'barplots_pvalues_L30_R70.xlsx');
title = sprintf('Mean gene expressions of 30 percent patients with low %s vs rest', ...
    gene_name);
getBarPlot(structNumericDataH.h70, structNumericDataL.l30, geneNames, title, ...
    filename, gene_name, fileOutput, 'expression', 70);
title = sprintf('Box plot of gene expressions of 30 percent patients with low %s vs rest',...
    gene_name);
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