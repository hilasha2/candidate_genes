function kaplan_meier_by_clustering(outputDir, numericData, geneNames, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff) 

for numGroups = 2:6
    clusteringPlotTitle = sprintf('Clustering by k-means, numGroups = %d', numGroups);
    kmPlotTitle = sprintf('Kaplan of patients clustered by K-mean, numGroups = %d', numGroups); 
    clustering_method = 'kmeans_only';
    fileName = sprintf('gene_averages_per_group_%d_groups_clustered_by_kmeans.xlsx', numGroups); 
    fileOutput = fullfile(outputDir,fileName); 
    kaplan_meier_by_clustering_method(numericData, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff, clusteringPlotTitle, ...
    kmPlotTitle, numGroups, clustering_method, geneNames, fileOutput);

    clusteringPlotTitle = sprintf('Clustering by k-means after pca, numGroups = %d', numGroups);
    kmPlotTitle = sprintf('Kaplan of patients clustered by K-mean after pca, numGroups = %d', numGroups); 
    clustering_method = 'kmeans_after_pca';
    fileName = sprintf('gene_averages_per_group_%d_groups_clustered_by_kmeans_after_pca.xlsx', numGroups);
    fileOutput = fullfile(outputDir,fileName); 
    kaplan_meier_by_clustering_method(numericData, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff, clusteringPlotTitle, ...
    kmPlotTitle, numGroups, clustering_method, geneNames, fileOutput);
end
end

%% Supporting functions

% We cluster the patients by k-means
function groupsIdx = getGroupsByKmeans(numericData, numGroups, plotTitle)
groupsIdx = kmeans(numericData', numGroups); 
plotScatteredGroups(numericData', groupsIdx, plotTitle)
end


% Clustering patients where numericData is spanned by pca space. 
function groupsIdx = getGroupsByKmeansAfterPca(numericData, numGroups, plotTitle)
% pca(groupsData) - rows, n - observations, columns, p - variables. 
% groupsData - n*p
% dataInPcaSpace - n*p, pcaEigenvalues - p*p
[~,dataInPcaSpace] = pca(numericData'); % zscore(data) is zscore(groupsData*coeff);  
groupsIdx = getGroupsByKmeans(dataInPcaSpace', numGroups, plotTitle);
end


% Plotting the patients by their groups in a scatter plot. 
% Since the patients have p dimensions (number of genes), we choose 
% to plot the patients by a different space spanned by pca components.
function plotScatteredGroups(groupsData, groupsIdx, plotTitle)
[~, dataInPcaSpace, pcaEigenvalues] = pca(groupsData);  
percEigenvalues = pcaEigenvalues./sum(pcaEigenvalues) * 100;
figure('Name', plotTitle, 'visible', 'off');
gscatter(dataInPcaSpace(:,1),dataInPcaSpace(:,2), groupsIdx);
xlabel(['First Principal Component ' num2str(percEigenvalues(1)) '%']);
ylabel(['Second Principal Component ' num2str(percEigenvalues(2)) '%']);
title(plotTitle);
end 

% Getting all required inputs for kaplan meier analysis. 
% The groups for the kaplan meier are the patients separated by our desired
% clustering method.
function [timeVar, censVar, groupVar] = ...
    get_time_cens_groups(patientsNames, patientsNamesKM, groupsIdx, timeData, cens)
groups = unique(groupsIdx);
numGroups = length(groups);
timeVar = [];
censVar = [];
groupVar = []; 
for idx = 1:numGroups
    patientsInGroups = patientsNames(groupsIdx == idx); 
    [timeGroup, censGroup] = KM_support_functions.get_time_and_cens(...
        patientsNamesKM, timeData, cens, patientsInGroups);
    timeVar = [timeVar; timeGroup];
    censVar = [censVar; censGroup];
    group = repmat(idx, length(patientsInGroups), 1);
    groupVar = [groupVar; group]; 
end
% event = 1, cens = 0 - so the cens' matrices have to be reversed.
censVar = ~censVar; 
% Apparently, we can get more than 2 groups only if groups are cell array
% so we have to do the following: 
categoricalGroups = categorical(groupVar);
groupVar = cellstr(categoricalGroups);
end


function stats = kaplan_meier(TimeVar, EventVar, GroupVar, timeCutOff, titlePlot)
[~, fh, stats] = MatSurv(TimeVar, EventVar, GroupVar, 'XLim', timeCutOff,...
    'Title', titlePlot, 'NoRiskTable', true, 'PairwiseP', true, 'TitleOptions', {'FontSize', 11});
fh.Name = titlePlot;
end 


function kaplan_meier_by_clustering_method(numericData, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff, clusteringPlotTitle, ...
    kmPlotTitle, numGroups, clustering_method, geneNames, outputFile)
switch clustering_method
    case 'kmeans_only'
        groupsIdx = getGroupsByKmeans(numericData, numGroups, clusteringPlotTitle);
    case 'kmeans_after_pca'
        groupsIdx = getGroupsByKmeansAfterPca(numericData, numGroups, clusteringPlotTitle);
    otherwise
        disp('No clustering method was procided')
        return
end 
[timeVar, censVar, groupVar] = ...
    get_time_cens_groups(patientsNames, patientsNamesKM, groupsIdx, timeData, cens);
kaplan_meier(timeVar, censVar, groupVar, timeCutOff, kmPlotTitle)
silhouetteTitle = sprintf('Silohouette of - %s', clusteringPlotTitle); 
figure('Name', silhouetteTitle, 'visible', off);
silhouette(numericData', groupsIdx);
title(silhouetteTitle);
if(numGroups == 2)
    averageGeneExpressionForEachGroup(geneNames, numericData, groupsIdx, outputFile)
end
end 


function averageGeneExpressionForEachGroup(geneNames, numericData, groupsIdx, outputFile)
groups = unique(groupsIdx);
numGroups = length(groups);
numericDataTransposed = numericData';
averages_per_gene = [];

tbl1 = ["Gene Names"; geneNames];
tbl2 = ["Average of gene expression for group"; "Standard Deviation for group"];
for idx = 1:numGroups
    % Gene expression table for the patients in group #idx
    group = numericDataTransposed(idx == groupsIdx, :); 
    gene_average_per_group = mean(group, 1);
    averages_per_gene = [averages_per_gene, gene_average_per_group'];
    average_all_genes_per_group = mean(gene_average_per_group);
    std_genes_per_group = std(gene_average_per_group);
    header_for_group = sprintf("group %d", idx);
    tmp = [header_for_group; gene_average_per_group'];
    tbl1 = [tbl1, tmp];
    tbl2 = [tbl2, [average_all_genes_per_group; std_genes_per_group]];
end
table = [tbl1;tbl2];
if numGroups == 2
    [~, pvalue] = ttest(averages_per_gene(:,1), averages_per_gene(:,2));
    empty_spaces = strings(size(table, 1) - 2, 1);
    tmp = ["P-value of t test"; pvalue; empty_spaces];
    table = [table, tmp];
end
xlswrite(outputFile, table);
end