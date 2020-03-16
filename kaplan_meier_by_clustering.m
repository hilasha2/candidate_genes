function kaplan_meier_by_clustering(outputDir, numericData, geneNames, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff) 

for numGroups = 2:6
    clusteringPlotTitle = sprintf('Clustering by k-means, numGroups = %d', numGroups);
    kmPlotTitle = sprintf('Kaplan of patients clustered by K-mean, numGroups = %d', numGroups); 
    clustering_method = 'kmeans only';
    fileName = sprintf('gene_averages_per_group_%d_groups_clustered_by_kmeans.xlsx', numGroups); 
    fileOutput = fullfile(outputDir,fileName); 
    kaplan_meier_by_clustering_method(numericData, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff, clusteringPlotTitle, ...
    kmPlotTitle, numGroups, clustering_method, geneNames, fileOutput);

    clusteringPlotTitle = sprintf('Clustering by k-means after pca, numGroups = %d', numGroups);
    kmPlotTitle = sprintf('Kaplan of patients clustered by K-mean after pca, numGroups = %d', numGroups); 
    clustering_method = 'kmeans after pca';
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
    'Title', titlePlot, 'NoRiskTable', true, 'PairwiseP', true, 'NoPlot', true, ...
    'TitleOptions', {'FontSize', 11});
fh.Name = titlePlot;
end 


function kaplan_meier_by_clustering_method(numericData, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff, clusteringPlotTitle, ...
    kmPlotTitle, numGroups, clustering_method, geneNames, outputFile)
switch clustering_method
    case 'kmeans only'
        groupsIdx = getGroupsByKmeans(numericData, numGroups, clusteringPlotTitle);
    case 'kmeans after pca'
        groupsIdx = getGroupsByKmeansAfterPca(numericData, numGroups, clusteringPlotTitle);
    otherwise
        disp('No clustering method was procided')
        return
end 
[timeVar, censVar, groupVar] = ...
    get_time_cens_groups(patientsNames, patientsNamesKM, groupsIdx, timeData, cens);
kaplan_meier(timeVar, censVar, groupVar, timeCutOff, kmPlotTitle)
silhouetteTitle = sprintf('Silohouette of - %s', clusteringPlotTitle); 
figure('Name', silhouetteTitle, 'visible', 'off');
silhouette(numericData', groupsIdx);
title(silhouetteTitle);
if(numGroups == 2)
    [averages_per_group, ~, stde_per_group] =...
        descriptiveStatisticsForGroups(geneNames, numericData, groupsIdx, outputFile);
    patientsGroups(patientsNames, groupsIdx, outputFile);
    meanGeneExpressionBarPlot(geneNames, averages_per_group, ....
    stde_per_group, numGroups, clustering_method)
end
end 

% Calculating average, standard deviation, standard error, t-test (for 2
% groups), and anova for each gene in groups.
function [averages_per_group, std_per_group, stde_per_group, pvalues_anova]...
    = descriptiveStatisticsForGroups(geneNames, numericData, groupsIdx, outputFile)
groups = unique(groupsIdx);
numGroups = length(groups);
numPatients = size(numericData, 2);
range_groups = 1:numGroups;
numericDataTransposed = numericData'; % rows - patients, % cols - genes
% Average, standard deviation and error of gene expression per gene per
% group
averages_per_group = [];
std_per_group = [];
stde_per_group = [];
 
% Handling the table for the results 
tbl1 = ["Gene Names"; geneNames];
template_headers = ["Average of gene expression for group", ...
    "Standard Deviation for group", "Standard Error for group"];
tmp_headers_repeated = repelem(template_headers, numGroups);
tmp_groups_repeated = repmat(range_groups, 1, length(template_headers));
% headers should look like: ["hello 1", "hello 2", "world 1", "world 2"] 
headers = tmp_headers_repeated + " " + tmp_groups_repeated;

% Calculating average, standard deviation and standard error
for idx = range_groups
    % group - Gene expression table for the patients in group #idx
    group = numericDataTransposed(idx == groupsIdx, :); 
    gene_average_per_group = mean(group, 1);
    gene_std_per_group = std(group, 1);
    gene_stde_per_group = gene_std_per_group/sqrt(numPatients);
    
    averages_per_group = [averages_per_group, gene_average_per_group'];
    std_per_group = [std_per_group, gene_std_per_group'];
    stde_per_group = [stde_per_group, gene_stde_per_group'];
end

% Calculating anova 
pvalues_anova = anovaPerGene(numericDataTransposed, groupsIdx);
headers = [headers, "Pvalue ANOVA"];
tbl = [tbl1, [headers; averages_per_group, ...
    std_per_group, stde_per_group, pvalues_anova]];

% Calculating t-test
if numGroups == 2
    pvalues_ttest = ttestPerGene(numericDataTransposed, groupsIdx);
    tmp = ["P-value of t test"; pvalues_ttest];
    tbl = [tbl, tmp];
end

%empty_spaces = strings(size(table, 1) - 2, 1);
writematrix(tbl, outputFile, 'Sheet', 1);
end

% Calculating one way anova for each gene while the groups are based on the
% patients
function pvalues = anovaPerGene(data, groupsIdx)
numGenes = size(data, 2);
pvalues = ones(numGenes, 1);
for i = 1:numGenes
   pvalues(i) = anova1(data(:, i), groupsIdx, 'off');
end
end

% Calculating t-test for each gene while the groups are based on the
% patients
function pvalues = ttestPerGene(data, groupsIdx)
numGenes = size(data, 2);
pvalues = ones(numGenes, 1);
group1 = data(groupsIdx == 1, :);
group2 = data(groupsIdx == 2, :);
for i = 1:numGenes
    [~,pvalues(i)] = ttest2(group1(:,i), group2(:,i));
end
end


function patientsGroups(patientsNames, groupsIdx, outputFile)
groups = unique(groupsIdx);
numGroups = length(groups);
headers = ["Patients", "Groups"];
patients = [];
groups = [];
for i = 1:numGroups 
    patients_in_group = patientsNames(groupsIdx == i);
    patients = [patients; patientsNames(groupsIdx == i)'];
    group_num_repeated = repmat(i, length(patients_in_group), 1);
    groups = [groups; group_num_repeated];
end
tbl = [headers; [patients, num2cell(groups)]];
writematrix(tbl, outputFile, 'Sheet', 2);
end


function meanGeneExpressionBarPlot(geneNames, meanExpClustered, ....
    errorClustered, numGroups, clusteringMethod)
numBars = length(geneNames);
figName = sprintf("Average gene expression per cluster - %d clusters %s",...
    numGroups, clusteringMethod);
figure('Name', figName, 'Visible', 'off');
hBar = bar(1:numBars, meanExpClustered, 'grouped');
hold on
for i = 1:numGroups
    barCoor = bsxfun(@plus, hBar(i).XData, [hBar(i).XOffset]');
    errorbar(barCoor, meanExpClustered(:,i), errorClustered(:,i), '.');
end
hold off
clear title xlabel ylabel; 
xticks(1:numBars);
xticklabels(geneNames);
str = sprintf("Average Gene Expression Per Cluster\nClustering method: %s\nNumber of groups: %d",...
    clusteringMethod, numGroups);
title(str);
xlabel('Gene Names');
ylabel('Average Expression of Genes');
end