function kaplan_meier_by_clustering(outputDir, numericData, geneNames, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff, printAllFigures) 

for numGroups = 2:6
    clusteringPlotTitle = sprintf('Clustering by k-means, numGroups = %d', numGroups);
    kmPlotTitle = sprintf('Kaplan Meier of patients clustered by K-mean, numGroups = %d', numGroups); 
    clustering_method = 'kmeans only';
    fileName = sprintf('gene_exp_averages_per_group_%d_groups_clustered_by_kmeans.xlsx', numGroups); 
    fileOutput = fullfile(outputDir,fileName); 
    kaplan_meier_by_clustering_method(numericData, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff, clusteringPlotTitle, ...
    kmPlotTitle, numGroups, clustering_method, geneNames, fileOutput, ...
    printAllFigures);

    clusteringPlotTitle = sprintf('Clustering by k-means after pca, numGroups = %d', numGroups);
    kmPlotTitle = sprintf('Kaplan of patients clustered by K-mean after pca, numGroups = %d', numGroups); 
    clustering_method = 'kmeans after pca';
    fileName = sprintf('gene_exp_averages_per_group_%d_groups_clustered_by_kmeans_after_pca.xlsx', numGroups);
    fileOutput = fullfile(outputDir,fileName); 
    kaplan_meier_by_clustering_method(numericData, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff, clusteringPlotTitle, ...
    kmPlotTitle, numGroups, clustering_method, geneNames, fileOutput, ...
    printAllFigures);
end
end

%% Supporting functions

% ---
% We cluster the patients by k-means
function groupsIdx = getGroupsByKmeans(numericData, numGroups, ...
    plotTitle, printAllFigures)
groupsIdx = kmeans(numericData', numGroups); 
if printAllFigures
    plotScatteredGroups(numericData', groupsIdx, plotTitle)
end
end

% ---
% Clustering patients where numericData is spanned by pca space. 
function groupsIdx = getGroupsByKmeansAfterPca(numericData, numGroups, ...
    plotTitle, printAllFigures)
% pca(groupsData) - rows, n - observations (patients), columns, p - variables (genes). 
% groupsData = numericData' - n*p
% dataInPcaSpace - n*p, pcaEigenvalues (coeff) - p*p
[~,dataInPcaSpace] = pca(numericData'); % zscore(dataInPcaSpace) is zscore(groupsData*coeff);  
groupsIdx = getGroupsByKmeans(dataInPcaSpace', numGroups, plotTitle, ...
    printAllFigures);
end

% ---
% Plotting the patients by their groups in a scatter plot. 
% The patients have p dimensions (number of genes).
% We choose to plot the patients by a different space spanned by pca components.
% We plot the data by the first 2 pca components (eigen-values). 
function plotScatteredGroups(groupsData, groupsIdx, plotTitle)
[~, dataInPcaSpace, pcaEigenvalues] = pca(groupsData);  
percEigenvalues = pcaEigenvalues./sum(pcaEigenvalues) * 100;
figure('Name', plotTitle, 'visible', 'off');
gscatter(dataInPcaSpace(:,1),dataInPcaSpace(:,2), groupsIdx);
xlabel(['First Principal Component ' num2str(percEigenvalues(1)) '%']);
ylabel(['Second Principal Component ' num2str(percEigenvalues(2)) '%']);
title(plotTitle);
end 

% ---
% Getting all required inputs for kaplan meier analysis. 
% The groups for the kaplan meier are the patients separated by our desired
% clustering method.
function [timeVar, censVar, groupVar] = ...
    get_time_cens_groups(patientsNames, patientsNamesKM, groupsIdx, timeData, cens)
groups = unique(groupsIdx); % unique orders the result.
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
    group = repmat(idx, length(timeGroup), 1);
    groupVar = [groupVar; group]; 
end
% event = 1, cens = 0 - so the cens' matrices have to be reversed.
censVar = ~censVar; 
% Apparently, we can get more than 2 groups only if groups are cell array
% so we have to do the following: 
categoricalGroups = categorical(groupVar);
groupVar = cellstr(categoricalGroups);
end

% ---
% Kaplan Meier function 
function stats = kaplan_meier(TimeVar, EventVar, GroupVar, timeCutOff, ...
    titlePlot, printAllFigures)
if printAllFigures
    [~, fh, stats] = MatSurv(TimeVar, EventVar, GroupVar, 'XLim', timeCutOff,...
        'Title', titlePlot, 'NoRiskTable', false, 'PairwiseP', true, ...
        'TitleOptions', {'FontSize', 11});
    fh.Name = titlePlot;
    fh.Visible = 'off';
else
    [~, ~, stats] = MatSurv(TimeVar, EventVar, GroupVar, 'XLim', timeCutOff,...
        'NoPlot', true,'PairwiseP', true);
end
end 

% ---
% Main function
function kaplan_meier_by_clustering_method(numericData, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff, clusteringPlotTitle, ...
    kmPlotTitle, numGroups, clustering_method, geneNames, outputFile, ...
    printAllFigures)
switch clustering_method
    case 'kmeans only'
        groupsIdx = getGroupsByKmeans(numericData, numGroups, ...
            clusteringPlotTitle, printAllFigures);
    case 'kmeans after pca'
        groupsIdx = getGroupsByKmeansAfterPca(numericData, numGroups, ...
            clusteringPlotTitle, printAllFigures);
    otherwise
        disp('No clustering method was procided')
        return
end 
[timeVar, censVar, groupVar] = ...
    get_time_cens_groups(patientsNames, patientsNamesKM, groupsIdx, timeData, cens);
kaplan_meier(timeVar, censVar, groupVar, timeCutOff, kmPlotTitle, ...
    printAllFigures);
if printAllFigures
    silhouetteTitle = sprintf('Silohouette of - %s', clusteringPlotTitle); 
    figure('Name', silhouetteTitle, 'visible', 'off');
    silhouette(numericData', groupsIdx);
    title(silhouetteTitle);
end
if(numGroups == 2)
    [averages_per_group, ~, stde_per_group] =...
        descriptiveStatisticsForGroups(geneNames, numericData, groupsIdx, outputFile);
    patientsGroups(patientsNames, groupsIdx, outputFile);
    divideGenesAndAvereageGeneExpressionIntoGroups ...
        (geneNames, averages_per_group, numGroups, clustering_method, ...
        outputFile);
    if printAllFigures
        meanGeneExpressionBarPlot(geneNames, averages_per_group, ....
        stde_per_group, numGroups, clustering_method);
        pearsonCorr(groupsIdx, numericData, geneNames, numGroups, ...
            clustering_method);
    end
end
end 
%% ---- Functions for analysing the data after the clustering.

% ---
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

% ---
% Calculating one way anova for each gene while the groups are based on the
% patients
function pvalues = anovaPerGene(data, groupsIdx)
numGenes = size(data, 2);
pvalues = ones(numGenes, 1);
for i = 1:numGenes
   pvalues(i) = anova1(data(:, i), groupsIdx, 'off');
end
end

% ---
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

% ---
% Writing to an excel sheet the patients and their division into groups.
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

% ---
% Plotting a bar plot of average gene expression per gene. Number of bar
% per gene is the number of groups clustered. Includes standard errors. 
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

% ---
function divideGenesAndAvereageGeneExpressionIntoGroups ...
    (geneNames, meanExpClustered, numGroups, clusteringMethod, ...
    outputFile)
if numGroups ~= 2
    disp('Number of patients clusters must be 2!');
    return
end
% Dividing genes according to their expression: above/bellow average in group 1 and
% above/bellow average in group 2. Average is 0.
idx_pos_pos = meanExpClustered(:,1) > 0 & meanExpClustered(:,2) > 0;
idx_pos_neg = meanExpClustered(:,1) > 0 & meanExpClustered(:,2) < 0;
idx_neg_pos = meanExpClustered(:,1) < 0 & meanExpClustered(:,2) > 0;
idx_neg_neg = meanExpClustered(:,1) < 0 & meanExpClustered(:,2) < 0;

scatterPlotPosNeg(idx_pos_pos, geneNames, meanExpClustered,...
    "above", "above", clusteringMethod);
scatterPlotPosNeg(idx_pos_neg, geneNames, meanExpClustered,...
    "above", "bellow", clusteringMethod);
scatterPlotPosNeg(idx_neg_pos, geneNames, meanExpClustered,...
    "bellow", "above", clusteringMethod);
scatterPlotPosNeg(idx_neg_neg, geneNames, meanExpClustered,...
    "bellow", "bellow", clusteringMethod);


tbl1 = createTablePosNeg(idx_pos_pos, geneNames, meanExpClustered, ...
    "above", "above");
tbl2 = createTablePosNeg(idx_pos_neg, geneNames, meanExpClustered, ...
    "above", "bellow");
tbl3 = createTablePosNeg(idx_neg_pos, geneNames, meanExpClustered, ...
    "bellow", "above");
tbl4 = createTablePosNeg(idx_neg_neg, geneNames, meanExpClustered, ...
    "bellow", "bellow");
tbl12 = outerjoin(tbl1, tbl2);
tbl34 = outerjoin(tbl3, tbl4);
tbl = outerjoin(tbl12, tbl34);
writetable(tbl, outputFile, 'Sheet', 3);
end

% ---
% Scatter plot of average gene expression - patients divided by 2 clusters, 
% and by a specific 'group'. 'group' is based on the value of expression
% above/bellow average expression. Average expression is 0.
function scatterPlotPosNeg(idx, geneNames, meanExpClustered,...
    group1avg, group2avg, clusteringMethod)
figName = sprintf('Scatter plot of average gene expression per cluster %s grp1 %s grp 2 %s', ...
    clusteringMethod, group1avg, group2avg);
figure('Name', figName, 'Visible', 'off'); 
geneNamesInGroups = geneNames(idx);
meanExpClusteredInGroup = meanExpClustered(idx, :);
cat_geneNamesInGroups = categorical(geneNamesInGroups);
sh1 = scatter(cat_geneNamesInGroups, meanExpClusteredInGroup(:, 1), 4);  
hold on
sh2 = scatter(cat_geneNamesInGroups, meanExpClusteredInGroup(:, 2), 4); 
hold off
clear title xlabel ylabel; 
titleStr = sprintf("Average Gene Expression Per Cluster\nClustering method: %s",...
    clusteringMethod);
annotationStr = sprintf("Genes in group 1: have %s avg exp\nGenes in group 2: have %s avg exp",...
                group1avg, group2avg);
title(titleStr);
annotation('textbox', [.2 .5 .3 .3], 'String', annotationStr,...
    'FitBoxToText', 'on');
legend([sh1 sh2], 'group 1', 'group 2'); 
xlabel('Gene Names');
ylabel('Average Expression of Genes');
end

% ---
function tbl = createTablePosNeg(idx, geneNames, meanExpClustered, ...
    grp1avg, grp2avg)
grp1Title = sprintf('group 1 (%s average)', grp1avg);
grp2Title = sprintf('group 2 (%s average)', grp2avg);
variableNames = {'Genes', grp1Title, grp2Title}; 
tbl = table(geneNames(idx), meanExpClustered(idx, 1), meanExpClustered(idx, 2), ...
    'VariableNames', variableNames);
end

% ---
% Pearson correlation between the genes in each group
function pearsonCorr(idx, numericData, geneNames, numGroups, clusteringMethod)
transposedNumericData = numericData'; % transposedNum.. row - patients, cols - genes.
for i = 1:numGroups
    [correlation_mat, pvals] = corrcoef(transposedNumericData(idx == i, :));
    coef_str = "Pearson Correlation Coefficients";
    pval_str = "Pearson Correlation Pvalues";
    pearsonFig(coef_str, correlation_mat, geneNames, clusteringMethod, numGroups, i);
    pearsonFig(pval_str, pvals, geneNames, clusteringMethod, numGroups, i);
end
end

% ---
% Creating figures for pearson correlation
function pearsonFig(plotString, data, geneNames, clusteringMethod, numGroups, groupID)
numGenes = length(geneNames);
figNameCorr = sprintf("%s - %d clusters %s group number %d",...
        plotString, numGroups, clusteringMethod, groupID);
figure('Name', figNameCorr, 'Visible', 'off');
imagesc(data);
colorbar;
set(gca,'XTick', 1:numGenes,'XTickLabel', geneNames, 'XTickLabelRotation', 90, ...
    'YTick', 1:numGenes, 'YTickLabel', geneNames);
titlePlot = sprintf("%s\nnumber of clusters: %d, group: %d \nclustering method: %s", ...
   plotString, numGroups, groupID, clusteringMethod); 
title(titlePlot);
end