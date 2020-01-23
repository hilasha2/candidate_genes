function mutation_analysis(filename_mut, numericData, geneNames, geneNamesMut, patientsNames, ...
    patientsNamesMut, mainGeneName, mainGeneIdx, outputDir)

    fileOutput = fullfile(outputDir, 'mutation_and_expression_80_20.xlsx');
    barPlotMutatedGenesSplitByMainGeneExpression(filename_mut, numericData, geneNames, ... 
    geneNamesMut, patientsNames, patientsNamesMut, mainGeneName, mainGeneIdx, 0.2, 0.8, fileOutput);

    fileOutput = fullfile(outputDir, 'mutation_in_main_gene_vs_candidate_genes.xlsx');
    barPlotMutatedGenesVsMutatedMainGene(filename_mut, geneNames, geneNamesMut, ...
    patientsNames, patientsNamesMut, mainGeneName, mainGeneIdx, fileOutput);
    

    barPlotAvgExpressionOfMainGeneMutatedVsNot(numericData, geneNamesMut, ...
    patientsNames, patientsNamesMut, mainGeneName, mainGeneIdx)
end


%% Supporting functions

% Splits patients into topPerc and lowPerc. Then the patients are split
% again according to whether they have a mutation in a certain gene or not.
function barPlotMutatedGenesSplitByMainGeneExpression(filename_mut, numericData, geneNames, ...
    geneNamesMut, patientsNames, patientsNamesMut, mainGeneName, mainGeneIdx, topPerc, lowPerc, fileOutput)
 [~, ~, patientsNamesTop, patientsNamesLow]...
    = split_gene_data_by_percentage(numericData, patientsNames, topPerc, lowPerc, mainGeneIdx);
tableTitle = {filename_mut, '', '', '', '', ''; ...
    'Patients expression was split by the following main gene:', ...
    'top perc.', 'low perc.', '', '', ''; ...
    mainGeneName, num2str(topPerc*100), num2str(lowPerc*100), '', '', ''; ...
    'Candidate gene name', 'value', ...
    'High Expression in main gene - mutation in candidate gene', ...
    'High Expression in main gene - no mutation in candidate gene', ...
    'Low Expression in main gene - mutation in candidate gene', ...
    'Low Expression in main gene - no mutation in candidate gene'};
table = {};
for geneIdx = 1:length(geneNames)
    candidateGene = geneNames{geneIdx};
    [patientsWithMutatedGene, patientsWithNonMutatedGene] = ...
        splitPatientsByGeneMutation(candidateGene, geneNamesMut, patientsNames, patientsNamesMut);
    patientsHM = intersect(patientsNamesTop, patientsWithMutatedGene);
    patientsHNM = intersect(patientsNamesTop, patientsWithNonMutatedGene);
    patientsLM = intersect(patientsNamesLow, patientsWithMutatedGene);
    patientsLNM = intersect(patientsNamesLow, patientsWithNonMutatedGene);
    splitPatientsCell = {patientsHM, patientsHNM, patientsLM, patientsLNM};
    histData = percentage(splitPatientsCell, patientsNames);
    pvalues = pvaluesChi2(splitPatientsCell);
    table1 = {candidateGene, 'no. of patients'; '', 'percantage of patients'; '', 'p-value'}; 
    table2 = [doubleArrayToCellArray(cellfun(@length, splitPatientsCell));
        doubleArrayToCellArray(histData);
        doubleArrayToCellArray(pvalues)];
    table = [table; table1, table2];
    
    % ploting the histogram
    plotTitle = sprintf('Perc. of patients split by %s %d%% vs. %d%% & mutation vs. non in %s',...
        mainGeneName, lowPerc*100, topPerc*100, candidateGene);
    figure('Name', plotTitle, 'visible', 'off');
    bar(histData);
    clear title xlabel ylabel; 
    xLabels = {'High-Mutated', 'High-Non Mutated', 'Low-Mutated', 'Low-Non Mutated'};
    xticklabels(xLabels);
    xtickangle(45);
    ylabel('Percentage of patients'); 
    title(plotTitle, 'FontSize', 16);
    for i = 1:length(xLabels)
        text(i, histData(i),num2str(histData(i),2),'vert','bottom');
    end
end
table = [tableTitle; table];
xlswrite(fileOutput, table);
end

% Splitting the patients according to mutation/non-mutation in main gene
% and later the patients are split again according to mutation/non-mutation
% in another gene. 
function barPlotMutatedGenesVsMutatedMainGene(filename_mut, geneNames, geneNamesMut, ...
    patientsNames, patientsNamesMut, mainGeneName, mainGeneIdx, fileOutput)
    [patientsWithMutatedGene1, patientsWithNonMutatedGene1] = ...
        splitPatientsByGeneMutation(mainGeneName, geneNamesMut, patientsNames, patientsNamesMut);
   tableTitle = {filename_mut, '', '', '', '', ''; ...
    'Patients were split first by the following main gene:', ...
    mainGeneName, '', '', '', ''; ...
    'Candidate gene name', 'value', ...
    'Mutation in main gene - mutation in candidate gene', ...
    'Mutation in main gene - no mutation in candidate gene', ...
    'No mutation in main gene - mutation in candidate gene', ...
    'No mutation in main gene - no mutation in candidate gene'};
    table = {};
   for geneIdx = 1:length(geneNames)
       if geneIdx ~= mainGeneIdx
           candidateGene = geneNames{geneIdx};
           [patientsWithMutatedGene2, patientsWithNonMutatedGene2] = ...
               splitPatientsByGeneMutation(candidateGene, geneNamesMut, patientsNames, patientsNamesMut);
               patientsM1_M2 = intersect(patientsWithMutatedGene1, patientsWithMutatedGene2);
               patientsM1_NM2 = intersect(patientsWithMutatedGene1, patientsWithNonMutatedGene2);
               patientsNM1_M2 = intersect(patientsWithNonMutatedGene1, patientsWithMutatedGene2);
               patientsNM1_NM2 = intersect(patientsWithNonMutatedGene1, patientsWithNonMutatedGene2);
               splitPatientsCell = {patientsM1_M2, patientsM1_NM2, patientsNM1_M2, patientsNM1_NM2};
               histData = percentage(splitPatientsCell, patientsNames);
               pvalues = pvaluesChi2(splitPatientsCell);
               table1 = {candidateGene, 'no. of patients'; '', 'percantage of patients'; '', 'p-value'}; 
               table2 = [doubleArrayToCellArray(cellfun(@length, splitPatientsCell));
                   doubleArrayToCellArray(histData);
                   doubleArrayToCellArray(pvalues)];
               table = [table; table1, table2];
               
               % ploting the histogram
               plotTitle = sprintf('Perc. of patients with mutation vs. non mutation in genes %s (1) & %s (2)',...
                    mainGeneName, candidateGene);
                figure('Name', plotTitle, 'visible', 'off');
                bar(histData);
                clear title xlabel ylabel; 
                xLabels = {'Mut1-Mut2', 'Mut1-NonMut2', 'NonMut1-Mut2', 'NonMut1-NonMut2'};
                xticklabels(xLabels);
                xtickangle(45);
                ylabel('Percentage of patients'); 
                title(plotTitle, 'FontSize', 16);
                for i = 1:length(xLabels)
                    text(i, histData(i),num2str(histData(i),2),'vert','bottom');
                end
       end
   end
   table = [tableTitle; table];
   xlswrite(fileOutput, table);
end

% Comparing average of expression of patients with mutations in main gene
% vs. no mutation.
function barPlotAvgExpressionOfMainGeneMutatedVsNot(numericData, geneNamesMut, ...
    patientsNames, patientsNamesMut, mainGeneName, mainGeneIdx)
[patientsWithMutatedGene, patientsWithNonMutatedGene] = ...
    splitPatientsByGeneMutation(mainGeneName, geneNamesMut, patientsNames, patientsNamesMut);
patientsIdxMut= contains(patientsNames, patientsWithMutatedGene);
patientsIdxNonMut= contains(patientsNames, patientsWithNonMutatedGene);
expressionMut = numericData(mainGeneIdx, patientsIdxMut);
expressionNonMut = numericData(mainGeneIdx, patientsIdxNonMut);
avgMut = mean(expressionMut);
avgNonMut = mean(expressionNonMut);
histData = [avgMut, avgNonMut]; 
[~, pvalue] = ttest2(expressionMut, expressionNonMut);

% ploting the histogram
plotTitle = sprintf('Average of %s expression of patients with mutation vs. non mutation',...
                    mainGeneName);
                figure('Name', plotTitle, 'visible', 'off');
                bar(histData);
                clear title xlabel ylabel; 
                xLabels = {'Mutated gene', 'Non mutated gene'};
                xticklabels(xLabels);
                xtickangle(45);
                ylabel('Mean expression of gene'); 
                title(plotTitle, 'FontSize', 16);
                for i = 1:length(xLabels)
                    text(i, histData(i),num2str(histData(i),2),'vert','bottom');
                end
                dim = [.7 .3 .4 .5];
                str = sprintf('p-value: %f', pvalue);
                annotation('textbox', dim, 'String', str, 'FitBoxToText','on');

end

function [patientsWithMutatedGene, patientsWithNonMutatedGene] = ...
    splitPatientsByGeneMutation(candidateGene, geneNamesMut, patientsNames, patientsNamesMut) % TODO - there's a problem here.
    geneIndices = strcmp(geneNamesMut, candidateGene);
    patientsWithMutatedGene = unique(patientsNamesMut(geneIndices)); 
    [~, patientsIdx] = intersect(patientsNames, patientsWithMutatedGene);
    patientsWithMutatedGene = patientsNames(patientsIdx)'; % col vector 
    patientsNames(patientsIdx) = [];
    patientsWithNonMutatedGene = patientsNames'; % col vector
end

function percentageOfPatientsHist = percentage(splitPatientsCell, patientsNames)
    splitPatientsHist = cellfun(@length, splitPatientsCell);
    percentageOfPatientsHist = splitPatientsHist./length(patientsNames)*100;
end


function pvalue = pvaluesChi2(splitPatientsCell)
sumOfCategories = cellfun(@length, splitPatientsCell);
cat1a = sumOfCategories(1); 
cat1b = sumOfCategories(2);
cat2a = sumOfCategories(3);
cat2b = sumOfCategories(4);
sumCat1 = cat1a + cat1b;
sumCat2 = cat2a + cat2b;
sumCatA = cat1a + cat2a;
sumCatB = cat1b + cat2b;
sumAll = sum(sumOfCategories);

expectedValue1a = expectedValue(sumCat1, sumCatA, sumAll);
expectedValue1b = expectedValue(sumCat1, sumCatB, sumAll);
expectedValue2a = expectedValue(sumCat2, sumCatA, sumAll);
expectedValue2b = expectedValue(sumCat2, sumCatB, sumAll);
expectedValues = [expectedValue1a, expectedValue1b, expectedValue2a, expectedValue2b];
chi2 = (sumOfCategories - expectedValues).^2./expectedValues;
pvalue = chi2cdf(chi2,3, 'upper');
end

function expectedVal = expectedValue(sumCat1, sumCat2, sumAll)
expectedVal = sumCat1 .* sumCat2 ./ sumAll;
end


function cellArray = doubleArrayToCellArray(rowDoubleArray)
tmp = sprintf('%f*', rowDoubleArray);
tmp(end) = [];
cellArray = regexp(tmp, '*', 'split');
end