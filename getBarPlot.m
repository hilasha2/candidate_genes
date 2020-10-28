% Comparing numericDataA vs. numericDataB.
% For each gene, we compare the gene expression of patients with 
% high gene expression/mutation (numericDataA) vs. patients with 
% low gene expression/non-mutation (numericDataB).
% P-values are calculated using two-sample t-test.
% Barplots are also plotted if printAllFigures is 'true'.
% 'analysisType' ('mutation' or 'expression')
% determinces whether we split the patients by their 
% mutation/non-mutation in the driver gene (geneName) or by 
% its expression high/low. 
function tbl = getBarPlot(numericDataA, numericDataB, geneNames, ...
    titlePlot, inputFile, geneName, outputFile, analysisType, ...
    printAllFigures, varargin)
dataMeanA = mean(numericDataA, 2);
dataMeanB = mean(numericDataB, 2);
[~, pvalues] = ttest2(numericDataA', numericDataB');

if printAllFigures == true
    barPlot(dataMeanA, dataMeanB, numericDataA, numericDataB, ...
        geneNames, titlePlot, geneName, analysisType);
    bool = (pvalues' < 0.05);

    % Plotting bar plot with significant genes only 
    numericDataSignificantA = numericDataA(bool, :);
    numericDataSignificantB = numericDataB(bool, :);
    dataMeanSignificantA = dataMeanA(bool);
    dataMeanSignificantB = dataMeanB(bool);
    geneNamesSignificant = geneNames(bool);
    titlePlotSignficant = sprintf('%s - significant genes only', titlePlot);
    barPlot(dataMeanSignificantA, dataMeanSignificantB, numericDataSignificantA, ...
        numericDataSignificantB, geneNamesSignificant, titlePlotSignficant, ...
        geneName, analysisType);
end

switch analysisType
    case 'expression'
        tbl = createTable(geneName, geneNames, inputFile, numericDataA, ...
            numericDataB, dataMeanA, dataMeanB, pvalues, analysisType, varargin);
    case 'mutation'
        tbl = createTable(geneName, geneNames, inputFile, numericDataA, ...
            numericDataB, dataMeanA, dataMeanB, pvalues, analysisType);
end

xlswrite(outputFile, tbl);
end

% Plotting a bar graph of two groups of patients. 
% One bar for every gene. Groups are split based on mutation/non-mutation
% or expression high/low in main gene.
% dataA is always high/mutated.
% dataB is low/non-mutated
function barPlot(dataMeanA, dataMeanB, numericDataA, numericDataB, ...
    geneNames, titlePlot, geneName, analysisType)
barData = cat(2, dataMeanA, dataMeanB);
numBars = length(geneNames);
figure('Name', titlePlot, 'visible', 'off');
hBar = bar(1:numBars, barData);
bar1Coor = bsxfun(@plus, hBar(1).XData, [hBar(1).XOffset]'); 
bar2Coor = bsxfun(@plus, hBar(2).XData, [hBar(2).XOffset]');
hold on
stdA = std(numericDataA, 0, 2)./sqrt(length(dataMeanA));
stdB = std(numericDataB, 0, 2)./sqrt(length(dataMeanB));
errorbar(bar1Coor, dataMeanA, stdA, '.');
errorbar(bar2Coor, dataMeanB, stdB, '.');
hold off
clear title xlabel ylabel; 
xticks(1:numBars);
xticklabels(geneNames);
xtickangle(90);
switch analysisType
    case 'expression'
        groupA = sprintf('High %s', geneName);
        groupB = sprintf('Low %s', geneName);
    case 'mutation' 
        groupA = sprintf('Mutated %s', geneName);
        groupB = sprintf('Non-Mutated %s', geneName);
end

legend(groupA, groupB, 'Location', 'northeastoutside');
xlabel('Gene Names');
ylabel('Expression of Genes');
title(titlePlot);
end

% Constructing a table for excel summarising the results of t-test2.
function tbl = createTable(geneName, geneNames, inputFile, numericDataA, ...
    numericDataB, dataMeanA, dataMeanB, pvalues, analysisType, varargin)
geneNameChar = sprintf('%s', geneName);
% First 2 rows:
tbl = {inputFile, '', '', '', '', '', ''; ...
    'data split according to', geneNameChar, '', '', '', '', ''};
% 3rd & 4th rows
switch analysisType
    case 'expression'
        numHighPatientsStr = sprintf('no. of high-%s patients', geneName);
        numLowPatientsStr = sprintf('no. of low-%s patients', geneName);
        tbl = [tbl; {'top percentage', numHighPatientsStr, ...
            'low percentage', numLowPatientsStr, '', '', ''}];
        topPerc = varargin{1};
        tbl = [tbl; {topPerc{1}, size(numericDataA, 2), 100-topPerc{1}, ...
            size(numericDataB, 2), '', '', ''}];
    case 'mutation'
        numMutatedPatientsStr = sprintf('no. of mutated-%s patients', geneName);
        numNonMutatedPatientsStr = sprintf('no. of non mutated-%s patients', geneName);
        tbl = [tbl; {numMutatedPatientsStr, numNonMutatedPatientsStr, ...
            '', '', '', '', ''}];
        tbl = [tbl; {size(numericDataA, 2), size(numericDataB, 2), ...
            '', '', '', '', ''}];
end

% 5th row
switch analysisType
    case 'expression'
        tbl = [tbl; {'Gene Names', 'mean high expression in main gene', ...
            'standard error for patients with high expression', ...
            'mean low expression in main gene', ...
            'standard error for patients with low expression', ...
            'P-value', 'Significant'}];
    case 'mutation'
        tbl = [tbl; {'Gene Names', ...
            'mean expression of group with mutation in main gene', ...
            'standard error for patients with mutation', ...
            'mean expression of group without mutaion in main gene', ...
            'standard error for patients without mutation', ...
            'P-value', 'Significant'}];
end
% The rest of the rows
stdA = std(numericDataA, 0, 2)./sqrt(length(dataMeanA));
stdB = std(numericDataB, 0, 2)./sqrt(length(dataMeanB));
% bool = (pvalues < 0.05 & pvalues <= pvalues(geneIdx));
bool = (pvalues < 0.05);
tmp = [geneNames, num2cell(dataMeanA), num2cell(stdA), ...
    num2cell(dataMeanB), num2cell(stdB), num2cell(pvalues'), num2cell(bool')];
tmp = sortrows(tmp, 6); % sort by p-values
tbl = [tbl; tmp];
end