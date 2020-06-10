classdef KM_support_functions
methods(Static)

    
% ---------- KM - Splitting cohort into 4 groups   
    
% kaplan_meier_gene1_vs_gene2_encapsulated - Top level function for KM analysis. 
% * type_of_splitting - Can be either 'mut_exp' or 'exp_exp'. 
% 'mut_exp' - Calculates KM of mutated (or not) gene 1 vs expression level of gene 2.
% 'exp_exp' - Calculates KM of expression level of gene 1 vs expression level of gene 2.
% *  patientsNamesGene1GA - Names of the patients in group 1 split by gene 1. If 'mut_exp' then it's patientsNamesMut, 
% if 'exp_exp' then it's patientsNamesTopGene1. 
% * patientsNamesGene1GB - Names of the patients in group 2 split by gene 1. If 'mut_exp' then it's patientsNamesNotMut, 
% if 'exp_exp' then it's patientsNamesLowGene1. 
% * percTop - The percentage of the top x population. Thus far percLow = 1 - percTop, but it is not necessary. 
function kaplan_meier_gene1_vs_gene2_encapsulated(filename_km, fileOutput, outputDir, ...
	geneName1, numericData, patientsNames, patientsNamesGene1GA, patientsNamesGene1GB, patientsNamesKM,...
    geneNames, timeData, cens, timeCutOff, percTop, percLow, type_of_splitting, extra_data_splitting)
fileOutput = fullfile(outputDir, fileOutput); 

% Initiating table to add to an excel file
geneName1 = sprintf('%s', geneName1);
tmpStr = sprintf('Significant in comparison to %ss p-value', geneName1);
if strcmp('mut_exp', type_of_splitting)
	filename_mut = extra_data_splitting{1}; 
	table = {'clinical data file:', filename_km, 'mutations file:', ...
     filename_mut, 'Main gene',geneName1,'','','',''; 'Gene split by mutation',...
     'Gene split by expression', 'High perc. gene 2', 'Low perc. gene 2', ...
     'Initial no. of patients MH', 'Initial no. of patients ML',...
     'Initial no. of patients NMH', 'Initial no. of patients NML',...
     'p-value', tmpStr};
elseif strcmp('exp_exp', type_of_splitting)
	pvalueGene1 = extra_data_splitting{1};
	table = {filename_km, '','','','','','','','','','','', '', '';...
    'Gene Name', 'High perc.', 'Low perc.', 'Initial no. of patients - High',...
    'Initial no. of patients - Low', 'p-value', '', '', '', '', '', '', '', '';
    geneName1, percTop * 100, round(percLow * 100), ...
    length(patientsNamesGene1GA), length(patientsNamesGene1GB), pvalueGene1,...
    '', '', '', '', '', '', '', '';
    'Gene 1 Name', 'Gene 2 Name', 'High perc. gene 1', 'Low perc. gene 1',...
    'High perc. gene 2', 'Low perc. gene 2',...
    'Initial no. of patients HH', 'Initial no. of patients HL',...
    'Initial no. of patients LH', 'Initial no. of patients LL',...
    'p-value of co-expression', tmpStr, 'p-value of expression', ...
    'Significant (p-value expression)'};
end


KM_support_functions.kaplan_meier_gene1_vs_gene2(fileOutput, geneName1, ...
    numericData, patientsNames, patientsNamesKM, patientsNamesGene1GA, ...
    patientsNamesGene1GB, geneNames, timeData, cens, timeCutOff,...
    percTop, percLow, percTop, percLow, table, type_of_splitting, ...
    extra_data_splitting);
% Switching splitting percentage: percTop = 1 - percTop, percLow = 1 - percLow.
KM_support_functions.kaplan_meier_gene1_vs_gene2(fileOutput, geneName1,...
    numericData, patientsNames, patientsNamesKM, patientsNamesGene1GA, ...
    patientsNamesGene1GB, geneNames, timeData, cens, timeCutOff, ...
    percTop, percLow, 1 - percTop, 1 - percLow, table, ...
    type_of_splitting, extra_data_splitting);

end


function kaplan_meier_gene1_vs_gene2(fileOutput, geneName1, numericData, patientsNames, ...
	patientsNamesKM, patientsNamesGene1GA, patientsNamesGene1GB, ...
    geneNames, timeData, cens, timeCutOff, percTopG1, percLowG1, percTopG2, percLowG2, ...
	table, type_of_splitting, extra_data_splitting)
if isempty(patientsNamesGene1GA) || isempty(patientsNamesGene1GB) 
    fprintf('Spliting cohort by gene did not succeed: %s\n', geneName1)
    if strcmp('mut_exp', type_of_splitting)
        fprintf('Most probably there was not any mutation in said gene in cohort\n');
    end
    return;
end
if strcmp('mut_exp', type_of_splitting)
    % P-value of gene2 is compared to gene1's p-value. 
    % Since there's no relevance for KM in 'mut_exp' for gene1
    % it is denoted by 2 an impossible value for p-value
	pvalueGene1 = 2; 
elseif strcmp('exp_exp', type_of_splitting)
	pvalueGene1 = extra_data_splitting{1}; 
end
	
tableData = {};
% Going over all the genes besides the main gene (geneName1)
for geneIdx = 1:length(geneNames)
    geneName2 = geneNames{geneIdx};
	if ~isequal(geneName2, geneName1) || strcmp('mut_exp', type_of_splitting)
	
		if strcmp('mut_exp', type_of_splitting)
		titlePlot = sprintf('(G1) %s (mutated or not) vs. (G2) %s (H%d %dL) expression level - survial of cohort', ...
				geneName1, geneName2, round(percTopG2 * 100), round(percLowG2 * 100));
		elseif strcmp('exp_exp', type_of_splitting)
		titlePlot = sprintf('(G1) %s (H%d L%d) vs. (G2) %s (H%d %dL) expression - survial of cohort', ...
				geneName1, round(percTopG1 * 100), round(percLowG1 * 100), geneName2, ...
                round(percTopG2 * 100), round(percLowG2 * 100));
		end
		
		[~, ~, patientsNamesTopGene2, patientsNamesLowGene2] ...
			= split_gene_data_by_percentage(numericData, patientsNames, percTopG2, percLowG2, geneIdx);
        % KM just for geneName2 - data split by geneName2's expression -
        % high vs. low (we're comparing these results with the
        % co-expression KM),
        if strcmp('exp_exp', type_of_splitting)
            titlePlot2 = sprintf('Kaplan Meier by gene expression of %s (H%d%% vs. L%d%%)',...
                geneName2, round(percTopG2 * 100), round(percLowG2 * 100));
            statsGene2Only = KM_support_functions.get_log_rank_two_grps(...
                patientsNamesKM, timeData, cens, patientsNamesTopGene2,...
                patientsNamesLowGene2, timeCutOff, titlePlot2);
            KM_support_functions.plot_kaplan_meier_two_groups(...
                patientsNamesKM, timeData, cens, patientsNamesTopGene2,...
                patientsNamesLowGene2, titlePlot2, statsGene2Only, geneName2);
            extra_data = {statsGene2Only.p}; 
        elseif strcmp('mut_exp', type_of_splitting)
            extra_data = {};
        end
        
		% Getting 4 groups split by intersection of all the 4 groups.   
		[patientsGroup1A2Top, patientsGroup1A2Low, patientsGroup1B2Top, patientsGroup1B2Low] ...
			= KM_support_functions.interesect_patients_by_mut_and_exp(patientsNamesGene1GA, ...
			patientsNamesGene1GB, patientsNamesTopGene2, patientsNamesLowGene2); 
		% If one of the groups is empty then we don't bother to calculate 
		if isempty(patientsGroup1A2Top) || isempty(patientsGroup1A2Low) ||...
                isempty(patientsGroup1B2Top) || isempty(patientsGroup1B2Low)
			fprintf('One of the groups for kaplan meier analysis is empty when split by gene: %s', geneName2); 
			return;
			% Previously:  
			% If the patients' names of the high/low gene expression values are
			% the same both in gene 1 and the other gene, then we compare only 2
			% groups and not 4:
			%if isempty(highG1_lowG2) || isempty(lowG1_highG2) 
			%   stats = get_log_rank_two_grps(patientsNamesKM, timeData, cens, highG1_highG2,...
			%       lowG1_lowG2, timeCutOff, titlePlot);
			%   plot_kaplan_meier_two_groups(patientsNamesKM, timeData, cens, highG1_highG2,...
			%       lowG1_lowG2, titlePlot, stats, geneName1);
		else
			% Getting log rank p-value info
			stats = KM_support_functions.get_log_rank_four_grps(patientsNamesKM, timeData, ...
						cens, patientsGroup1A2Top, patientsGroup1A2Low, patientsGroup1B2Top, ...
						patientsGroup1B2Low, timeCutOff, titlePlot, type_of_splitting);
						
			% Plot the Kaplan Meier graph for the 4 groups
			KM_support_functions.plot_kaplan_meier_four_groups(patientsNamesKM, timeData, cens, ...
						patientsGroup1A2Top, patientsGroup1A2Low, patientsGroup1B2Top, ...
						patientsGroup1B2Low, titlePlot, stats, geneName1, geneName2, type_of_splitting);
			tableRow = KM_support_functions.createRowTable(geneName1, geneName2,  percTopG1, ...
						percLowG1, percTopG2, percLowG2, patientsGroup1A2Top, patientsGroup1A2Low, ...
						patientsGroup1B2Top, patientsGroup1B2Low, pvalueGene1, stats.p,...
                        type_of_splitting, extra_data); 
			tableData = [tableData; tableRow];
		end
	end
end

if strcmp('mut_exp', type_of_splitting)
	pvalueCol = 9; 
elseif strcmp('exp_exp', type_of_splitting)
	pvalueCol = 11;  
end
tableData = sortrows(tableData, pvalueCol);
table = [table; tableData];
tmpStr = sprintf('_G2_H%d_L%d.xlsx', round(percTopG2 * 100), round(percLowG2 * 100)); 
fileOutput = strcat(fileOutput, tmpStr); 
xlswrite(fileOutput, table);
end 

% get_log_rank_four_grps - Calculates the Kaplan Meier estimator for 4 groups. 
% * 'patientsNamesHighHigh' is basically the group of patients who have both high expression in gene 1 
% and in gene 2. If mutation is involved instead of expression, then 'High' would be equivalent to 
% 'mutation', while 'Low' would be equivalent to 'not mutated'.  
% * type_of_splitting - values are 'true'/'false'. If 'true' it means that we do kaplan meier calculation 
% that includes a group that was split by mutation, if 'false' then groups were split only by expression level.
function [stats] = get_log_rank_four_grps(patientsNamesKM, timeData, cens, patientsNamesHighHigh,...
    patientsNamesHighLow, patientsNamesLowHigh, patientsNamesLowLow, timeCutOff, titlePlot, type_of_splitting)
[timeDataHH,censHH] = KM_support_functions.get_time_and_cens(patientsNamesKM, timeData, cens, patientsNamesHighHigh);
[timeDataHL, censHL]= KM_support_functions.get_time_and_cens(patientsNamesKM, timeData, cens,  patientsNamesHighLow);
[timeDataLH,censLH] = KM_support_functions.get_time_and_cens(patientsNamesKM, timeData, cens, patientsNamesLowHigh);
[timeDataLL, censLL]= KM_support_functions.get_time_and_cens(patientsNamesKM, timeData, cens,  patientsNamesLowLow);

TimeVar = [timeDataHH; timeDataHL; timeDataLH; timeDataLL];
EventVar = [~censHH; ~censHL; ~censLH; ~censLL];
g1 = cell(size(timeDataHH));
g2 = cell(size(timeDataHL));
g3 = cell(size(timeDataLH));
g4 = cell(size(timeDataLL));

if strcmp('mut_exp', type_of_splitting)
	g1(:) = {'Mut-High'};
	g2(:) = {'Mut-Low'};
	g3(:) = {'NotMut-High'};
	g4(:) = {'NotMut-High'};
elseif strcmp('exp_exp', type_of_splitting) 
	g1(:) = {'High-High'};
	g2(:) = {'High-Low'};
	g3(:) = {'Low-High'};
	g4(:) = {'Low-High'};
end 
GroupVar = [g1; g2; g3; g4];
[~, ~, stats] = MatSurv(TimeVar, EventVar, GroupVar, 'XLim', timeCutOff, 'Title', titlePlot, 'NoPlot', true);
end


function plot_kaplan_meier_four_groups(patientsNamesKM, timeData, cens, patientsNamesHighHigh,...
    patientsNamesHighLow, patientsNamesLowHigh, patientsNamesLowLow, titlePlot, stats, geneName1,...
	geneName2, type_of_splitting)
[timeDataHH,censHH] = KM_support_functions.get_time_and_cens(patientsNamesKM, timeData, cens, patientsNamesHighHigh);
[timeDataHL, censHL]= KM_support_functions.get_time_and_cens(patientsNamesKM, timeData, cens,  patientsNamesHighLow);
[timeDataLH,censLH] = KM_support_functions.get_time_and_cens(patientsNamesKM, timeData, cens, patientsNamesLowHigh);
[timeDataLL, censLL]= KM_support_functions.get_time_and_cens(patientsNamesKM, timeData, cens,  patientsNamesLowLow);
figure('Name', titlePlot, 'visible', 'off');
clear title xlabel ylabel;

hold on
labels = {};
if strcmp('mut_exp', type_of_splitting)
	labels = KM_support_functions.drawIfLenCensNotOne(censHH, timeDataHH, labels, 'Mut G1 - High G2');
	labels = KM_support_functions.drawIfLenCensNotOne(censHL, timeDataHL, labels, 'Mut G1 - Low G2');
	labels = KM_support_functions.drawIfLenCensNotOne(censLH, timeDataLH, labels, 'NotMut G1 - High G2');
	labels = KM_support_functions.drawIfLenCensNotOne(censLL, timeDataLL, labels, 'NotMut G1 - Low G2');
elseif strcmp('exp_exp', type_of_splitting)
	labels = KM_support_functions.drawIfLenCensNotOne(censHH, timeDataHH, labels, 'High G1 - High G2');
	labels = KM_support_functions.drawIfLenCensNotOne(censHL, timeDataHL, labels, 'High G1 - Low G2');
	labels = KM_support_functions.drawIfLenCensNotOne(censLH, timeDataLH, labels, 'Low G1 - High G2');
	labels = KM_support_functions.drawIfLenCensNotOne(censLL, timeDataLL, labels, 'Low G1 - Low G2');
end

hold off
title(titlePlot);
xlabel('Time(Month)');
ylabel('Survival Probability');
ylim([0,1]); 
legend(labels);
pvalue = stats.p;
dim = [.2 .2 .3 .3];
initialNumHH = length(patientsNamesHighHigh);
initialNumHL = length(patientsNamesHighLow);
initialNumLH = length(patientsNamesLowHigh);
initialNumLL = length(patientsNamesLowLow);
% str = sprintf("Initial no. of patients with:\nhigh %s and high %s expression: %d\nhigh %s and low %s expression: %d\nlow %s and high %s expression: %d\nlow %s and low %s expression: %d\np-value: %.8f", ...
%     geneName1, geneName2, initialNumHH, geneName1, geneName2, initialNumHL, ...
%     geneName1, geneName2, initialNumLH, geneName1, geneName2, initialNumLL, pvalue);
str = sprintf("p-value: %.8f", pvalue); 
annotation('textbox', dim, 'String', str, 'FitBoxToText','on');
end


function [patientsGroup1A2A, patientsGroup1A2B, patientsGroup1B2A, patientsGroup1B2B] ...
    = interesect_patients_by_mut_and_exp(patientsNamesGroup1A, patientsNamesGroup1B, ...
    patientsNamesGroup2A, patientsNamesGroup2B)
patientsGroup1A2A = intersect(patientsNamesGroup1A, patientsNamesGroup2A); 
patientsGroup1A2B = intersect(patientsNamesGroup1A, patientsNamesGroup2B);
patientsGroup1B2A = intersect(patientsNamesGroup1B, patientsNamesGroup2A);
patientsGroup1B2B = intersect(patientsNamesGroup1B, patientsNamesGroup2B);
end


function tableRow = createRowTable(GeneName1, GeneName2,  percTopG1, ...
percLowG1, percTopG2, percLowG2, patientsGroup1A2A, patientsGroup1A2B, ...
patientsGroup1B2A, patientsGroup1B2B, PValueGene1, PValueGene2, ...
type_of_splitting, extra_data) 
tableRow = {};
tableRow{end + 1} = GeneName1;
tableRow{end + 1} = GeneName2;
if strcmp('exp_exp', type_of_splitting)
    tableRow{end + 1} = round(percTopG1 * 100);
    tableRow{end + 1} = round(percLowG1 * 100);
end 
tableRow{end + 1} = round(percTopG2 * 100);
tableRow{end + 1} = round(percLowG2 * 100);
tableRow{end + 1} = length(patientsGroup1A2A);
tableRow{end + 1} = length(patientsGroup1A2B);
tableRow{end + 1} = length(patientsGroup1B2A);
tableRow{end + 1} = length(patientsGroup1B2B);
tableRow{end + 1} = PValueGene2;
% Whether p-value of gene2 (co-expression) is significant. 
% True if gene 2 is statistically significant and smaller than p-value of
% gene 1.
if PValueGene2 < PValueGene1 && PValueGene2 < 0.05
    tableRow{end + 1} = 'TRUE';
else
    tableRow{end + 1} = 'FALSE';
end
if strcmp('exp_exp', type_of_splitting)
    PValueGene2_NotCoexpression = extra_data{1};
    tableRow{end + 1} = PValueGene2_NotCoexpression; 
    if PValueGene2_NotCoexpression < 0.05
        tableRow{end + 1} = 'TRUE';
    else
        tableRow{end + 1} = 'FALSE';
    end
end
end

function legend = drawIfLenCensNotOne(cens, timeData, legend, label)
if length(cens) > 1
ecdf(timeData, 'censoring', cens, 'function', 'survivor');
legend = [legend, label];
end
end


function [timeDataNew, censNew] = get_time_and_cens(patientsNamesKM, timeData, cens, patientsNames)
[~, idx] = intersect(patientsNamesKM, patientsNames);
timeDataNew = timeData(idx);
censNew = cens(idx);
end



%% ------------ KM functions for splitting a cohort into two groups

% * 'patientsNamesTop' - A group of patients with high gene expression
% according to a certain gene. 
% * 'patientsNamesLow' - patients with low gene expression in said gene. 
function [stats] = get_log_rank_two_grps(patientsNamesKM, timeData, cens, patientsNamesTop,...
    patientsNamesLow, timeCutOff, titlePlot)
[timeDataTop, censTop] = KM_support_functions.get_time_and_cens(patientsNamesKM, timeData, cens,...
						patientsNamesTop);
[timeDataLow, censLow]= KM_support_functions.get_time_and_cens(patientsNamesKM, timeData, cens, ...
						patientsNamesLow);

TimeVar = [timeDataTop; timeDataLow];
EventVar = [~censTop; ~censLow]; % event = 1, cens = 0 - so the cens' matrices have to be reversed.
g1 = cell(size(timeDataTop));
g1(:) = {'Group 1'};
g2 = cell(size(timeDataLow));
g2(:) = {'Group 2'};
GroupVar = [g1; g2];
[~, ~, stats] = MatSurv(TimeVar, EventVar, GroupVar, 'XLim', timeCutOff, 'Title', titlePlot, 'NoPlot', true);
end


function plot_kaplan_meier_two_groups(patientsNamesKM, timeData, cens, patientsNamesTop,...
    patientsNamesLow, titlePlot, stats, geneName)
[timeDataTop,censTop] = KM_support_functions.get_time_and_cens(patientsNamesKM, timeData, cens, ... 
						patientsNamesTop);
[timeDataLow, censLow]= KM_support_functions.get_time_and_cens(patientsNamesKM, timeData, cens, ...
						patientsNamesLow);
figure('Name', titlePlot, 'visible', 'off');
clear title xlabel ylabel;
ecdf(timeDataTop, 'censoring', censTop, 'function', 'survivor');
hold on
ecdf(timeDataLow, 'censoring', censLow, 'function', 'survivor');
hold off
title(titlePlot);
xlabel('Time(Month)');
ylabel('Survival Probability');
ylim([0,1]); 
legend('high', 'low');
pvalue = stats.p;
HR = stats.HR;
low95 = stats.up95;
up95 = stats.low95;
dim = [.2 .2 .3 .3];
initialNumTop = length(patientsNamesTop);
initialNumLow = length(patientsNamesLow);
str = sprintf("Initial no. of patients with high %s: %d\nInitial no. of patients with low %s: %d\np-value: %f\nHR: %f\nlow95: %f\nup95: %f",...
    geneName, initialNumTop, geneName, initialNumLow, pvalue, HR, low95, up95);
annotation('textbox', dim, 'String', str, 'FitBoxToText','on');
end
end
end 