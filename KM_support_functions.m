 classdef KM_support_functions
methods(Static)

    
% ---------- KM - Splitting cohort into 4 groups   
    
% kaplan_meier_gene1_vs_gene2_encapsulated - Top level function for KM analysis. 
% * type_of_splitting - Can be either 'mut_exp' or 'exp_exp'. 
%    'mut_exp' - Calculates KM of mutated (or not) gene 1 vs expression 
%     level of gene 2.
%    'exp_exp' - Calculates KM of expression level of gene 1 vs 
%      expression level of gene 2.
% * patientsNamesGene1GA - Names of the patients in group 1 split by gene 1. 
%    If 'mut_exp' then it's patientsNamesMut, 
%     i.e. patiens with mutations in gene 1. 
%    If 'exp_exp' then it's patientsNamesTopGene1, 
%     i.e. patients with high expression in gene 1. 
% * patientsNamesGene1GB - Names of the patients in group 2 split by gene 1. 
%    If 'mut_exp' then it's patientsNamesNotMut,
%     i.e. patients with no mutation in gene 1.
%    If 'exp_exp' then it's patientsNamesLowGene1,
%     i.e. patients with low gene expression in gene 1.
% *  Later, each group, group 1 and 2 are split into 2 groups themselves
%     based on gene 2's expression. 
% * percTop - The percentage of the top x population. 
%    Thus far percLow = 1 - percTop, but it is not necessary. 
function kaplan_meier_gene1_vs_gene2_encapsulated(filename_km, fileOutput, outputDir, ...
	geneName1, numericData, patientsNames, patientsNamesGene1GA, ...
    patientsNamesGene1GB, patientsNamesKM,geneNames, timeData, cens, ...
    timeCutOff, percTop, percLow, type_of_splitting, ...
    extra_data_splitting, print_config)
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
    extra_data_splitting, print_config);
% Switching splitting percentage: percTop = 1 - percTop, percLow = 1 - percLow.
KM_support_functions.kaplan_meier_gene1_vs_gene2(fileOutput, geneName1,...
    numericData, patientsNames, patientsNamesKM, patientsNamesGene1GA, ...
    patientsNamesGene1GB, geneNames, timeData, cens, timeCutOff, ...
    percTop, percLow, 1 - percTop, 1 - percLow, table, ...
    type_of_splitting, extra_data_splitting, print_config);

end


function kaplan_meier_gene1_vs_gene2(fileOutput, geneName1, numericData, patientsNames, ...
	patientsNamesKM, patientsNamesGene1GA, patientsNamesGene1GB, ...
    geneNames, timeData, cens, timeCutOff, percTopG1, percLowG1, percTopG2, percLowG2, ...
	table, type_of_splitting, extra_data_splitting, print_config)

% Checking that the groups patientsNamesGene1GA & patientsNamesGene1GB are
% not empty.
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
% Going over all the genes (geneName2) [besides the main gene (geneName1) if
% in exp_exp].
% Splitting data now further according to geneName2's expression and doing
% kaplan meier analyis.
for geneIdx = 1:length(geneNames)
    geneName2 = geneNames{geneIdx};
    bool_create_plot = KM_support_functions.create_plot(geneName2, print_config);

    % If we are not geneName1 OR we are in 'mut_exp' calculations.
    % in 'exp_exp' we skip geneName1 because data was already split by its
    % genge expresion. However, in 'mut_exp' data was split first by
    % mutation/no-mutation, so we can split by its expression.
	if ~isequal(geneName2, geneName1) || strcmp('mut_exp', type_of_splitting)

        if bool_create_plot
            titlePlot = KM_support_functions.create_title_plot_km_four_groups(...
                geneName1, geneName2, percTopG1, percLowG1, percTopG2,...
                percLowG2, type_of_splitting);
        end
		[~, ~, patientsNamesTopGene2, patientsNamesLowGene2] ...
			= split_gene_data_by_percentage(numericData, patientsNames, ...
            percTopG2, percLowG2, geneIdx);

        % KM just for geneName2 - data split by geneName2's expression -
        % high vs. low (we're comparing these results with the
        % co-expression KM).
        if strcmp('exp_exp', type_of_splitting)
            statsGene2Only = KM_support_functions.get_log_rank_two_grps(...
                patientsNamesKM, timeData, cens, patientsNamesTopGene2,...
                patientsNamesLowGene2, timeCutOff);
            if bool_create_plot
                titlePlot2 = KM_support_functions.create_title_plot_km_two_groups(...
                    geneName2, percTopG2, percLowG2);
                KM_support_functions.plot_kaplan_meier_two_groups(...
                    patientsNamesKM, timeData, cens, patientsNamesTopGene2,...
                    patientsNamesLowGene2, titlePlot2, statsGene2Only, geneName2);
            end
            extra_data = {statsGene2Only.p}; 
        elseif strcmp('mut_exp', type_of_splitting)
            extra_data = {};
        end
        
		% Getting 4 groups split by intersection of all the 4 groups.   
		[patientsGroup1A2Top, patientsGroup1A2Low, ...
            patientsGroup1B2Top, patientsGroup1B2Low] = ...
            KM_support_functions.interesect_patients_four_groups(...
            patientsNamesGene1GA, patientsNamesGene1GB, ...
            patientsNamesTopGene2, patientsNamesLowGene2); 
		% If one of the groups is empty then we don't bother to calculate 
		if isempty(patientsGroup1A2Top) || isempty(patientsGroup1A2Low) ||...
                isempty(patientsGroup1B2Top) || isempty(patientsGroup1B2Low)
            fprintf('One of the groups for kaplan meier analysis is empty when split by gene: %s', ...
                geneName2);
            stats.p = NaN;
		else
			% Getting log rank p-value info
			stats = KM_support_functions.get_log_rank_four_grps(...
                patientsNamesKM, timeData, cens, patientsGroup1A2Top, ...
                patientsGroup1A2Low, patientsGroup1B2Top, ...
                patientsGroup1B2Low, timeCutOff, type_of_splitting);
						
			% Plot the Kaplan Meier graph for the 4 groups
            if bool_create_plot 
                [strG1H, strG1L] = KM_support_functions.get_hr_or_lr_strings(percTopG1, 1);
                [strG2H, strG2L] = KM_support_functions.get_hr_or_lr_strings(percTopG2, 1);
                strings_LHR = {strG2H, strG2L, strG1H, strG1L};
                KM_support_functions.plot_kaplan_meier_four_groups(...
                    patientsNamesKM, timeData, cens, ...
                    patientsGroup1A2Top, patientsGroup1A2Low, ...
                    patientsGroup1B2Top, patientsGroup1B2Low, ...
                    titlePlot, stats, strings_LHR, ...
                    type_of_splitting);
            end
        end
        tableRow = KM_support_functions.createRowTable(geneName1, ...
            geneName2,  percTopG1, percLowG1, percTopG2, percLowG2, ...
            patientsGroup1A2Top, patientsGroup1A2Low, ...
            patientsGroup1B2Top, patientsGroup1B2Low, pvalueGene1, ...
            stats.p, type_of_splitting, extra_data); 
        tableData = [tableData; tableRow];
	end
end

if strcmp('mut_exp', type_of_splitting)
	pvalueCol = 9; 
elseif strcmp('exp_exp', type_of_splitting)
	pvalueCol = 11;  
end
tableData = sortrows(tableData, pvalueCol);
table = [table; tableData];
[strG2H, strG2L] = KM_support_functions.get_hr_or_lr_strings(percTopG2, 0);
tmpStr = sprintf('_G2_%s%d_%s%d.xlsx', strG2H, round(percTopG2 * 100), ...
    strG2L, round(percLowG2 * 100)); 
fileOutput = strcat(fileOutput, tmpStr); 
xlswrite(fileOutput, table);
end 

% get_log_rank_four_grps - 
% Calculates the Kaplan Meier estimator for 4 groups. 
% * 'patientsNamesHighHigh' is basically the group of patients 
% who have both high expression in gene 1 and in gene 2. 
% If mutation is involved instead of expression, then 'High' 
% would be equivalent to 'mutation', while 'Low' would be 
% equivalent to 'not mutated'.  
% * type_of_splitting - 'mut_exp' or 'exp_exp'.
function [stats] = get_log_rank_four_grps(patientsNamesKM, timeData, ...
        cens, patientsNamesHighHigh, patientsNamesHighLow, ...
        patientsNamesLowHigh, patientsNamesLowLow, timeCutOff, ...
        type_of_splitting)
    
[timeDataHH,censHH] = KM_support_functions.get_time_and_cens(...
    patientsNamesKM, timeData, cens, patientsNamesHighHigh);
[timeDataHL, censHL]= KM_support_functions.get_time_and_cens(...
    patientsNamesKM, timeData, cens,  patientsNamesHighLow);
[timeDataLH,censLH] = KM_support_functions.get_time_and_cens(...
    patientsNamesKM, timeData, cens, patientsNamesLowHigh);
[timeDataLL, censLL]= KM_support_functions.get_time_and_cens(...
    patientsNamesKM, timeData, cens,  patientsNamesLowLow);

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
[~, ~, stats] = MatSurv(TimeVar, EventVar, GroupVar, 'XLim', ...
    timeCutOff, 'NoPlot', true);
end


function plot_kaplan_meier_four_groups(patientsNamesKM, timeData, cens, ...
        patientsNamesHighHigh, patientsNamesHighLow, ...
        patientsNamesLowHigh, patientsNamesLowLow, ...
        titlePlot, stats, stringsLHR, type_of_splitting)
    
[timeDataHH,censHH] = KM_support_functions.get_time_and_cens(...
    patientsNamesKM, timeData, cens, patientsNamesHighHigh);
[timeDataHL, censHL]= KM_support_functions.get_time_and_cens(...
    patientsNamesKM, timeData, cens,  patientsNamesHighLow);
[timeDataLH,censLH] = KM_support_functions.get_time_and_cens(...
    patientsNamesKM, timeData, cens, patientsNamesLowHigh);
[timeDataLL, censLL]= KM_support_functions.get_time_and_cens(...
    patientsNamesKM, timeData, cens,  patientsNamesLowLow);

figure('Name', titlePlot, 'visible', 'off');
clear title xlabel ylabel;

hold on
labels = {};

if strcmp('mut_exp', type_of_splitting)
    str = sprintf('Mut G1 - %s G2', stringsLHR{1});
	labels = KM_support_functions.drawIfLenCensNotOne(censHH, ...
        timeDataHH, labels, str);
    str = sprintf('Mut G1 - %s G2', stringsLHR{2});
	labels = KM_support_functions.drawIfLenCensNotOne(censHL, ...
        timeDataHL, labels, str);
    str = sprintf('NotMut G1 - %s G2', stringsLHR{1});
	labels = KM_support_functions.drawIfLenCensNotOne(censLH, ...
        timeDataLH, labels, str);
    str = sprintf('NotMut G1 - %s G2', stringsLHR{2});
	labels = KM_support_functions.drawIfLenCensNotOne(censLL, ...
        timeDataLL, labels, str);
elseif strcmp('exp_exp', type_of_splitting)
    str = sprintf('%s G1 - %s G2', stringsLHR{3}, stringsLHR{1});
	labels = KM_support_functions.drawIfLenCensNotOne(censHH, ...
        timeDataHH, labels, str);
    str = sprintf('%s G1 - %s G2', stringsLHR{3}, stringsLHR{2});
	labels = KM_support_functions.drawIfLenCensNotOne(censHL, ...
        timeDataHL, labels, str);
    str = sprintf('%s G1 - %s G2', stringsLHR{4}, stringsLHR{1});
	labels = KM_support_functions.drawIfLenCensNotOne(censLH, ...
        timeDataLH, labels, str);
    str = sprintf('%s G1 - %s G2', stringsLHR{4}, stringsLHR{2});
	labels = KM_support_functions.drawIfLenCensNotOne(censLL, ...
        timeDataLL, labels, str);
end

hold off
title(titlePlot);
xlabel('Time(Month)');
ylabel('Survival Probability');
ylim([0,1]); 
legend(labels);
pvalue = stats.p;
dim = [.2 .2 .3 .3];
str = sprintf("p-value: %.8f", pvalue); 
annotation('textbox', dim, 'String', str, 'FitBoxToText','on');
end


function [patientsGroup1A2A, patientsGroup1A2B, patientsGroup1B2A, ...
        patientsGroup1B2B] = interesect_patients_four_groups(...
        patientsNamesGroup1A, patientsNamesGroup1B, ...
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
if isnan(PValueGene2)
    tableRow{end + 1} = NaN;
else
    % Whether p-value of gene2 (co-expression) is significant. 
    % True if gene 2 is statistically significant and smaller than
    % p-value of gene 1.
    if PValueGene2 < PValueGene1 && PValueGene2 < 0.05
        tableRow{end + 1} = 'TRUE';
    else
        tableRow{end + 1} = 'FALSE';
    end
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


function [timeDataNew, censNew] = get_time_and_cens(patientsNamesKM, ...
        timeData, cens, patientsNames)
[~, idx] = intersect(patientsNamesKM, patientsNames);
timeDataNew = timeData(idx);
censNew = cens(idx);
end

function titlePlot = create_title_plot_km_four_groups(geneName1, ...
        geneName2, percTopG1, percLowG1, percTopG2, percLowG2, ...
        type_of_splitting)
[strG1H, strG1L] = KM_support_functions.get_hr_or_lr_strings(percTopG1, 0);
[strG2H, strG2L] = KM_support_functions.get_hr_or_lr_strings(percTopG2, 0);

if strcmp('mut_exp', type_of_splitting)
    titlePlot = sprintf('Kaplan Meier of (G1) %s (mutated or not) vs. expression of (G2) %s (%s%d %s%d)', ...
				geneName1, geneName2, strG2H, round(percTopG2 * 100), ...
                strG2L, round(percLowG2 * 100));
elseif strcmp('exp_exp', type_of_splitting)
    titlePlot = sprintf('Kaplan Meier of (G1) %s (%s%d %s%d) vs. (G2) %s (%s%d %s%d)', ...
        geneName1, strG1H, round(percTopG1 * 100), strG1L, ...
        round(percLowG1 * 100), geneName2, strG2H, ...
        round(percTopG2 * 100), strG2L, round(percLowG2 * 100));
end
end

% Print if getGeneImages ~= {'none'}
% and if geneName2 is in the cell array getGeneImages.
function tf = create_plot(geneName2, print_config)
print_gene = 0;
if ~isempty(print_config.genes2print)
    print_gene = ismember(1,(strcmp(geneName2, print_config.genes2print)));
end    
tf = print_gene || ~print_config.print_none;
end

% In order to correctly define which type of expression-split we did we
% look at the percentage by which the patients were split according to. For
% example, if patients were split into 30% patients with high expression
% vs. 70% patients with "low" expression we represent that with:
% strH = 'H' (or 'High' if outputFullStr)
% strL = 'R' (or 'Rest').
% If the patients were split into 70% patients with "high" expression vs.
% 30% patients with low expression we represent that with:
% strH = 'R' (or 'Rest')
% strL = 'L' (or 'Low')
function [strH, strL] = get_hr_or_lr_strings(percTop, outputFullStr)
if (percTop < 0.5) && (outputFullStr)
    strH = 'High';
    strL = 'Rest';
elseif (percTop < 0.5) && (~outputFullStr)
    strH = 'H';
    strL = 'R';
elseif (percTop >= 0.5) && (outputFullStr) 
    strH = 'Rest';
    strL = 'Low';
elseif (percTop >=0.5) && (~outputFullStr)
    strH = 'R';
    strL = 'L';
end
end
%% ------------ KM functions for splitting a cohort into two groups

% * 'patientsNamesTop' - A group of patients with high gene expression
% according to a certain gene. 
% * 'patientsNamesLow' - patients with low gene expression in said gene. 
function [stats] = get_log_rank_two_grps(patientsNamesKM, timeData, cens, ...
        patientsNamesTop, patientsNamesLow, timeCutOff)
[timeDataTop, censTop] = KM_support_functions.get_time_and_cens(...
    patientsNamesKM, timeData, cens, patientsNamesTop);
[timeDataLow, censLow]= KM_support_functions.get_time_and_cens(...
    patientsNamesKM, timeData, cens, patientsNamesLow);

TimeVar = [timeDataTop; timeDataLow];
% event = 1, cens = 0 - so the cens' matrices have to be reversed.
EventVar = [~censTop; ~censLow]; 
g1 = cell(size(timeDataTop));
g1(:) = {'Group 1'};
g2 = cell(size(timeDataLow));
g2(:) = {'Group 2'};
GroupVar = [g1; g2];
[~, ~, stats] = MatSurv(TimeVar, EventVar, GroupVar, 'XLim', ...
    timeCutOff, 'NoPlot', true);
end


function plot_kaplan_meier_two_groups(patientsNamesKM, timeData, cens, ...
        patientsNamesTop, patientsNamesLow, titlePlot, stats, geneName)
    
[timeDataTop, censTop] = KM_support_functions.get_time_and_cens(...
    patientsNamesKM, timeData, cens, patientsNamesTop);
[timeDataLow, censLow]= KM_support_functions.get_time_and_cens(...
    patientsNamesKM, timeData, cens, patientsNamesLow);
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


function titlePlot = create_title_plot_km_two_groups(geneName, ...
        percTop, percLow)
[strH, strL] = KM_support_functions.get_hr_or_lr_strings(percTop, 0);
titlePlot = sprintf('Kaplan Meier by gene expression of %s (%s%d%% vs. %s%d%%)', ...
                geneName, strH, round(percTop * 100), strL, ...
                round(percLow * 100));
end
end 
end