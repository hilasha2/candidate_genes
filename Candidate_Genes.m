% Some info before using this file:
% * TODO - calculations in the code that need to be changed.
% * predicate - Places in the code where we assume some conditions the input
% provides, such as, we assume the patients names in the expression data
% file are placed in the 1st column.


%% Input Parameters

% * filename - fullpath of the gene expressions excel file 
% (I call it gene expression, but it can be rna expression, cna expression etc. as long as it 
% has the same table structure as the gene expression table).
% * sheet_name - current sheet to work on in filename. Either the exact name
% or number.
% * filename_km - fullpath of the clinical data with survival analysis.
% * do_genes_expression_calculations - Whether to do gene expressions
% calculations from filename.
% * do_km_analysis - Whether to do kaplan meier analysis based on driver gene and
% candidate genes expressions according to survival data found in filename_km,
% * do_subtype_histograms - whether to do cancer's subtypes histogram calculations. 
% * use_days_not_months - whether to use the 'days' column (OS_DAYS)  in
% filname_km instead of the months column (OS_MONTHS). Only reason to use - more data
% in the 'days' column.
% * colTimeData - Column in filename_km where the overall survival in months
% is (e.g OS_MONTHS).
% * colSurvivalStatus - Column in filename_km where the overall survival
% status is (e.g OS_STATUS or VITAL_STATUS).
% * colCancerSubtypes - Column in filename_km where the cancer subtypes is
% located, (e.g. CLAUDIN_SUBTYPE). It's relevant if do_subtype_histograms
% is true. 
% * colPatientsNamesKM -  Column in filename_km where the patients' names /
% IDs is (PATIENT_ID).
% * livingStatusStr - How the "living/alive/etc." is written in the
%  survival status column. E.g - "LIVING", "Alive", etc. 
% * deceasedButNotCancerStr - Add if there's a status in the survival
% status column which indicates that the patient died, but not due to
% cancer. For example: "Died of Other Causes".
% * timeCutOff - overall period in months for Kaplan Meier analysis. Use
% 120 if you don't know what to set.
% * outputDir - Directory where all the outputs are saved.
% * extraClinicalDataCellMatrix - A n*2 matrix, where the first column
% represents the name of the specific clinical data, like ETHNICITY, 
% and the 2nd column represents the column number in filename_km. E.g
% -{'sex', 2; 'cns status', 8; 'testicular involvement', 9}; -  not yet
% applicable
% * do_mutation_analysis - whether to do mutation analysis.
% * colMutationsTypes - column in filename_mut with the different mutations
% classifications (e.g. - Variant_Classification)
% * colPatientsNamesMutation - column in filename_mut with patients' names
% (pick the shorter names, e.g. choose the col called Tumor_Sample_Barcode
% and not Matched_Norm_Sample_Barcode).
% * colGeneNamesMut - column in filename_mut of the genes' names
% (e.g. - Hugo_Symbol). 
% * filename_mut - full path of the mutation data excel file. 
% * gene_nomenclature - Either 'hugo', 'entrez', or 'both'. Basically means
% which nomenclature is used for the gene names, i.e. Hugo symbol will have
% 'MET' while 'Entrez_gene_id' will have 4233. If both is passed then that
% means the first column is 'Hugo_symbol' and the second is
% 'Entrez_gene_id'.
% * gene_name - Indicates which gene the analysis should focus on. Default
% is 'MET'. Either 'string' or 'char'. 
% * use_cbioportal_not_lab_std - Whether to use the standard of data
% sources as in cbiportal or the standard that will be decided in the lab.
% 1 is cbioportal, 0 is lab standard.
% * rowMutationsData - The row in filename_mut where the data begins.
% * getGeneImages - Char cell array of gene names whose output images are desired, i.e.
% {'ODC1', 'BCAT1', 'COL29A1'} (do not include the driver gene). 
% If all gene images are desired then write - 
% {'all'}, and if no image is desired write {'none'}. 
% Predicate - no gene is called 'none' or 'all'.
% * do_km_by_mutation_and_expression - boolean. Whether to calculate KM of the data
% where expression levels and mutations of the genes are combined.
% * do_km_by_clustering - Whether to do Kaplan Meier calculations based on
% clustering of filename. Needs filename_km for survival data. 

function Candidate_Genes(filename, filename_km, sheet_name, ...
    do_genes_expression_calculations, do_km_analysis, ... 
    do_subtype_histograms, use_cbioportal_not_lab_std,...
    gene_nomenclature, gene_name, colPatientsNamesKM,...
    colTimeData, colSurvivalStatus, colCancerSubtypes, livingStatusStr, ...
    deceasedButNotCancerStr, timeCutOff, use_days_not_months, outputDir,...
    extraClinicalDataCellMatrix, do_mutation_analysis, colMutationsTypes,...
    colPatientsNamesMutation, colGeneNamesMut, filename_mut, rowMutationsData, ...
    getGeneImages, do_km_by_mutation_and_expression, do_km_by_clustering)
disp('Code running will print "Done!" when finished');
% Checking some variables 
assert_msg = sprintf("Length of variable 'getGeneImages' is 0." + ...
    "\nShould be {'none'}, {'all'}, or cell array with gene names.");
assert((iscell(getGeneImages))&&(~isempty(getGeneImages)), ...
    assert_msg);
% Close any opened excel file.
!taskkill -f -im EXCEL.exe
%% -------- Reading the data

%  -------- Reading candidate genes data.
    % Obeservations - patients. 
    % Variables - genes.
    % numericData - A numeric matrix of gene expressions data.
    % patientsNames - A character cell array of the patients IDs.
    % geneNames - A character cell array with all the genes.

    first_row = readmatrix(filename, 'FileType', 'spreadsheet', 'Range', ...
        '1:1', 'OutputType', 'char', 'Sheet', 1, 'UseExcel', true);
    first_col = readmatrix(filename, 'FileType', 'spreadsheet', 'Range', ...
        'A:A', 'OutputType', 'char', 'Sheet', 1, 'UseExcel', true);
    
    
    if strcmp(gene_nomenclature, 'entrez') 
        numericData = readmatrix(filename, 'FileType', 'spreadsheet', ...
            'Range', 'B2', 'OutputType', 'double', 'Sheet', 1, 'UseExcel', true); 
        patientsNames = first_row(2:end);
        geneNames = first_col;
    elseif strcmp(gene_nomenclature, 'hugo')
        numericData = readmatrix(filename, 'FileType', 'spreadsheet', ...
            'Range', 'B2', 'OutputType', 'double', 'Sheet', 1, 'UseExcel', true);
        patientsNames = first_row(2:end);
        geneNames = first_col(2:end);
    elseif strcmp(gene_nomenclature, 'both')
        numericData = readmatrix(filename, 'FileType', 'spreadsheet', ...
            'Range', 'C2', 'OutputType', 'double', 'Sheet', 1, 'UseExcel', true);
        patientsNames = first_row(3:end);
        geneNames = first_col(2:end);
    end


    % If the numeric data has rows with NaN then:
    % Calculate without the rows (genes) which have NaN.
    nanIndices = isnan(numericData); 
    if ismember(1, nanIndices)
        [rowsWithNaN, ~] = find(nanIndices);
        rowsWithNaN = unique(rowsWithNaN);

        numericData(rowsWithNaN, :) = [];
        geneNames(rowsWithNaN) = [];
    end

    % main gene's index in geneNames
    [~, geneIdx] = ismember(gene_name ,geneNames); % predicate - gene_name is correct 
    if geneIdx == 0
        disp("The gene name that was given cannot be found in database");
        return
    end
    gene_name = char(string(gene_name));
    
    % Normalization of gene expression data - mean 0, stdv - 1 
    % Zscore (normalization) on rows (per gene) - 
    % Checked how cbioportal.org normalize their gene expression table and
    % it's per gene. Gene expressions should be the same if you got data
    % from cbioportal which was already zscored.
    % TODO - check whether data from NIH should be treated differently. 
    numericData = zscore(numericData, 0, 1); % Might contain 'NaN' data. 

% -------- Reading clinical data
if do_km_analysis || do_km_by_mutation_and_expression || ...
        do_km_by_clustering || do_subtype_histograms

    [numericDataKM, txtDataKM] = xlsread(filename_km, 1);
    
    % timeData - Time failure in Kaplan Meier analysis, in our case:
    % Overall survival in months since the initial diagnosis.
    % survivalStatus - Survival status - e.g - Living, Deceased - due to cancer or not. 
    % patientsNamesKM - Patients' IDs in filename_km.

    if use_cbioportal_not_lab_std    
        timeData = numericDataKM(3:end, colTimeData - 1);
        survivalStatus = txtDataKM(6:end, colSurvivalStatus);
        patientsNamesKM = txtDataKM(6:end, colPatientsNamesKM); 
    else
        timeData = numericDataKM(3:end, colTimeData - 1);
        survivalStatus = txtDataKM(4:end, colSurvivalStatus);
        patientsNamesKM = txtDataKM(4:end, colPatientsNamesKM); 
    end
    
    notEmptyIndicesTime = ~isnan(timeData);
    notEmptyIndicesSurvival = ~findEmptyString(survivalStatus);

    bool = notEmptyIndicesTime & notEmptyIndicesSurvival;
    timeData = timeData(bool);
    if use_days_not_months
        % months are rounded up from days survival data.
        timeData = ceil((timeData./365)*12);
    end
    survivalStatus = survivalStatus(bool); 
    patientsNamesKM = patientsNamesKM(bool);
    
    living = repmat(livingStatusStr, size(survivalStatus, 1), 1);
    died_not_of_breast_cancer = repmat(deceasedButNotCancerStr, size(survivalStatus, 1), 1);
    boolLiving = cellfun(@isequal, survivalStatus, living);
    boolDied = cellfun(@isequal, survivalStatus, died_not_of_breast_cancer);
    % censoring data.
    cens = (boolDied == 1) | (boolLiving == 1) | (timeData > timeCutOff); 
    
    patientsNames = patientsNamesFromDiffDatasets(patientsNames, patientsNamesKM, use_cbioportal_not_lab_std);
end

% -------- Reading Mutation Data

if do_mutation_analysis || do_km_by_mutation_and_expression
    [numericDataMut, txtDataMut] = xlsread(filename_mut); 
     
    % predicate
    if ~use_cbioportal_not_lab_std
        if strcmp(gene_nomenclature, 'entrez') 
            geneNamesMut = numericDataMut((rowMutationsData - 1):end, colGeneNamesMut);
        elseif strcmp(gene_nomenclature, 'hugo') || strcmp(gene_nomenclature, 'both')
            geneNamesMut = txtDataMut(rowMutationsData:end, colGeneNamesMut);
        end
    else
        geneNamesMut = txtDataMut(rowMutationsData:end, colGeneNamesMut);
    end
    
    mutationsTypes = txtDataMut(rowMutationsData:end, colMutationsTypes); 
    patientsNamesMut = txtDataMut(rowMutationsData:end, colPatientsNamesMutation);     
    
    if strcmp(gene_nomenclature, 'entrez') 
        emptyIdxGenes = isnan(geneNamesMut);
    else
        emptyIdxGenes = geneNamesMut == "";
    end
    emptyIdxPatients = findEmptyString(patientsNamesMut);
    emptyMutTypes = findEmptyString(mutationsTypes);
    nanIndices = emptyIdxPatients | emptyMutTypes | emptyIdxGenes;
    
    geneNamesMut(nanIndices) = [];
    mutationsTypes(nanIndices) = [];
    patientsNamesMut(nanIndices) = [];
    
    patientsNamesMut = patientsNamesFromDiffDatasets(patientsNamesMut, patientsNames, use_cbioportal_not_lab_std);
end

close all hidden;

[structNumericDataH, structNumericDataL, ...
    structPatientsH, structPatientsL] = getSplitData(numericData, patientsNames, geneIdx);
%% -------- Handling figure printing configurations
% Do we print all figures? none? or just some?
print_all = unique(eq(getGeneImages{1}, "all"));
print_none =  unique(eq(getGeneImages{1}, "none"));

%% -------- Genes expressions calculations
if do_genes_expression_calculations
    expression_data_analysis(filename, outputDir, numericData, geneNames, patientsNames, ...
    geneIdx, gene_name, structNumericDataH, structNumericDataL, structPatientsH, print_all);
end

%% --------- Kaplan Meier analysis by expression
if do_km_analysis
    kaplan_meier_data_analysis(filename_km, outputDir, numericData, geneNames, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff, structPatientsH, structPatientsL, ...
    gene_name);
end 

%% --------- Subtype calculations
if do_subtype_histograms
     cancer_subtypes_analysis(filename_km, outputDir, numericData, patientsNames, geneNames,...
     txtDataKM, patientsNamesKM, geneIdx, colCancerSubtypes, use_cbioportal_not_lab_std, gene_name);
end

%% --------- Mutation calculations
if do_mutation_analysis
    mutation_analysis(filename_mut, numericData, geneNames, geneNamesMut, ...
        patientsNames, patientsNamesMut, gene_name, geneIdx, outputDir);
end

%% ---------- Mutation and Kaplan Meier calculations
if do_km_by_mutation_and_expression
   kaplan_meier_and_mutation_analysis(filename_km, filename_mut, outputDir, gene_name, numericData, ...
    patientsNames, patientsNamesMut, patientsNamesKM, geneNames, geneNamesMut, timeData, cens, timeCutOff); 
end

%% ---------- Kaplan Meier analysis by clustering
if do_km_by_clustering
       kaplan_meier_by_clustering(outputDir, numericData, geneNames, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff)
end

%% ---------- Saving all opened figures as images

FolderName = strcat(outputDir,'\graphs');   % Destination folder for plots
mkdir(FolderName); 
% Getting all the opened figures, from last created to first.
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
FigList = flip(FigList); % Fliping the order of the figures.

if unique(eq(getGeneImages{1}, "none"))
    disp('No images saved');
else
    if unique(eq(getGeneImages{1}, "all"))
        disp("Saving all images");
    end
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName = get(FigHandle, 'Name');
        % prerequisite - no gene is called "all". 
        if contains(FigName, getGeneImages) || unique(eq(getGeneImages{1}, "all"))
            disp(FigName);
            save_figure(FigHandle, FigName, FolderName)     
        end
    end
end

disp('Done!');
end

% ---------- Helper functions

% Saving a figure by its handle, FigHandle,  and name, FigName. Image is
% saved in the folder FolderName.
function save_figure(FigHandle, FigName, FolderName)
    % FigHandle.Position = get(0, 'ScreenSize');
    FigHandle.PaperPositionMode = 'auto';
    str = sprintf('%s.png', FigName);
    file = fullfile(FolderName, str);
    saveas(FigHandle, file);
    if unique(class(FigHandle) == 'matlab.ui.Figure')
        str = sprintf('%s.fig', FigName);  
        file = fullfile(FolderName, str); 
        savefig(FigHandle, file, 'compact'); 
    end
end

% Checking whether PATIENT_ID's length in patientsNames and
% patientsNamesKM is equal.
% predicate - patientsNames look like 'TCGA-3B-A9HI-01', while
% patientsNamesKM looks like 'TCGA-3B-A9HI'  
function patientsNames = patientsNamesFromDiffDatasets(patientsNames, patientsNamesOtherDataset, use_cbioportal_not_lab_std)
if (use_cbioportal_not_lab_std)
    if (min(cellfun(@(x)length(x),patientsNames)) > min(cellfun(@(x)length(x),patientsNamesOtherDataset)))
      tmp = split(patientsNames, '-');
      if isrow(patientsNames)
          patientsNames = join(tmp(:,:,1:(end-1)), '-');
      else
          patientsNames = join(tmp(:,1:(end-1)), '-');
      end          
    end
end
end

% Finds all instances of non-existing info in the string vector 'str'.
% Yes, that's a black list that can only grow, it's wrong, but this is what
% works thus far. 
function bool = findEmptyString(str)
    bool = (str == "NA") | (str == "Nan") | (str == "NaN") | (str == "[Not Available]")...
        | (str == "") | (str == "Not Available");
end