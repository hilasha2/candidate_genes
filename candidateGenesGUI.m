classdef candidateGenesGUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        FirstStepTab                   matlab.ui.container.Tab
        CANDIDATEGENESGUILabel         matlab.ui.control.Label
        ChooseDatasetDropDown          matlab.ui.control.DropDown
        ChooseadrivergeneDropDown      matlab.ui.control.DropDown
        SpecifyDriverGeneEditField     matlab.ui.control.EditField
        OutputImagesDropDown           matlab.ui.control.DropDown
        OutputDirectoryButton          matlab.ui.control.Button
        ChooseDatasetDropDownLabel     matlab.ui.control.Label
        ChooseDriverGeneDropDownLabel  matlab.ui.control.Label
        SpecifcyDriverGeneEditFieldLabel  matlab.ui.control.Label
        ChooseOutputDirectoryLabel     matlab.ui.control.Label
        OutputImagesDropDownLabel      matlab.ui.control.Label
        UploadFileButton               matlab.ui.control.Button
        UploadGeneExpressionFileLabel  matlab.ui.control.Label
        GeneExpressionFileTextArea     matlab.ui.control.TextArea
        OutputDirectoryTextArea        matlab.ui.control.TextArea
        TextAreaLabel                  matlab.ui.control.Label
        TextArea                       matlab.ui.control.TextArea
        SecondStepTab                  matlab.ui.container.Tab
        CANDIDATEGENESGUILabel_2       matlab.ui.control.Label
        GeneExpressionCheckBox         matlab.ui.control.CheckBox
        KMbyExpressionCheckBox         matlab.ui.control.CheckBox
        KMbyClusteringCheckBox         matlab.ui.control.CheckBox
        KMbyMutationCheckBox           matlab.ui.control.CheckBox
        CancerSubtypeCheckBox          matlab.ui.control.CheckBox
        MutationsTypeCheckBox          matlab.ui.control.CheckBox
        DoneButton                     matlab.ui.control.Button
        ChoosewhichanalysestodoLabel   matlab.ui.control.Label
    end

    
    properties (Access = private)
        funcArgs
    end
    
    methods (Access = private)
        
        function turnOffDriverGene(app)
            app.SpecifcyDriverGeneEditFieldLabel.Enable = 'off';
            app.SpecifyDriverGeneEditField.Enable = 'off';
            app.SpecifyDriverGeneEditField.Editable = 'off';
            app.SpecifyDriverGeneEditField.Value = '';
        end
        
        % checking whether the driver gene name is within Hugo Symbol's
        % guidelines:
        % Examples:
        % Correct:
        % JOI0-35345-3453
        % DJKFJGLK-3453-ET53
        % C3orf2B4BBB123C-C23orf123B34 % not correct irl but can't be helped.
        % A342-C123orf234-B4234
        % A342-C123orf21B34-B4234
        % 
        % Wrong:
        % S-ERWRW-3213
        % A-
        % QEQ-
        % 243DFG
        % AC23orf23
        % WEP304-2C23orf23
        
        % geneSymbol - a character vector if one gene inserted 
        % or a cell array with characters if multiple genes. 
        function  isValid = checkValidGeneSymbol(~, geneSymbol)
            
            expression = '^(([A-Z][A-Z0-9]+)|(C\d{1,4}orf\d{1,4}[A-Z0-9]{0,3}))(([-]?[A-Z0-9]+)|(-C\d{1,4}orf\d{1,4}[A-Z0-9]{0,3})){0,2}';
            
            [~, endIdx] = regexp(geneSymbol,expression);
            if iscellstr(geneSymbol)
                len_arr = cellfun(@length, geneSymbol);
                lenAsCell = num2cell(len_arr);
                isSameLength = cellfun(@isequal, lenAsCell, endIdx);
                isValid = isSameLength & ~isempty(geneSymbol);
            elseif ischar(geneSymbol)
                len = length(geneSymbol);
                isValid = isequal(len, endIdx) & ~isempty(geneSymbol);
            end
        end
        
        % Pop an error message if genes are not according to hugo symbol
        % nomenclature.
        function faultyGenesMsg(~, isValid, genes)
            if iscellstr(genes)
                if length(isValid) > 1
                    genes = genes(~isValid);
                end
                invalid_genes = strjoin(genes,' ');
            else
                invalid_genes = genes;
            end
            msg = sprintf('Invalid gene name(s): %s', invalid_genes);
            errordlg(msg, 'cannot start calculation',"replace");
        end
        
        function turnOffImagesText(app)
            app.TextAreaLabel.Enable = 'off';
            app.TextArea.Editable = 'off';
            app.TextArea.Enable = 'off';
        end
        
        % Since Matlab has some problems saving files within different 
        % computers unless the given path is in its IP form - I change some
        % computer names into their IP address. 
        function newPath = changePathToIP(~, path)
            if contains(path, '\\metlab24')
                newPath = strrep(path, '\\metlab24', '\\132.66.95.239');
            elseif contains(path, '\\metlab25')
                newPath = strrep(path, '\\metlab25', '\\132.66.95.241');
            elseif contains(path, '\\metlab26')
                newPath = strrep(path, '\\metlab26', '\\132.66.207.18');
            else
                newPath = path;
            end
            
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)

            % Initializing accodring to germ_cell_1 dataset.
            % clinical data file:
            app.funcArgs.filename_km = '\\metlab25\G\TsarfatyLabUsers2020\Nizan\Hila\DB for human survival analysis\germ cell\Testicular Germ Cell Tumors (TCGA, PanCancer Atlas)\data_clinical_patient.xlsx';
            % mutations file:
            app.funcArgs.filename_mut = '\\metlab25\G\TsarfatyLabUsers2020\Nizan\Hila\DB for human survival analysis\germ cell\Testicular Germ Cell Tumors (TCGA, PanCancer Atlas)\data_mutations_extended.xlsx';
            app.funcArgs.sheet_name = 1; % sheet index in gene expression file.
            app.funcArgs.use_days_not_months = 0; % whether to use days vs. months as the time unit for survival analysis.
            app.funcArgs.colTimeData = 32; % overall survival in months column in clinical file - OS_MONTHS.
            app.funcArgs.colSurvivalStatus = 31; % OS_STATUS (LIVING/DECEASED) column in clinical file.
            app.funcArgs.colCancerSubtypes = 0; % usually not needed if do_subtype_histograms == 0.
            app.funcArgs.colPatientsNamesKM = 1; % PATIENT_ID column in clinical data.
            app.funcArgs.livingStatusStr = "LIVING";
            app.funcArgs.deceasedButNotCancerStr = ""; % usually not needed. Only in breast cancer. 
            app.funcArgs.timeCutOff = 120; % Stop calculating the survival analysis after 120 months. 
            app.funcArgs.extraClinicalDataCellMatrix = []; % redundant currently.
            app.funcArgs.colMutationsTypes = 10; % Where the types of mutations are located in mutation file.
            app.funcArgs.colPatientsNamesMutation = 17; % patients id column in mutations table (tumor_sample_barcode). 
            app.funcArgs.colGeneNamesMut = 1; % the column with gene names in mutations file. 
            app.funcArgs.rowMutationsData = 2; % Where the data begins in mutations table.
            app.funcArgs.use_cbioportal_not_lab_std = 1; % by which standard the dataset is presented, can be the lab's standard 0, or cbioportal 1.
            app.funcArgs.gene_nomenclature = "both"; % nomenclatures of genes can be 'hugo', 'entrez', or 'both'.
            
            % Chosen values by the user:
            app.funcArgs.filename = 0; % gene expression excel file.
            app.funcArgs.gene_name = 'MET'; % driver gene name.
            app.funcArgs.getGeneImages = {'none'}; % which images to get. can be either 'all', 'none' or specific by gene names.
            app.funcArgs.outputDir = 0; % output directory.
            % Analyses
            app.funcArgs.do_genes_expression_calculations = 0; 
            app.funcArgs.do_km_analysis = 0; 
            app.funcArgs.do_km_by_mutation_and_expression = 0; 
            app.funcArgs.do_km_by_clustering = 0; 
            app.funcArgs.do_mutation_analysis = 0; 
            app.funcArgs.do_subtype_histograms = 0;
        end

        % Value changed function: ChooseDatasetDropDown
        function chooseDataset(app, event)
            value = app.ChooseDatasetDropDown.Value;
            switch value
                case 'germ_cell_1'
                    app.funcArgs.filename_km = '\\metlab26\d\Users\Hila Shacham\Candidate Genes\DB for human survival analysis\germ cell\Testicular Germ Cell Tumors (TCGA, PanCancer Atlas)\data_clinical_patient.xlsx';
                    app.funcArgs.filename_mut = '\\metlab26\d\Users\Hila Shacham\Candidate Genes\DB for human survival analysis\germ cell\Testicular Germ Cell Tumors (TCGA, PanCancer Atlas)\data_mutations_extended.xlsx';
                    app.funcArgs.colTimeData = 32; 
                    app.funcArgs.colSurvivalStatus = 31; 
                    app.funcArgs.colCancerSubtypes = 0;
                    app.funcArgs.colPatientsNamesKM = 1; 
                    app.funcArgs.livingStatusStr = "LIVING";
                    app.funcArgs.deceasedButNotCancerStr = ""; 
                    app.funcArgs.colMutationsTypes = 10; 
                    app.funcArgs.colPatientsNamesMutation = 17; 
                    app.funcArgs.colGeneNamesMut = 1;
                    app.funcArgs.rowMutationsData = 2; 
                case 'lymphoma_1'
                    app.funcArgs.filename_km = '\\metlab26\d\Users\Hila Shacham\Candidate Genes\DB for human survival analysis\lymphoma\Pediatric Acute Lymphoid Leukemia - Phase II (TARGET, 2018)\data_clinical_patient.xlsx';
                    app.funcArgs.filename_mut = '\\metlab26\d\Users\Hila Shacham\Candidate Genes\DB for human survival analysis\lymphoma\Pediatric Acute Lymphoid Leukemia - Phase II (TARGET, 2018)\data_mutations_extended.xlsx';
                    app.funcArgs.colTimeData = 22; 
                    app.funcArgs.colSurvivalStatus = 20; 
                    app.funcArgs.colCancerSubtypes = 0; 
                    app.funcArgs.colPatientsNamesKM = 1; 
                    app.funcArgs.livingStatusStr = "LIVING";
                    app.funcArgs.deceasedButNotCancerStr = "";
                    app.funcArgs.colMutationsTypes = 10; 
                    app.funcArgs.colPatientsNamesMutation = 17; 
                    app.funcArgs.colGeneNamesMut = 1;
                    app.funcArgs.rowMutationsData = 2; 
                case 'sarcoma_1'
                    app.funcArgs.filename_km = '\\metlab26\d\Users\Hila Shacham\Candidate Genes\DB for human survival analysis\sarcoma\Adult Soft Tissue Sarcomas TCGA, Cell 2017\data_clinical_patient.xlsx';
                    app.funcArgs.filename_mut = '\\metlab26\d\Users\Hila Shacham\Candidate Genes\DB for human survival analysis\sarcoma\Adult Soft Tissue Sarcomas TCGA, Cell 2017\data_mutations_extended.xlsx';
                    app.funcArgs.colTimeData = 13; 
                    app.funcArgs.colSurvivalStatus = 11;
                    app.funcArgs.colCancerSubtypes = 0; 
                    app.funcArgs.colPatientsNamesKM = 1; 
                    app.funcArgs.livingStatusStr = "LIVING";
                    app.funcArgs.deceasedButNotCancerStr = "";  
                    app.funcArgs.colMutationsTypes = 10; 
                    app.funcArgs.colPatientsNamesMutation = 17;  
                    app.funcArgs.colGeneNamesMut = 1;
                    app.funcArgs.rowMutationsData = 2;
                case 'breast_cancer_1'
                    app.funcArgs.filename_km = '\\metlab26\d\Users\Hila Shacham\breast cancer sources\brca_metabric\data_clinical_patient.xlsx';
                    app.funcArgs.filename_mut = '\\metlab26\d\Users\Hila Shacham\breast cancer sources\brca_metabric\data_mutations_extended.xlsx';
                    app.funcArgs.colTimeData = 13; 
                    app.funcArgs.colSurvivalStatus = 17; 
                    app.funcArgs.colCancerSubtypes = 15; 
                    app.funcArgs.colPatientsNamesKM = 1; 
                    app.funcArgs.livingStatusStr = "Living";
                    app.funcArgs.deceasedButNotCancerStr = "Died of Other Causes"; 
                    app.funcArgs.colMutationsTypes = 9; 
                    app.funcArgs.colPatientsNamesMutation = 17; 
                    app.funcArgs.colGeneNamesMut = 1;
                    app.funcArgs.rowMutationsData = 3;  
                case 'ovary_cancer_1'
                    app.funcArgs.filename_km = '\\metlab26\d\Users\Hila Shacham\Candidate Genes\DB for human survival analysis\ovary cancer sources\Ovarian Serous Cystadenocarcinoma (TCGA, Provisional)\data_bcr_clinical_data_patient.xlsx';
                    app.funcArgs.filename_mut = '\\metlab26\d\Users\Hila Shacham\Candidate Genes\DB for human survival analysis\ovary cancer sources\Ovarian Serous Cystadenocarcinoma (TCGA, Provisional)\data_mutations_extended.xlsx';
                    app.funcArgs.colTimeData = 57; 
                    app.funcArgs.colSurvivalStatus = 56; 
                    app.funcArgs.colCancerSubtypes = 0; % no subtypes 
                    app.funcArgs.colPatientsNamesKM = 1; 
                    app.funcArgs.livingStatusStr = "LIVING";
                    app.funcArgs.deceasedButNotCancerStr = ""; 
                    app.funcArgs.colMutationsTypes = 9; 
                    app.funcArgs.colPatientsNamesMutation = 16; 
                    app.funcArgs.colGeneNamesMut = 1;
                    app.funcArgs.rowMutationsData = 2; 
            end
            
        end

        % Button pushed function: UploadFileButton
        function UploadFileButtonPushed(app, event)
            [file, path] = uigetfile('*.xlsx', 'Select gene expression file');
            if isequal(file, 0) && isequal(app.funcArgs.filename, 0)
                app.funcArgs.filename = 0;
            elseif ~isequal(file, 0)
                newPath = changePathToIP(app, path);
                app.funcArgs.filename = fullfile(newPath, file);
                app.GeneExpressionFileTextArea.Value = app.funcArgs.filename;
            end
        end

        % Value changed function: ChooseadrivergeneDropDown
        function ChooseadrivergeneDropDownValueChanged(app, event)
            value = app.ChooseadrivergeneDropDown.Value;
            switch value
                case "MET"
                    app.funcArgs.gene_name = 'MET';
                    app.turnOffDriverGene();
                case "TP53"
                    app.funcArgs.gene_name = 'TP53';
                    app.turnOffDriverGene();
                case "BRCA1"
                    app.funcArgs.gene_name = 'BRCA1';
                    app.turnOffDriverGene();
                case "Other"
                    app.funcArgs.gene_name = '';
                    app.SpecifcyDriverGeneEditFieldLabel.Enable = 'on';
                    app.SpecifyDriverGeneEditField.Enable = 'on';
                    app.SpecifyDriverGeneEditField.Editable = 'on';
            end
            
        end

        % Value changed function: SpecifyDriverGeneEditField
        function SpecifyDriverGeneEditFieldValueChanged(app, event)
            value = app.SpecifyDriverGeneEditField.Value;
            app.funcArgs.gene_name = value;
        end

        % Value changed function: OutputImagesDropDown
        function OutputImagesDropDownValueChanged(app, event)
            value = app.OutputImagesDropDown.Value;
            switch value
                case "none"
                    app.funcArgs.getGeneImages = {'none'};
                    app.turnOffImagesText();
                case "all"
                    app.funcArgs.getGeneImages = {'all'};
                    app.turnOffImagesText();
                case "other"
                    app.TextAreaLabel.Enable = 'on';
                    app.TextArea.Editable = 'on';
                    app.TextArea.Enable = 'on';
            end
        end

        % Value changed function: TextArea
        function TextAreaValueChanged(app, event)
            value = app.TextArea.Value;
            if iscellstr(value)
                value = strjoin(value,' ');
            end
            genes = split(value, {',', ', ', ' ,', ' '});
            genes(cellfun(@isempty, genes)) = [];
            app.funcArgs.getGeneImages = genes;
        end

        % Button pushed function: OutputDirectoryButton
        function OutputDirectoryButtonPushed(app, event)
            dirPath= uigetdir;
            if ~isequal(dirPath, 0)
                newDirPath = changePathToIP(app, dirPath);
                app.funcArgs.outputDir = newDirPath;
                app.OutputDirectoryTextArea.Value = app.funcArgs.outputDir;
            end
        end

        % Value changed function: GeneExpressionCheckBox
        function GeneExpressionCheckBoxValueChanged(app, event)
            value = app.GeneExpressionCheckBox.Value;
            app.funcArgs.do_genes_expression_calculations = value;
        end

        % Value changed function: KMbyExpressionCheckBox
        function KMbyExpressionCheckBoxValueChanged(app, event)
            value = app.KMbyExpressionCheckBox.Value;
            app.funcArgs.do_km_analysis = value;
        end

        % Value changed function: KMbyMutationCheckBox
        function KMbyMutationCheckBoxValueChanged(app, event)
            value = app.KMbyMutationCheckBox.Value;
            app.funcArgs.do_km_by_mutation_and_expression = value;
        end

        % Value changed function: KMbyClusteringCheckBox
        function KMbyClusteringCheckBoxValueChanged(app, event)
            value = app.KMbyClusteringCheckBox.Value;
            app.funcArgs.do_km_by_clustering = value; 
        end

        % Button pushed function: DoneButton
        function DoneButtonPushed(app, event)
            % Checking whether gene expression table is OK:
            if isequal(app.funcArgs.filename, 0)
                errordlg('No gene expression file chosen', 'cannot start calculation',"replace");
                return
            end
            [~,~,ext] = fileparts(app.funcArgs.filename);
            if ~isequal(ext,'.xlsx')
                 errordlg('Wrong file extension for gene expression file', 'cannot start calculation',"replace");
                 return
            end
            
            % Check whether the driver gene was correctly inserted:
            if isempty(app.funcArgs.gene_name)
                errordlg('No driver gene chosen', 'cannot start calculation',"replace");
                return
            end
            validDriverGene = checkValidGeneSymbol(app, app.funcArgs.gene_name);
            if ~validDriverGene
                faultyGenesMsg(app, validDriverGene, app.funcArgs.gene_name);
                return
            end
            
            % Check whether the gene names for which to output images were
            % inserted correctly:
            if isequal(app.OutputImagesDropDown.Value, 'other')
                if isempty(app.funcArgs.getGeneImages)
                   errordlg('No gene names were inserted for output images', 'cannot start calculation',"replace");
                   return
                end
                validGeneNames = checkValidGeneSymbol(app, app.funcArgs.getGeneImages);
                if sum(~validGeneNames) > 0
                    faultyGenesMsg(app, validGeneNames, app.funcArgs.getGeneImages);
                    return
                end
            end
            
            % Check whether an output directory was given:
            if isequal(app.funcArgs.outputDir, 0)
                errordlg('No output directory chosen', 'cannot start calculation',"replace");
                return
            end
            
            % Check that at least one analysis was chosen:
            anyAnalysis = app.funcArgs.do_genes_expression_calculations |...
            app.funcArgs.do_km_analysis | ...
            app.funcArgs.do_km_by_mutation_and_expression | ...
            app.funcArgs.do_km_by_clustering | ...
            app.funcArgs.do_mutation_analysis | ... 
            app.funcArgs.do_subtype_histograms;
            
            if isequal(anyAnalysis, 0)
                errordlg('Choose at least one analysis', 'cannot start calculation',"replace");
                return
            end
            
            progress = uiprogressdlg(app.UIFigure, ...
                'Title', 'Code Running', ...
                'Message', 'Your code is running, will close when done', ...
                'Indeterminate', 'on');
            
            Candidate_Genes(app.funcArgs.filename, ...
                app.funcArgs.filename_km, ...
                app.funcArgs.sheet_name,...
                app.funcArgs.do_genes_expression_calculations, ...
                app.funcArgs.do_km_analysis, ... 
                app.funcArgs.do_subtype_histograms, ...
                app.funcArgs.use_cbioportal_not_lab_std,...
                app.funcArgs.gene_nomenclature, ...
                app.funcArgs.gene_name, ...
                app.funcArgs.colPatientsNamesKM,...
                app.funcArgs.colTimeData, ...
                app.funcArgs.colSurvivalStatus, ...
                app.funcArgs.colCancerSubtypes, ...
                app.funcArgs.livingStatusStr, ...
                app.funcArgs.deceasedButNotCancerStr, ...
                app.funcArgs.timeCutOff, ...
                app.funcArgs.use_days_not_months, ...
                app.funcArgs.outputDir,...
                app.funcArgs.extraClinicalDataCellMatrix, ...
                app.funcArgs.do_mutation_analysis, ...
                app.funcArgs.colMutationsTypes,...
                app.funcArgs.colPatientsNamesMutation, ...
                app.funcArgs.colGeneNamesMut, ...
                app.funcArgs.filename_mut, ...
                app.funcArgs.rowMutationsData, ...
                app.funcArgs.getGeneImages, ...
                app.funcArgs.do_km_by_mutation_and_expression, ...
                app.funcArgs.do_km_by_clustering);
            close(progress);

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 599 581];
            app.UIFigure.Name = 'UI Figure';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [2 1 599 581];

            % Create FirstStepTab
            app.FirstStepTab = uitab(app.TabGroup);
            app.FirstStepTab.Title = 'First Step';

            % Create CANDIDATEGENESGUILabel
            app.CANDIDATEGENESGUILabel = uilabel(app.FirstStepTab);
            app.CANDIDATEGENESGUILabel.HorizontalAlignment = 'center';
            app.CANDIDATEGENESGUILabel.FontName = 'Yu Gothic UI Semibold';
            app.CANDIDATEGENESGUILabel.FontSize = 18;
            app.CANDIDATEGENESGUILabel.FontWeight = 'bold';
            app.CANDIDATEGENESGUILabel.Position = [177 508 200 27];
            app.CANDIDATEGENESGUILabel.Text = 'CANDIDATE GENES GUI';

            % Create ChooseDatasetDropDown
            app.ChooseDatasetDropDown = uidropdown(app.FirstStepTab);
            app.ChooseDatasetDropDown.Items = {'Testicular germ cell tumor - TCGA PanCancer Atlas', 'Limphoid Leukemia - TARGET 2018', 'Sarcoma soft tissue - TCGA 2017', 'Breast cancer - METABRIC Nature 2016', 'Ovarian Serous Cystadenocarcinoma (TCGA, Provisional)'};
            app.ChooseDatasetDropDown.ItemsData = {'germ_cell_1', 'lymphoma_1', 'sarcoma_1', 'breast_cancer_1', 'ovary_cancer_1'};
            app.ChooseDatasetDropDown.ValueChangedFcn = createCallbackFcn(app, @chooseDataset, true);
            app.ChooseDatasetDropDown.FontName = 'Yu Gothic UI';
            app.ChooseDatasetDropDown.FontSize = 14;
            app.ChooseDatasetDropDown.Position = [131 460 418 22];
            app.ChooseDatasetDropDown.Value = 'germ_cell_1';

            % Create ChooseadrivergeneDropDown
            app.ChooseadrivergeneDropDown = uidropdown(app.FirstStepTab);
            app.ChooseadrivergeneDropDown.Items = {'MET', 'TP53', 'BRCA1', 'Other'};
            app.ChooseadrivergeneDropDown.ValueChangedFcn = createCallbackFcn(app, @ChooseadrivergeneDropDownValueChanged, true);
            app.ChooseadrivergeneDropDown.FontName = 'Yu Gothic UI';
            app.ChooseadrivergeneDropDown.FontSize = 14;
            app.ChooseadrivergeneDropDown.Position = [171 325 100 22];
            app.ChooseadrivergeneDropDown.Value = 'MET';

            % Create SpecifyDriverGeneEditField
            app.SpecifyDriverGeneEditField = uieditfield(app.FirstStepTab, 'text');
            app.SpecifyDriverGeneEditField.ValueChangedFcn = createCallbackFcn(app, @SpecifyDriverGeneEditFieldValueChanged, true);
            app.SpecifyDriverGeneEditField.Editable = 'off';
            app.SpecifyDriverGeneEditField.FontName = 'Yu Gothic UI';
            app.SpecifyDriverGeneEditField.Enable = 'off';
            app.SpecifyDriverGeneEditField.Position = [365 289 100 22];

            % Create OutputImagesDropDown
            app.OutputImagesDropDown = uidropdown(app.FirstStepTab);
            app.OutputImagesDropDown.Items = {'None', 'All', 'Specific IDMGs (excluding driver gene)'};
            app.OutputImagesDropDown.ItemsData = {'none', 'all', 'other'};
            app.OutputImagesDropDown.ValueChangedFcn = createCallbackFcn(app, @OutputImagesDropDownValueChanged, true);
            app.OutputImagesDropDown.FontName = 'Yu Gothic UI';
            app.OutputImagesDropDown.FontSize = 14;
            app.OutputImagesDropDown.Position = [441 247 100 22];
            app.OutputImagesDropDown.Value = 'none';

            % Create OutputDirectoryButton
            app.OutputDirectoryButton = uibutton(app.FirstStepTab, 'push');
            app.OutputDirectoryButton.ButtonPushedFcn = createCallbackFcn(app, @OutputDirectoryButtonPushed, true);
            app.OutputDirectoryButton.FontName = 'Yu Gothic UI';
            app.OutputDirectoryButton.FontSize = 14;
            app.OutputDirectoryButton.Position = [177 66 120 28];
            app.OutputDirectoryButton.Text = 'Output Directory';

            % Create ChooseDatasetDropDownLabel
            app.ChooseDatasetDropDownLabel = uilabel(app.FirstStepTab);
            app.ChooseDatasetDropDownLabel.FontName = 'Yu Gothic UI';
            app.ChooseDatasetDropDownLabel.FontSize = 14;
            app.ChooseDatasetDropDownLabel.Position = [20 460 102 22];
            app.ChooseDatasetDropDownLabel.Text = 'Choose Dataset';

            % Create ChooseDriverGeneDropDownLabel
            app.ChooseDriverGeneDropDownLabel = uilabel(app.FirstStepTab);
            app.ChooseDriverGeneDropDownLabel.FontName = 'Yu Gothic UI';
            app.ChooseDriverGeneDropDownLabel.FontSize = 14;
            app.ChooseDriverGeneDropDownLabel.Position = [20 325 136 22];
            app.ChooseDriverGeneDropDownLabel.Text = 'Choose a driver gene';

            % Create SpecifcyDriverGeneEditFieldLabel
            app.SpecifcyDriverGeneEditFieldLabel = uilabel(app.FirstStepTab);
            app.SpecifcyDriverGeneEditFieldLabel.FontName = 'Yu Gothic UI';
            app.SpecifcyDriverGeneEditFieldLabel.Enable = 'off';
            app.SpecifcyDriverGeneEditFieldLabel.Position = [20 289 337 22];
            app.SpecifcyDriverGeneEditFieldLabel.Text = 'Specify other driver gene name (Hugo symbol, case sensitive):';

            % Create ChooseOutputDirectoryLabel
            app.ChooseOutputDirectoryLabel = uilabel(app.FirstStepTab);
            app.ChooseOutputDirectoryLabel.FontName = 'Yu Gothic UI';
            app.ChooseOutputDirectoryLabel.FontSize = 14;
            app.ChooseOutputDirectoryLabel.Position = [20 69 156 22];
            app.ChooseOutputDirectoryLabel.Text = 'Choose output directory';

            % Create OutputImagesDropDownLabel
            app.OutputImagesDropDownLabel = uilabel(app.FirstStepTab);
            app.OutputImagesDropDownLabel.FontName = 'Yu Gothic UI';
            app.OutputImagesDropDownLabel.FontSize = 14;
            app.OutputImagesDropDownLabel.Position = [21 247 410 22];
            app.OutputImagesDropDownLabel.Text = 'Output images? (If you have too many genes, don''t choose ''All''):';

            % Create UploadFileButton
            app.UploadFileButton = uibutton(app.FirstStepTab, 'push');
            app.UploadFileButton.ButtonPushedFcn = createCallbackFcn(app, @UploadFileButtonPushed, true);
            app.UploadFileButton.FontName = 'Yu Gothic UI';
            app.UploadFileButton.FontSize = 14;
            app.UploadFileButton.Position = [239 400 100 28];
            app.UploadFileButton.Text = 'Upload File';

            % Create UploadGeneExpressionFileLabel
            app.UploadGeneExpressionFileLabel = uilabel(app.FirstStepTab);
            app.UploadGeneExpressionFileLabel.FontName = 'Yu Gothic UI';
            app.UploadGeneExpressionFileLabel.FontSize = 14;
            app.UploadGeneExpressionFileLabel.Position = [20 403 213 22];
            app.UploadGeneExpressionFileLabel.Text = 'Upload gene expression excel file';

            % Create GeneExpressionFileTextArea
            app.GeneExpressionFileTextArea = uitextarea(app.FirstStepTab);
            app.GeneExpressionFileTextArea.Editable = 'off';
            app.GeneExpressionFileTextArea.FontName = 'Arial';
            app.GeneExpressionFileTextArea.FontColor = [0.502 0.502 0.502];
            app.GeneExpressionFileTextArea.Position = [20 363 521 28];

            % Create OutputDirectoryTextArea
            app.OutputDirectoryTextArea = uitextarea(app.FirstStepTab);
            app.OutputDirectoryTextArea.Editable = 'off';
            app.OutputDirectoryTextArea.FontName = 'Arial';
            app.OutputDirectoryTextArea.FontColor = [0.502 0.502 0.502];
            app.OutputDirectoryTextArea.Position = [20 24 521 28];

            % Create TextAreaLabel
            app.TextAreaLabel = uilabel(app.FirstStepTab);
            app.TextAreaLabel.FontName = 'Yu Gothic UI';
            app.TextAreaLabel.Enable = 'off';
            app.TextAreaLabel.Position = [20 205 510 32];
            app.TextAreaLabel.Text = {'Images for the following IDMGs (Hugo symbol, case sensitive, separated by space or comma):'; 'Example 1: BAX, ATM, PDZD2; or Example 2: BAX ATM PDZD2'};

            % Create TextArea
            app.TextArea = uitextarea(app.FirstStepTab);
            app.TextArea.ValueChangedFcn = createCallbackFcn(app, @TextAreaValueChanged, true);
            app.TextArea.Editable = 'off';
            app.TextArea.Enable = 'off';
            app.TextArea.Position = [21 112 509 84];

            % Create SecondStepTab
            app.SecondStepTab = uitab(app.TabGroup);
            app.SecondStepTab.Title = 'Second Step';

            % Create CANDIDATEGENESGUILabel_2
            app.CANDIDATEGENESGUILabel_2 = uilabel(app.SecondStepTab);
            app.CANDIDATEGENESGUILabel_2.HorizontalAlignment = 'center';
            app.CANDIDATEGENESGUILabel_2.FontName = 'Yu Gothic UI Semibold';
            app.CANDIDATEGENESGUILabel_2.FontSize = 18;
            app.CANDIDATEGENESGUILabel_2.FontWeight = 'bold';
            app.CANDIDATEGENESGUILabel_2.Position = [177 508 200 27];
            app.CANDIDATEGENESGUILabel_2.Text = 'CANDIDATE GENES GUI';

            % Create GeneExpressionCheckBox
            app.GeneExpressionCheckBox = uicheckbox(app.SecondStepTab);
            app.GeneExpressionCheckBox.ValueChangedFcn = createCallbackFcn(app, @GeneExpressionCheckBoxValueChanged, true);
            app.GeneExpressionCheckBox.Text = 'Gene Expression Calculations';
            app.GeneExpressionCheckBox.FontName = 'Yu Gothic UI';
            app.GeneExpressionCheckBox.FontSize = 14;
            app.GeneExpressionCheckBox.Position = [35 403 203 22];

            % Create KMbyExpressionCheckBox
            app.KMbyExpressionCheckBox = uicheckbox(app.SecondStepTab);
            app.KMbyExpressionCheckBox.ValueChangedFcn = createCallbackFcn(app, @KMbyExpressionCheckBoxValueChanged, true);
            app.KMbyExpressionCheckBox.Text = 'Kaplan Meier by expression and co-expression';
            app.KMbyExpressionCheckBox.FontName = 'Yu Gothic UI';
            app.KMbyExpressionCheckBox.FontSize = 14;
            app.KMbyExpressionCheckBox.Position = [34 368 312 22];

            % Create KMbyClusteringCheckBox
            app.KMbyClusteringCheckBox = uicheckbox(app.SecondStepTab);
            app.KMbyClusteringCheckBox.ValueChangedFcn = createCallbackFcn(app, @KMbyClusteringCheckBoxValueChanged, true);
            app.KMbyClusteringCheckBox.Text = 'Kaplan Meier after clustering methods on expression';
            app.KMbyClusteringCheckBox.FontName = 'Yu Gothic UI';
            app.KMbyClusteringCheckBox.FontSize = 14;
            app.KMbyClusteringCheckBox.Position = [34 294 352 22];

            % Create KMbyMutationCheckBox
            app.KMbyMutationCheckBox = uicheckbox(app.SecondStepTab);
            app.KMbyMutationCheckBox.ValueChangedFcn = createCallbackFcn(app, @KMbyMutationCheckBoxValueChanged, true);
            app.KMbyMutationCheckBox.Text = 'Kaplan Meier by mutation and expression';
            app.KMbyMutationCheckBox.FontName = 'Yu Gothic UI';
            app.KMbyMutationCheckBox.FontSize = 14;
            app.KMbyMutationCheckBox.Position = [34 331 281 22];

            % Create CancerSubtypeCheckBox
            app.CancerSubtypeCheckBox = uicheckbox(app.SecondStepTab);
            app.CancerSubtypeCheckBox.Enable = 'off';
            app.CancerSubtypeCheckBox.Text = 'Cancer subtype analysis';
            app.CancerSubtypeCheckBox.FontName = 'Yu Gothic UI';
            app.CancerSubtypeCheckBox.FontSize = 14;
            app.CancerSubtypeCheckBox.Position = [34 259 169 22];

            % Create MutationsTypeCheckBox
            app.MutationsTypeCheckBox = uicheckbox(app.SecondStepTab);
            app.MutationsTypeCheckBox.Enable = 'off';
            app.MutationsTypeCheckBox.Text = 'Mutations type analysis';
            app.MutationsTypeCheckBox.FontName = 'Yu Gothic UI';
            app.MutationsTypeCheckBox.FontSize = 14;
            app.MutationsTypeCheckBox.Position = [34 224 167 22];

            % Create DoneButton
            app.DoneButton = uibutton(app.SecondStepTab, 'push');
            app.DoneButton.ButtonPushedFcn = createCallbackFcn(app, @DoneButtonPushed, true);
            app.DoneButton.FontName = 'Yu Gothic UI';
            app.DoneButton.FontSize = 14;
            app.DoneButton.Position = [250 132 100 28];
            app.DoneButton.Text = 'Done';

            % Create ChoosewhichanalysestodoLabel
            app.ChoosewhichanalysestodoLabel = uilabel(app.SecondStepTab);
            app.ChoosewhichanalysestodoLabel.FontName = 'Yu Gothic UI';
            app.ChoosewhichanalysestodoLabel.FontSize = 14;
            app.ChoosewhichanalysestodoLabel.FontWeight = 'bold';
            app.ChoosewhichanalysestodoLabel.Position = [35 448 193 22];
            app.ChoosewhichanalysestodoLabel.Text = 'Choose which analyses to do:';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = candidateGenesGUI

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end