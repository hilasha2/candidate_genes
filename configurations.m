%% Configurations
% breast cancer:
filename = '\\metlab24\d\InPut\Hila\breast cancer sources\Michal_gene_expression_summary_table.xlsx';
filename_km = '\\132.66.34.9\F\ilants2018\Grants\Current\1BCRF\Report0618Final\ZReport0618Data\Users\Michal\PatientsData\BC\BCP_data_clinical_patient_v2.xlsx';
sheet_name = 1;
do_genes_expression_calculations = 0;
do_clinical_calculations = 1;
do_subtype_histograms = 0;
use_days_not_months = 0;
colTimeData = 2;
colSurvivalStatus = 4;
colCancerSubtypes = 16;
colPatientsNamesKM = 1;
livingStatusStr = "Living";
deceasedButNotCancerStr = "Died of Other Causes"; 
timeCutOff = 120;
outputDir = 'E:\Hila\results michal\MET';
extraClinicalDataCellMatrix = [];
do_mutation_analysis = 0;
colMutationsTypes = 9;
colPatientsNamesMutation = 17;
colGeneNamesMut = 1;
filename_mut = '\\metlab24\d\InPut\Hila\ovary cancer sources\ov_tcga\data_mutations_extended.xlsx';
ppt_template = '\\metlab24\d\InPut\Hila\Candidate Genes\BRC_Candidate_Genes.pptx';
gene_nomenclature = "both";
gene_name = "MET";
use_cbioportal_not_lab_std = 1;
rowMutationsData = 3;
getGeneImages = {'none'};
do_mutation_and_clinical_analysis = 1;
%%
% breast cancer new:
filename = '\\132.66.95.239\D\InPut\Hila\breast cancer sources\bc_expression_summary.xlsx';
filename_km = '\\132.66.95.239\D\InPut\Hila\breast cancer sources\brca_metabric\data_clinical_patient.xlsx';
sheet_name = 'expression 24';
do_genes_expression_calculations = 0;
do_clinical_calculations = 1;
do_subtype_histograms = 1;
use_days_not_months = 0;
colTimeData = 13;
colSurvivalStatus = 14;
colCancerSubtypes = 15;
colPatientsNamesKM = 1;
livingStatusStr = "Living";
deceasedButNotCancerStr = "Died of Other Causes"; 
timeCutOff = 120;
outputDir = '\\132.66.95.239\d\InPut\Hila\Candidate Genes\BC_new\exp_24';
extraClinicalDataCellMatrix = [];
do_mutation_analysis = 0;
colMutationsTypes = 10;
colPatientsNamesMutation = 17;
colGeneNamesMut = 1;
filename_mut = '\\132.66.95.239\D\InPut\Hila\breast cancer sources\brca_metabric\data_mutations_extended.xlsx';
ppt_template = '132.66.95.239\D\InPut\Hila\Candidate Genes\BRC_Candidate_Genes.pptx';
gene_nomenclature = "both";
gene_name = "IGF1R";
use_cbioportal_not_lab_std = 1;
rowMutationsData = 3;
do_mutation_and_clinical_analysis = 0;
%%
sheet_name = 'expression 48';
outputDir = '\\132.66.95.239\d\InPut\Hila\Candidate Genes\BC_new\exp_48';

%%
sheet_name = 'phospho 24';
outputDir = '\\132.66.95.239\d\InPut\Hila\Candidate Genes\BC_new\phospho_24';

%%
sheet_name = 'phospho 48';
outputDir = '\\132.66.95.239\d\InPut\Hila\Candidate Genes\BC_new\phospho_48';

%% Breast Cancer Ori

% breast cancer new:
filename = '\\132.66.95.239\E\18\MetlabUsers2017\Ori\Candidate Genese Analysis\breast_cancer_gene_expression_summary.xlsx';
filename_km = '\\metlab26\D\Users\Hila Shacham\breast cancer sources\brca_metabric\data_clinical_patient.xlsx';
sheet_name = 'mimp';
do_genes_expression_calculations = 1;
do_clinical_calculations = 0;
do_subtype_histograms = 1;
use_days_not_months = 0;
colTimeData = 13;
colSurvivalStatus = 14;
colCancerSubtypes = 15;
colPatientsNamesKM = 1;
livingStatusStr = "Living";
deceasedButNotCancerStr = "Died of Other Causes"; 
timeCutOff = 120;
outputDir = '\\metlab26\D\Users\Hila Shacham\Candidate Genes\Ori\BC\29.10.2019\MTCH2';
extraClinicalDataCellMatrix = [];
do_mutation_analysis = 0;
colMutationsTypes = 9;
colPatientsNamesMutation = 17;
colGeneNamesMut = 1;
filename_mut = '\\metlab26\D\Users\Hila Shacham\breast cancer sources\brca_metabric\data_mutations_extended.xlsx';
ppt_template = '\\metlab26\D\Users\Hila Shacham\Candidate Genes\BRC_Candidate_Genes.pptx';
gene_nomenclature = "both";
gene_name = "MTCH2";
use_cbioportal_not_lab_std = 1;
rowMutationsData = 3;
getGeneImages = {'all'};
do_mutation_and_clinical_analysis = 1;

%% Breast Cancer Ophir


filename = '\\metlab25\E\Hila\Ophir_gene_expression_table_v4.xlsx';
filename_km = '\\132.66.95.239\D\InPut\Hila\breast cancer sources\brca_metabric\data_clinical_patient.xlsx';
sheet_name = 1;
do_genes_expression_calculations = 0;
do_clinical_calculations = 1;
do_subtype_histograms = 0;
use_days_not_months = 0;
colTimeData = 13;
colSurvivalStatus = 14;
colCancerSubtypes = 15;
colPatientsNamesKM = 1;
livingStatusStr = "Living";
deceasedButNotCancerStr = "Died of Other Causes"; 
timeCutOff = 120;
outputDir = '\\metlab24\d\InPut\Hila\Candidate Genes\Ophir output\results May 27 2019\STAT3';
extraClinicalDataCellMatrix = [];
do_mutation_analysis = 0;
colMutationsTypes = 10;
colPatientsNamesMutation = 17;
colGeneNamesMut = 1;
filename_mut = '\\132.66.95.239\D\InPut\Hila\breast cancer sources\brca_metabric\data_mutations_extended.xlsx';
ppt_template = '132.66.95.239\D\InPut\Hila\Candidate Genes\BRC_Candidate_Genes.pptx';
gene_nomenclature = "both";
gene_name = "STAT3";
use_cbioportal_not_lab_std = 1;
rowMutationsData = 3;
%getGeneImages = {'Kaplan' 'ADAMTS15' 'ADAMTS9' 'ADAMTSL2' 'ALDH1A2' 'ALPK3' 'AMIGO3' 'ANKRD39' 'ANKRD9' 'ARHGAP18' 'ARHGAP23' 'ARNTL' 'AVIL' 'BASP1' 'BCL2L1' 'BDKRB1' 'BLOC1S2' 'C1orf116' 'CACYBP' 'CCDC86' 'CD248' 'CD55'};
%getGeneImages =  {'CDH17' 'CDK18' 'CENPF' 'CHAF1A' 'CHST2' 'CLDN15' 'CLU' 'CNTFR' 'CORO2B' 'COTL1' 'CYP26B1' 'CYP2S1' 'DACT2' 'DHRS9' 'DIAPH3' 'ENTPD6' 'EXOG' 'FAT3' 'GGCT' 'GJA5' 'GMPPB' 'GNGT2' 'GPM6A' 'GRIN2A' 'HIST1H2BC' 'HOMER1'};
% getGeneImages =  {'HYOU1' 'ISLR2' 'ITGA11' 'ITPKC' 'KLHL13' 'LEO1' 'LIPG' 'LPIN1' 'LRRN2' 'LRRN4' 'LYPD1' 'MANF' 'MAST1' 'MCPH1' 'MSLN' 'MTSS1L' 'MUC16' 'NCOR2' 'NLE1' 'NOTCH3' 'NPAS2' 'NPR3' 'NR2F6' 'NT5E' 'NTRK1' 'PDIA5' 'PDZD2' 'PER2'};
% getGeneImages = {'PIK3C2B' 'PLEKHG2' 'PPA1' 'PTHLH' 'RAB20' 'RGS2' 'RHOC' 'RPS27L' 'RRP1B' 'RRP9' 'RYR3' 'SDAD1' 'SDF2L1' 'SEC61B' 'SFRP4' 'SLC1A2' 'SLPI' 'SPEN' 'STMN4' 'TAC1' 'TAF4B' 'TBRG4' 'TFPI2' 'THY1' 'TICAM1' 'TOP1MT' 'TRPC6' 'TSKU' 'TSR1' 'TTN' 'TUBB2B' 'TWIST1' 'UGP2' 'UNC13D' 'UNC5B' 'UPK1B' 'USP2' 'VSNL1' 'WDR77' 'WDR90' 'WNT10B' 'ZNF593'};
%getGeneImages = {'ODC1', 'BCAT1', 'COL29A1', 'LEPREL1', 'CXCR7'};
%getGeneImages = {'C14orf169' 'ADAM19' 'AHSA2' 'JUB' 'C14orf104' 'C21orf70' 'MLL2' 'MOBKL2A' 'FAM47E-STBD1' 'C4orf43' 'ADAMTS15' 'ADAMTS9' 'ADAMTSL2' 'ALDH1A2' 'ALPK3' 'AMIGO3' 'ANKRD39' 'ANKRD9' 'ARHGAP18' 'ARHGAP23' 'ARNTL' 'AVIL' 'BASP1' 'BCL2L1' 'BDKRB1' 'BLOC1S2' 'C1orf116' 'CACYBP' 'CCDC86' 'CD248' 'CD55'};
% getGeneImages = textscan(txt, '%s', 'TextType', 'string'); 
getGeneImages = {'Kaplan'};
do_mutation_and_clinical_analysis = 0;

%%
outputDir = '\\metlab25\e\Hila\results May 27 2019\MYC 2';
gene_name = "MYC";

%%
outputDir = '\\metlab25\e\Hila\results May 27 2019\NFKB1';
gene_name = "NFKB1";

%%
% lymphoma data set 3:
filename = '\\132.66.34.9\F\ilants2018\Grants\Current\3DotanlHematoOncology2016\Dotan2018Report\Data\OG_3HumanDatabase\RawSummary\Lek_Lym_DBPatients Summary.xlsx';
sheet_name = 2;
do_genes_expression_calculations = true;
do_clinical_calculations = true;
do_subtype_histograms = false;
use_days_not_months = false;
filename_km = '\\132.66.34.9\F\ilants2018\Grants\Current\3DotanlHematoOncology2016\Dotan2018Report\Data\OG_3HumanDatabase\RawDATALymphoma\Ly3\data_bcr_clinical_data_patient.xlsx';
colTimeData = 83; 
colSurvivalStatus = 28; 
colPatientsNamesKM = 2;
colCancerSubtypes = 0;
livingStatusStr = "Alive";
deceasedButNotCancerStr = ""; 
timeCutOff = 120;
extraClinicalDataCellMatrix = [];
outputDir = '\\132.66.95.239\D\InPut\Hila\Candidate Genes\ly3';
do_mutation_analysis = 0;
colMutationsTypes = 0;
colPatientsNamesMutation = 0;
colGeneNamesMut = 0;
filename_mut = '';
ppt_template = 'D:\InPut\Hila\Candidate Genes\BRC_Candidate_Genes.pptx';
gene_nomenclature = "both";
gene_name = "MET";
use_cbioportal_not_lab_std = 1;
rowMutationsData = 0;

Candidate_Genes(filename, filename_km, sheet_name,ppt_template, ...
    do_genes_expression_calculations, do_clinical_calculations, ... 
    do_subtype_histograms, use_cbioportal_not_lab_std,...
    gene_nomenclature, gene_name, colPatientsNamesKM,...
    colTimeData, colSurvivalStatus, colCancerSubtypes, livingStatusStr, ...
    deceasedButNotCancerStr, timeCutOff, use_days_not_months, outputDir,...
    extraClinicalDataCellMatrix, do_mutation_analysis, colMutationsTypes,...
    colPatientsNamesMutation, colGeneNamesMut, filename_mut, rowMutationsData)

% leukemia data set 4:
sheet_name = 3;
use_days_not_months = true;
filename_km = '\\132.66.34.9\F\ilants2018\Grants\Current\3DotanlHematoOncology2016\Dotan2018Report\Data\OG_3HumanDatabase\RawDATA_Leukemia\Le4\data_clinical_patient.xlsx';
colTimeData = 21; 
colSurvivalStatus = 20; 
colPatientsNamesKM = 1;
livingStatusStr = "LIVING";
%extraClinicalDataCellMatrix = {'sex', 2;'cns status', 8; 'testicular involvement', 9};
outputDir = '\\132.66.95.239\D\InPut\Hila\Candidate Genes\le4';
do_mutation_analysis = true;
colMutationsTypes = 10;
colPatientsNamesMutation = 17;
colGeneNamesMut = 1;
filename_mut = '\\132.66.34.9\f\ilants2018\Grants\Current\3DotanlHematoOncology2016\Dotan2018Report\Data\OG_3HumanDatabase\RawDATA_Leukemia\Le4\data_mutations_extended.xlsx';
rowMutationsData = 2;

Candidate_Genes(filename, filename_km, sheet_name,ppt_template, ...
    do_genes_expression_calculations, do_clinical_calculations, ... 
    do_subtype_histograms, use_cbioportal_not_lab_std,...
    gene_nomenclature, gene_name, colPatientsNamesKM,...
    colTimeData, colSurvivalStatus, colCancerSubtypes, livingStatusStr, ...
    deceasedButNotCancerStr, timeCutOff, use_days_not_months, outputDir,...
    extraClinicalDataCellMatrix, do_mutation_analysis, colMutationsTypes,...
    colPatientsNamesMutation, colGeneNamesMut, filename_mut, rowMutationsData) 
% leukemia data set 7:
 sheet_name = 5;
 use_days_not_months = false;
 filename_km = '\\132.66.34.9\f\ilants2018\Grants\Current\3DotanlHematoOncology2016\Dotan2018Report\Data\OG_3HumanDatabase\RawDATA_Leukemia\Le7\data_bcr_clinical_data_patient.xlsx';
 colTimeData = 38; 
 colSurvivalStatus = 37; 
 colPatientsNamesKM = 2;
 extraClinicalDataCellMatrix = [];
 outputDir = '\\132.66.95.239\D\InPut\Hila\Candidate Genes\le7';
do_mutation_analysis = 0;

Candidate_Genes(filename, filename_km, sheet_name,ppt_template, ...
    do_genes_expression_calculations, do_clinical_calculations, ... 
    do_subtype_histograms, use_cbioportal_not_lab_std,...
    gene_nomenclature, gene_name, colPatientsNamesKM,...
    colTimeData, colSurvivalStatus, colCancerSubtypes, livingStatusStr, ...
    deceasedButNotCancerStr, timeCutOff, use_days_not_months, outputDir,...
    extraClinicalDataCellMatrix, do_mutation_analysis, colMutationsTypes,...
    colPatientsNamesMutation, colGeneNamesMut, filename_mut, rowMutationsData)

% leukemia data set 8: 
sheet_name = 6;
filename_km = '\\132.66.34.9\F\ilants2018\Grants\Current\3DotanlHematoOncology2016\Dotan2018Report\Data\OG_3HumanDatabase\RawDATA_Leukemia\Le8\data_clinical_patient.xlsx';
colTimeData = 32; 
colSurvivalStatus = 30; 
colPatientsNamesKM = 1;
outputDir = '\\132.66.95.239\D\InPut\Hila\Candidate Genes\le8';

Candidate_Genes(filename, filename_km, sheet_name,ppt_template, ...
    do_genes_expression_calculations, do_clinical_calculations, ... 
    do_subtype_histograms, use_cbioportal_not_lab_std,...
    gene_nomenclature, gene_name, colPatientsNamesKM,...
    colTimeData, colSurvivalStatus, colCancerSubtypes, livingStatusStr, ...
    deceasedButNotCancerStr, timeCutOff, use_days_not_months, outputDir,...
    extraClinicalDataCellMatrix, do_mutation_analysis, colMutationsTypes,...
    colPatientsNamesMutation, colGeneNamesMut, filename_mut, rowMutationsData)

% lymphoma data set 2:
sheet_name = 1;
do_clinical_calculations = false;
outputDir = '\\132.66.95.239\D\InPut\Hila\Candidate Genes\ly2';

Candidate_Genes(filename, filename_km, sheet_name,ppt_template, ...
    do_genes_expression_calculations, do_clinical_calculations, ... 
    do_subtype_histograms, use_cbioportal_not_lab_std,...
    gene_nomenclature, gene_name, colPatientsNamesKM,...
    colTimeData, colSurvivalStatus, colCancerSubtypes, livingStatusStr, ...
    deceasedButNotCancerStr, timeCutOff, use_days_not_months, outputDir,...
    extraClinicalDataCellMatrix, do_mutation_analysis, colMutationsTypes,...
    colPatientsNamesMutation, colGeneNamesMut, filename_mut, rowMutationsData)

% leukemia data set 6:
sheet_name = 4;
outputDir = '\\132.66.95.239\D\InPut\Hila\Candidate Genes\le6';

Candidate_Genes(filename, filename_km, sheet_name,ppt_template, ...
    do_genes_expression_calculations, do_clinical_calculations, ... 
    do_subtype_histograms, use_cbioportal_not_lab_std,...
    gene_nomenclature, gene_name, colPatientsNamesKM,...
    colTimeData, colSurvivalStatus, colCancerSubtypes, livingStatusStr, ...
    deceasedButNotCancerStr, timeCutOff, use_days_not_months, outputDir,...
    extraClinicalDataCellMatrix, do_mutation_analysis, colMutationsTypes,...
    colPatientsNamesMutation, colGeneNamesMut, filename_mut, rowMutationsData)


%% Sarcoma

%filename = '\\metlab22\F\ilants2018\Grants\Current\2ISF2017\Sarcoma human database\sarcoma_rna_expression_summary.xlsx';
%filename = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V3 Updated gene list\Sarcoma\sarcoma_summary_table.xlsx';
filename = '\\metlab24\D\InPut\Hila\Candidate Genes\Nizan\July 2019\Sarcoma\sarcoma gene expression table tcga.xlsx';
sheet_name = 1;
filename_km = '\\metlab22\F\ilants2018\Grants\Current\2ISF2017\Sarcoma human database\sarc_tcga\data_bcr_clinical_data_patient.xlsx';
do_genes_expression_calculations = 1;
do_clinical_calculations = 1; 
do_subtype_histograms = 0;
use_days_not_months = 0;
colTimeData = 69;
colSurvivalStatus = 68;
colCancerSubtypes = 18;
colPatientsNamesKM = 2; 
livingStatusStr = "LIVING"; 
deceasedButNotCancerStr = "";
timeCutOff = 120; 
%outputDir = '\\metlab25\E\Hila\Nizan\sarc'; 
outputDir = '\\metlab24\D\InPut\Hila\Candidate Genes\Nizan\July 2019\Sarcoma';
extraClinicalDataCellMatrix = {}; 
do_mutation_analysis = 1;
colMutationsTypes = 10; 
colPatientsNamesMutation = 17;
colGeneNamesMut = 1;
filename_mut = "\\metlab22\F\ilants2018\Grants\Current\2ISF2017\Sarcoma human database\sarc_tcga\data_mutations_extended.xlsx";
gene_nomenclature = "both";
gene_name = "TP53";
use_cbioportal_not_lab_std = 1;
ppt_template = "D:\InPut\Hila\Candidate Genes\BRC_Candidate_Genes.pptx";
rowMutationsData = 2;
use_cbioportal_not_lab_std = 1;
getGeneImages = {'all'};
do_mutation_and_clinical_analysis = 1;
%% Sarcoma - Nizan's experiment tumors cc.line

filename = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V2 Updated gene list\sarcoma\expression_summary_sarc.xlsx';
sheet_name = 1;
filename_km = '\\metlab22\F\ilants2018\Grants\Current\2ISF2017\Sarcoma human database\sarc_tcga\data_bcr_clinical_data_patient.xlsx';
do_genes_expression_calculations = 1;
do_clinical_calculations = 1; 
do_subtype_histograms = 0;
use_days_not_months = 0;
colTimeData = 69;
colSurvivalStatus = 68;
colCancerSubtypes = 18;
colPatientsNamesKM = 2; 
livingStatusStr = "LIVING"; 
deceasedButNotCancerStr = "";
timeCutOff = 120; 
outputDir = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V2 Updated gene list\sarcoma\TP53'; 
extraClinicalDataCellMatrix = {}; 
do_mutation_analysis = 0;
colMutationsTypes = 10; 
colPatientsNamesMutation = 17;
colGeneNamesMut = 1;
filename_mut = "\\metlab22\F\ilants2018\Grants\Current\2ISF2017\Sarcoma human database\sarc_tcga\data_mutations_extended.xlsx";
gene_nomenclature = "both";
gene_name = "TP53";
use_cbioportal_not_lab_std = 1;
ppt_template = "\\metlab24\d\InPut\Hila\Candidate Genes\BRC_Candidate_Genes.pptx";
rowMutationsData = 2;

%% Lymphoma - Nizan's experiment tumors cc.line 

filename = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V2 Updated gene list\lymphoma\expression_summary_ly.xlsx';
sheet_name = 1;
filename_km = '\\metlab22\F\ilants2018\Grants\Current\3DotanlHematoOncology2016\Dotan2018Report\Data\Nizan\Lymphoma\Ly3\data_bcr_clinical_data_patient.xlsx';
do_genes_expression_calculations = 1;
do_clinical_calculations = 1; 
do_subtype_histograms = 0;
use_days_not_months = 0;
colTimeData = 83;
colSurvivalStatus = 82;
colCancerSubtypes = 0;
colPatientsNamesKM = 2; 
livingStatusStr = "LIVING"; 
deceasedButNotCancerStr = "";
timeCutOff = 120; 
outputDir = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V2 Updated gene list\lymphoma\TP53'; 
extraClinicalDataCellMatrix = {}; 
do_mutation_analysis = 0;
colMutationsTypes = 0; 
colPatientsNamesMutation = 0;
colGeneNamesMut = 0;
filename_mut = "";
gene_nomenclature = "both";
gene_name = "TP53"; % No TP53 in gene expressions
use_cbioportal_not_lab_std = 1;
ppt_template = "\\metlab24\d\InPut\Hila\Candidate Genes\BRC_Candidate_Genes.pptx";
rowMutationsData = 0;

%% Leukemia - Nizan's experiment tumors cc.lines
%% le7
filename = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V2 Updated gene list\leukemia\le7\expression_summary_le7.xlsx';
sheet_name = 1;
filename_km = '\\132.66.34.9\f\ilants2018\Grants\Current\3DotanlHematoOncology2016\Dotan2018Report\Data\OG_3HumanDatabase\RawDATA_Leukemia\Le7\data_bcr_clinical_data_patient.xlsx';
do_genes_expression_calculations = 1;
do_clinical_calculations = 1; 
do_subtype_histograms = 0;
use_days_not_months = 0;
colTimeData = 38;
colSurvivalStatus = 37;
colCancerSubtypes = 0;
colPatientsNamesKM = 2; 
livingStatusStr = "LIVING"; 
deceasedButNotCancerStr = "";
timeCutOff = 120; 
outputDir = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V2 Updated gene list\leukemia\le7\TP53'; 
extraClinicalDataCellMatrix = {}; 
do_mutation_analysis = 0;
colMutationsTypes = 0; 
colPatientsNamesMutation = 0;
colGeneNamesMut = 0;
filename_mut = "";
gene_nomenclature = "both";
gene_name = "TP53";
use_cbioportal_not_lab_std = 1;
ppt_template = "\\metlab24\d\InPut\Hila\Candidate Genes\BRC_Candidate_Genes.pptx";
rowMutationsData = 0;

%% le8 
filename = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V2 Updated gene list\leukemia\le8\expression_summary_le8.xlsx';
filename_km = '\\132.66.34.9\f\ilants2018\Grants\Current\3DotanlHematoOncology2016\Dotan2018Report\Data\OG_3HumanDatabase\RawDATA_Leukemia\Le8\data_clinical_patient.xlsx';
colTimeData = 32;
colSurvivalStatus = 30;
colPatientsNamesKM = 1; 
outputDir = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V2 Updated gene list\leukemia\le8\TP53';
gene_name = "TP53";

%% le4
filename = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V2 Updated gene list\leukemia\le4\expression_summary_le4.xlsx';
filename_km = '\\metlab22\F\ilants2018\Grants\Current\3DotanlHematoOncology2016\Dotan2018Report\Data\Nizan\Leukemia\Le4\data_clinical_patient.xlsx';
colTimeData = 22;
colSurvivalStatus = 20;
colPatientsNamesKM = 1; 
outputDir = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V2 Updated gene list\leukemia\le4\TP53';
gene_name = "TP53";

%% breast cancer - Nizan's experiment tumors cc.lines

filename = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V2 Updated gene list\breast cancer\expression_summary_bc.xlsx';
filename_km = '\\132.66.95.239\D\InPut\Hila\breast cancer sources\brca_metabric\data_clinical_patient.xlsx';
sheet_name = 1;
do_genes_expression_calculations = 1;
do_clinical_calculations = 1;
do_subtype_histograms = 1;
use_days_not_months = 0;
colTimeData = 13;
colSurvivalStatus = 14;
colCancerSubtypes = 15;
colPatientsNamesKM = 1;
livingStatusStr = "Living";
deceasedButNotCancerStr = "Died of Other Causes"; 
timeCutOff = 120;
outputDir = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\V2 Updated gene list\breast cancer\TP53';
extraClinicalDataCellMatrix = [];
do_mutation_analysis = 0;
colMutationsTypes = 10;
colPatientsNamesMutation = 17;
colGeneNamesMut = 1;
filename_mut = '\\132.66.95.239\D\InPut\Hila\breast cancer sources\brca_metabric\data_mutations_extended.xlsx';
ppt_template = '\\132.66.95.239\D\InPut\Hila\Candidate Genes\BRC_Candidate_Genes.pptx';
gene_nomenclature = "both";
gene_name = "TP53";
use_cbioportal_not_lab_std = 1;
rowMutationsData = 3;

%% Ovary cancer

% data set 1
filename = '\\metlab24\d\InPut\Hila\ovary cancer sources\ovary_cancer_expression_summary.xlsx';
sheet_name = 1;
filename_km = '\\metlab24\d\InPut\Hila\ovary cancer sources\ov_tcga_pub\data_clinical_patient.xlsx';
do_genes_expression_calculations = 1;
do_clinical_calculations = 1; 
do_subtype_histograms = 0;
use_days_not_months = 0;
colTimeData = 5;
colSurvivalStatus = 4;
colCancerSubtypes = 0;
colPatientsNamesKM = 1; 
livingStatusStr = "LIVING"; 
deceasedButNotCancerStr = "";
timeCutOff = 120; 
outputDir = '\\metlab24\d\InPut\Hila\Candidate Genes\output_ovary\ov_tcga_pub\HGF'; 
extraClinicalDataCellMatrix = {}; 
do_mutation_analysis = 1;
colMutationsTypes = 9; 
colPatientsNamesMutation = 16;
colGeneNamesMut = 1;
filename_mut = '\\metlab24\d\InPut\Hila\ovary cancer sources\ov_tcga_pub\data_mutations_extended.xlsx';
gene_nomenclature = "both";
gene_name = "HGF";
use_cbioportal_not_lab_std = 1;
ppt_template = '\\metlab24\d\InPut\Hila\Candidate Genes\BRC_Candidate_Genes.pptx';
rowMutationsData = 3;

Candidate_Genes(filename, filename_km, sheet_name,ppt_template, ...
    do_genes_expression_calculations, do_clinical_calculations, ... 
    do_subtype_histograms, use_cbioportal_not_lab_std,...
    gene_nomenclature, gene_name, colPatientsNamesKM,...
    colTimeData, colSurvivalStatus, colCancerSubtypes, livingStatusStr, ...
    deceasedButNotCancerStr, timeCutOff, use_days_not_months, outputDir,...
    extraClinicalDataCellMatrix, do_mutation_analysis, colMutationsTypes,...
    colPatientsNamesMutation, colGeneNamesMut, filename_mut, rowMutationsData)

%% data set 2

sheet_name = 2;
filename_km = '\\metlab24\d\InPut\Hila\ovary cancer sources\ov_tcga\data_bcr_clinical_data_patient.xlsx';
colTimeData = 57;
colSurvivalStatus = 56;
colPatientsNamesKM = 2; 
livingStatusStr = "LIVING"; 
outputDir = '\\metlab24\d\InPut\Hila\Candidate Genes\output_ovary\ov_tcga\HGF'; 
do_mutation_analysis = 1;
colMutationsTypes = 9; 
colPatientsNamesMutation = 16;
colGeneNamesMut = 1;
filename_mut = '\\metlab24\d\InPut\Hila\ovary cancer sources\ov_tcga\data_mutations_extended.xlsx';
rowMutationsData = 2;

Candidate_Genes(filename, filename_km, sheet_name,ppt_template, ...
    do_genes_expression_calculations, do_clinical_calculations, ... 
    do_subtype_histograms, use_cbioportal_not_lab_std,...
    gene_nomenclature, gene_name, colPatientsNamesKM,...
    colTimeData, colSurvivalStatus, colCancerSubtypes, livingStatusStr, ...
    deceasedButNotCancerStr, timeCutOff, use_days_not_months, outputDir,...
    extraClinicalDataCellMatrix, do_mutation_analysis, colMutationsTypes,...
    colPatientsNamesMutation, colGeneNamesMut, filename_mut, rowMutationsData)

%% Hala 
% BC 

filename = '\\metlab26\D\Users\Hila Shacham\breast cancer sources\Hala bc expression genes summary v3.xlsx';
filename_km = '\\metlab26\D\Users\Hila Shacham\breast cancer sources\brca_metabric\data_clinical_patient.xlsx';
sheet_name = 1;
do_genes_expression_calculations = 1;
do_clinical_calculations = 1;
do_subtype_histograms = 10;
use_days_not_months = 0;
colTimeData = 2;
colSurvivalStatus = 4;
colCancerSubtypes = 16;
colPatientsNamesKM = 1;
livingStatusStr = "Living";
deceasedButNotCancerStr = "Died of Other Causes"; 
timeCutOff = 120;
outputDir = '\\metlab26\D\Users\Hila Shacham\Candidate Genes\Hala\bc';
extraClinicalDataCellMatrix = [];
do_mutation_analysis = 1;
colMutationsTypes = 9;
colPatientsNamesMutation = 17;
colGeneNamesMut = 1;
filename_mut = '\\metlab26\D\Users\Hila Shacham\breast cancer sources\brca_metabric\data_mutations_extended.xlsx';
ppt_template = '\\metlab26\D\Users\Hila Shacham\Candidate Genes\BRC_Candidate_Genes.pptx';
gene_nomenclature = "both";
gene_name = "BRCA1";
use_cbioportal_not_lab_std = 1;
rowMutationsData = 3;
getGeneImages = {'none'};
do_mutation_and_clinical_analysis = 1;

%% Ovary 

filename = '\\metlab24\d\InPut\Hila\ovary cancer sources\Hala_ovary_cancer_expression_tcga.xlsx';
sheet_name = 1;
filename_km = '\\metlab24\d\InPut\Hila\ovary cancer sources\ov_tcga\data_bcr_clinical_data_patient.xlsx';
do_genes_expression_calculations = 1;
do_clinical_calculations = 1; 
do_subtype_histograms = 0;
use_days_not_months = 0;
colTimeData = 57;
colSurvivalStatus = 56;
colCancerSubtypes = 0;
colPatientsNamesKM = 2; 
livingStatusStr = "LIVING"; 
deceasedButNotCancerStr = "";
timeCutOff = 120; 
outputDir = '\\metlab25\e\Hila\bc\Hala\ovary'; 
extraClinicalDataCellMatrix = {}; 
do_mutation_analysis = 1;
colMutationsTypes = 11; 
colPatientsNamesMutation = 17;
colGeneNamesMut = 1;
filename_mut = '\\metlab24\d\InPut\Hila\ovary cancer sources\ov_tcga\data_mutations_extended.xlsx';
gene_nomenclature = "both";
gene_name = "BRCA1";
use_cbioportal_not_lab_std = 1;
ppt_template = '\\metlab24\d\InPut\Hila\Candidate Genes\BRC_Candidate_Genes.pptx';
rowMutationsData = 3;
getGeneImages = {'none'};
do_mutation_and_clinical_analysis = 1;

%% Lymphoma Nizan

filename = 'D:\Users\Hila Shacham\Candidate Genes\Nizan\July 2019\Lymphoid leukemia\all_phase2_target_2018_pub\lymphoma gene expression table all_phase2_target_2018_pub.xlsx';
sheet_name = 1;
filename_km = '\\metlab24\d\InPut\Hila\lymphoma cancer sources\all_phase2_target_2018_pub\data_clinical_patient.xlsx';
do_genes_expression_calculations = 1;
do_clinical_calculations = 0; 
do_subtype_histograms = 0;
use_days_not_months = 0;
colTimeData = 22;
colSurvivalStatus = 20;
colCancerSubtypes = 0;
colPatientsNamesKM = 1; 
livingStatusStr = "LIVING"; 
deceasedButNotCancerStr = "";
timeCutOff = 120; 
outputDir = '\\metlab25\e\Hila\Nizan\Lymphoma\allphase_no';
%outputDir = '\\metlab24\d\InPut\Hila\Candidate Genes\Nizan\July 2019\all_phase2_target_2018_pub'; 
extraClinicalDataCellMatrix = {}; 
do_mutation_analysis = 0;
colMutationsTypes = 10; 
colPatientsNamesMutation = 17;
colGeneNamesMut = 1;
filename_mut = '\\metlab24\d\InPut\Hila\lymphoma cancer sources\all_phase2_target_2018_pub\data_mutations_extended.xlsx';
gene_nomenclature = "both";
gene_name = "TP53";
use_cbioportal_not_lab_std = 1;
ppt_template = '\\metlab24\d\InPut\Hila\Candidate Genes\BRC_Candidate_Genes.pptx';
rowMutationsData = 2;
getGeneImages = {'all'};
do_mutation_and_clinical_analysis = 0;
%% bc - sheba
filename = 'D:\Users\Hila Shacham\Candidate Genes\Sheba\Output\list 4 short\list 4 short gene expression table.xlsx';
filename_km = '\\metlab26\D\Users\Hila Shacham\breast cancer sources\brca_metabric\data_clinical_patient.xlsx';
sheet_name = 1;
do_genes_expression_calculations = 1;
do_clinical_calculations = 1;
do_subtype_histograms = 10;
use_days_not_months = 0;
colTimeData = 2;
colSurvivalStatus = 4;
colCancerSubtypes = 16;
colPatientsNamesKM = 1;
livingStatusStr = "Living";
deceasedButNotCancerStr = "Died of Other Causes"; 
timeCutOff = 120;
outputDir = 'D:\Users\Hila Shacham\Candidate Genes\Sheba\Output\list 4 short';
extraClinicalDataCellMatrix = [];
do_mutation_analysis = 1;
colMutationsTypes = 9;
colPatientsNamesMutation = 17;
colGeneNamesMut = 1;
filename_mut = '\\metlab26\D\Users\Hila Shacham\breast cancer sources\brca_metabric\data_mutations_extended.xlsx';
ppt_template = '\\metlab26\D\Users\Hila Shacham\Candidate Genes\BRC_Candidate_Genes.pptx';
gene_nomenclature = "both";
gene_name = "STAT3";
use_cbioportal_not_lab_std = 1;
rowMutationsData = 3;
getGeneImages = {'none'};
do_mutation_and_clinical_analysis = 1;


%%
%txt = 'ABTB2 ADAMTS12 ADORA2A ALPL ANLN ANXA3 ANXA8 APLP1 ARHGDIB ASPM ATP6V0A4 BCAT1 BST1 BYSL CCL2 CCL7 CCRL2 CD14 CD44 CDC45 CDH12 CDH17 CDH2 CDK18 CDK5R1 CDKN1A CENPE CENPF CENPL CFI CHAF1A CHD7 CHORDC1 CHST2 CIAPIN1 CILP CLCF1 CLDN15 CLU CNTFR COL11A2 COL3A1 COL5A3 COL8A1 CORO2B COTL1 CREB3L1 CRELD2 CRLF1 CXCL1 CXCL10 CXCL2 CYP1B1 CYP26B1 CYP2S1 DACT2 DBN1 DGKH DHRS9 DIAPH3 DNAJA1 DNAJA4 DNAJB1 DNAJB2 DPF1 EDN3 EEF1E1 ELN ENTPD6 EXOG FAT3 FBN1 FERMT3 FGF18 FGFR4 FKBP11 FOSL1 FSCN1 FXYD5 GFRA1 GGCT GJA5 GMPPB GNGT2 GPM6A GRIN2A GRWD1 HAS2 HIF3A HIST1H2BC HOMER1 HSP90AA1 HSPA4L HSPA5 HSPA8 HSPB1 HSPB8 HSPD1 HSPH1 HYOU1 IGF1 IL6 IRF5 ISLR2 ITGA11 ITPKC KDM6B KLHL13 LAMC2 LBP LEO1 LIF LIPG LOXL2 LPIN1 LPL LRRN2 LRRN4 LYPD1 MANF MAPT MAST1 MCPH1 MFAP5 MLKL MMP14 MMP3 MSLN MTSS1L MUC16 NCOR2 NFKB1 NFKB2 NFKBIA NGF NGFR NIP7 NLE1 NOLC1 NOTCH3 NPAS2 NPR3 NR2F6 NT5E NTRK1 ODC1 P4HA3 PANX1 PCDH10 PDE8B PDIA5 PDIA6 PDLIM4 PDZD2 PER2 PIK3C2B PIM1 PLAT PLAU PLAUR PLEKHG2 PLK4 PPA1 PPRC1 PRC1 PRKG2 PTHLH PTK2B PVR RAB20 REL RELN RGS2 RHOC RIPK3 RND1 RPS19BP1 RPS27L RRP1B RRP9 RYR3 S100A9 SDAD1 SDF2L1 SEC61B SERPINE1 SFRP4 SLC15A3 SLC1A2 SLPI SOCS1 SORCS2 SPEN SPON2 STIP1 STMN4 SYK SYVN1 TAC1 TAF4B TBRG4 TFPI2 THBS1 THBS2 THY1 TICAM1 TIMP1 TLR2 TMSB10 TNFAIP3 TNFRSF11B TNFRSF12A TOP1MT TOP2A TPX2 TRIM69 TRPC6 TSKU TSR1 TTN TUBB2B TWIST1 TWIST2 UCHL3 UCK2 UGDH UGP2 UNC13D UNC5B UPK1B UPK3B USP2 VSNL1 WDR77 WDR90 WNT10B ZC3H12A ZNF593 ZSWIM6 ACKR3 COL6A5 P3H2 FGF16 MYC';
%gene_names = textscan(txt, '%s', 'TextType', 'string'); 
%for i = 1:length(gene_names{1})
%    gene_name = gene_names{1}{i};
Candidate_Genes(filename, filename_km, sheet_name,ppt_template, ...
    do_genes_expression_calculations, do_clinical_calculations, ... 
    do_subtype_histograms, use_cbioportal_not_lab_std,...
    gene_nomenclature, gene_name, colPatientsNamesKM,...
    colTimeData, colSurvivalStatus, colCancerSubtypes, livingStatusStr, ...
    deceasedButNotCancerStr, timeCutOff, use_days_not_months, outputDir,...
    extraClinicalDataCellMatrix, do_mutation_analysis, colMutationsTypes,...
    colPatientsNamesMutation, colGeneNamesMut, filename_mut, rowMutationsData, ...
    getGeneImages, do_mutation_and_clinical_analysis)
%end