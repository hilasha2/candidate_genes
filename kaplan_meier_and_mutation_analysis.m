function kaplan_meier_and_mutation_analysis(filename_km, filename_mut, ...
    outputDir, gene_name, numericData, patientsNames, ...
    patientsNamesMutOtherDB, patientsNamesKM, geneNames, ...
    geneNamesMut, timeData, cens, timeCutOff, print_config)

% Spliting patients according to whether they have a mutation ..
% in driver gene, or not. 
[~, ~, patientsNamesMut, patientsNamesNotMut]...
    = split_gene_data_by_patients_w_mutation(numericData, patientsNames, ...
    patientsNamesMutOtherDB, gene_name, geneNamesMut);

% Kaplan Meier Analysis:
% Mutated 1st Gene (MG1), NotMutated 1st Gene (NMG1), 
% High 2nd Gene expression (HG2), Low 2nd Gene expression (LG2) 
% with the following combinations: 
% MG1 + HG2, MG1 + LG2, NMG1 + HG2, NMG1 + LG2.
fileOutput = sprintf('KM_CandidateGenesVsMutated_%s_', gene_name);
  KM_support_functions.kaplan_meier_gene1_vs_gene2_encapsulated(...
      filename_km, fileOutput, outputDir, ...
      gene_name, numericData, patientsNames, patientsNamesMut, ...
      patientsNamesNotMut, patientsNamesKM,...
      geneNames, timeData, cens, timeCutOff, 0.2, 0.8, 'mut_exp', ...
      {filename_mut}, print_config);

KM_support_functions.kaplan_meier_gene1_vs_gene2_encapsulated(...
    filename_km, fileOutput, outputDir, ...
	gene_name, numericData, patientsNames, patientsNamesMut, ...
    patientsNamesNotMut, patientsNamesKM,...
    geneNames, timeData, cens, timeCutOff, 0.8, 0.2, 'mut_exp', ...
    {filename_mut}, print_config);

KM_support_functions.kaplan_meier_gene1_vs_gene2_encapsulated(...
    filename_km, fileOutput, outputDir, ...
	gene_name, numericData, patientsNames, patientsNamesMut, ...
    patientsNamesNotMut, patientsNamesKM,...
    geneNames, timeData, cens, timeCutOff, 0.3, 0.7, 'mut_exp', ...
    {filename_mut}, print_config);

KM_support_functions.kaplan_meier_gene1_vs_gene2_encapsulated(...
    filename_km, fileOutput, outputDir, ...
	gene_name, numericData, patientsNames, patientsNamesMut, ...
    patientsNamesNotMut, patientsNamesKM,...
    geneNames, timeData, cens, timeCutOff, 0.7, 0.3, 'mut_exp', ...
    {filename_mut}, print_config);

end

