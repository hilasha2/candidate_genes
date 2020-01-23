function [mutData, nonMutData, patientsNamesMut, patientsNamesNotMut]...
    = split_gene_data_by_patients_w_mutation(numericData, patientsNames, ...
    patientsNamesMutOtherDB, gene_name, geneNamesMut)

idx_gene_mut = cellfun(@(x) strcmp(gene_name, x), geneNamesMut);
patientsNamesMutMainGene = unique(patientsNamesMutOtherDB(idx_gene_mut));
idx_patients_mut = ismember(patientsNames, patientsNamesMutMainGene); 
patientsNamesMut = patientsNames(idx_patients_mut); 
patientsNamesNotMut = patientsNames(~idx_patients_mut);
mutData = numericData(:, idx_patients_mut);
nonMutData = numericData(:, ~idx_patients_mut);
end