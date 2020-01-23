function kaplan_meier_by_clustering(filename_km, outputDir, numericData, geneNames, patientsNames, ....
    patientsNamesKM, timeData, cens, timeCutOff) 



end

function calculate_clustergram(numericData, geneNames, patientsNames, title)

CGobj = clustergram(numericData, 'ColumnLabels', patientsNames, 'RowLabels', geneNames,... 
        'ImputeFun', @knnimpute, 'Standardize', 'row', 'ColumnLabelsRotate', 90);
addTitle(CGobj, title);
clusterGroup(cgobj, groupIndex, 'row'. tf_InfoOnly, 'true')
%divide_patients_by_clustering(numberOfGroups, 
end

