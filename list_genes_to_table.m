[file_exp, path_exp] = uigetfile('*.xlsx', 'Expression table');
if file_exp == 0
    disp("no expression file chosen");
    return;
end
[file_gene_lst, path_gene_lst] = uigetfile('*.xlsx', 'Gene list table');
if file_gene_lst == 0
    disp("no gene list file chosen");
    return;
end
dir = uigetdir();
if dir == 0
    disp("no output directory chosen");
    return;
end
xls_name = fullfile(dir, 'result.xlsx');
% v1 - [~,~,expression_table] = xlsread(fullfile(path_exp, file_exp));
expression_table = readtable(fullfile(path_exp, file_exp), 'FileType', 'spreadsheet'); 
fullfile_gene_lst = fullfile(path_gene_lst, file_gene_lst);
[~,sheets] = xlsfinfo(fullfile_gene_lst);
% v1 - gene_names = expression_table(:, 1); 
gene_names = expression_table.Hugo_Symbol; 

% v1 - table1 = expression_table(1,:);

%%
for  sheet = sheets
   [~,gene_lst_txt] = xlsread(fullfile_gene_lst, sheet{1});
   % v1 - gene_lst_txt = unique(gene_lst_txt);
   [~,idx] = intersect(gene_names, gene_lst_txt);
   % table2 = expression_table(idx,:);
   table = expression_table(idx, :);
   table.Properties.VariableNames = expression_table.Properties.VariableNames;
   % v1 - table = [table1; table2];
   % v1 - writetable(cell2table(table), xls_name, 'Sheet', sheet{1}, 'FileType', 'spreadsheet', 'WriteVariableNames', false);
   writetable(table, xls_name, 'Sheet', sheet{1}, 'FileType', 'spreadsheet', 'WriteVariableNames', true);
end
