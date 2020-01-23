function barplot_pvalue_graph(titlePlot, upperCutoff, lowerCutoff, xtickAngle)
    close all;
    [file, path] = uigetfile('*.xlsx', 'Choose bar plot excel');
    file_path = fullfile(path, file);
    table = readtable(file_path, 'FileType', 'spreadsheet', 'TextType', 'string');
    genes = table2array(table(5:end, 1));
    ratio = str2double(table2array(table(5:end, 2)));
    pvalue = str2double(table2array(table(5:end, 3)));
    idx1 = ((abs(ratio) > upperCutoff) ) & (pvalue < 0.05);
    x1 = categorical(genes(idx1));
    y1 = abs(ratio(idx1));
    bulletsize1 = -log(pvalue(idx1))*12;
    idx2 = ((abs(ratio) < lowerCutoff) ) & (pvalue < 0.05);
    x2 = categorical(genes(idx2));
    y2 = abs(ratio(idx2));
    bulletsize2 = -log(pvalue(idx2))*12;
    plotGraph(x1, y1, bulletsize1, upperCutoff, 'upper', titlePlot, xtickAngle);
    plotGraph(x2, y2, bulletsize2, lowerCutoff, 'lower', titlePlot, xtickAngle);


    function plotGraph(x,y,bulletsize,cutOff, lowerOrUpper,titlePlot, xtickAngle)
        if ~isempty(y) 
            figure('visible', 'on')
            if strcmp(lowerOrUpper, 'lower') 
                scatter(x, y, bulletsize, 'g', 'o', 'filled', 'MarkerEdgeColor', 'k');
                miny = (min(y) - 1);
                ylim([miny 2]); 
            elseif strcmp(lowerOrUpper, 'upper')
                scatter(x, y, bulletsize, 'r', 'o', 'filled', 'MarkerEdgeColor', 'k');
                maxy = (max(y) + 1);
                ylim([0 maxy]);
            end    
            yline(1); 
            title(titlePlot);
            xlabel('Genes');
            ylabel('Expression Ratio');
            xtickangle(xtickAngle);
            str = sprintf('cutoff: %0.1f', cutOff);
            annotation('textbox', 'String', str, 'Position', [0.75 0.85 0.13 0.075]); 
        end
    end

end