function colNum = ExcelStringColConverterToNum(colString)
if ~isstring(colString) && ~ischar(colString) 
    disp("The input enetered is not a string or a char");
    return
elseif isempty(colString) || colString == ""
    disp("Input string is empty");
    return
elseif ~all(isletter(char(colString)))
    disp("Input has non letters characters");
    return
end

colString = char(upper(colString));

colNum = 0;
for ch = colString
        colNum = colNum.*26;
        colNum = colNum + (ch - 'A' + 1);
end

end 
