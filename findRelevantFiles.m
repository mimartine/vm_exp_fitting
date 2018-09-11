function fileList = findRelevantFiles(inputList,targetStr)

listLength = size(inputList,1);
lv = nan(listLength,1);

for i = 1:listLength
    if isempty(strfind(inputList(i,:),targetStr))
        lv(i) = 1;
    end
    
end % for i ...

fileList = inputList(lv ~= 1,:);

end % of function...