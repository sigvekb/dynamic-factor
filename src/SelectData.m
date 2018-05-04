function [selectedData, newStruct, varsFound] = SelectData(data, structure)

selectedData = data;
selectedData(selectedData == 0) = NaN;

% Define block structure in terms of new ordering
[d,b]=size(structure);
newStruct = zeros(d,b);
varsFound = [];
for i=1:b
    block = structure(~isnan(structure(:,i)),i);
    len = length(block);
    for j=1:len
        var = structure(j,i);
        found = find(varsFound==var);
        if isempty(found)
            newStruct(j,i) = length(varsFound)+1;
            varsFound = [varsFound var];
        else
            newStruct(j,i) = found;
        end
    end
end

selectedData = selectedData(:, varsFound);