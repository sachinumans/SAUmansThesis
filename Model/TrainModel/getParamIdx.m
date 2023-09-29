function paramIdx = getParamIdx(paramNameStr,paramList)
% getParamIdx(paramNameStr,paramList) returns a column vector of indexes at
% which the elements of paramNameStr are found in paramList.

paramIdx = find(ismember(paramList,paramNameStr))';

if length(paramIdx) ~= length(paramNameStr)
    missingParams = [];
    for k=1:length(paramNameStr)
        if ~ismember(paramNameStr{k},paramList)
            missingParams = [missingParams paramNameStr{k} ', '];
            % fprintf('\nParameter ''%s'' is not in the parameter list.',paramNameStr{k});
        end
    end

    answer = questdlg([missingParams(1:end-2) ' not in the parameter list.'], ...
        'Missing parameters', ...
        'Ignore missing parameters', ...
        'Set breakpoint',...
        'Stop execution','Ignore missing parameters');

    switch answer
        case 'Ignore missing parameters'
            warning('Returned only the indexes of found parameters.')
        case 'Set breakpoint'
            dbstop in getParamIdx;
        case 'Stop execution'
            error(['Parameters ' missingParams(1:end-2) ' not found.'])
    end
end

end

