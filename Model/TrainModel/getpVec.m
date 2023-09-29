function pVec = getpVec(p_nonOpt, pOpt_list, paramList)
% getpVec(p_nonOpt, pOpt_list, paramList) constructs an anonymous function
% that creates a parameter vector for optimization. The function preserves 
% the value of non-optimized parameters and replaces the elements
% corresponding to to-be-optimized parameters with 'wildcards'.
%
% Inputs:
% p_nonOpt  : nParams-by-1 vector containing the values of the non-optimized 
%             parameters, where nParams is the total number of parameters. 
%             The order must correspond to the one of the paramList.
% pOpt_list : list of parameters to be optimized.
% paramList : ordered list of all parameter names.
%
% Outputs:
% pVec : function handle, expects input in column vector form (nOpt-by-1).

nPar = length(paramList); % total number of parameters
nOpt = length(pOpt_list); % number of optimized parameters
pOpt_idx = getParamIdx(pOpt_list,paramList); % get indexes of optimized params

% Get indexes of non-optimized parameters in the parameter list
p_nonOpt_select = double(~ismember(paramList,pOpt_list))'; 

% Vector whose entries are 0 if the corresponding parameter is optimized, or
% the non-optimized value otherwise
p_nonOpt_vec = (p_nonOpt_select'*diag(p_nonOpt))';

% Column k corresponds to the k-th optimized parameter in pOpt_list and
% consists of the vector indicating the logical position of the optimized
% parameter in paramList.
p_Opt_select_mat = zeros(nPar,nOpt); % nParams-by-nOpt
for k=1:nOpt
    p_Opt_select_mat(pOpt_idx(k),k) = 1;
end

% Second term: ector whose entries are 0 if the corresponding parameter is 
% non-optimized, or to-be-optimized parameter otherwise
pVec = @(p) (p_nonOpt_vec +  (p_Opt_select_mat*p));

end