function display_errors(err,grad_err,max_err)
% DISPLAY_ERRORS - Create error table given errors for continuous mesh
% levels
%
% Syntax:
%     display_errors(err,grad_err,max_err)
%
% Inputs:
%     err - array of L2 errors for mesh levels corresponding to indices
%     grad_err - array of L2 gradient errors for mesh levels corresponding 
%         to indices
%     max_err - array of max errors for mesh levels corresponding to
%         indicies
%
% Outputs:
%
% Author: Nicole Stock
% Date: Spring 2020

[~,n] = size(err);

fprintf(' Mesh\t| Error\t\t| Convergence\t|');
if exist('grad_err','var')
    fprintf(' Grad. Error\t| Convergence\t|');
end
if exist('max_err','var')
    fprintf(' Max Error\t| Convergence');
end
fprintf('\n------------------------------------------------------------');
fprintf('-------------------------------------------\n');
for i = 1:n
    fprintf(' %d\t| %2.2e\t|', i, err(i));
    if i >= 2
        conv = log2(err(i-1)/err(i));
        fprintf(' %2.2f\t\t|', conv);
    else
        fprintf('\t\t|')
    end
    if exist('grad_err','var')

        fprintf(' %2.2e\t|', grad_err(i));
        if i >= 2
            conv = log2(grad_err(i-1)/grad_err(i));
            fprintf(' %2.2f\t\t|', conv);
        else
            fprintf('\t\t|')
        end
    end
    if exist('max_err','var')
        fprintf(' %2.2e\t|', max_err(i));
        if i >= 2
            conv = log2(max_err(i-1)/max_err(i));
            fprintf(' %2.2f\t\t', conv);
        else
            fprintf('\t\t')
        end
    end
    fprintf('\n');
end