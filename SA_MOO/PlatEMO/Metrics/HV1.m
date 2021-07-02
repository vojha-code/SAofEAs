function hv = HV1(Population)
% <max>
% Hypervolume
% Implementation: https://uk.mathworks.com/matlabcentral/fileexchange/30785-hypervolume-computation
    hv = 0;
    F = Population.best.objs;
    [M, l] = size(F);
	samples = 100000;
    
    ub = (max(F,[],2));
	lb = min(F')';
    
    F_samples = repmat(lb,1,samples) + rand(M,samples) .* repmat((ub - lb),1,samples);
	is_dominated_count = 0;
    for i = 1:samples
		for j = 1:l
			if (dominates(F(:,j), F_samples(:,i)))
				is_dominated_count = is_dominated_count + 1;
				break;
			end
		end
	end
	hv = prod(ub - lb) * (is_dominated_count / samples);
end

function d = dominates(fA, fB)
% [d] = dominates(fA, fB)
%
% Compares two solutions A and B given their objective function
% values fA and fB. Returns whether A dominates B.
%
% Input:
% - fA					- The objective function values of solution A
% - fB					- The objective function values of solution B
%
% Output:
% - d					- d is 1 if fA dominates fB, otherwise d is 0 
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011
	% Elegant, but not very efficient
	%d = (all(fA <= fB) && any(fA < fB));
	% Not so elegant, but more efficient
	d = false;
	for i = 1:length(fA)
		if (fA(i) > fB(i))
			d = false;
			return
		elseif (fA(i) < fB(i))
			d = true;
		end
	end
end