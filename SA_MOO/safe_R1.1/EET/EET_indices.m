function  [ mi, sigma, EE, mi_sd, sigma_sd, mi_lb, sigma_lb, mi_ub, sigma_ub, mi_all, sigma_all ] = EET_indices(r,xmin,xmax,X,Y,design_type,varargin)
%
% Compute the sensitivity indices according to the Elementary Effects Test
% (Saltelli, 2008) or 'method of Morris' (Morris, 1991).
% These are: the mean (mi) of the EEs associated to input 'i',
% which measures the input influence; and the standard deviation (sigma)
% of the EEs, which measures its level of interactions with other inputs.
%
% Basic usage:
% [mi, sigma, EE] = EET_indices(r,xmin,xmax,X,Y,design_type)
%
% Input:
%           r = number of sampling point                           - scalar
%        xmin = lower bounds of input ranges                 - vector (1,M)
%        xmax = upper bounds of input ranges                 - vector (1,M)
%           X = matrix of sampling datapoints where EE must be computed
%                                                      - matrix (r*(M+1),M)
%           Y = associated output values               - vector (r*(M+1),1)
% design_type = design type (string)
%               [options: 'radial','trajectory']
% Output:
%          mi = mean of the elementary effects               - vector (1,M)
%       sigma = standard deviation of the elementary effects - vector (1,M)
%          EE = matrix of elementary effects                 - matrix (r,M)
%
%
% Advanced usage:
%
% [mi, sigma, EE] = EET_indices(r,xmin,xmax,X,Y,design_type,Nboot)
% [mi, sigma, EE] = EET_indices(r,xmin,xmax,X,Y,design_type,Nboot,alfa)
%
% Optional input:
%       Nboot = number of resamples used for boostrapping (default:0)
%        alfa = significance level for the confidence intervals estimated
%               by bootstrapping (default: 0.05)
% In this case, the output 'mi' and 'sigma' are the mean and standard
% deviation of the EEs averaged over Nboot resamples.
%
% Advanced usage/2:
%
% [mi, sigma, EE, mi_sd, sigma_sd, mi_lb, sigma_lb, mi_ub, sigma_ub] = ... 
%             EET_indices(r,xmin,xmax,X,Y,design_type,Nboot)
%
% Optional output:
%    mi_sd = standard deviation of 'mi' across Nboot resamples
% sigma_sd = standard deviation of 'sigma' across Nboot resamples
%    mi_lb = lower bound of 'mi' (at level alfa) across Nboot resamples
% sigma_lb = lower bound of 'sigma' across Nboot resamples
%    mi_ub = upper bound of 'mi' (at level alfa) across Nboot resamples
% sigma_ub = upper bound of 'sigma' across Nboot resamples
%                                                - all the above are 
%                                                 vector (1,M) if Nboot>1
%                                                 (empty vector otherwise)
% Or:
%
% [mi, sigma, EE, mi_sd, sigma_sd, mi_lb, sigma_lb, mi_ub, sigma_ub, ...
%      mi_all, sigma_all ] = EET_indices(r,xmin,xmax,X,Y,design_type,Nboot)
%
% Optional output:
%    mi_all = Nboot estimates of 'mi'                    - matrix (Nboot,M)
% sigma_all = Nboot estimates of 'sigma'                 - matrix (Nboot,M)
%
% REFERENCES:
% 
% Morris, M.D. (1991), Factorial sampling plans for preliminary          
% computational experiments, Technometrics, 33(2).
%
% Saltelli, A., et al. (2008) Global Sensitivity Analysis, The Primer,
% Wiley.

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isscalar(r); error('''r'' must be a scalar'); end
if r<=0; error('''r'' must be positive' ); end
if abs(r-round(r)); error('''r'' must be integer'); end
[N,M] = size(xmin)   ;
[n,m] = size(xmax) ;
if N~=1 ;error('''xmin'' must be a row vector'); end
if n~=1 ;error('''xmax'' must be a row vector'); end
if M~=m ;error('''xmin'' and ''xmax'' must be the same size'); end
Dr = xmax - xmin ;
if any(Dr<=0)
    error('all components of ''xmax'' must be higher than the corresponding ones in ''xmin''')
end
[n,m] = size(X) ;
if n~=r*(M+1) ;error('''X'' must have r*(M+1) rows'); end
if m~=M ;error('''X'' must have M columns'); end
[n,m] = size(Y) ;
if n~=r*(M+1) ;error('''Y'' must have r*(M+1) rows'); end
if m~=1 ;error('''Y'' must be a column vector'); end
if ~ischar(design_type); error('''design_type'' must be a string'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<7
    Nboot=1;
else
    Nboot=varargin{1};
    if ~isscalar(Nboot); error('''Nboot'' must be scalar'); end
    if Nboot<0; error('''Nboot'' must be nonnegative (if 0, bootstrapping is not used)' ); end
    if abs(Nboot-round(Nboot)); error('''Nboot'' must be an integer'); end
end
if nargin<8
    alfa=0.05;
else
    alfa=varargin{2};
    if ~isscalar(alfa); error('''alfa'' must be scalar'); end
    if any([alfa<0,alfa>1]); error('''alfa'' must be in [0,1]' ); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Elementary Effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EE = nan(r,M) ; % matrix of elementary effects
k  = 1 ;
ki = 1 ;
for i=1:r
    for j=1:M
        if strcmp(design_type,'radial') % radial design: EE is the difference 
            % between output at one point in the i-th block and output at
            % the 1st point in the block
            EE(i,j) = abs( Y(k+1) - Y(ki) ) / abs( X(k+1,j)-X(ki,j) )*Dr(j) ;
        elseif strcmp(design_type,'trajectory') % trajectory design: EE is the difference 
            % between output at one point in the i-th block and output at
            % the previous point in the block (the "block" is indeed a
            % trajectory in the input space composed of points that
            % differ in one component at the time)
            idx = find( abs( X(k+1,:)-X(k,:) ) > 0 );  % if using 'morris' 
            % sampling, the points in the block may not
            % be in the proper order, i.e. each point in the block differs
            % from the previous/next one by one component but we don't know
            % which one; this is here computed and saved in 'idx'
            if isempty(idx)  ; error('X(%d,:) and X(%d,:) are equal',[k,k+1]); end
            if length(idx)>1 ; error('X(%d,:) and X(%d,:) differ in more than one component',[k,k+1]); end
            EE(i,idx) = abs( Y(k+1) - Y(k) ) / abs( X(k+1,idx)-X(k,idx) )*Dr(idx) ;
        else
            error('''design_type'' must be one among {''radial'',''trajectory''}')
        end
        k=k+1 ;
    end
    k=k+1;
    ki=k ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Mean and Standard deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Nboot>1
    
    bootsize=r;
    B   = floor((rand(bootsize,Nboot)*r+1));
    mi_all=nan(Nboot,M) ;
    sigma_all=nan(Nboot,M);

    for n=1:Nboot
        mi_all(n,:)= mean(EE(B(:,n),:));
        sigma_all(n,:)=std(EE(B(:,n),:));
     end

    mi = mean(mi_all) ;
    mi_sd = std(mi_all)  ;
    mi_lb = sort(mi_all) ; mi_lb = mi_lb(max(1,round(Nboot*alfa/2)),:)  ;
    mi_ub = sort(mi_all) ; mi_ub  = mi_ub (round(Nboot*(1-alfa/2)),:)     ;

    sigma = mean(sigma_all) ;
    sigma_sd = std(sigma_all)  ;
    sigma_lb = sort(sigma_all) ; sigma_lb = sigma_lb(max(1,round(Nboot*alfa/2)),:)  ;
    sigma_ub = sort(sigma_all) ; sigma_ub  = sigma_ub (round(Nboot*(1-alfa/2)),:)     ;

    fprintf('\n\t mean(EE) std(EE)\n');
    fprintf('X%d:\t %2.3f\t %2.3f\n',[ 1:M; mi; sigma ]);
else
    mi       = mean(EE);
    sigma    = std(EE);
    mi_sd    = [] ;
    sigma_sd = [] ; 
    mi_lb    = [] ;
    sigma_lb = [] ;
    mi_ub    = [] ;
    sigma_ub = [] ;
    mi_all   = [] ;
    sigma_all= [] ;
    fprintf('\n\t mean(EE) std(EE)\n');
    fprintf('X%d:\t %2.3f\t %2.3f\n',[ 1:M; mi; sigma ]);
end