
function [W,H] = nmf(X,K,alg,maxiter,speak)
%
% NMF wrapper function
% function [W,H] = nmf(X,K,alg[,maxiter,speak])
%
% INPUT:
%           'X'     Inputmatrix
%           'K'     Number of components
%           'alg'   Algorithm to use: 
%                   'mm'     multiplicative updates using euclidean
%                            distance. Lee, D..D., and Seung, H.S., (2001)
%                   'cjlin'  alternative non-negative least squares using 
%                            projected gradients, author: Chih-Jen Lin, 
%                            National Taiwan University.
%                   'prob'   probabilistic NFM interpretating X as samples
%                            from a multinomial, author: Lars Kai Hansen,
%                            Technical University of Denmark
%                   'als'    Alternating Least Squares. Set negative
%                            elements to zero. 
%                   'alsobs' Alternating Least Squares. Set negative elements
%                            to zero and adjusts the other elements acording
%                            to Optimal Brain Surgeon. 
%           'maxiter'   Maximum number of iterations, default = 1000.
%           'speak'     Print information to screen unless speak = 0,
%                       default = 0
%
% OUTPUT:
% W       : N x K matrix
% H       : K x M matrix
%
% Kasper Winther Joergensen
% Informatics and Mathematical Modelling
% Technical University of Denmark
% kwj@imm.dtu.dk
% 2006/12/15

switch(nargin)
    case {0,1,2}
        error('Missing parameter. Type "help nmf" for usage.');
        
    case 3
        maxiter = 1000;
        speak = 0;
    case 4
        speak = 0;
    case 5
        % empty
    otherwise
        error('Too many parameters. Type "help nmf" for usage.');
        
end

% find dimensionallity of X
[D,N] = size(X);

% switch algorithm 
switch lower(alg)
    case 'mm'
        if speak, disp('Using mm algorithm'),end
        [W,H]=nmf_mm(X,K,maxiter,speak);
    case 'prob' 
        if speak, disp('Using prob algorithm'),end
        [W,H]=nmf_prob(X,K,maxiter,speak);
    case 'cjlin'
        if speak, disp('Using cjlin algorithm'),end
        [W,H]=nmf_cjlin(X,rand(D,K),rand(K,N),0.000001,10000,maxiter);
    case 'als'
        if speak, disp('Using als algorithm'),end
        [W,H]=nmf_als(X,K,maxiter,speak);
    case 'alsobs'
        if speak, disp('Using alsobs algorithm'),end
        [W,H]=nmf_alsobs(X,K,maxiter,speak);
    otherwise
        error('Unknown method. Type "help nmf" for usage.');
        
end

[W,H,nrgy] = order_comp(W,H);

function [W,H,nrgy]=order_comp(W,H)
%
% Order components according to "energy"
%
[D,K]=size(W);
[K,N]=size(H);
%
nrgy=zeros(K,1);
wsum=sum(W,1);
hsum=sum(H,2)';
nrgy=wsum.*hsum;
[nrgy,index]=sort(-nrgy);
nrgy=-nrgy;
W=W(:,index);
H=H(index,:);
function [W,H]=nmf_als(X,K,Nitsmax,speak)
%
% Truncated linear solution to LS -NMF
%
% INPUT:
% X (N,M) : N (dimensionallity) x M (samples) non negative input matrix
% K       : Number of components
% maxiter : Maximum number of iterations to run
% speak   : prints iteration count and changes in connectivity matrix 
%           elements unless speak is 0
%
% OUTPUT:
% W       : N x K matrix
% H       : K x M matrix
%
% Lars Kai Hansen, IMM-DTU (c) October 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print_iter = 50; % iterations between print on screen

[D,N]=size(X);
Xscale=sum(sum(X));
%INIT
W=rand(D,K);
H=rand(K,N);
Rscale=sum(sum(W*H));
sqrnorm=sqrt(Rscale/Xscale);
H=H/sqrnorm;
W=W/sqrnorm;

Xr_old = W*H;

%ITERATE
for n=1:Nitsmax,
    %    W=X*(pinv(H*H')*H)'; % old updates
    W = ((pinv(H*H')*H)*X')';
    W=(W>0).*W;
    W=W./(repmat(sum(W),D,1)+eps); % normalize columns to unit length

    H=(W*pinv(W'*W))'*X;
    H=H.*(H>0);

    % print to screen
    if (rem(n,print_iter)==0) & speak,
        Xr = W*H;
        diff = sum(sum(abs(Xr_old-Xr)));
        Xr_old = Xr;
        eucl_dist  = nmf_euclidean_dist(X,W*H);
        errorx=mean(mean(abs(X-W*H)))/mean(mean(X));
        disp(['Iter = ',int2str(n),...
            ', relative error = ',num2str(errorx),...
            ', diff = ', num2str(diff),...
            ', eucl dist ' num2str(eucl_dist)])
        if errorx < 10^(-5), break, end
    end
end
