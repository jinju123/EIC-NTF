function [t]=mvntftv(X,R,L,A,B,C,options,Atrue)
% x: the HSI in XJK where the I J K denote the width height and band number
% R: the number of endmembers
% L: Rank of each abudance map
% A: initialization of A
% B: initial value of B
% C: initial value of C
% options:
%    derta: the sum to one penalty
%        
% maxiters=options.maxiters;%number of maximum iteration
convergeNum=options.convergeNum;%converge number
derta=options.derta;%sum to one penalty
mu=options.mu;%regular parameter
lambda=options.lambda;
[I,J,K]=size(X);     % Size of the problem
% Matrix Unfoldings of X
X1=tens2mat(X,1);  % X1 is IxJK 
X2=tens2mat(X,2);  % X2 is JxIK (or smaller dimensions if compression was done)
X3=tens2mat(X,3);  % X3 is KxJI (idem)
Ps=kron(eye(R),ones(L,1)); % Pattern matrix
% LOOP for alternating updates
objNew=norm(X1-A*myKr(C,B,ones(1,R),ones(1,R)*L)','fro');%min_{A}||X_{JK¡ÁI}?(B@C) ¡¤ AT||2 F
btdRes{1}=A;
btdRes{2}=B;
btdRes{3}=C;
abundanceRes=normAbundance(btdRes,R,L);%abundances
U=A;
V=B;
for i=1:R
    W=btdRes{1}(:,(L*(i-1))+1:(L*i));
    H=btdRes{2}(:,(L*(i-1))+1:(L*i));
    E{i}=W*H';%A_r*B_r
end
times=0;
objhistory=[];
objhistory=[objhistory objNew];
sad=0;
sadhis=[];
rmse=0;
[sparsity]=calSparsity(abundanceRes');
sumToOne=ones(I,J);
iter=1;
while 1
    oldEndmember=C;
    oldAbundance=abundanceRes;
    objOld=objNew;
    oldSparsity=sparsity;
    oldSad=sad;
    oldRmse=rmse;  
    %Aold=A;
    M1=myKr(C,B,ones(1,R),ones(1,R)*L);%B@C
    A=A.*((X1*M1+derta*sumToOne*B+mu*U)./(A*M1'*M1+1e-4+derta*A*B'*B+mu*A));
%   A=(X1*M1+derta*sumToOne*B+mu*U)/(M1'*M1+1e-4+derta*B'*B+mu);
    M2=myKr(C,A,ones(1,R),ones(1,R)*L);
     %B=(X2*M2+derta*sumToOne'*Aold+mu*V)/(M2'*M2+1e-4+derta*Aold'*Aold+mu);
    B=B.*((X2*M2+derta*sumToOne'*A+mu*V)./(B*M2'*M2+1e-4+derta*B*A'*A+mu*B));
    
    M=kr(B,A)*Ps;
     %C=(X3*M)/(M'*M+eps);
    C=C.*((X3*M)./(C*M'*M));
   
    for i=1:R
    W=A(:,(L*(i-1))+1:(L*i));
    uu=U(:,(L*(i-1))+1:(L*i));
    vv=V(:,(L*(i-1))+1:(L*i));
    ee=E{i};
     %U(:,(L*(i-1))+1:(L*i))=(ee*vv+W)/(1+vv'*vv);
    U(:,(L*(i-1))+1:(L*i))=uu.*(ee*vv+W)./(uu+uu*vv'*vv);
    end
    %updateU
    for i=1:R  
    H=B(:,(L*(i-1))+1:(L*i));
    uu=U(:,(L*(i-1))+1:(L*i));
    vv=V(:,(L*(i-1))+1:(L*i));
    ee=E{i};
    % V(:,(L*(i-1))+1:(L*i))=(ee'*uu+H)/(1+uu'*uu);
    V(:,(L*(i-1))+1:(L*i))=vv.*(ee'*uu+H)./(vv+vv*uu'*uu);
    end
%     %updateV
    for i=1:R
    uu=U(:,(L*(i-1))+1:(L*i));
    vv=V(:,(L*(i-1))+1:(L*i));
    Ehat{i}=uu*vv';                  
    %E{i}=fgp_denoise(Ehat{i},17,0.1,lambda/mu,200);
    E{i}=denoise_bound(Ehat{i},lambda/mu,min(min(Ehat{i})),max(max(Ehat{i})),200,0.001,'iso');
    %E{i}=(mu*uu*vv'+lambda*eehat)/(lambda+mu);
    end
    
%     for i=1:R
%     ee=E{i};
%     ds=3;
%     Ds=7;
%     h=10;
%     Ehat{i}=fastNLmeans1(ee,ds,Ds,h);  
%     end
%     %compute Ehat
%     for i=1:R
%     rankVector=ones(1,R)*L;
%     uu=U(:,sumRank(rankVector,i-1)+1:sumRank(rankVector,i));
%     vv=V(:,sumRank(rankVector,i-1)+1:sumRank(rankVector,i));
%     eehat=Ehat{i};
%     E{i}=(mu*uu*vv'+lambda*eehat)/(lambda+mu);
%     end
%     ssum=0;
%     for i=1:R
%       ssum=ssum+norm(E{i}-Ehat{i},'fro');  
%     end
    objNew=0.5*norm(X1-A*myKr(C,B,ones(1,R),ones(1,R)*L)','fro')+0.5*derta*norm(sumToOne-A*B','fro');
    objhistory=[objhistory objNew];
    btdRes{1}=A;
    btdRes{2}=B;
    btdRes{3}=C;
    %abundanceRes=normAbundance(btdRes,R,L);
      for i=1:R
     e=reshape(E{i},I*J,1);
     abundanceRes(:,i)=e;
      end
    toltemp = abs(objOld - objNew)/objOld;
    fprintf('iter [%d]: obj [%d],sad [%d],rmse [%d],sparsity [%d], relative error [%d], C[%d]\n ', iter,objNew,sad,rmse,sparsity,toltemp,norm(C(:),2));
    oldsad=sad;
    [sad,allSadDistance,sor]=cosDistance(C,Atrue);
     sadhis=[sadhis oldsad];
    [rmse]=HyperRmse(oldAbundance,abundanceRes,sor);
    [sparsity]=calSparsity(abundanceRes);
    %     ||(abs(oldSparsity-sparsity)/oldSparsity<2*1e-4&&abs(oldSad-sad)/oldSad<2*(1e-4)&&abs(oldRmse-rmse)/oldRmse<2*(1e-4))
    %
%     if (abs(oldSparsity-sparsity)/oldSparsity<2e-4&&abs(oldSad-sad)/oldSad<(2e-4)&&abs(oldRmse-rmse)/oldRmse<(2e-4))
%         times = times + 1;
%     else
%         times=0;
%     end
%     if toltemp<1e-4
%         times = times + 1;
%     else
%         times=0;
%     end
%     
%      times=0;
%      if iter>300&&oldsad<sad
%          times=times+1;
     %end
if toltemp<1e-4||iter==convergeNum
    
        t{1}=C;
        t{2}=abundanceRes;
        t{3}=objhistory;
        t{4}=btdRes;
        t{5}=sadhis;
        break;
end
    iter = iter+1;
%     if iter==convergeNum
%         t{1}=C;
%         t{2}=abundanceRes;
%         t{3}=objhistory;
%         t{4}=btdRes;
%         break;
%     end
    % Define a stop criterion
    
end

%    STEP3: COME BACK to ORIGINAL SPACE and PERFORM A FEW MORE ITERATIONS IN ORIGINAL SPACE

end

%*******************************************************************************
function [U1,U2,U3,S,S1,S2,S3] = mlsvd3(X,size_core)
%MLSVD3 Multilinear singular value decomposition of a third-order tensor.
[I1,I2,I3]=size(X);
[U1,S1,temp]=svd(reshape(X,I1,I3*I2),'econ'); S1=diag(S1);
[U2,S2,temp]=svd(reshape(permute(X,[2 3 1]),I2,I1*I3),'econ'); S2=diag(S2);
[U3,S3,temp]=svd(reshape(permute(X,[3 1 2]),I3,I2*I1),'econ'); S3=diag(S3);
if nargin==2
    U1=U1(:,1:min(size_core(1),I2*I3));
    U2=U2(:,1:min(size_core(2),I1*I3));
    U3=U3(:,1:min(size_core(3),I1*I2));
end
S=tmprod(tmprod(tmprod(X,U1',1),U2',2),U3',3);
end

%*******************************************************************************
function X_out = tmprod(X,U,mode)
%TMPROD mode-n tensor-matrix product.
[I,J,K]=size(X);
[M,N]=size(U);
if (mode~=1) && (mode~=2) && (mode~=3)
    error('The input variable mode should be 1, 2 or 3')
end
if N~=size(X,mode)
    error(['The number of columns of the input matrix should be equal to dimension ',int2str(mode),' of the input tensor'])
end
if mode==1
    X_out = reshape(U*reshape(X,I,J*K) ,M,J,K);
elseif mode==2
    X_out = permute(reshape (U*reshape(permute(X,[2 1 3]),J,I*K), M,I,K),[2 1 3]);
elseif mode==3
    X_out = permute(reshape (U*reshape(permute(X,[3 1 2]),K,I*J), M,I,J),[2 3 1]);
end
end

% %*******************************************************************************
% function [X_mat]=tens2mat(X,mode)
% %TENS2MAT Matrix unfoldings of a 3rd order tensor X along a given mode
% % INPUTS: - X : tensor of size (IxJxK)
% %         - mode = 1 or 2 or 3
% % OUTPUTS: X_mat is the matrix unfolding representation of X
% % if mode==1:  X_mat is IKxJ  (i=1,...,I, is the slowly varying index)
% % if mode==2:  X_mat is JIxK  (j=1,...,J, is the slowly varying index)
% % if mode==3   X_mat is KJxI  (k=1,...,K  is the slowly varying index)
% [I,J,K]=size(X);
% if mode==1
%     X_mat=reshape(permute(X,[3 1 2]),I*K,J);
% elseif mode==2
%     X_mat=reshape(X,J*I,K);
% elseif mode==3
%     X_mat=reshape(permute(X,[2 3 1]),K*J,I);
% else
%     error('Input argument mode must be 1, 2 or 3');
% end
% end

%*******************************************************************************
function C = kr(A,B)
%KR Khatri-Rao product.
[I R1]=size(A); J=size(B,1);
C=zeros(I*J,R1);
for j=1:R1
    C(:,j)=reshape(B(:,j)*A(:,j).',I*J,1);
end
end

%*******************************************************************************
function Mat = kr_part(B,C,partB,partC)
%KR_PART Partition-Wise Kronecker product
[J M]=size(B);
[K N]=size(C);
if (sum(partB)~=M)
    error(['Error: a matrix with ',int2str(M),' columns can not be partitioned in such a way'])
end
if (sum(partC)~=N)
    error(['Error: a matrix with ',int2str(N),' columns can not be partitioned in such a way'])
end
if length(partB)~=length(partC)
    error('Error: the 2 input matrices do not have the same number of blocks')
end

indB=[0 cumsum(partB)];
indC=[0 cumsum(partC)];
indMat=[0 cumsum(partB.*partC)];

Mat=zeros(J*K,sum(partB.*partC));
for i=1:length(partC)
    Mat(:,indMat(i)+1:indMat(i+1))=fast_kron( B(:,indB(i)+1:indB(i+1)) , C(:,indC(i)+1:indC(i+1)));
end
end






