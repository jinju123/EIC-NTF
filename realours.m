addpath('supplement');
addpath('data');
addpath initial
%addpath fast_strategy
%load jjdata1.mat
%load dc1Initrank40.mat
%load xiong_25dB.mat;
load JasperRidge.mat;
%load syntheticDataIter1rank19.mat;
load jasperinitsum1_95;
%load jasinittv;
%load jasvca;
realEndmember=M;
imgsize=100;
%realAbundance=X';
%Y=x_n;
%  [w Rw] = estNoise(Y,'additive');
%  [~, E]=hysime(Y,w,Rw);
%  [E,~,~]= svd(Y,'econ');
  %[E,~,~]= svd(Y,'econ');
 Y=(Y-min(min(Y)))/(max(max(Y))-min(min(Y)));
 %Y(Y<0.38)=0.38;
tensorData=reshape( Y',100,100,198);
% btdInit=temp.btdInit; 
% btdInit{1}=max(Jasinit30{1},1e-4);
% btdInit{2}=max(Jasinit30{2},1e-4);
% btdInit{3}=max(Jasinit30{3},1e-4);
EndNum=4;
rankNumber=95;
% [A_VCA,~] = hyperVca(Y,EndNum);
% S_init = fcls(A_VCA,Y);
%S=reshape((Y-A_VCA*S_init)',100,100,198);
%initialation
btdInit{1}=max(Jasinit30{1},1e-4);
btdInit{2}=max(Jasinit30{2},1e-4);
%btdInit{2}=rand(100,EndNum*rankNumber);
 A_rand=rand(198,EndNum);
for jjj=1:4
    A_smooth(:,jjj)=bfilt_gray(A_rand(:,jjj),3,1,1);
end
btdInit{3}=A_smooth;
allSad=zeros(1,10);
allRmse=zeros(1,10);


%our mvntf
% ouroptions.convergeNum=3000;
% ouroptions.derta=0.4;
% ouroptions.lambda=1;
% tours=new_mvntf(tensorData,EndNum,rankNumber,btdInit{1},btdInit{2},btdInit{3},ouroptions,realEndmember);
% ccour=tours{1};
for i=1:4
     ww=btdInit{1}(:,(i-1)*95+1:i*95)*btdInit{2}(:,(i-1)*95+1:i*95)';
      WI(i,:)=reshape(ww,1,10000); 
end
Ar=reshape(WI',100,100,4);
for i=1:4
      A(:,(i-1)*100+1:i*100)=Ar(:,:,i);
end

% ours
ouroptions.convergeNum=1000;
ouroptions.derta=5;
ouroptions.lambda1=3;
ouroptions.lambda2=1;
ouroptions.mu=0.01;
tours=EIC_NTF(tensorData,EndNum,A,btdInit{3},ouroptions,realEndmember);
ccour=tours{1};
 

%mvntf
options.convergeNum=1000;
options.derta=3;
t=mvntf(tensorData,EndNum,rankNumber,btdInit{1},btdInit{2},btdInit{3},options,realEndmember);
cc=t{1};

%mvntf_tv
options.convergeNum=1000;
options.derta=3;
options.lambda=1;
options.mu=0.01;
ttv=mvntftv(tensorData,EndNum,rankNumber,btdInit{1},btdInit{2},btdInit{3},options,M);
cctv=ttv{1};
              
           
 %ourresult
  fprintf('our result is: ');
[oursad,ourDistance,oursor]=cosDistance(ccour,realEndmember);
fprintf('the sad is :[%d]\n',oursad);  

 %mvntf result   
fprintf('mvntf result is:' );
[sadmvntf,mvntfDistance,mvntfsor]=cosDistance(cc,realEndmember);
fprintf('the sad is :[%d]\n',sadmvntf);

%mvntf_tv result
fprintf('mvntf_tv result is:' );
[sadtv,tvDistance,tvsor]=cosDistance(cctv,realEndmember);
fprintf('the sad is :[%d] \n',sadtv);

  
    