clear
addpath('Generate_datacube');
load A.mat;
[Y_true, X] =  getSynData(A,17,0,0.8);
% SNR = 30; % [10 20 30 40]
% [h,w,L]=size(Y_true);
% N=h*w;
% std_noise = sqrt(sum(sum((A*X).^2))/N/L/10^(SNR/10));
% noise = std_noise*randn(L,N);
% Y_noise = A*X + noise;
% xx=reshape(X',64,64,6);
EndNum=6;
ranknumber=20;
 S=reshape(X',64,64,EndNum);
 for i=1:6
      [W,H]=nmf(S(:,:,i),ranknumber,'verbose',1,'method','hals');
      data2init{1}(:,(i-1)*ranknumber+1:i*ranknumber)=W;
      data2init{2}(:,(i-1)*ranknumber+1:i*ranknumber)=H';
 end
% save data2theta05 A X Y_true Y_noise177
 save data2init1 A X Y_true 