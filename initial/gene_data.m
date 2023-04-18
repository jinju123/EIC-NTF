[Y_true, X] =  getSynData(M,1,0.8,17);
SNR = 30; % [10 20 30 40]
[h,w,L]=size(Y_true);
N=h*w;
std_noise = sqrt(sum(sum((A*X).^2))/N/L/10^(SNR/10));
noise = std_noise*randn(L,N);
Y_noise = A*X + noise;
 X(:,93026:307*307)=X(:,34:(307*307-93025)+33);
 ranknumber=177;
 S=reshape(X',307,307,4);
 for i=1:4
      [W,H]=nmf(S(:,:,i),ranknumber,'verbose',1,'method','hals');
      urbaninit30{1}(:,(i-1)*ranknumber+1:i*ranknumber)=W;
      urbaninit30{2}(:,(i-1)*ranknumber+1:i*ranknumber)=H';
 end
 A_rand=rand(162,4);
 for jjj=1:4
 A_init(:,jjj)=bfilt_gray(A_rand(:,jjj),3,1,1);
 end
 save urban3 X S urbaninit30 A_init