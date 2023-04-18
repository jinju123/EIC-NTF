function [mixed, abf] = getSynData(A, pure,beta,win)

%
% Generate synthetic data.
% The spectra of simulated data is obtained from the USGS library "signatures"
%
% Input
%   - A: matrix of reflectances
%   - win: size of smoothing filter
%   - pure: 0 - no pure pixels, 1 - exist pure pixel
%
% Output
%   - mixed: generated synthetic mixed data
%   - abf: actual abundance fractions
%
% The pure pixels can be removed by adding the following two lines
%        ----Index = ceil(find(abf>0.8)/c);
%        ----abf(:,Index) = 1/c*ones(c,1)*ones(1,length(Index));
%


[band, c] = size(A);
dim = 100;

%label = ones((dim/10)^2,1);
%label = ones(dim,1);
label = ones(10*10,1);
num = floor(length(label)/c);

for i=1:c-1
    label((i-1)*num+1:i*num) = (i+1); 
end
        
ridx = randperm(length(label));%Ëæ»ú´òÂÒ
label = label(ridx)';
label = reshape(label,ceil(dim/10),ceil(dim/10));
abf = zeros(dim,dim,c);
img = zeros(dim,dim);
for i=1:dim
    for j=1:dim
        for cls = 1:c
            if label(floor((i-1)/10)+1,floor((j-1)/10)+1) == cls
                tmp = zeros(c,1);
                tmp(cls) = beta;
                gf=randperm(c);
                gf(gf==cls)=[];
                clss=gf(randi(c-1));
                tmp(clss) = 1-beta;
%                 ggf=randperm(c);
%                 ggf(ggf==cls)=[];
%                 ggf(ggf==clss)=[];
%                 clsss=ggf(randi(c-2));
%                 tmp(clsss) = 1-1.5*beta;
                abf(i,j,:) = tmp;
                img(i,j) = c;
            end
        end
    end
end

%low pass filter
%H = ones(win,win)/(win*win);
%H=fspecial('gaussian',[win,win],2);
H=fspecial('gaussian',2);
img_fil = filter2(H,img);
for i=1:c
    abf(:,:,i) = filter2(H,abf(:,:,i));
end
%abf = abf(ceil(win/2):end+floor(win/2),ceil(win/2):end+floor(win/2),:);


% generate mixtures
[M,N,c] = size(abf);
abf = reshape(abf,M*N,c)';

% remove pure pixels
if pure == 0
    Index = ceil(find(abf>0.8)/c);
    abf(:,Index) = 1/c*ones(c,1)*ones(1,length(Index));
end

for i=1:dim^2
    
    if sum(abf(:,i))<1
    abf(:,i) = ones(c,1)*((1-sum(abf(:,i)))/4)+abf(:,i);
    end
end
mixed = reshape((A*abf)',M,N,band);

 