function g = bfilt_gray(f,r,a,b)
% f灰度图；r滤波半径；a全局方差；b局部方差
[x,~] = meshgrid(-r:r);
w1 = exp(-(x.^2)/(2*a^2));
%f = tofloat(f);%f = im2double(f);
 
 
[m,~] = size(f);
f_temp = padarray(f,[r r],'symmetric');
%g = zeros(m,n);

    for j = r+1:m+r
        temp = f_temp(j-r:j+r);
        w2 = exp(-(temp-f(j-r)).^2/(2*b^2));
        w = w1(1,:).*w2;
        s = temp.*w;
        g(j-r) = sum(s(:))/sum(w(:));
    end
    

