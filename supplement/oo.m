
clc
f0=A(:,5);
%f0=f0(:,:,1);
[m,n]=size(f0);
f0=double(f0);

lamda=0.005; % smoothness paramter, the larger the smoother
tao=0.125; % fixed do not change it.

p1=zeros(m,n);
p2=zeros(m,n);

tic
for step=1:100
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    div_p=div(p1,p2);
    %cx=Fx(div_p-f0/lamda);
    cy=Fy(div_p-f0/lamda);
absc=sqrt(cy.^2);
%p1=(p1+tao*cx)./(1+tao*abs_c);
p2=(p2+tao*cy)./(1+tao*absc);

end

u=f0-lamda*div_p;
toc
% figure; imagesc(f0); colormap(gray); axis off; axis equal;
% figure; imagesc(u); colormap(gray); axis off; axis equal;

% Compute divergence using backward derivative
function f = div(a,b)
f = By(b);
end

% Forward derivative operator on x with boundary condition u(:,:,1)=u(:,:,1)
function Fxu = Fx(u)
[m,n] = size(u);
Fxu = circshift(u,[0 -1])-u;
Fxu(:,n) = zeros(m,1);
end

% Forward derivative operator on y with boundary condition u(1,:,:)=u(m,:,:)
function Fyu = Fy(u)
[m,n] = size(u);
Fyu = circshift(u,[-1 0])-u;
Fyu(m,:) = zeros(1,n);
end

% Backward derivative operator on x with boundary condition Bxu(:,1)=u(:,1)
function Bxu = Bx(u)
[~,n] = size(u);
Bxu = u - circshift(u,[0 1]);
Bxu(:,1) = u(:,1);
Bxu(:,n) = -u(:,n-1);
end

% Backward derivative operator on y with boundary condition Bxu(1,:)=u(1,:)
function Byu = By(u)
[m,~] = size(u);
Byu = u - circshift(u,[1 0]);
Byu(1,:) = u(1,:);
Byu(m,:) = -u(m-1,:);
end