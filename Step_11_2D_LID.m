clear; clc; clf;
tic
Ra=1000; Pr=0.71;
N = 40; M = 40;
x = linspace(0,1,N+1); dx = x(2)-x(1);
y = linspace(0,1,M+1); dy = y(2)-y(1);
[Y,X] = meshgrid(y,x);
u = zeros(N+1,M+1); v = zeros(N+1,M+1); p = zeros(N+1,M+1);t=zeros(N+1,M+1);
t(1,:)=1;% boundary condition for u, the other conditions are given in initial values


err = 1;
while err > 1e-3
    dt = min(0.5/nu/(1/dx^2+1/dy^2),2*nu/max(max(u.^2+v.^2)));
    un = u; vn = v; tn = t;
    for j = 2 : M
        for i = 2 : N
            u(i,j) = un(i,j) ...
                +dt*(-un(i,j)*(un(i+1,j)-un(i-1,j))/2/dx ...
                -vn(i,j)*(un(i,j+1)-un(i,j-1))/2/dy ...
                +Pr*(un(i-1,j)-2*un(i,j)+un(i+1,j))/dx^2 ...
                +Pr*(un(i,j-1)-2*un(i,j)+un(i,j+1))/dy^2); % this is u*
            
            v(i,j) = vn(i,j) ...
                +dt*(-un(i,j)*(vn(i+1,j)-vn(i-1,j))/2/dx ...
                -vn(i,j)*(vn(i,j+1)-vn(i,j-1))/2/dy ...
                +Pr*(vn(i-1,j)-2*vn(i,j)+vn(i+1,j))/dx^2 ...
                +Pr*(vn(i,j-1)-2*vn(i,j)+vn(i,j+1))/dy^2)+Ra*Pr*t; % this is v*
            
            t(i,j) = tn(i,j) ...
                +dt*(-un(i,j)*(tn(i+1,j)-tn(i-1,j))/2/dx ...
                -vn(i,j)*(tn(i,j+1)-tn(i,j-1))/2/dy ...
                +(tn(i-1,j)-2*tn(i,j)+tn(i+1,j))/dx^2 ...
                +(tn(i,j-1)-2*tn(i,j)+tn(i,j+1))/dy^2); % this is theta*
        end
    end
    t(:,1) = t(:,2); p(:,N+1) = p(:,N);
    % now solve PPE
    err_p = 1; 
    while err_p > 1e-3
        p_old = p;
        for j = 2 : M
            for i = 2 : N
                p(i,j) = (p(i-1,j)+p(i+1,j)+...
                    p(i,j-1)+p(i,j+1) ...
                    -(dx/2/dt)*(u(i+1,j)-u(i-1,j) ...
                    +v(i,j+1)-v(i,j-1)))/4;
            end
        end
        p(1,:) = p(2,:); p(:,1) = p(:,2);
        p(N+1,:)=p(N,:); p(:,M+1) = p(:,M);
        err_p = max(max(abs( (p-p_old)./p)));
    end
    % now I will correct my u and v
    for j = 2 : M
        for i = 2 : N
            u(i,j) = u(i,j) - dt*(p(i+1,j)-p(i-1,j))/2/dx;
            v(i,j) = v(i,j) - dt*(p(i,j+1)-p(i,j-1))/2/dy;
        end
    end
    err = max(max(abs( (u-un)./u )));
    %surf(X,Y,u), axis([0 1 0 1 0 0.08]); pause(0.02);
    quiver(X,Y,u,v,'r'); axis ([0 1 0 1]); pause(0.01);
    answer=max(0.5,j)
end
toc
hold on
plot(0.5+0.5*u(ceil(N/2),:),y)
plot([0.5 0.5],[0, 1]);
plot([0, 1],[0.5, 0.5])
plot(x,0.5*v(:,ceil(M/2))+0.5);
% Ghia Data Re = 100
yyy = [0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1];
xxx=[0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1];
uuu=0.01*[0,-3.717,-4.192,-4.775,-6.434,-10.15,-15.662,-21.09,-20.581,-13.641,0.332,23.151,68.717,73.722,78.871,84.123,100];
vvv=0.01*[0,9.233,10.091,10.89,12.317,16.077,17.507,17.527,5.454,-24.533,-22.445,-16.914,-10.313,-8.864,-7.391,-5.906,0];
plot(0.5+0.5*uuu,yyy,'o',xxx,0.5+0.5*vvv,'s')
axis equal
axis([0 1 0 1]);
hold off