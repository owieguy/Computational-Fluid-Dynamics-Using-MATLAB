clear; clc; clf;
tic
L=1; k=50; pi=20; pe=0
N = 40; M = 40;
x = linspace(0,1,N+1); dx = x(2)-x(1);
y = linspace(0,1,M+1); dy = y(2)-y(1);
A = linspace(0.5,0.1,N+1);
[Y,X] = meshgrid(y,x);
u = zeros(N+1,M+1); v = zeros(N+1,M+1); p = zeros(N+1,M+1);


err = 1;
while err > 1e-4
    dt = min(0.5/Pr/(1/dx^2+1/dy^2),2*Pr/max(max(u.^2+v.^2)));
    un = u; vn = v;
    for j = 2 : M
        for i = 2 : N
            u(i,j) = un(i,j)*(1-k*dt*un(i,j)); %This is u*
        end
    end
    u(1,:)=u(2,:); 
    u(N+1,:)=u(N,:);
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
        p(N+1,:)=0; p(:,M+1) = p(:,M);
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
    quiver(X,Y,u,v,'r'); axis ([0 1 0 1]); pause(0.02);
end
toc
hold on
plot(u(1,:),y,'o',y/2.*(1-y),y)
axis equal
axis([0 1 0 1]);
hold off