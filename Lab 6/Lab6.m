clear; clc; clf;
tic
Ra=1000; Pr=0.71;
N = 30; M = 30;
x = linspace(0,1,N+1); dx = x(2)-x(1);
y = linspace(0,1,M+1); dy = y(2)-y(1);
[Y,X] = meshgrid(y,x);
u = zeros(N+1,M+1); v = zeros(N+1,M+1); p = zeros(N+1,M+1); t=zeros(N+1,M+1);
t(1,:)=1;% boundary condition for u, the other conditions are given in initial values


err = 1;
while err > 1e-3
    dt = 0.5*min(0.5/Pr/(1/dx^2+1/dy^2),2*Pr/max(max(u.^2+v.^2)));
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
                +Pr*(vn(i,j-1)-2*vn(i,j)+vn(i,j+1))/dy^2....
                +Ra*Pr*tn(i,j)); % this is v*
            
            t(i,j) = tn(i,j) ...
                +dt*(-un(i,j)*(tn(i+1,j)-tn(i-1,j))/2/dx ...
                -vn(i,j)*(tn(i,j+1)-tn(i,j-1))/2/dy ...
                +(tn(i-1,j)-2*tn(i,j)+tn(i+1,j))/dx^2 ...
                +(tn(i,j-1)-2*tn(i,j)+tn(i,j+1))/dy^2); % this is theta*
        end
    end
    t(:,1) = t(:,2); t(:,N+1) = t(:,N);
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
    err = max(max(abs( (t-tn)./t )));
    %surf(X,Y,u), axis([0 1 0 1 0 0.08]); pause(0.02);
    quiver(X,Y,u,v,'r'); axis ([0 1 0 1]); %pause(0.01);
end
toc
answer=max(u(15,:))
answer1=max(u(16,:))
answer2=max(u(14,:))