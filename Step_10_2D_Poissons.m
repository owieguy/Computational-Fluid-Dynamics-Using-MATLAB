clear; clc;
tic
N=40; M=40;
x=linspace(0,1,N+1); dx=x(2)-x(1);
y=linspace(0,1,M+1); dy=y(2)-y(1);

[Y,X]=meshgrid(y,x);
u=zeros(N+1,M+1);

err = 1;
while err>1e-8
    u_old=u;
    for j=2:M
        for i=2:N
            u(i,j)=(dx^2.....
                + u(i-1,j) + u(i+1,j)...
                + u(i,j-1) + u(i,j+1))*0.25;
        end 
    end
    err=max(max(abs((u-u_old)./u)));
    %surf(X,Y,u),axis([0 1 0 1 0 0.08]);pause(0.02);
end
toc
surf(X,Y,u),axis([0 1 0 1 0 0.08]);pause(0.02);