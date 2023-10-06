clear; clc;

N=80; M=80; nu=1;
x=linspace(-2,4,N+1); dx=x(2)-x(1);
y=linspace(-2,4,M+1); dy=y(2)-y(1);

[Y,X]=meshgrid(y,x);
u=ones(N+1,M+1);

for j=1:M+1
    for i=1:N+1
        if x(i)>=0.5 && x(i)<=1 && y(j)>=0.5 && y(j)<=1
            u(i,j)=1000;
        end
    end 
end
dt=0.5/(nu*(1/dx^2+1/dy^2));

for n=1:ceil(0.5/dt)
    un=u; u(1,:)=1; u(:,1)=1;u(N+1,:)=1;u(:,M+1)=1;
    for j=2:M
        for i=2:N
            u(i,j)=un(i,j)....
                +(nu*dt/dx^2)*(un(i-1,j)-2*un(i,j)+un(i+1,j))...
                +(nu*dt/dy^2)*(un(i,j-1)-2*un(i,j)+un(i,j+1));
        end 
    end
    surf(X,Y,u),axis([-2 4 -2 4 1 2]);pause(0.2);
end