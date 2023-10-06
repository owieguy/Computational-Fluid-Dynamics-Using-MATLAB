clear; clc;

N=25; M=25; cx=1; cy=1;
x=linspace(0,2,N+1); dx=x(2)-x(1);
y=linspace(0,2,M+1); dy=y(2)-y(1);

[Y,X]=meshgrid(y,x);
u=ones(N+1,M+1);

for j=1:M+1
    for i=1:N+1
        if x(i)>=0.5 && x(i)<=1 && y(j)>=0.5 && y(j)<=1
            u(i,j)=2;
        end
    end 
end
dt=min(0.5*dx/cx,0.5*dy/cy);

for n=1:ceil(0.5/dt)
    un=u; u(1,:)=1; u(:,1)=1;
    for j=2:M+1
        for i=2:N+1
            u(i,j)=un(i,j)....
                -(cx*dt/dx)*(un(i,j)-un(i-1,j))...
                -(cy*dt/dy)*(un(i,j)-un(i,j-1));
        end 
    end
    surf(X,Y,u), pause(0.2)
end