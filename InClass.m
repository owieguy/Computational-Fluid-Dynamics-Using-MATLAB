clear; clc;
N = 141; nu=1;
x = linspace(0,2,N+1); 
dx= x(2)-x(1);
dt=1.5*dx/c;
tf=1.5;
M=tf/dt;
%++++++++++++++++++++++++++++++++++++++++++++++++
u = ones(N+1,1);
for i=1:N+1
	if x(i)>0.5 && x(i)<1
		u(i)=2;
	end
end
%++++++++++++++++++++++++++++++++++++++++++++++++
for n=1:M
	un=u;
    for i=2:N+1
		u(i)=un(i)-(c*dt/dx)*(un(i)-un(i-1));
    end
    u(1)=1;
    plot(x,u,'-o')
    pause(0.2)
end