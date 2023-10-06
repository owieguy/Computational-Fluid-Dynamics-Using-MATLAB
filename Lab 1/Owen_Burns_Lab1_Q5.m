clear; clc;
N = 41; 
x = linspace(0,2,N+1); 
u=sin(pi*x);
dx= x(2)-x(1);
dt=dx;
tf=30*dt;
M=tf/dt;
%++++++++++++++++++++++++++++++++++++++++++++++++
for n=1:M
    un=u;
    for i=2:N
        if un(i)>0
            u(i)=un(i)-(un(i)*dt/dx)*(un(i)-un(i-1));
        else
            u(i)=un(i)-(un(i)*dt/dx)*(un(i+1)-un(i));
        end
    end
    
    u(1)=0; 
    u(N)=0;
    plot(x,u,'-o')
    pause(0.2)
end