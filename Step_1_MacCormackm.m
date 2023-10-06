clear; clc;
N = 81; c=1;
x = linspace(0,2,N+1); 
dx= x(2)-x(1);
dt=1.*dx/c;; 
tf=3;
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
		ubar(i)=un(i)-(c*dt/dx)*(un(i)-un(i-1));

    end
    ubar(1)=ubar(N+1);
    for i=1:N
        ubarbar(i)=un(i)-(c*dt/dx)*(ubar(i+1)-ubar(i));
    end
    ubarbar(N+1)=ubarbar(1);
 
    u=0.5*(ubar+ubarbar);

    plot(x,u,'-o')
    pause(0.2)
end
