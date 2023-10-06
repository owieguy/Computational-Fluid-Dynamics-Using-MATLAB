R = 0.1;      %m
k = 20;       %W/mK
S = 1000;     %W/m^3
h = 50;       %W/m^2K
Tw = 20;      %C

N = 81;
r = linspace(0,R,N+1); 
dr = r(2) - r(1);
dt = 0.025*dr^2; 
tf = 5000000*dt;
M = tf/dt;

u = 20*ones(N+1,1);

for n=1:M
	un=u;
    for i=2:N
		u(i)=un(i)+dt*(k*(un(i+1)-un(i-1))/(2*r(i)*dr)+k*(un(i+1)-2*un(i)+un(i-1))/dr^2+S);
    end
    u(1)=u(2);
    %u(N+1)=(k*u(N)+Tw*dr*h)/(dr*h+k);
    u(N+1)=u(N)-(dr*h/k)*(u(N)-Tw);
    %u(N+1)=(Tw+(k/(h*dr)))/(1+k/(h*dr));
    %u(N+1)=(k*u(N)+Tw)/(dr*h+k);
end
plot(r,u,'-o')
Final=u(1)