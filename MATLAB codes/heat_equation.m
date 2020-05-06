function heat_equation(n,T)
% function Example4_1(n,T)
% solve Ut=Uxx, U(0,t)=U(1,t)=0, U(0,x)=sin(3*pi*x).
%    n: n+1 is the number of intervals in x-direction, 
%    T: the no. of intervals in t-dir.

clf,

t=0.1; dt=t/T; dx=1/(n+1); 
r=dt/(dx*dx);
x=[dx:dx:1-dx]; % all the inner points

%analytical solution when t=0.1:
ex=[0:0.005:1];
u_ex=exp(-4*(pi*pi)*t)*sin(2*pi*ex);

%up is the previous solution and un is the next solution.
up=zeros(size(x)); un=zeros(size(x));

%initial condition:
up=sin(2*pi*x);

%the calculating loop:
for m=1:T,
  for j=2:n-1,
    un(j) = r*up(j-1)+(1-2*r)*up(j)+r*up(j+1);
  end,
%take care of boundary conditions:
  un(1)=(1-2*r)*up(1)+r*up(2); 
  un(n)=r*up(n-1)+(1-2*r)*up(n);
%put un into up:
  up=un;
end,

p=plot(ex,u_ex,[0 x 1],[0 up 0],'--',[0 x 1],[0 up 0],'o');
h=text(0.01,-0.017,sprintf('dx=%6.3g, dt=%g, r=%g',dx,dt,r)); 
set(p,'LineWidth',1.2), set(h,'FontSize',18),

print -deps Example4_1.eps