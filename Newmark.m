function [A,B,D,E] = Newmark(dt,x,dx,ddx,beta,gamma)

A = 1/(beta*dt^2);
D = gamma/(beta*dt);
B = -A*(x+dt*dx+(0.5-beta)*dt^2*ddx);
E = dx+dt*((1-gamma)*ddx+gamma*B);