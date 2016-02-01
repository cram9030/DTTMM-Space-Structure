t = cputime;
[T,Y] = ode45(@GFEM,[0 1e3],zeros(length(A),1));
e = cputime-t