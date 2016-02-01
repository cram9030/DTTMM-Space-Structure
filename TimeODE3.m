t = cputime;
[YE3] = ode3(@GFEM,0:1e-5:1e-3,zeros(80,1));
e = cputime-t