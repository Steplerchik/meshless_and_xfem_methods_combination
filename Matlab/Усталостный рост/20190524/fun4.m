function w = fun4(nu1,nu2,eta,Kmat,crxy,step)

dz1=1/4*(1+nu1)*(1+nu2);
dz2=1/4*(1-nu1)*(1+nu2);

ks2=-(1+(1+eta(2))*(dz1+dz2-1));
ks1=1+(eta(1)-1)/(-eta(2)-1)*(-ks2-1)-2*dz1;

x=crxy(1) + ks1*step;
y=crxy(2) + ks2*step;

j1 = (1+nu2)/8;
j2 = 2*(1+eta(2));
j12= j1*j2;

plot(x,y,'.')

w = Kmat(x,y)*step^2*j12;

end