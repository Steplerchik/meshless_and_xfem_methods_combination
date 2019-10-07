function w = Iy_tip3(nn,eta,Kmat,crxy,step)

F= @(y) Ix_tip3(y,nn,eta,Kmat,crxy,step);

aa=-1+1e-10;
b=1-1e-10;
h=(b-aa)/(2*nn); %определяем шаг

integ = F(aa);

for i=1:1:(nn-1) %сам алгоритм Симпсона

y=aa+2*h*i;

integ=integ+2*F(y)+4*F(y-h);

end

integ=integ+4*F(aa+h*(2*nn-1))+F(aa+h*2*nn);

integ=h*integ/3;

w=integ;


end