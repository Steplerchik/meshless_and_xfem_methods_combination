function w = Ix_tip1(y_,nn,eta,Kmat,crxy,step)

F = @(x) fun1(x,y_,eta,Kmat,crxy,step);

aa=-1+1e-10;
a=1-1e-10;
h=(a-aa)/(2*nn); %определяем шаг

integ = F(aa);

for i=1:1:(nn-1) %сам алгоритм Симпсона

x=aa+2*h*i;

integ=integ+2*F(x)+4*F(x-h);

end

integ=integ+4*F(aa+h*(2*nn-1))+F(aa+h*2*nn);

integ=h*integ/3;
w=integ;
end