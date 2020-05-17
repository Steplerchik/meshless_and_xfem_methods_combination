function w = Ix_tip3(y_,nn,eta,Kmat,crxy,step)

F = @(x) fun3(x,y_,eta,Kmat,crxy,step);

aa=-1+1e-10;
a=1-1e-10;
h=(a-aa)/(2*nn); %определяем шаг

integ = F(aa);

% plot(aa,y_,'.')
for i=1:1:(nn-1) %сам алгоритм Симпсона

x=aa+2*h*i;

integ=integ+2*F(x)+4*F(x-h);

% plot(x-h,y_,'.')
% plot(x,y_,'.')
end

integ=integ+4*F(aa+h*(2*nn-1))+F(aa+h*2*nn);
% plot(aa+h*(2*nn-1),y_,'.')
% plot(aa+h*2*nn,y_,'.')
integ=h*integ/3;

w=integ;

end