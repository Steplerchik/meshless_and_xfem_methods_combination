function w = Iy_reg_dwn(step,ycen,n,P,X,Y,D,a,kol,b,mas,i_edge,i_tip1,i_tip2,Xc,Yc)

F= @(y) Ix_reg(n,y,P,X,Y,D,a,kol,step,i_edge,i_tip1,i_tip2,Xc,Yc);

nn=round((ycen-2.5*step)/step); %количество частей деления

aa=0;
h=(ycen-2.5*step-aa)/(2*nn); %определяем шаг

integ = F(aa);

for i=1:1:(nn-1+1e-10) %сам алгоритм Симпсона
i
y=aa+2*h*i;

integ=integ+2*F(y)+4*F(y-h);

end

integ=integ+4*F(aa+h*(2*nn-1))+F(aa+h*2*nn);

integ=h*integ/3;

w=integ;


end