function w = Iy_reg_add1(step,ycen1,ycen2,n,P,X,Y,D,a,kol,b,mas,i_edge,i_tip1,i_tip2,Xc,Yc)

F= @(y) Ix_reg_add1(n,y,P,X,Y,D,a,kol,step,i_edge,i_tip1,i_tip2,Xc,Yc);

if ycen1<=ycen2
   nn=round((ycen2+5*step-(ycen1+5*step))/step)*2; %количество частей деления
   aa=ycen1+5*step;
   h=(ycen2+5*step-aa)/(2*nn); %определяем шаг
else
   nn=round((ycen1-5*step-(ycen2-5*step))/step)*2; %количество частей деления
   aa=ycen2-5*step;
   h=(ycen1-5*step-aa)/(2*nn); %определяем шаг
end

integ = F(aa);

for i=1:1:(nn-1) %сам алгоритм Симпсона
i
y=aa+2*h*i;

integ=integ+2*F(y)+4*F(y-h);

end

integ=integ+4*F(aa+h*(2*nn-1))+F(aa+h*2*nn);

integ=h*integ/3;

w=integ;


end