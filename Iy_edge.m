function w = Iy_edge(step,ycen_d,ycen_u,n,P,X,Y,D,a,kol,b,mas,i_edge,i_tip1,i_tip2,Xc,Yc)

F= @(y) Ix_edge(n,y,P,X,Y,D,a,kol,step,i_edge,i_tip1,i_tip2,Xc,Yc);

nn=round((ycen_u-ycen_d+10*step)/step)*2; %количество частей деления

aa=ycen_d-5*step;
h=(ycen_u+5*step-aa)/(2*nn); %определяем шаг

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