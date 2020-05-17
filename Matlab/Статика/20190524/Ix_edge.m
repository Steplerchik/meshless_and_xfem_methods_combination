function w = Ix_edge(n,y_,P,X,Y,D,a,kol,step,i_edge,i_tip1,i_tip2,Xc,Yc)

F= @(x) K_mat(n,x,y_,P,X,Y,D,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc);

nn=1*2; %количество частей деления

aa=5*step;
h=(6*step-aa)/(2*nn); %определяем шаг

integ = F(aa);

plot(aa,y_,'.')
for i=1:1:(nn-1) %сам алгоритм Симпсона

x=aa+2*h*i;

integ=integ+2*F(x)+4*F(x-h);

plot(x-h,y_,'.')
plot(x,y_,'.')
end

integ=integ+4*F(aa+h*(2*nn-1))+F(aa+h*2*nn);
plot(aa+h*(2*nn-1),y_,'.')
plot(aa+h*2*nn,y_,'.')
integ=h*integ/3;

w=integ;

end