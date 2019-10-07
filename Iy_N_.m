function w = Iy_N(n,P,X,Y,kol,b,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)

for i=1:length(Y)-1
    if (Y(i)<Yc(1))&&(Y(i+1)>Yc(1))
        step=Y(i+1)-Y(i);
        ycen=Y(i)+step/2;
    end
end
x_=0;
F= @(y) N(n,x_,y,P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc);

nn=b/step; %количество частей деления

aa=0;
h=(b/a*kol*step-aa)/(2*nn); %определяем шаг

integ = F(aa);

for i=1:1:(nn-1) %сам алгоритм Симпсона

y=aa+2*h*i;

integ=integ+2*F(y)+4*F(y-h);

end

integ=integ+4*F(aa+h*(2*nn-1))+F(aa+h*2*nn);

integ=h*integ/3;

w=integ';

end