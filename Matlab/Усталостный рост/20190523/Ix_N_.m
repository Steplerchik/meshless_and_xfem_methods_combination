function w = Ix_N(n,P,X,Y,kol,b,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)

for i=1:length(Y)-1
    if (Y(i)<Yc(1))&&(Y(i+1)>Yc(1))
        step=Y(i+1)-Y(i);
        ycen=Y(i)+step/2;
    end
end
y_=0;
F= @(x) N(n,x,y_,P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc);

nn=kol; %количество частей деления

aa=0;
h=(kol*step-aa)/(2*nn); %определяем шаг

integ = F(aa);

for i=1:1:(nn-1) %сам алгоритм Симпсона

x=aa+2*h*i;

integ=integ+2*F(x)+4*F(x-h);

end

integ=integ+4*F(aa+h*(2*nn-1))+F(aa+h*2*nn);

integ=h*integ/3;

w=integ';

end