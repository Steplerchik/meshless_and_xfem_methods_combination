function w = weight_func(Xk,Yk,x,y,kol,a)

%support domain    
sd=domain(kol,a);

distance=sqrt((Xk-x)^2+(Yk-y)^2);

if (sd < distance)
    w(1,1)=0;
    w(2,1)=0;
    w(3,1)=0;
else
    
d=sd;
r=distance/d;
w(1,1)=(exp(-(r/0.3)^2)).*(r<=1);

if (r==0)
    w(2,1)=0;
    w(3,1)=0;
else
w(2,1)=((exp(-(r/0.3)^2))*(-2*(r/0.3))*(x-Xk)/(r*d^2)).*(r<=1);
w(3,1)=((exp(-(r/0.3)^2))*(-2*(r/0.3))*(y-Yk)/(r*d^2)).*(r<=1);
end

% w(1,1)=(2/3-4*r^2+4*r^3).*(r<=1/2)+(4/3-4*r+4*r^2-4/3*r^3).*(r>1/2).*(r<=1);
% if (r==0)
%     w(2,1)=0;
%     w(3,1)=0;
% else
% w(2,1)=((-8*r+12*r^2)*(x-Xk)/(r*d^2)).*(r<=1/2)+((-4+8*r-4*r^2)*(x-Xk)/(r*d^2)).*(r>1/2).*(r<=1);
% w(3,1)=((-8*r+12*r^2)*(y-Yk)/(r*d^2)).*(r<=1/2)+((-4+8*r-4*r^2)*(y-Yk)/(r*d^2)).*(r>1/2).*(r<=1);
% end


end
end