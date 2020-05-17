function w = N(n,x,y,P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc)

W=zeros(2*n,2*n);
k=0;
check=0;
for i=1:2*n 
     if (check==0)
        k=k+1;
        omega=weight_func(X(k),Y(k),x,y,kol,a);
      check=1;
     else
        check=0;
     end
     W(i,i)=omega(1,1);
end    
A=P'*W*P;
B=P'*W;
U=chol(A);
A_1=inv(U)*inv(U');

ppx=[1 x y 0 0 0];
ppy=[0 0 0 1 x y];

N(1,:)=ppx*A_1*B;
N(2,:)=ppy*A_1*B;

zn=crack_side(x,y,Xc,Yc);
n_=0;
N_e = [];
for i = i_edge
        n_=n_+1;
        N_e(1,n_)=N(1,2*i-1)*zn;
        N_e(2,n_)=N(2,2*i-1)*zn;
        n_=n_+1;
        N_e(1,n_)=N(1,2*i)*zn;
        N_e(2,n_)=N(2,2*i)*zn;
end

if length(Xc) > 2
    v11 = [Xc(1)-Xc(2) Yc(1)-Yc(2)];
    if Yc(1)-Yc(2) > 0
        alfa1 = acos((Xc(1)-Xc(2))/norm(v11));
    else
        alfa1 = -acos((Xc(1)-Xc(2))/norm(v11));
    end
    
    v12 = [Xc(3)-Xc(2) Yc(3)-Yc(2)];
    if Yc(3)-Yc(2) > 0
        beta1 = acos((Xc(3)-Xc(2))/norm(v12));
    else
        beta1 = -acos((Xc(3)-Xc(2))/norm(v12));
    end
    gamma1 = (alfa1+beta1)/2;
    c1 = tan(gamma1);
    d1 = Yc(2)-c1*Xc(2);
    if ((c1 < 0)&&(y < c1*x+d1))||((c1 > 0)&&(y > c1*x+d1))
        key1 = 0;
    else
        key1 = 1;
    end
else
    key1 = 1;
end

step = a/kol;
if key1 == 1
par = polar1(Xc(1),Xc(1)-Xc(2),Yc(1),Yc(1)-Yc(2),x,y);
fi=par(5);
n_=0;
else
%     xcf = Xc(2)+0.5*step*cos(atan((Yc(2)-Yc(3))/(Xc(2)-Xc(3))));
%     ycf = Yc(2)+0.5*step*sin(atan((Yc(2)-Yc(3))/(Xc(2)-Xc(3))));
% par = polar1(xcf,xcf-Xc(2),ycf,ycf-Yc(2),x,y);
% fi=par(5);
% n_=0;
par = polar1(Xc(1),Xc(1)-x,Yc(1),Yc(1)-y,x,y);
fi=par(5);
n_=0;
end
for i = i_tip1
        n_=n_+1;
        N_t1(1,n_)=N(1,2*i-1) * par(1)*cos(fi);
        N_t1(2,n_)=N(2,2*i) * par(1)*sin(fi);
        n_=n_+1;
        N_t1(1,n_)=N(1,2*i-1) * par(1)*(-sin(fi));
        N_t1(2,n_)=N(2,2*i) * par(1)*cos(fi);
        
        n_=n_+1;
        N_t1(1,n_)=N(1,2*i-1) * par(2)*cos(fi);
        N_t1(2,n_)=N(2,2*i) * par(2)*sin(fi);
        n_=n_+1;
        N_t1(1,n_)=N(1,2*i-1) * par(2)*(-sin(fi));
        N_t1(2,n_)=N(2,2*i) * par(2)*cos(fi);
        
        n_=n_+1;
        N_t1(1,n_)=N(1,2*i-1) * par(3)*cos(fi);
        N_t1(2,n_)=N(2,2*i) * par(3)*sin(fi);
        n_=n_+1;
        N_t1(1,n_)=N(1,2*i-1) * par(3)*(-sin(fi));
        N_t1(2,n_)=N(2,2*i) * par(3)*cos(fi);
        
        n_=n_+1;
        N_t1(1,n_)=N(1,2*i-1) * par(4)*cos(fi);
        N_t1(2,n_)=N(2,2*i) * par(4)*sin(fi);
        n_=n_+1;
        N_t1(1,n_)=N(1,2*i-1) * par(4)*(-sin(fi));
        N_t1(2,n_)=N(2,2*i) * par(4)*cos(fi);
end
    

lx = length(Xc);
if length(Xc) > 2
    v21 = [Xc(lx-2)-Xc(lx-1) Yc(lx-2)-Yc(lx-1)];
    if Yc(lx-2)-Yc(lx-1) > 0
        alfa2 = acos((Xc(lx-2)-Xc(lx-1))/norm(v21));
    else
        alfa2 = -acos((Xc(lx-2)-Xc(lx-1))/norm(v21));
    end
    
    v22 = [Xc(lx)-Xc(lx-1) Yc(lx)-Yc(lx-1)];
    if Yc(lx)-Yc(lx-1) > 0
        beta2 = acos((Xc(lx)-Xc(lx-1))/norm(v22));
    else
        beta2 = -acos((Xc(lx)-Xc(lx-1))/norm(v22));
    end
    gamma2 = (alfa2+beta2)/2;
    c2 = tan(gamma2);
    d2 = Yc(lx-1)-c2*Xc(lx-1);
    if ((c2 < 0)&&(y > c2*x+d2))||((c2 > 0)&&(y < c2*x+d2))
        key2 = 0;
    else
        key2 = 1;
    end
else
    key2 = 1;
end
% plot([Xc(length(Xc)-1) Xc(length(Xc)-1)-step],[c2*Xc(length(Xc)-1)+d2 c2*(Xc(length(Xc)-1)-step)+d2])
% plot([Xc(2) Xc(2)-step],[c1*Xc(2)+d1 c1*(Xc(2)-step)+d1])

if key2 == 1
par = polar1(Xc(length(Xc)),Xc(length(Xc))-Xc(length(Xc)-1),Yc(length(Yc)),Yc(length(Yc))-Yc(length(Xc)-1),x,y);
fi=par(5);
n_=0;
else
%     xcfe = Xc(lx-1)-0.5*step*cos(atan((Yc(lx-2)-Yc(lx-1))/(Xc(lx-2)-Xc(lx-1))));
%     ycfe = Yc(lx-1)-0.5*step*sin(atan((Yc(lx-2)-Yc(lx-1))/(Xc(lx-2)-Xc(lx-1))));
% par = polar1(xcfe,xcfe-Xc(lx-1),ycfe,ycfe-Yc(lx-1),x,y);
% fi=par(5);
% n_=0;
par = polar1(Xc(length(Xc)),Xc(length(Xc))-x,Yc(length(Yc)),Yc(length(Yc))-y,x,y);
fi=par(5);
n_=0;
end
for i = i_tip2
        n_=n_+1;
        N_t2(1,n_)=N(1,2*i-1) * par(1)*cos(fi);
        N_t2(2,n_)=N(2,2*i) * par(1)*sin(fi);
        n_=n_+1;
        N_t2(1,n_)=N(1,2*i-1) * par(1)*(-sin(fi));
        N_t2(2,n_)=N(2,2*i) * par(1)*cos(fi);
        
        n_=n_+1;
        N_t2(1,n_)=N(1,2*i-1) * par(2)*cos(fi);
        N_t2(2,n_)=N(2,2*i) * par(2)*sin(fi);
        n_=n_+1;
        N_t2(1,n_)=N(1,2*i-1) * par(2)*(-sin(fi));
        N_t2(2,n_)=N(2,2*i) * par(2)*cos(fi);
        
        n_=n_+1;
        N_t2(1,n_)=N(1,2*i-1) * par(3)*cos(fi);
        N_t2(2,n_)=N(2,2*i) * par(3)*sin(fi);
        n_=n_+1;
        N_t2(1,n_)=N(1,2*i-1) * par(3)*(-sin(fi));
        N_t2(2,n_)=N(2,2*i) * par(3)*cos(fi);
        
        n_=n_+1;
        N_t2(1,n_)=N(1,2*i-1) * par(4)*cos(fi);
        N_t2(2,n_)=N(2,2*i) * par(4)*sin(fi);
        n_=n_+1;
        N_t2(1,n_)=N(1,2*i-1) * par(4)*(-sin(fi));
        N_t2(2,n_)=N(2,2*i) * par(4)*cos(fi);
end

N_new=[N,N_e,N_t1,N_t2];

w=N_new;

end