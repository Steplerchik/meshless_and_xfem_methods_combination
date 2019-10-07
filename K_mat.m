function w = K_mat(n,x,y,P,X,Y,D,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc)

W=zeros(2*n,2*n);
Wdx=zeros(2*n,2*n);
Wdy=zeros(2*n,2*n);
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
     Wdx(i,i)=omega(2,1);
     Wdy(i,i)=omega(3,1);
end    
A=P'*W*P;
B=P'*W;
Adx=P'*Wdx*P;
Bdx=P'*Wdx;
Ady=P'*Wdy*P;
Bdy=P'*Wdy;

U=chol(A);
A_1=inv(U)*inv(U');
A_1dx=-A_1*Adx*A_1;
A_1dy=-A_1*Ady*A_1;

ppx=[1 x y 0 0 0 ];
ppdx=[0 1 0 0 0 0];
ppy=[0 0 0 1 x y];
ppdy=[0 0 0 0 0 1];

ppdy2=[0 0 1 0 0 0];
ppdx2=[0 0 0 0 1 0]; 

Bfem(1,:)=ppdx*A_1*B+ppx*A_1dx*B+ppx*A_1*Bdx;
Bfem(2,:)=ppdy*A_1*B+ppy*A_1dy*B+ppy*A_1*Bdy;
Bfem(3,:)=ppdy2*A_1*B+ppx*A_1dy*B+ppx*A_1*Bdy+ppdx2*A_1*B+ppy*A_1dx*B+ppy*A_1*Bdx;

N(1,:)=ppx*A_1*B;
N(2,:)=ppy*A_1*B;

zn=crack_side(x,y,Xc,Yc);
n_=0;
Bfem_e = [];
for i = i_edge
        n_=n_+1;
        Bfem_e(1,n_)=Bfem(1,2*i-1)*zn;
        Bfem_e(2,n_)=Bfem(2,2*i-1)*zn;
        Bfem_e(3,n_)=Bfem(3,2*i-1)*zn;
        n_=n_+1;
        Bfem_e(1,n_)=Bfem(1,2*i)*zn;
        Bfem_e(2,n_)=Bfem(2,2*i)*zn;
        Bfem_e(3,n_)=Bfem(3,2*i)*zn;
end

par = polar1(Xc(1),Xc(1)-Xc(2),Yc(1),Yc(1)-Yc(2),x,y);
fi=par(5);
par_ = polar2(Xc(1),Xc(1)-Xc(2),Yc(1),Yc(1)-Yc(2),x,y);
n_=0;
for i = i_tip1
        n_=n_+1;
        Bfem_t1(1,n_)=Bfem(1,2*i-1)*par(1)*cos(fi) + N(1,2*i-1) * par_(1)*cos(fi);
        Bfem_t1(2,n_)=Bfem(2,2*i)*par(1)*sin(fi) + N(2,2*i) * par_(2)*sin(fi);
        Bfem_t1(3,n_)=Bfem(3,2*i-1)*par(1)*cos(fi)+Bfem(3,2*i)*par(1)*sin(fi) + N(1,2*i-1) * par_(2)*cos(fi) + N(2,2*i) * par_(1)*sin(fi);
        n_=n_+1;
        Bfem_t1(1,n_)=Bfem(1,2*i-1)*par(1)*(-sin(fi)) + N(1,2*i-1) * par_(1)*(-sin(fi));
        Bfem_t1(2,n_)=Bfem(2,2*i)*par(1)*cos(fi) + N(2,2*i) * par_(2)*cos(fi);
        Bfem_t1(3,n_)=Bfem(3,2*i-1)*par(1)*(-sin(fi))+Bfem(3,2*i)*par(1)*cos(fi) + N(1,2*i-1) * par_(2)*(-sin(fi)) + N(2,2*i) * par_(1)*cos(fi);
        
        n_=n_+1;
        Bfem_t1(1,n_)=Bfem(1,2*i-1)*par(2)*cos(fi) + N(1,2*i-1) * par_(3)*cos(fi);
        Bfem_t1(2,n_)=Bfem(2,2*i)*par(2)*sin(fi) + N(2,2*i) * par_(4)*sin(fi);
        Bfem_t1(3,n_)=Bfem(3,2*i-1)*par(2)*cos(fi)+Bfem(3,2*i)*par(2)*sin(fi) + N(1,2*i-1) * par_(4)*cos(fi) + N(2,2*i) * par_(3)*sin(fi);
        n_=n_+1;
        Bfem_t1(1,n_)=Bfem(1,2*i-1)*par(2)*(-sin(fi)) + N(1,2*i-1) * par_(3)*(-sin(fi));
        Bfem_t1(2,n_)=Bfem(2,2*i)*par(2)*cos(fi) + N(2,2*i) * par_(4)*cos(fi);
        Bfem_t1(3,n_)=Bfem(3,2*i-1)*par(2)*(-sin(fi))+Bfem(3,2*i)*par(2)*cos(fi) + N(1,2*i-1) * par_(4)*(-sin(fi)) + N(2,2*i) * par_(3)*cos(fi);

        n_=n_+1;
        Bfem_t1(1,n_)=Bfem(1,2*i-1)*par(3)*cos(fi) + N(1,2*i-1) * par_(5)*cos(fi);
        Bfem_t1(2,n_)=Bfem(2,2*i)*par(3)*sin(fi) + N(2,2*i) * par_(6)*sin(fi);
        Bfem_t1(3,n_)=Bfem(3,2*i-1)*par(3)*cos(fi)+Bfem(3,2*i)*par(3)*sin(fi) + N(1,2*i-1) * par_(6)*cos(fi) + N(2,2*i) * par_(5)*sin(fi);
        n_=n_+1;
        Bfem_t1(1,n_)=Bfem(1,2*i-1)*par(3)*(-sin(fi)) + N(1,2*i-1) * par_(5)*(-sin(fi));
        Bfem_t1(2,n_)=Bfem(2,2*i)*par(3)*cos(fi) + N(2,2*i) * par_(6)*cos(fi);
        Bfem_t1(3,n_)=Bfem(3,2*i-1)*par(3)*(-sin(fi))+Bfem(3,2*i)*par(3)*cos(fi) + N(1,2*i-1) * par_(6)*(-sin(fi)) + N(2,2*i) * par_(5)*cos(fi);

        n_=n_+1;
        Bfem_t1(1,n_)=Bfem(1,2*i-1)*par(4)*cos(fi) + N(1,2*i-1) * par_(7)*cos(fi);
        Bfem_t1(2,n_)=Bfem(2,2*i)*par(4)*sin(fi) + N(2,2*i) * par_(8)*sin(fi);
        Bfem_t1(3,n_)=Bfem(3,2*i-1)*par(4)*cos(fi)+Bfem(3,2*i)*par(4)*sin(fi) + N(1,2*i-1) * par_(8)*cos(fi) + N(2,2*i) * par_(7)*sin(fi);
        n_=n_+1;
        Bfem_t1(1,n_)=Bfem(1,2*i-1)*par(4)*(-sin(fi)) + N(1,2*i-1) * par_(7)*(-sin(fi));
        Bfem_t1(2,n_)=Bfem(2,2*i)*par(4)*cos(fi) + N(2,2*i) * par_(8)*cos(fi);
        Bfem_t1(3,n_)=Bfem(3,2*i-1)*par(4)*(-sin(fi))+Bfem(3,2*i)*par(4)*cos(fi) + N(1,2*i-1) * par_(8)*(-sin(fi)) + N(2,2*i) * par_(7)*cos(fi);
        
end

par = polar1(Xc(length(Xc)),Xc(length(Xc))-Xc(length(Xc)-1),Yc(length(Yc)),Yc(length(Yc))-Yc(length(Xc)-1),x,y);
fi=par(5);
par_ =polar2(Xc(length(Xc)),Xc(length(Xc))-Xc(length(Xc)-1),Yc(length(Yc)),Yc(length(Yc))-Yc(length(Xc)-1),x,y);
n_=0;
for i = i_tip2
        n_=n_+1;
        Bfem_t2(1,n_)=Bfem(1,2*i-1)*par(1)*cos(fi) + N(1,2*i-1) * par_(1)*cos(fi);
        Bfem_t2(2,n_)=Bfem(2,2*i)*par(1)*sin(fi) + N(2,2*i) * par_(2)*sin(fi);
        Bfem_t2(3,n_)=Bfem(3,2*i-1)*par(1)*cos(fi)+Bfem(3,2*i)*par(1)*sin(fi) + N(1,2*i-1) * par_(2)*cos(fi) + N(2,2*i) * par_(1)*sin(fi);
        n_=n_+1;
        Bfem_t2(1,n_)=Bfem(1,2*i-1)*par(1)*(-sin(fi)) + N(1,2*i-1) * par_(1)*(-sin(fi));
        Bfem_t2(2,n_)=Bfem(2,2*i)*par(1)*cos(fi) + N(2,2*i) * par_(2)*cos(fi);
        Bfem_t2(3,n_)=Bfem(3,2*i-1)*par(1)*(-sin(fi))+Bfem(3,2*i)*par(1)*cos(fi) + N(1,2*i-1) * par_(2)*(-sin(fi)) + N(2,2*i) * par_(1)*cos(fi);
        
        n_=n_+1;
        Bfem_t2(1,n_)=Bfem(1,2*i-1)*par(2)*cos(fi) + N(1,2*i-1) * par_(3)*cos(fi);
        Bfem_t2(2,n_)=Bfem(2,2*i)*par(2)*sin(fi) + N(2,2*i) * par_(4)*sin(fi);
        Bfem_t2(3,n_)=Bfem(3,2*i-1)*par(2)*cos(fi)+Bfem(3,2*i)*par(2)*sin(fi) + N(1,2*i-1) * par_(4)*cos(fi) + N(2,2*i) * par_(3)*sin(fi);
        n_=n_+1;
        Bfem_t2(1,n_)=Bfem(1,2*i-1)*par(2)*(-sin(fi)) + N(1,2*i-1) * par_(3)*(-sin(fi));
        Bfem_t2(2,n_)=Bfem(2,2*i)*par(2)*cos(fi) + N(2,2*i) * par_(4)*cos(fi);
        Bfem_t2(3,n_)=Bfem(3,2*i-1)*par(2)*(-sin(fi))+Bfem(3,2*i)*par(2)*cos(fi) + N(1,2*i-1) * par_(4)*(-sin(fi)) + N(2,2*i) * par_(3)*cos(fi);

        n_=n_+1;
        Bfem_t2(1,n_)=Bfem(1,2*i-1)*par(3)*cos(fi) + N(1,2*i-1) * par_(5)*cos(fi);
        Bfem_t2(2,n_)=Bfem(2,2*i)*par(3)*sin(fi) + N(2,2*i) * par_(6)*sin(fi);
        Bfem_t2(3,n_)=Bfem(3,2*i-1)*par(3)*cos(fi)+Bfem(3,2*i)*par(3)*sin(fi) + N(1,2*i-1) * par_(6)*cos(fi) + N(2,2*i) * par_(5)*sin(fi);
        n_=n_+1;
        Bfem_t2(1,n_)=Bfem(1,2*i-1)*par(3)*(-sin(fi)) + N(1,2*i-1) * par_(5)*(-sin(fi));
        Bfem_t2(2,n_)=Bfem(2,2*i)*par(3)*cos(fi) + N(2,2*i) * par_(6)*cos(fi);
        Bfem_t2(3,n_)=Bfem(3,2*i-1)*par(3)*(-sin(fi))+Bfem(3,2*i)*par(3)*cos(fi) + N(1,2*i-1) * par_(6)*(-sin(fi)) + N(2,2*i) * par_(5)*cos(fi);

        n_=n_+1;
        Bfem_t2(1,n_)=Bfem(1,2*i-1)*par(4)*cos(fi) + N(1,2*i-1) * par_(7)*cos(fi);
        Bfem_t2(2,n_)=Bfem(2,2*i)*par(4)*sin(fi) + N(2,2*i) * par_(8)*sin(fi);
        Bfem_t2(3,n_)=Bfem(3,2*i-1)*par(4)*cos(fi)+Bfem(3,2*i)*par(4)*sin(fi) + N(1,2*i-1) * par_(8)*cos(fi) + N(2,2*i) * par_(7)*sin(fi);
        n_=n_+1;
        Bfem_t2(1,n_)=Bfem(1,2*i-1)*par(4)*(-sin(fi)) + N(1,2*i-1) * par_(7)*(-sin(fi));
        Bfem_t2(2,n_)=Bfem(2,2*i)*par(4)*cos(fi) + N(2,2*i) * par_(8)*cos(fi);
        Bfem_t2(3,n_)=Bfem(3,2*i-1)*par(4)*(-sin(fi))+Bfem(3,2*i)*par(4)*cos(fi) + N(1,2*i-1) * par_(8)*(-sin(fi)) + N(2,2*i) * par_(7)*cos(fi);
end

Bfem_new=[Bfem,Bfem_e,Bfem_t1,Bfem_t2];

w=Bfem_new'*D*Bfem_new;

end