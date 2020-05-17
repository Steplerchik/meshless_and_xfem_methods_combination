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

par = polar1(Xc(1),Xc(1)-Xc(2),Yc(1),Yc(1)-Yc(2),x,y);
fi=par(5);
n_=0;
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

par = polar1(Xc(length(Xc)),Xc(length(Xc))-Xc(length(Xc)-1),Yc(length(Yc)),Yc(length(Yc))-Yc(length(Xc)-1),x,y);
fi=par(5);
n_=0;
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