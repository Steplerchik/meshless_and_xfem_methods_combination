function w = Stat_BC(n,P,X,Y,kol,b,sigma,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)

F_new=sigma;

F_real=Ix_N(n,P,X,Y,kol,b,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)*F_new;

for uzli=1:n

if (X(uzli)~=0)&&(Y(uzli)>b)
    hola = 5
    F_real(2*uzli,1)=0;
end
end

w=F_real;
end