function w = Stat_BC(n,P,X,Y,kol,b,sigma,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)

F_new=sigma;
F_neww(1,1) = sigma(2,1);
F_neww(2,1) = sigma(1,1);
F_neww = 5*F_neww;

% F_real=Ix_N(n,P,X,Y,kol,b,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)*F_new + Ix_N_(n,P,X,Y,kol,b,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)*((-1)*F_new) + Iy_N(n,P,X,Y,kol,b,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)*F_neww + Iy_N_(n,P,X,Y,kol,b,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)*((-1)*F_neww);
F_real=Ix_N(n,P,X,Y,kol,b,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)*F_new;

for uzli=1:n

if (Y(uzli)==0)
    F_real(2*uzli,1)=0;
end

if (X(uzli)==0)&&(Y(uzli)==0)
    F_real(2*uzli-1,1)=0;
end

% if (X(uzli)==0)&&(Y(uzli)>b)
%     F_real(2*uzli-1,1)=0;
% end

% if (Y(uzli)>b)
%     F_real(2*uzli-1,1)=0;
% end

end

w=F_real;
end