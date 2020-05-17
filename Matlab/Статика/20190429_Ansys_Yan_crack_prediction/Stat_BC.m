function w = Stat_BC(n,P,X,Y,kol,b,sigma,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)

F_new=sigma;
F_neww(1,1) = sigma(2,1);
F_neww(2,1) = sigma(1,1);
F_neww = 5*F_neww;

F_real=Ix_N(n,P,X,Y,kol,b,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)*F_new + Iy_N(n,P,X,Y,kol,b,a,mas,i_edge,i_tip1,i_tip2,Xc,Yc)*F_neww;

w=F_real;
end