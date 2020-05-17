function w = predict(kin,C,m,step)

if kin(1,1) > kin(2,1)
    K1max = kin(1,1);
    K1min = kin(2,1);
else
    K1max = kin(2,1);
    K1min = kin(1,1);    
end

if K1min < 0
    del_K1 = K1max;
else
    del_K1 = K1max - K1min;
end
% del_K1 = K1max - K1min;

if abs(kin(1,2)) > abs(kin(2,2))
    K2max = kin(1,2);
    K2min = kin(2,2);
else
    K2max = kin(2,2);
    K2min = kin(1,2);
end

del_K2 = K2max - K2min;

if (K2max > 0)&&(K2min > 0)
    tet = -acos((3*K2max^2+K1max*sqrt(K1max^2+8*K2max^2))/(K1max^2+9*K2max^2));
else
    tet = acos((3*K2max^2+K1max*sqrt(K1max^2+8*K2max^2))/(K1max^2+9*K2max^2));
end

del_Keq = 1/2*cos(tet/2)*(del_K1*(1+cos(tet))-3*del_K2*sin(tet));

if del_Keq < 0
    del_l = 0;
else
%     N=100000;
%     del_l = C*(del_Keq^m)*N;
    del_l = 1*step;
end

w = [del_K1, del_K2,tet/pi*180,del_l];

end
    
