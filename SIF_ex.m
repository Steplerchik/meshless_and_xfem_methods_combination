function w = SIF_ex(tip, n,P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc,G,k,u,r1,r2)

if tip == 1
    xt = Xc(1);
    xd = Xc(1)-Xc(2);
    yt = Yc(1);
    yd = Yc(1)-Yc(2);
    
    Np=N(n,xt+0.00000000001,yt,P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc);
    u_t=Np*u;
    
elseif tip == 2
    xt = Xc(length(Xc));
    xd = Xc(length(Xc))-Xc(length(Xc)-1);
    yt = Yc(length(Yc));
    yd = Yc(length(Yc))-Yc(length(Xc)-1);
    
    Np=N(n,xt-0.00000000001,yt,P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc);
    u_t=Np*u;
    
end

% l = sqrt(xd^2+yd^2);
% tet = -pi+pi/360:(2*pi-2*pi/360)/35:pi-pi/360;
% 
% for i = 1:length(tet)
% %     r(i) = 0.03*l + rand(1)*(0.2*l - 0.03*l);
%      r(i) = 0.03*l + (0.2*l - 0.03*l);
% end

l = sqrt(xd^2+yd^2);
tet_ = -pi+pi/360:(2*pi-2*pi/360)/35:pi-pi/360;
r_ = r1*l*ones(1,length(tet_));
tet = [];
r = [];
count = 5;
for j = 1:count+1
    tet = [tet,tet_];
    r = [r,r_];
    r_ = r_ + (r2*l - r1*l)/count*ones(1,length(tet_));
end

if xd >= 0
    fi = atan(yd/xd);
else
    if yd >= 0
    fi = pi + atan(yd/xd);
    else
    fi = -pi + atan(yd/xd);    
    end
end

alf=tet+fi;

X_=xt+r.*cos(alf);
Y_=yt+r.*sin(alf);
plot(X_,Y_)
for point=1:length(X_)
Np=N(n,X_(point),Y_(point),P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc);
u_(:,point)=Np*u-u_t;
end

sif = SIF_uv(fi,u_,tet,r,G,k);


w = [sif(1) sif(2)];
end
    
