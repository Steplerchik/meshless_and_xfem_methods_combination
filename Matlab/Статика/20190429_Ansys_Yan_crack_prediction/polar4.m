function w = polar(xt,xd,yt,yd,x,y)
mu = 0.3;
E=2.1*10^5;
G=E/(2*(1+mu));
kap = (3 - mu)/(1 + mu);

r = sqrt((x-xt)^2+(y-yt)^2);
    
if x >= xt
    alf = atan((y-yt)/(x-xt));
else
    if y >= yt
    alf = pi + atan((y-yt)/(x-xt));
    else
    alf = -pi + atan((y-yt)/(x-xt));    
    end
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

tet=alf-fi;
if tet > pi
    tet=tet-2*pi;
end
if tet < -pi
    tet=tet+2*pi;
end

    u = 1/(2*G)*sqrt(r/(2*pi))*cos(tet/2)*(kap - 1 + 2*(sin(tet/2))^2);
    v = 1/(2*G)*sqrt(r/(2*pi))*sin(tet/2)*(kap + 1 - 2*(cos(tet/2))^2);
    
    V=v*cos(fi)+u*sin(fi);
    U=u*cos(fi)-v*sin(fi);
    
    dudx = 1/(2*G*sqrt(2*pi*r))*cos(tet/2)*((1-mu)/(1+mu)-sin(tet/2)*sin(3*tet/2));
    dvdy = 1/(2*G*sqrt(2*pi*r))*cos(tet/2)*((1-mu)/(1+mu)+sin(tet/2)*sin(3*tet/2));

    dudy = 1/(2*G*sqrt(2*pi*r))*(1/2*(kap - 1 + 2*(sin(tet/2))^2)*sin(tet/2)+cos(tet/2)*sin(tet)*cos(tet));
    dvdx = 1/(2*G*sqrt(2*pi*r))*(1/2*(kap + 1 - 2*(cos(tet/2))^2)*(-sin(tet/2))-sin(tet/2)*(sin(tet))^2);
    
    dxdX=cos(fi);
    dydX=-sin(fi);
    dxdY=sin(fi);
    dydY=cos(fi);

    dudX=dudx*dxdX+dudy*dydX;
    dvdX=dvdx*dxdX+dvdy*dydX;
    dudY=dudx*dxdY+dudy*dydY;
    dvdY=dvdx*dxdY+dvdy*dydY;
        
    dUdX=dudX*cos(fi)-dvdX*sin(fi);
    dUdY=dudY*cos(fi)-dvdY*sin(fi);
    dVdX=dvdX*cos(fi)+dudX*sin(fi);
    dVdY=dvdY*cos(fi)+dudY*sin(fi);
    
    
% plot(x+u*10^10,y+v*10^10,'.')
% w=[tet/pi*180 r u v];
w = [U dUdX dUdY V dVdX dVdY];
end
    
