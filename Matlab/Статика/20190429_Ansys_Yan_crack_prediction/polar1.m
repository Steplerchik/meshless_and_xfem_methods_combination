function w = polar1(xt,xd,yt,yd,x,y)

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

    u1 = sqrt(r)*cos(tet/2);
    u2 = sqrt(r)*sin(tet/2);
    u3 = sqrt(r)*sin(tet/2)*sin(tet);
    u4 = sqrt(r)*cos(tet/2)*sin(tet);
        
% plot(x+u*10^10,y+v*10^10,'.')
% w=[tet/pi*180 r u v];
w = [u1 u2 u3 u4 fi];
end
    
