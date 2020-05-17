function w = polar2(xt,xd,yt,yd,x,y)

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

    du1dx = 1/(2*sqrt(r))*cos(tet/2);
    du1dy = 1/(2*sqrt(r))*sin(tet/2);
    
    du2dx = -1/(2*sqrt(r))*sin(tet/2);
    du2dy = 1/(2*sqrt(r))*cos(tet/2);
    
    du3dx = -1/(2*sqrt(r))*sin(tet)*sin(3*tet/2);
    du3dy = 1/(2*sqrt(r))*sin(tet)*cos(tet/2)+1/(sqrt(r))*sin(tet/2)*cos(tet)*cos(tet);
    
    du4dx = -1/(2*sqrt(r))*sin(tet)*cos(3*tet/2);
    du4dy = 1/(2*sqrt(r))*sin(tet)*sin(tet/2)+1/(sqrt(r))*cos(tet/2)*cos(tet)*cos(tet);
    
    
    dxdX=cos(fi);
    dydX=-sin(fi);
    dxdY=sin(fi);
    dydY=cos(fi);

    du1dX=du1dx*dxdX+du1dy*dydX;
    du1dY=du1dx*dxdY+du1dy*dydY;
    
    du2dX=du2dx*dxdX+du2dy*dydX;
    du2dY=du2dx*dxdY+du2dy*dydY; 
    
    du3dX=du3dx*dxdX+du3dy*dydX;
    du3dY=du3dx*dxdY+du3dy*dydY; 
    
    du4dX=du4dx*dxdX+du4dy*dydX;
    du4dY=du4dx*dxdY+du4dy*dydY;    
    
% plot(x+u*10^10,y+v*10^10,'.')
% w=[tet/pi*180 r u v];
w = [du1dX du1dY du2dX du2dY du3dX du3dY du4dX du4dY];
end
    
