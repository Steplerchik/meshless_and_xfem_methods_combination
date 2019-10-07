function w = growth(xt,xd,yt,yd,tet_c,l)

tet_c = tet_c/180*pi;

if xd >= 0
    fi = atan(yd/xd);
else
    if yd >= 0
    fi = pi + atan(yd/xd);
    else
    fi = -pi + atan(yd/xd);    
    end
end

alf = fi + tet_c;

x = xt + l*cos(alf);
y = yt + l*sin(alf);

w = [x, y];
end
    
