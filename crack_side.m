function w = crack_side(x,y,Xc,Yc)
h=1;
v = [x-Xc(h) y-Yc(h)];
r = norm(v);
for i=1:length(Xc)
    vi = [x-Xc(i) y-Yc(i)];
    if norm(vi) < r
        h = i;
        v = [x-Xc(h) y-Yc(h)];
        r = norm(v);
    end
end


if h < length(Xc)
    v1 = [Xc(h+1)-Xc(h) Yc(h+1)-Yc(h)];
    k=v1(2)/v1(1);
    v2=[-k 1];
    offset=r*abs((v*v2')/(norm(v)*norm(v2))); % нашли расстояние от точки до трещины
    
    if h>1
    v1_ = [Xc(h-1)-Xc(h) Yc(h-1)-Yc(h)];
    k_ = v1_(2)/v1_(1);
    v2_=[-k_ 1];
    offset_=r*abs((v*v2_')/(norm(v)*norm(v2_))); % нашли расстояние от точки до трещины (соседняя ломаная)
    
    if (offset <= offset_)
        if (acos((v*v2')/(norm(v)*norm(v2)))/pi*180 <= 90)
            qui = 1;
        else
            qui = -1;
        end
    else
        if (acos((v*v2_')/(norm(v)*norm(v2_)))/pi*180 <= 90)
            qui = 1;
        else
            qui = -1;
        end
    end
    
    else
        if (acos((v*v2')/(norm(v)*norm(v2)))/pi*180 <= 90)
            qui = 1;
        else
            qui = -1;
        end
    end
  
else    
    v1 = [Xc(h-1)-Xc(h) Yc(h-1)-Yc(h)];
    k = v1(2)/v1(1);
    v2=[-k 1];

    if (acos((v*v2')/(norm(v)*norm(v2)))/pi*180 <= 90)
            qui = 1;
    else
            qui = -1;
    end
    
end
w = qui;
end