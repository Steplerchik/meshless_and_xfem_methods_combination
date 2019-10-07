function w = crack_dist(x,y,Xc,Yc,bound,step)
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

if (norm(v)==0)&&(h~=1)&&(h~=length(Xc))
    qui=1;
% elseif (norm([x-Xc(1) y-Yc(1)])<0.6*step)||(norm([x-Xc(length(Xc)) y-Yc(length(Yc))])<0.6*step)
% qui = 0;
else
    
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
    
    if (offset <= bound)&&(offset_ <= bound)
        qui = 1; % сверили к какой ломаной точка ближе
    elseif (offset <= bound)&&((v*v1')/(norm(v)*norm(v1)) > 0)
        qui = 1;
    elseif (offset_ <= bound)&&((v*v1_')/(norm(v)*norm(v1_)) > 0)
        qui = 1;
    else
        qui = 0;
    end
    
    elseif (offset <= bound)&&((v*v1')/(norm(v)*norm(v1)) > 0)&&(r*(v*v1')/(norm(v)*norm(v1)) > step+1^(-10))
        qui = 1;
    else
        qui = 0;
    end  
   
else    
    v1 = [Xc(h-1)-Xc(h) Yc(h-1)-Yc(h)];
    k = v1(2)/v1(1);
    v2=[-k 1];
    offset=r*abs((v*v2')/(norm(v)*norm(v2)));
    
    if (offset <= bound)&&((v*v1')/(norm(v)*norm(v1)) > 0)&&(r*(v*v1')/(norm(v)*norm(v1)) > step+1^(-10))
        qui = 1;
    else
        qui = 0;
    end  
end

end
w=qui;

end