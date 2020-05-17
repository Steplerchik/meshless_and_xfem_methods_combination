clc
clear all
close all
% x=1;
% y=0;
% a=[1 0];
% b=[x y];
% 
% if (y==0)&&(x<0)
%     pi/pi*180
% else
% acos(sum(a.*b)/norm(a)/norm(b))/pi*180*sign(y)
% end
xd=-1;
yd=0;

if xd >= 0
    fi = atan(yd/xd)/pi*180
else
    if yd >= 0
    fi = (pi + atan(yd/xd))/pi*180
    else
    fi = (-pi + atan(yd/xd))/pi*180 
    end
end