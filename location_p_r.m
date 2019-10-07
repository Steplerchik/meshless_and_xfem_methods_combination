function w = location_p_r(a,b,klu,kol)

n_x=zeros(klu+klu+1,kol+1);
n_y=zeros(klu+klu+1,kol+1);

for fi=1:(klu+1)
    sred=(b-a*cos(pi/4/klu*(fi-1)))/kol;
for i=1:(kol+1)
    n_x(fi,i)=a*cos(pi/4/klu*(fi-1))+sred*(i-1);
    n_y(fi,i)=a*sin(pi/4/klu*(fi-1))+sred*(i-1)*tan(pi/4/klu*(fi-1));
end
end

for fi=(klu+2):(klu+klu+1)
    sred=(b-a*sin(pi/4/klu*(fi-1)))/kol;
for i=1:(kol+1)
    n_x(fi,i)=a*cos(pi/4/klu*(fi-1))+sred*(i-1)/tan(pi/4/klu*(fi-1));
    n_y(fi,i)=a*sin(pi/4/klu*(fi-1))+sred*(i-1);
end
end

n=(klu+klu+1)*(kol+1);
w=zeros(2,n);
n_=0;
for fi=1:(klu+klu+1)
for i=1:(kol+1)
    n_=n_+1;
    w(1,n_)=n_x(fi,i);
    w(2,n_)=n_y(fi,i);
end
end

end