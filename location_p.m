function w = location_p(a,b,klu,kol)
sc=4;
n_x=zeros(klu+klu+1,kol+1);
n_y=zeros(klu+klu+1,kol+1);
sum=0;
for c=1:(kol-1)
        sum=sum+c;
end

for fi=1:(klu+1)
    sred=(b-a*cos(pi/4/klu*(fi-1)))/kol;
    step=kol*sred/(kol+(sc-1)*sum/(kol-1));
    delta=(sc-1)*step/(kol-1);
    cc=0;
for i=1:(kol+1)
    if (i>2)
        cc=cc+(i-2);
    end
    n_x(fi,i)=a*cos(pi/4/klu*(fi-1))+step*(i-1)+cc*delta;
    n_y(fi,i)=a*sin(pi/4/klu*(fi-1))+(step*(i-1)+cc*delta)*tan(pi/4/klu*(fi-1));
end
end

for fi=(klu+2):(klu+klu+1)
    sred=(b-a*sin(pi/4/klu*(fi-1)))/kol;
    step=kol*sred/(kol+(sc-1)*sum/(kol-1));
    delta=(sc-1)*step/(kol-1);
    cc=0;
for i=1:(kol+1)
    if (i>2)
        cc=cc+(i-2);
    end
    n_x(fi,i)=a*cos(pi/4/klu*(fi-1))+(step*(i-1)+cc*delta)/tan(pi/4/klu*(fi-1));
    n_y(fi,i)=a*sin(pi/4/klu*(fi-1))+step*(i-1)+cc*delta;
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