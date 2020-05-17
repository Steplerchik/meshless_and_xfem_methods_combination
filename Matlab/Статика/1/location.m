function w = location(a,b,kol)

kol1=kol;
kol2=kol*b/a;
n_x=zeros(kol1+1,kol1+1);
n_y=zeros(kol2+1,kol2+1);
h=a/kol;

for i=1:(kol1+1)
    for j=1:(kol2+1)
  
    n_x(i,j)=h*(i-1);

    n_y(i,j)=h*(j-1);

    end
end

n_=0;
for i=1:(kol1+1)
    for j=1:(kol2+1)
    n_=n_+1;
    w(1,n_)=n_x(i,j);
    w(2,n_)=n_y(i,j);
    end
end

end