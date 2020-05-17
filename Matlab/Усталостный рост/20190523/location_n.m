function w = location_n(a,b,kol)
sc=2;
n_x=zeros(kol+1,kol+1);
n_y=zeros(kol+1,kol+1);
sum=0;
h=b/kol;
for c=1:(kol-1)
        sum=sum+c;
end

step=kol*h/(kol+(sc-1)*sum/(kol-1));
delta=(sc-1)*step/(kol-1);
ccc=0;
for i=1:(kol+1)
     if (i>2)
        ccc=ccc+(i-2);
     end
      cc=0;
    for j=1:(kol+1)
        n_y(i,j)=b-(step*(i-1)+ccc*delta);
        if (j>2)
        cc=cc+(j-2);
        end

    if (n_y(i,j)<=a)
        hh=(b-sqrt(a^2-(n_y(i,j))^2))/kol;
        step2=kol*hh/(kol+(sc-1)*sum/(kol-1));
        delta2=(sc-1)*step2/(kol-1);
        n_x(i,j)=b-(step2*(j-1)+cc*delta2);
    else
    n_x(i,j)=b-(step*(j-1)+cc*delta);
    end
    end
end
n=(kol+1)*(kol+1);
w=zeros(2,n);
n_=0;
for i=1:(kol+1)
    for j=1:(kol+1)
    n_=n_+1;
    w(1,n_)=n_x(i,j);
    w(2,n_)=n_y(i,j);
end
end
    
 

end