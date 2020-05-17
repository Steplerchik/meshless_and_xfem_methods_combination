function w = P_mat(n,m,X,Y)

w=zeros(2*n,m);

k=0;
check=0;
for i=1:2*n   
  
    if (check==0)
        key=1;
        k=k+1;
      p=[1; X(k); Y(k)];
      check=1;
    else
         key=4;
        check=0;
    end
    
key1=0;
for j=1:m
    if (j>=key)&&(j<=key+2)
        key1=key1+1;
    w(i,j)=p(key1);
    end
end

end

end