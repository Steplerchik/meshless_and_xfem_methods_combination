function w = SIF_uv(fi,u_,tet,r,G,k)

U=u_(1,:);
V=u_(2,:);

v=V*cos(fi)-U*sin(fi);
u=U*cos(fi)+V*sin(fi);
n_=0;
for i=1:length(u)
    n_=n_+1;
    Uv(n_)=u(i);
    n_=n_+1;
    Uv(n_)=v(i);
end

Uv=Uv';

n_=0;
for i=1:length(u)
    n_=n_+1;
    H(n_,:)=[sqrt(r(i))*(k-cos(tet(i)))*cos(tet(i)/2) sqrt(r(i))*sin(tet(i)/2)*(k+2+cos(tet(i))) r(i) 0];
    n_=n_+1;
    H(n_,:)=[sqrt(r(i))*(k-cos(tet(i)))*sin(tet(i)/2) sqrt(r(i))*(-cos(tet(i)/2))*(k-2+cos(tet(i))) 0 r(i)];
end

SIF=2*G*sqrt(2*pi)*(inv(H'*H)*H')*Uv;
w = SIF;
end
    
