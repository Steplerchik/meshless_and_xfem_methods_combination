clc
clear all
close all
Start = cputime;
                            % Размеры пластины

b=100;        %мм
a=100;        %мм

                            % Cвойства материала
                            
mu=0.25;
E=2.0*10^5;  %МПа
G=E/(2*(1+mu));
C_Paris = 7e-8;
m_Paris = 2.1;

                            % Тип плоского состояния
                            
% Рассматривается случай ПНС
D=E/(1-mu^2)*[1 mu 0; mu 1 0; 0 0 (1-mu)/2];
k = (3 - mu)/(1 + mu);

                            % Параметры сетки узлов
setka=1;
kol=11*setka; % количество участков по х
step = a / kol;

                            % Расположение узлов
N_L=location(a,b,kol);
X=N_L(1,:);
Y=N_L(2,:);
n=length(X);

                            % Изображение сетки узлов
figure
hold on
plot(X,Y,'k.')

                            % Расположение трещины
% Начальная трещина
lc = a/3;

yc1 = b/2;
yc2 = b/2;
xc1 = a/2 + lc/2;
xc2 = a/2 - lc/2;

Yc = [yc1 yc2];
Xc = [xc1 xc2];

plot(Xc, Yc,'k')
%%
for steplerchik=1:8
                            % Назначение узлам около трещины специальных
                            % свойств
clear ('i_edge', 'i_tip1', 'i_tip2', 'K1', 'Kmat', 'K', 'F', 'K_new', 'u')

n_=0;
n_1=0;
n_2=0;
chack=0;

for i=1:n
    c_d=crack_dist(X(i),Y(i),Xc,Yc,1.2*step,step);
    ni_1=norm([X(i)-Xc(1) Y(i)-Yc(1)]); % расстояние
    ni_2=norm([X(i)-Xc(length(Xc)) Y(i)-Yc(length(Xc))]); % расстояние

    if (c_d == 1)
    plot(X(i),Y(i),'ro')
        n_=n_+1;
        i_edge(n_) = i;
        chack=0;
    end
    if chack==0
        if  (ni_1 < 2.5*step-1^(-10))
        plot(X(i),Y(i),'ks')
        n_1=n_1+1;
        i_tip1(n_1) = i;
        end
        
        if (ni_2 < 2.5*step-1^(-10))
        plot(X(i),Y(i),'bs')
        n_2=n_2+1;
        i_tip2(n_2) = i;
        end
    end
    chack=0;       
end

                            % Создание матрицы P
                            
m=6; % Столбцов в матрице P
P=P_mat(n,m,X,Y);

                            % Формирование матрицы К
                            
figure
hold on
mas=kol; % Размерность точек интегрирования

%Задание y-координаты центра квадрата интегрирования
for i=1:length(Y)-1
    if (Y(i)<=Yc(1))&&(Y(i+1)>Yc(1))
        step=Y(i+1)-Y(i);
        ycen1=Y(i)+step/2;
    end
     if (Y(i)<=Yc(length(Yc)))&&(Y(i+1)>Yc(length(Yc)))
        step=Y(i+1)-Y(i);
        ycen2=Y(i)+step/2;
     end
end

if ycen1<=ycen2
    ycen_d=ycen1;
    ycen_u=ycen2;
else
    ycen_d=ycen2;
    ycen_u=ycen1;
end

% Интегрирование и формирование глобальной матрицы жесткости
K1(:,:,1)=Iy_edge(step,ycen_d,ycen_u,n,P,X,Y,D,a,kol,b,mas,i_edge,i_tip1,i_tip2,Xc,Yc);
K1(:,:,2)=Iy_reg_up(step,ycen_u,n,P,X,Y,D,a,kol,b,mas,i_edge,i_tip1,i_tip2,Xc,Yc);
K1(:,:,3)=Iy_reg_dwn(step,ycen_d,n,P,X,Y,D,a,kol,b,mas,i_edge,i_tip1,i_tip2,Xc,Yc);
if ycen1~=ycen2
K1(:,:,4)=Iy_reg_add2(step,ycen1,ycen2,n,P,X,Y,D,a,kol,b,mas,i_edge,i_tip1,i_tip2,Xc,Yc);
K1(:,:,5)=Iy_reg_add1(step,ycen1,ycen2,n,P,X,Y,D,a,kol,b,mas,i_edge,i_tip1,i_tip2,Xc,Yc);

siz=2.5*step;

Kmat = @(x,y) K_mat(n,x,y,P,X,Y,D,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc);

crxy=[siz, ycen2];
eta=([Xc(length(Xc)), Yc(length(Yc))]-crxy)./siz;
ras=10;
K1(:,:,6)=Iy_tip1(ras,eta,Kmat,crxy,siz)+Iy_tip2(ras,eta,Kmat,crxy,siz)+Iy_tip3(ras,eta,Kmat,crxy,siz)+Iy_tip4(ras,eta,Kmat,crxy,siz);

crxy=[kol*step-siz, ycen1];
eta=([Xc(1), Yc(1)]-crxy)./siz;
ras=10;
K1(:,:,7)=Iy_tip1(ras,eta,Kmat,crxy,siz)+Iy_tip2(ras,eta,Kmat,crxy,siz)+Iy_tip3(ras,eta,Kmat,crxy,siz)+Iy_tip4(ras,eta,Kmat,crxy,siz);

else
siz=2.5*step;

Kmat = @(x,y) K_mat(n,x,y,P,X,Y,D,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc);

crxy=[siz, ycen2];
eta=([Xc(length(Xc)), Yc(length(Yc))]-crxy)./siz;
ras=10;
K1(:,:,4)=Iy_tip1(ras,eta,Kmat,crxy,siz)+Iy_tip2(ras,eta,Kmat,crxy,siz)+Iy_tip3(ras,eta,Kmat,crxy,siz)+Iy_tip4(ras,eta,Kmat,crxy,siz);

crxy=[kol*step-siz, ycen1];
eta=([Xc(1), Yc(1)]-crxy)./siz;
ras=10;
K1(:,:,5)=Iy_tip1(ras,eta,Kmat,crxy,siz)+Iy_tip2(ras,eta,Kmat,crxy,siz)+Iy_tip3(ras,eta,Kmat,crxy,siz)+Iy_tip4(ras,eta,Kmat,crxy,siz);
end    

K=sum(K1,3);

                            % Нагрузка амплитудными значениями циклов
                            % sigma_x и sigma_y
                            
sigma(:,1)=[0; 1]; %min
sigma(:,2)=[0; 2]; %max

figure
plot(X,Y,'k.')
hold on
plot(Xc,Yc,'k')

for i = 1:2
    
% Учет статических ГУ
    %МПа
F=Stat_BC(n,P,X,Y,kol,b,sigma(:,i),a,mas,i_edge,i_tip1,i_tip2,Xc,Yc);

% Учет кинематических ГУ
K_new=K_Kin_BC(X,Y,n,K,b);

% Вычисление вектора узловых перемещений
u=K_new\F;

kin_tip1(i,:) = SIF(1,n,P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc,G,k,u);
kin_tip2(i,:) = SIF(2,n,P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc,G,k,u);

end

crack1 = predict(kin_tip1,C_Paris,m_Paris);
crack2 = predict(kin_tip2,C_Paris,m_Paris); 

grow1 = growth(Xc(1),Xc(1)-Xc(2),Yc(1),Yc(1)-Yc(2),crack1(3),crack1(4));
grow2 = growth(Xc(length(Xc)),Xc(length(Xc))-Xc(length(Xc)-1),Yc(length(Xc)),Yc(length(Xc))-Yc(length(Xc)-1),crack2(3),crack2(4));

figure
hold on
for uzli=1:n
Nuz=N(n,X(uzli),Y(uzli),P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc);
uuz(:,uzli)=Nuz*u;
X_def(uzli)=X(uzli)+uuz(1,uzli)*10^4*3;
Y_def(uzli)=Y(uzli)+uuz(2,uzli)*10^4*3;
plot(X_def(uzli),Y_def(uzli),'k.')
end

Xc = [grow1(1), Xc, grow2(1)];
Yc = [grow1(2), Yc, grow2(2)];

figure
plot(X,Y,'k.')
hold on
plot(Xc,Yc,'k')
plot(Xc,Yc,'ko')

history1(steplerchik,:) = [crack1(3), crack1(4)];
history2(steplerchik,:) = [crack2(3), crack2(4)];
sifs1(steplerchik,:) = kin_tip1(2,:);
sifs2(steplerchik,:) = kin_tip2(2,:);
end

% al=1/2;
% F1=sqrt(sec(al*pi/2));
% F1_z=(1-0.025*al^2+0.06*al^4)*F1;
% Hura=sqrt(pi*a/4)*F1_z;

Elapsed = (cputime - Start) / 60

%%
                 % Изображение поля перемещений узлов
% figure
% hold on
% for uzli=1:n
% % X_def(uzli)=X(uzli)+u(2*uzli-1)*10^4*3;
% % Y_def(uzli)=Y(uzli)+u(2*uzli)*10^4*3;
% X_def(uzli)=X(uzli)+u(2*uzli-1)*10^4*3;
% Y_def(uzli)=Y(uzli)+u(2*uzli)*10^4*3;
% 
% plot(X_def(uzli),Y_def(uzli),'k.')
% end
%%
Xc2 = Xc;
Yc2 = Yc;
Yc = [Yc(2), Yc(3), Yc(4), Yc(5)];
Xc = [Xc(2), Xc(3), Xc(4), Xc(5)];
            % Вычисление интересующего поля
kolp=kol;

% Расположение точек
P_L=location_tr(a,b,kolp);
Xp=P_L(1,:);
Yp=P_L(2,:);
% Addition
ch=kolp*20;
%Xp_add=Xc(length(Xc))-step:(Xc(1)-Xc(length(Xc))+2*step)/ch:Xc(1)+step;
Xp_add=Xc(length(Xc)):(Xc(1)-Xc(length(Xc)))/ch:Xc(1);
%Xp_add=0:a/ch:a;
for i = 1:length(Xc)-1
    for j = 1:length(Xp_add)
        if (Xp_add(j) >= Xc(i+1))&&(Xp_add(j) <= Xc(i))
Yp_add1(j)=(Yc(i)-Yc(i+1))/((Xc(i)-Xc(i+1))).*(Xp_add(j)-Xc(i+1))+Yc(i+1)-0.001;
Yp_add2(j)=(Yc(i)-Yc(i+1))/((Xc(i)-Xc(i+1))).*(Xp_add(j)-Xc(i+1))+Yc(i+1)+0.001;
        end
    end
end
Xp=[Xp,Xp_add,Xp_add];
Yp=[Yp,Yp_add1,Yp_add2];
np=length(Xp);

Xp_def=zeros(1,np);
Yp_def=zeros(1,np);
SXp=zeros(1,np);
SYp=zeros(1,np);
SXYp=zeros(1,np);

up=zeros(2,np);
sp=zeros(3,np);

% Вычисление поля перемещений и поля напряжений
for point=1:np
Np=N(n,Xp(point),Yp(point),P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc);
up(:,point)=Np*u;

% Bp=B(n,Xp(point),Yp(point),P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc);
% sp(:,point)=D*Bp*u;
end
% 
% % 
% % Перезадание поля перемещений и поля напряжений
% for point=1:np
% % Xp_def(1,point)=Xp(point)+up(1,point)*10^10;
% % Yp_def(1,point)=Yp(point)+up(2,point)*10^10;
% 
% SXp(1,point)=sp(1,point);
% SYp(1,point)=sp(2,point);
% SXYp(1,point)=sp(3,point);
% 
% end
% % 
% % % Изображение аналитического и мешлесс вычисленных напряжений
% % figure
% % hold on
% % title('Sy (y=0)');
% % n_=0;
% % for point=1:np
% %     if (Yp(point)==0)
% %       plot(Xp(point),SYp(point),'b.')
% %       n_=n_+1;
% %       sssy(n_)=SYp(point);
% %     end
% % end
% % 
% % 
% % 
% % % Изображение поля перемещений
% % figure('Position',[400 80 600 600]);
% % set(axes, 'Position', [0.1 0.1 0.8 0.8]);
% % hold on
% % for point=1:np
% % Xp_def(1,point)=Xp(point)+up(1,point)*10^4*3;
% % Yp_def(1,point)=Yp(point)+up(2,point)*10^4*3;
% % xxx=[Xp(point),Xp_def(1,point)];
% % yyy=[Yp(point),Yp_def(1,point)];
% % 
% % plot(Xp_def(1,point),Yp_def(1,point),'c.')
% % plot(xxx,yyy,'b')
% % plot(Xp(point),Yp(point),'k.')
% % end
% 
figure('Position',[400 80 600 600]);
set(axes, 'Position', [0.1 0.1 0.8 0.8]);
hold on

for point=1:np
% Xp_def(1,point)=Xp(point)+up(1,point)*10^4*3;
% Yp_def(1,point)=Yp(point)+up(2,point)*10^4*3;
Xp_def(1,point)=Xp(point)+up(1,point)*10^4*3;
Yp_def(1,point)=Yp(point)+up(2,point)*10^4*3;
plot(Xp_def(1,point),Yp_def(1,point),'k.')
end

n_=0;
for point=np-ch:np
    n_=n_+1;
    p_crack_u(n_)=point;
end
n_=0;
for point=np-ch-1-ch:np-ch-1
    n_=n_+1;
    p_crack_d(n_)=point;
end
% for point=1:np-1
%     if (Yp(point)<yc)&&(Yp(point+1)>yc)
%         n_=n_+1;
%         p_crack_d(n_)=point;
%         p_crack_u(n_)=point+1;
%     end
% end
    
figure('Position',[400 80 600 600]);
set(axes, 'Position', [0.1 0.1 0.8 0.8]);
hold on

for point=1:np
plot(Xp(point),Yp(point),'k.')
end
for point=p_crack_d
plot(Xp(point),Yp(point),'ko')
end
for point=p_crack_u
plot(Xp(point),Yp(point),'ks')
end
%%
figure('Position',[400 80 600 600]);
set(axes, 'Position', [0.1 0.1 0.8 0.8]);
hold on
mnogh = 10^4;
for point=1:np
plot(Xp(point)+up(1,point)*mnogh,Yp(point)+up(2,point)*mnogh,'k.')
end
for point=p_crack_d
plot(Xp(point)+up(1,point)*mnogh,Yp(point)+up(2,point)*mnogh,'ko')
end
for point=p_crack_u
plot(Xp(point)+up(1,point)*mnogh,Yp(point)+up(2,point)*mnogh,'ks')
end
%%
figure('Position',[400 80 600 600]);
set(axes, 'Position', [0.1 0.1 0.8 0.8]);
hold on

% for point=1:np
% plot(Xp(point),Yp(point),'k.')
% end
for point=p_crack_d
    if crack_side(Xp(point),Yp(point),Xc,Yc) == -1
        plot(Xp(point),Yp(point),'bo')
    else
        plot(Xp(point),Yp(point),'ko')
    end
end
for point=p_crack_u
    if crack_side(Xp(point),Yp(point),Xc,Yc) == 1
        plot(Xp(point),Yp(point),'bs')
    else
        plot(Xp(point),Yp(point),'ks')
    end
end


% figure
% hold on
% title('Sy (p_crack_d)');
% n_=0;
% for point=p_crack_d
% plot(Xp(point),SYp(point),'bo')
% end
% 
% figure
% hold on
% title('Sy (p_crack_u)');
% n_=0;
% for point=p_crack_d
% plot(Xp(point),SYp(point),'bs')
% end
% 
% figure('Position',[400 80 600 600]);
% set(axes, 'Position', [0.1 0.1 0.8 0.8]);
% hold on
% 
% for point=1:np
% plot(Xp(point),Yp(point),'k.')
% end
% for point=p_crack_d
% plot(Xp(point),Yp(point),'ko')
% end
% for point=p_crack_u
% plot(Xp(point),Yp(point),'ks')
% end

%%
% if hui(2) > 0
%     tetc_1 = 2*atan(1/4*hui(1)/hui(2)-1/4*sqrt((hui(1)/hui(2))^2+8));
% else
%     tetc_1 = 2*atan(1/4*hui(1)/hui(2)+1/4*sqrt((hui(1)/hui(2))^2+8));
% end
% xd = Xc(1)-Xc(length(Xc));
% yd = Yc(1)-Yc(length(Yc));
% 
% if xd >= 0
%     fi = atan(yd/xd);
% else
%     if yd >= 0
%     fi = pi + atan(yd/xd);
%     else
%     fi = -pi + atan(yd/xd);    
%     end
% end
% 
% alf = fi + tetc_1;