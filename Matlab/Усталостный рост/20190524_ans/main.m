clc
clear all
close all
Start = cputime;
                            % Размеры пластины

b=200;        %мм
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
X_ini=N_L(1,:);
Y_ini=N_L(2,:);

% X = [X,a/2];
% Y = [Y,0];
% 
% X = [X,0];
% Y = [Y,b/2];


                            % Изображение сетки узлов
%%
                            % Расположение трещины
% Начальная трещина
lc = a/11*3;
teta = 15;  % По Цельсию
% Kin1 = 0.2896*sqrt(pi*lc/2);
% Kin2 = 0.4660*sqrt(pi*lc/2);
lc = lc/cos(teta/180*pi)

yc1 = b/2 + lc/2*sin(teta/180*pi);
yc2 = b/2 - lc/2*sin(teta/180*pi);
xc1 = a/2 + lc/2*cos(teta/180*pi);
xc2 = a/2 - lc/2*cos(teta/180*pi);
% xc1 = a-a/11*4;
% xc2 = a/11*2;


% Yc = [yc1 yc2 99.7974];
% Xc = [xc1 xc2 31.8682];
Yc = [yc1 yc1-1*step*sin(teta/180*pi) yc2+1*step*sin(teta/180*pi) yc2];
Xc = [xc1 xc1-1*step*cos(teta/180*pi) xc2+1*step*cos(teta/180*pi) xc2];
% Yc = [yc1 yc2 91.5762];
% Xc = [xc1 xc2 29.9577];
% Yc = [yc1 yc2 91.1687];
% Xc = [xc1 xc2 29.9140];
% % Yc = [yc1 yc2 91.2337];
% % Xc = [xc1 xc2 29.7202];
% Yc = [yc1 yc2 91.4854];
% Xc = [xc1 xc2 29.7453];
% Yc = [yc1 yc2];
% Xc = [xc1 xc2];

%%
for steplerchik=1:4
                            % Назначение узлам около трещины специальных
                            % свойств
clear ('i_edge', 'i_tip1', 'i_tip2', 'K1', 'Kmat', 'K', 'F', 'K_new', 'u','Xp_add','Yp_add1','Yp_add2')

figure
hold on
plot(Xc, Yc,'k')
X = X_ini;
Y = Y_ini;
for i = 1:length(Xc)
%     if (i==1)||(i==length(Xc))
        X = [X,Xc(i)];
        Y = [Y,Yc(i)];
%     end
end
X = [X,Xc(1)+1e-10,Xc(length(Xc))+1e-10];
Y = [Y,Yc(1)+1e-10,Yc(length(Yc))+1e-10];
n=length(X);
plot(X,Y,'k.')

n_=0;
n_1=0;
n_2=0;
chack=0;
i_edge = [];
for i=1:n
    c_d=crack_dist(X(i),Y(i),Xc,Yc,1.2*step,step);
    ni_1=norm([X(i)-Xc(1) Y(i)-Yc(1)]); % расстояние
    ni_2=norm([X(i)-Xc(length(Xc)) Y(i)-Yc(length(Xc))]); % расстояние

    if (c_d == 1)
    plot(X(i),Y(i),'ro')
        n_=n_+1;
        i_edge(n_) = i;
        chack=1;
    end
    if chack==0
%         if  (ni_1 < 2.5*step-1^(-10))
        if  (ni_1 < 2.5*step)
        plot(X(i),Y(i),'ks')
        n_1=n_1+1;
        i_tip1(n_1) = i;
        end
        
        if (ni_2 < 2.5*step)
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
key1 = 0;
key2 = 0;
for i=1:length(Y)-1
    if (key1==0)&&((Y(i)<=Yc(1))&&(Y(i+1)>Yc(1)))
        step=Y(i+1)-Y(i);
        ycen1=Y(i)+step/2;
        key1 = 1;
    end
     if (key2==0)&&((Y(i)<=Yc(length(Yc)))&&(Y(i+1)>Yc(length(Yc))))
        step=Y(i+1)-Y(i);
        ycen2=Y(i)+step/2;
        key2 = 1;
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
                            
sigma(:,2)=[0; 1]; %min
sigma(:,1)=[0; 2]; %max

figure
plot(X,Y,'k.')
hold on
plot(Xc,Yc,'k')

for i = 1:2
    
% Учет статических ГУ
    %МПа
F=Stat_BC(n,P,X,Y,kol,b,sigma(:,i),a,mas,i_edge,i_tip1,i_tip2,Xc,Yc);

% Учет кинематических ГУ
K_new=K_Kin_BC(X,Y,n,K,b,a);

% Вычисление вектора узловых перемещений
u=K_new\F;

r1 = 0.03;
r2 = 0.5;
% kin_tip1(i,:) = SIF_ex(1,n,P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc,G,k,u,r1,r2);
kin_tip2(i,:) = SIF_ex(2,n,P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc,G,k,u,r1,r2);

end

% crack1 = predict(kin_tip1,C_Paris,m_Paris);
crack2 = predict(kin_tip2,C_Paris,m_Paris,step); 

% grow1 = growth(Xc(1),Xc(1)-Xc(2),Yc(1),Yc(1)-Yc(2),crack1(3),crack1(4));
grow2 = growth(Xc(length(Xc)),Xc(length(Xc))-Xc(length(Xc)-1),Yc(length(Xc)),Yc(length(Xc))-Yc(length(Xc)-1),crack2(3),crack2(4));

ch=kol*20;
Xp_add=Xc(length(Xc)):(Xc(1)-Xc(length(Xc)))/ch:Xc(1);
for i = 1:length(Xc)-1
    for j = 1:length(Xp_add)
        if (Xp_add(j) >= Xc(i+1))&&(Xp_add(j) <= Xc(i))
Yp_add1(j)=(Yc(i)-Yc(i+1))/((Xc(i)-Xc(i+1))).*(Xp_add(j)-Xc(i+1))+Yc(i+1)-0.00001;
Yp_add2(j)=(Yc(i)-Yc(i+1))/((Xc(i)-Xc(i+1))).*(Xp_add(j)-Xc(i+1))+Yc(i+1)+0.00001;
        end
    end
end
Xp=[X,Xp_add,Xp_add];
Yp=[Y,Yp_add1,Yp_add2];
np=length(Xp);

figure
plot(Xp,Yp,'k.')

figure
hold on
for uzli=1:np
Nuz=N(n,Xp(uzli),Yp(uzli),P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc);
uuz(:,uzli)=Nuz*u;
X_def(uzli)=Xp(uzli)+uuz(1,uzli)*10^4*3;
Y_def(uzli)=Yp(uzli)+uuz(2,uzli)*10^4*3;
plot(X_def(uzli),Y_def(uzli),'k.')
end

% Xc = [grow1(1), Xc, grow2(1)];
% Yc = [grow1(2), Yc, grow2(2)];
Xc = [Xc, grow2(1)];
Yc = [Yc, grow2(2)];

figure
plot(X,Y,'k.')
hold on
plot(Xc,Yc,'k')
plot(Xc,Yc,'ko')

% history1(steplerchik,:) = [crack1(3), crack1(4)];
history2(steplerchik,:) = [crack2(3), crack2(4)];
% sifs1(steplerchik,:) = kin_tip1(2,:);
sifs2(steplerchik,:) = kin_tip2(2,:);
end

% al=1/2;
% F1=sqrt(sec(al*pi/2));
% F1_z=(1-0.025*al^2+0.06*al^4)*F1;
% Hura=sqrt(pi*a/4)*F1_z;
Elapsed = (cputime - Start) / 60
% %%
% r1 = 0.03;
% r2 = 0.5;
% Kin1_tip1_(1,:) = SIF_ex(2,n,P,X,Y,kol,a,i_edge,i_tip1,i_tip2,Xc,Yc,G,k,u,r1,r2);
% Kin_tip1_percent = [abs(Kin1_tip1_(1,1)-Kin1)/Kin1 abs(Kin1_tip1_(1,2)-Kin2)/Kin2]
