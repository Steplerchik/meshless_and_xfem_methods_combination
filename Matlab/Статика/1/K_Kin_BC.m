function w = K_Kin_BC(X,Y,n,K,b,a)

for uzli=1:n
if (Y(uzli)==0)
    for gr=1:length(K(1,:))
        K(gr,2*uzli)=0;
        K(2*uzli,gr)=0;
    end
K(2*uzli,2*uzli)=1;
end

if (X(uzli)==0)&&(Y(uzli)==0)
    hola = 1
    uzli
    for gr=1:length(K(1,:))
        K(gr,2*uzli-1)=0;
        K(2*uzli-1,gr)=0;
    end
K(2*uzli-1,2*uzli-1)=1;
end

if (X(uzli)==0)&&(Y(uzli)>b)
     hola = 2
    uz = uzli;
end

if (X(uzli)~=0)&&(Y(uzli)>b)
    hola = 3
    for gr=1:length(K(1,:))
        K(2*uzli,gr)=0;
    end
K(2*uzli,2*uzli)=1;
K(2*uzli,2*uz)=-1;
end




% 
% for uzli=1:n
% if (X(uzli)==2*a/11)&&(Y(uzli)==2*a/11)
%     for gr=1:length(K(1,:))
%         K(gr,2*uzli)=0;
%         K(2*uzli,gr)=0;
%     end
% K(2*uzli,2*uzli)=1;
% end
% 
% % if (X(uzli)==a/2)&&(Y(uzli)==0)
% if (X(uzli)==2*a/11)&&(Y(uzli)==2*a/11)
%     hola = 1
%     uzli
%     for gr=1:length(K(1,:))
%         K(gr,2*uzli-1)=0;
%         K(2*uzli-1,gr)=0;
%     end
% K(2*uzli-1,2*uzli-1)=1;
% end
% 
% if (X(uzli)==0)&&(Y(uzli)>b)
%     hola = 2
%     uzli
%     for gr=1:length(K(1,:))
%         K(gr,2*uzli-1)=0;
%         K(2*uzli-1,gr)=0;
%     end
% K(2*uzli-1,2*uzli-1)=1;
% end

% if (Y(uzli)>b)
%     for gr=1:length(K(1,:))
%         K(gr,2*uzli-1)=0;
%         K(2*uzli-1,gr)=0;
%     end
% K(2*uzli-1,2*uzli-1)=1;
% end


end
w=K;
end