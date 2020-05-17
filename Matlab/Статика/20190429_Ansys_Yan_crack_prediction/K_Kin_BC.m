function w = K_Kin_BC(X,Y,n,K,b)

for uzli=1:n
if (Y(uzli)==0)
    for gr=1:length(K(1,:))
        K(gr,2*uzli)=0;
        K(2*uzli,gr)=0;
    end
K(2*uzli,2*uzli)=1;
end

if (X(uzli)==0)
    for gr=1:length(K(1,:))
        K(gr,2*uzli-1)=0;
        K(2*uzli-1,gr)=0;
    end
K(2*uzli-1,2*uzli-1)=1;
end

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