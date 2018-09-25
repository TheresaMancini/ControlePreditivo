function [E,F,H]=diofantina(A_til, B,N,M)

na=length(A_til);
E(1,1)=1;

for i=2:na
    F(1,i-1)=-A_til(i);
end

H(1,:)=conv(E(1,:),B);
for i=2:N
    Rj=F(i-1,1);
    for j=2:na
        if(j<=na-1)
            F(i,j-1)=F(i-1,j)-Rj*A_til(j);
        else 
            F(i,j-1)=-Rj*A_til(j);
        end 
    end
    for j=1:i
        if(j==i)
            E(i,j)=Rj;
        else
            E(i,j)=E(i-1,j);
        end
    end
    %Calculando H 
   H(i-1,length(H(1,:))+1)=0;
   H(i,:)=conv(E(i,:),B);
    
end
