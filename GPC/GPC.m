function [E,F,H,G]=GPC(A,B,M,N,lambda)

R = ones(N,1);
A_til = [A 0] - [0 A];
Amostras=15;
[E,F,H]=diofantina(A_til,B,N,M);

for i=1:N
    for j=1:M
        if(j==i)
            G(i,j)=H(j,1);
        elseif(j<i)
            G(i,j)=H(i,i-j+1);
        else
            G(i,j)=0;
        end
    end
end

for k=1:Amostras 
    
    %---- Ler da Planta 
    if (k-2)>0
        ym(k) = -A(2)*ym(k-1) + B(1)*uctrl(k-1) + B(2)*uctrl(k-2);
    else
        if (k-1)>0
            ym(k) = -A(2)*ym(k-1) + B(1)*uctrl(k-1);
        else 
            ym(k)=0;
        end
    end
    
    for i=1:N
        if(k-1>0)
            F_livre(i,1)=H(i,i+1)*du(k-1)+F(i,1)*ym(k)+F(i,2)*ym(k-1);
        else
            F_livre(i,1)=0;
        end
    end

    deltau=inv(G'*G+lambda*eye(M,M))*G'*(R-F_livre);
    du(k) = deltau(1);
    
    if (k-1)>0
        uctrl(k)=uctrl(k-1)+deltau(1);
    else
        uctrl(k)=deltau(1);
    end
end

figure
plot(ym,'b');
grid on
figure
plot(uctrl,'r');
grid on
