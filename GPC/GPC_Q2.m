function [E,F,F_livre,H,G]=GPC_Q2(A,B,M,N,lambda,d)

A_til = [A 0] - [0 A];

R = ones(N,1);
Amostras=200;
[E,F,H]=diofantina(A_til,B,N,M);
%calculo do G
for i=1:N
    for j=1:M
        if(i==j)
            G(i,j)=0;
        elseif(j<i)
            G(i,j)=H(i,i-j);
        else
            G(i,j)=0;
        end
    end
end


for k=1:Amostras 
    
    %---- Ler da Planta 
    if (k-2)>0
        ym(k) = -A(2)*ym(k-1) -A(3)*ym(k-2) + B(1)*uctrl(k-1) + B(2)*uctrl(k-2);
    else
        if (k-1)>0
            ym(k) = -A(2)*ym(k-1) + B(1)*uctrl(k-1);
        else 
            ym(k)=0;
        end
    end
    
     for i=1:N            
        if(k-1>1)
            F_livre(i,1)=H(i,i)*du(k-1)+H(i,i+1)*du(k-2)+F(i,1)*ym(k)+F(i,2)*ym(k-1)+F(i,3)*ym(k-2); %olhar se isso ta certo.
        else
            if(k-1>0)
                 F_livre(i,1)=H(i,i+1)*du(k-1)+F(i,1)*ym(k)+F(i,2)*ym(k-1);
            else %precisa ajustar o H e os F 
                F_livre(i,1)=F(i,1)*ym(k);
            end
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


subplot(2,1,1);
plot(ym,'b');
title('Saida da Planta')

subplot(2,1,2);
plot(uctrl,'r');
title('Entrada de Controle')

grid on
