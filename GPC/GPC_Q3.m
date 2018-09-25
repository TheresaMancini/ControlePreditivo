function [E,F,F_livre,H,G]=GPC_Q3(A,B,M,N,lambda,d)

A_til = [A 0] - [0 A];

R = ones(N,1);
Amostras=150;
[E,F,H]=diofantina(A_til,B,N,M);

%calculo do G
for i=1:N
    for j=1:M
        if(i<(d+1)||i==(length(A)-length(B)))
            G(i,j)=0;
        elseif(j<=i-2)
            G(i,j)=H(i,i-j-1);
        else
            G(i,j)=0;
        end
    end
end

for k=1:Amostras 
    
    %---- Ler da Planta 
    if(k-3)>0
        ym(k) = -A(2)*ym(k-1) -A(3)*ym(k-2) -A(4)*ym(k-3)+ B(1)*uctrl(k-1) + B(2)*uctrl(k-2)+B(3)*uctrl(k-3);
    else
        if (k-2)>0
            ym(k) = -A(2)*ym(k-1) -A(3)*ym(k-2) + B(1)*uctrl(k-1) + B(2)*uctrl(k-2);
        else
            if (k-1)>0
                ym(k) = -A(2)*ym(k-1) + B(1)*uctrl(k-1);
            else 
                ym(k)=0;
            end
        end
    end   
    
    for i=1:N
        if(i==1)
            if(k-4>1)
                F_livre(i,1)=H(i,i+2)*du(k-2)+H(i,i)*du(k-3)+H(i,i+1)*du(k-4)+F(i,1)*ym(k)+F(i,2)*ym(k-1)+F(i,3)*ym(k-2)+F(i,4)*ym(k-3); %olhar se isso ta certo.
            else
                if(k-3>1)
                     F_livre(i,1)=H(i,i+2)*du(k-2)+H(i,i)*du(k-3)+F(i,1)*ym(k)+F(i,2)*ym(k-1)+F(i,3)*ym(k-2)+F(i,4)*ym(k-3);
                else
                    if(k-2>1)
                         F_livre(i,1)=H(i,i+2)*du(k-2)+H(i,i)*du(k-3)+F(i,1)*ym(k)+F(i,2)*ym(k-1)+F(i,3)*ym(k-2);
                    else
                        if(k-1>1)
                            F_livre(i,1)=F(i,1)*ym(k)+F(i,2)*ym(k-1);
                        else
                            F_livre(i,1)=F(i,1)*ym(k);
                        end
                    end
                end
            end
        
        else
            if(k-4>1)
                F_livre(i,1)=H(i,i-1)*du(k-1)+H(i,i+2)*du(k-2)+H(i,i)*du(k-3)+H(i,i+1)*du(k-4)+F(i,1)*ym(k)+F(i,2)*ym(k-1)+F(i,3)*ym(k-2)+F(i,4)*ym(k-3); %olhar se isso ta certo.
            else
                if(k-3>1)
                     F_livre(i,1)=H(i,i-1)*du(k-1)+H(i,i+2)*du(k-2)+H(i,i)*du(k-3)+F(i,1)*ym(k)+F(i,2)*ym(k-1)+F(i,3)*ym(k-2)+F(i,4)*ym(k-3);
                else
                    if(k-2>1)
                         F_livre(i,1)=H(i,i-1)*du(k-1)+H(i,i+2)*du(k-2)+H(i,i)*du(k-3)+F(i,1)*ym(k)+F(i,2)*ym(k-1)+F(i,3)*ym(k-2);
                    else
                        if(k-1>1)
                            F_livre(i,1)=H(i,i-1)*du(k-1)+F(i,1)*ym(k)+F(i,2)*ym(k-1);
                        else
                            F_livre(i,1)=F(i,1)*ym(k);
                        end
                    end
                end
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
