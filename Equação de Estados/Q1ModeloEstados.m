%Parâmetros de Sintonia

N=50;
lambda=10;

%Numero de Amostras
Amostras=500;

%Referência
R = 15*ones(N,1);

%Definição da Equação de Estado da Planta
% A=[0.7326 -0.1722; 0.08611 0.9909];
% B=[0.08611; 0.004528];
% C=[0 3];

A=[1.9354 -0.9365;1.0000 0];
B=[1;0];
C=[0.0024 0.0023];


P=zeros(N,N);
for i=1:N
    for j=1:N
        if(j<=i)
            P(i,j)=C*A^(i-j)*B;
        end
    end
end
[r,c]=size(C);
[rA,cA]=size(A);
Q=zeros(r,cA);
for i=1:N
    Q(i,:)=C*A^(i);
end



T_N=tril(ones(N,N));

%Calculo da Resposta Forçada
G=T_N*P;

x=[15; 15];
uctrl=0;
%Simulação da Plana ss2tf(A,B,C,0)
for k=1:Amostras
    if(k-2)>0
        ym(k) = 1.9354*ym(k-1) -0.9365*ym(k-2) + 0.00238*uctrl(k-1) + 0.002327*uctrl(k-2);
        x(:,k+1)=A*x(:,k)+B*uctrl(k);
    else
        if(k-1)>0
            ym(k)=1.9354*ym(k-1)+ 0.00238*uctrl(k-1);
            x(:,k+1)=A*x(:,k)+B*uctrl(k);
        else
            ym(k)=15;
            x(:,k+1)=A*x(:,k)+B*0;
        end
    end
    
    %Calculo da Resposta Livre
    if(k-1)>0
        delta_x=x(:,k)-x(:,k-1);
    else
        delta_x=x(:,k);
    end
    F=T_N*Q*delta_x+ones(N,1)*ym(k);
    
    %calculo da entrada de controle
    deltau=inv(G'*G+lambda*eye(N,N))*G'*(R-F);
    du(k)  = deltau(1);
    if (k-1)>0
        uctrl(k+1)=uctrl(k)+deltau(1);
    else
        uctrl(k+1)=deltau(1);
    end
end

figure
plot(ym,'b'); 
title('saída');
grid on;
figure
plot(uctrl,'r');
title('sinal de controle');
grid on;