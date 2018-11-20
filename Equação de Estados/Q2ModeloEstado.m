%Parâmetros de Sintonia

N=50;
lambda=10;

%Numero de Amostras
Amostras=200;

%Referência
R = ones(N,1);

%tf2ss(B,A);
%discretização fazendo sys=ss(A,B,C,D) 
%d_sys = c2d(sys,dt)
A= [0.1977 -2.539 -2.593;0.05403 0.8461 -0.1614;0.003363 0.09438 0.994];
B=[0.05403;0.003363;0.0001241];
C=[ 0  40  40];

x=[0; 0; 0];
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

uctrl=0;
%Simulação da Plana ss2tf(A,B,C,0)
for k=1:Amostras
    if(k-2)>0
        ym(k) = 1.7235*ym(k-1) -0.7403*ym(k-2) + 0.0136*uctrl(k-1) +0.0123*uctrl(k-2);
        x(:,k+1)=A*x(:,k)+B*uctrl(k);
    else
        if(k-1)>0
            ym(k)=1.7235*ym(k-1)+ 0.0136*uctrl(k-1);
            x(:,k+1)=A*x(:,k)+B*uctrl(k);
        else
            ym(k)=0;
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
grid on;
figure
plot(uctrl,'r');
grid on;