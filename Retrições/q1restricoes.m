clear all;
clc;

%Parametros do controlador

N = 15;  %Horizonte de predição 
M = 10;  %Ações de Controle 
lambda = 10; %Controla a agressvidade, menor lambida mais agressivo
Ns = 50; %Periodo de amostragem -> Truncamentos no momento em que a resposta ao degrau estabiliza
R = ones(N,1);

%------- Dado da questão
B =  [0.6445 -0.7176 -0.09059]; %Numérador da função de tranferência do domínio Z
A =  [1 -0.9183 0.084 -0.002029]; %Denominador da função de tranferência do domínio Z
PeriodoAmostragem = 0.2; 

%--- Dados da Simulação 
TempoSimulacao = 20;
k=1;
Amostras = TempoSimulacao/PeriodoAmostragem;


%---- Restrições

dumax=0;
dumin=-0.1;

ymax=2;
ymin=3;

umax=-0.4;
umin=-1;


%---- Aplicação do Degrau Unitário na Planta 

udegrau = ones(1,N+Ns);
g=ones(1,N+Ns);
for k=1:N+Ns
    if(k-3)>0
        g(k) = -A(2)*g(k-1) -A(3)*g(k-2) -A(4)*g(k-3) + B(1)*udegrau(k-1) +B(2)*udegrau(k-2) + B(3)*udegrau(k-3);
    else 
        if (k-2)>0
            g(k) = -A(2)*g(k-1) -A(3)*g(k-2) + B(1)*udegrau(k-1) + B(2)*udegrau(k-2);
        else
            if (k-1)>0
                g(k) = -A(2)*g(k-1) + B(1)*udegrau(k-1);
            else 
                g(k)=0;
            end
        end
    end
end


%----- Calcular a Matriz de Resposta Forçada 
G = zeros(N,M);

for i=1:N
    for j=1:M
        if(j<=i)
            G(i,j) = g(i-j+1);
        end
    end
end


%----- Executar o Controlador 
ym=zeros(1,Amostras); %escrevi depois
for k=1:Amostras % até aqui foi copiado do prof . %Quando k=67 problema. Excede tamanho do vetor g; 
    
    %---- Ler da Planta 
    if(k-3)>0
        ym(k) = -A(2)*ym(k-1) -A(3)*ym(k-2) -A(4)*ym(k-3) + B(1)*uctrl(k-1) +B(2)*uctrl(k-2) + B(3)*uctrl(k-3);
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
    
    %---- Calcular Vetor de Resposta Livres (F) vai de 1 até N
    F=zeros(N,1);
    soma=0;
    for j=1:N 
        soma=0;
        for i=1:Ns
            if(k-i)>0
                soma=soma+(g(i+j)-g(i))*du(k-i); %quando k>11 k-i>10
                %soma=0;
            end
        end
        F(j,1)=soma+ym(k);
    end
    
   % deltau=inv(G'*G+lambda*eye(M,M))*G'*(R-F);
    gg=2*(G'*G+lambda*eye(M,M));
    f=2*(F-R)'*G;
   
    %S=[tril(ones(M,M));-tril(ones(M,M))];
    S=[eye(M,M);-eye(M,M)];
    a=[ones(M,1)*dumax;-ones(M,1)*dumin];
%     if(k-1>0)
%         a=[ones(M,1)*(umax-uctrl(k-1));ones(M,1)*(uctrl(k-1)-umin)];
%     else
%         a=[ones(M,1)*(umax);ones(M,1)*(umin)];
%     end
        
    deltau=quadprog(gg,f',S,a);
    
    du(k)  = deltau(1);
    if (k-1)>0
        uctrl(k)=uctrl(k-1)+deltau(1);
    else
        uctrl(k)=deltau(1);
    end

end

for i=1:Amostras
    x(i)=i*PeriodoAmostragem;
end


figure
plot(x,ym,'b'); 
title("saida");
figure
plot(x,uctrl,'r');
title("entrada");
figure
plot(x,du);
title("varriação de controle");