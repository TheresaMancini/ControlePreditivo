% Planta de fase não mínima -> Tem um zero no semiplano direito. 
%Consequência a planta tem um undershoot (começa ao contrário) e só depois
%vai pra referência
clear all;
clc;

%Parametros do controlador

N = 15;  %Horizonte de predição 
M = 10;  %Ações de Controle 
lambda = 10; %Controla a agressvidade, menor lambida mais agressivo
Ns = 50; %Periodo de amostragem -> Truncamentos no momento em que a resposta ao degrau estabiliza
R = ones(N,1);


%Dados da planta simulada
% Código para discretizar a partir de função de transferência 
% c2d(G,0.2) -> Parametros tf e periodo de amostragem 

% Dado da questão/situação 
B =  [0.01036 0.0107]; %Numérador da função de tranferência do domínio Z
A =  [1 -2.126 1.105]; %Denominador da função de tranferência do domínio Z
PeriodoAmostragem = 0.1; 

%--- Dados da Simulação 
TempoSimulacao = 10;
k=1;
Amostras = TempoSimulacao/PeriodoAmostragem;


%---- Aplicação do Degrau Unitário na Planta 

udegrau = ones(1,N+Ns);
g=ones(1,N+Ns);
for k=1:N+Ns
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
    if (k-2)>0
        ym(k) = -A(2)*ym(k-1) -A(3)*ym(k-2) + B(1)*uctrl(k-1) + B(2)*uctrl(k-2);
    else
        if (k-1)>0
            ym(k) = -A(2)*ym(k-1) + B(1)*uctrl(k-1);
        else 
            ym(k)=0;
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
    
    deltau=inv(G'*G+lambda*eye(M,M))*G'*(R-F);
    du(k)  = deltau(1);
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