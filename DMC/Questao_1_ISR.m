%Modelo Gz 
N=15;
Ns=15;

B=[0.36 -0.2588];
A=[1 -0.4973 0.03668];
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