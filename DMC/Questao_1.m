%10x

N=16;
Ns=15;

B=[0.36 -0.2588];
A=[1 -0.4973 0.03668];

t = 0:0.3:9;
u=sin(2*t);
for i=1:30
    if(i==1)
        deltau(i)=u(i);
    else
        deltau(i)=u(i)-u(i-1);
    end
end
y=zeros(1,N);
for k=1:N
    for i=1:N
        if(k-i>1) 
            y(k)=y(k)+g(i)*deltau(k-i);
        end
    end
end

plot(y);