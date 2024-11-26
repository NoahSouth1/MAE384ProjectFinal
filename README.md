h=0.1;
T=30;
c=0.1;
t=1:h:T;
b0=0.3;
%Variable Beta
b=zeros(1,length(t));
k=1;
while k<=length(t)
b=b0*(1+5*sin(2*pi*(365/365)*t));
    k=k+1;
end
%
S=zeros(1,length(t));
I=zeros(1,length(t));
R=zeros(1,length(t));
S(1)=990;
I(1)=10;
N=S(1)+I(1)+R(1);
k=2;
while k<=length(t)
    S(k)=S(k-1)-((b(k-1)/N)*S(k-1)*I(k-1))*h;
    I(k)=I(k-1)+((b(k-1)/N)*S(k-1)*I(k-1)-c*I(k-1))*h;
    R(k)=R(k-1)+(c*I(k-1))*h;
    k=k+1;
end
figure(1)
plot(t,S,'g')
hold on
plot(t,I,'r')
plot(t,R,'b')
hold off
line=zeros(1,length(t));
k=1;
while k<=length(t)
    line(k)=0.3;
    k=k+1;
end
figure(2)
plot(t,b)
hold on
plot(t,line)
hold off

