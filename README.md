h=1;
T=100;
b=0.3;
c=0.1;
t=1:h:T;
S=zeros(1,length(t));
I=zeros(1,length(t));
R=zeros(1,length(t));
S(1)=990;
I(1)=10;
N=S(1)+I(1)+R(1);
k=2;
while k<=T
    S(k)=S(k-1)+((b/N)*-S(k-1)*I(k-1));
    I(k)=I(k-1)+((b/N)*S(k-1)*I(k-1)-c*I(k-1));
    R(k)=R(k-1)+(c*I(k-1));
    k=k+1;
end
figure(1)
plot(t,S,'g')
hold on
plot(t,I,'r')
plot(t,R,'b')
hold off