

%% CETVRTE VEZBE INICIJALIZACIJA

clear all; close all; clc;
N=1000;

%% GENERISANJE PODATAKA

%PRVA KLASA
M1=[2;2];
R1=2;
Tetax=rand(1,N)*2*pi;      %random vrednost ugla
Rx=rand(1,N)*R1;           %random Rx
X=[Rx.*cos(Tetax);Rx.*sin(Tetax)]+M1*ones(1,N);
%DRUGA KLASA
M2=[8;8];
R2=3; d=2;
Tetay=rand(1,N)*2*pi;
Ry=rand(1,N)*d+R2;
Y=[Ry.*cos(Tetay);Ry.*sin(Tetay)]+M2*ones(1,N);
figure; plot(X(1,:),X(2,:),'ro');
hold on;
plot(Y(1,:),Y(2,:),'bv');
legend('Klasa 1','Klasa 2'); hold on;

%% Linearni klasifikator na bazi zeljenog izlaza

Gama=[ones(N,1); 2*ones(N,1)];   %Gama zeljeni izlazi,mogu se menjati
Z=[ones(1,N) -ones(1,N);X -Y];   %definicija Z,1 i -1,X i -X
W=(Z*Z')^(-1)*Z*Gama; %W=pinv(Z')*Gama
%W(1)=v0; [W(2);W(3)]=v;
v0=W(1);
v=W(2:3);
x=[0,12];
%V(1)x+V(2)Y+VO
y=-v0/v(2)-v(1)/v(2).*x;   %Zasto je ovo y?
figure(1); plot(x,y,'k'); 

%% kvadratni klasifikator
% x'*Q*x+v'*x+v0=0
clear all; close all; clc;
%PRVA KLASA
N=1000;
M1=[8;8];
R1=2;
Tetax=rand(1,N)*2*pi;
Rx=rand(1,N)*R1;
X=[Rx.*cos(Tetax);Rx.*sin(Tetax)]+M1*ones(1,N);
%DRUGA KLASA
M2=[8;8];
R2=3; d=2;
Tetay=rand(1,N)*2*pi;
Ry=rand(1,N)*d+R2;
Y=[Ry.*cos(Tetay);Ry.*sin(Tetay)]+M2*ones(1,N);
figure; plot(X(1,:),X(2,:),'ro');
hold on;
plot(Y(1,:),Y(2,:),'bv'); hold on;
legend('Klasa 1','Klasa 2');
%%
Gama=[ones(N,1); 2*ones(N,1)];    %Gama zeljeni izlazi,mogu se menjati
Z=[ones(1,N) -ones(1,N); X -Y; X(1,:).^2 -Y(1,:).^2; 2.*X(1,:).*X(2,:) -2.*Y(1,:).*Y(2,:); X(2,:).^2 -Y(2,:).^2];
W=(Z*Z')^(-1)*Z*Gama;
v0=W(1); v=W(2:3); Q=[W(4) W(5);W(5) W(6)];

syms xp yp
[xp,yp,~,~]=solve(v0+v(1)*xp+v(2)*yp+xp^2*Q(1)+xp*yp*2*Q(2)+yp^2*Q(4),xp,yp,'returnconditions',true);
z=4:0.001:12;
xp=eval(xp); xplot=[xp(1,:) fliplr(xp(2,:))]; %da nastavi crta ceo krug npr
yp=eval(yp); yplot=[yp(1,:) fliplr(yp(2,:))];
xplot1=xplot(imag(xplot)==0);
yplot1=yplot(imag(xplot)==0);%za kmpleksne
plot(xplot1,yplot1,'k');


