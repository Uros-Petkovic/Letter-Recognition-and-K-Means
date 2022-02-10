
%% Generisanje klasa

clc; clear all; close all;

M1=[2;2];
S1=[1.5 0.5;0.5 1.5];
[Fi1,L1]=eig(S1);
T1=Fi1*L1^0.5;
N=500;

X1=zeros(2,N);
for i=1:N
    X1(:,i)=T1*randn(2,1)+M1; 
end

M2=[9;9];
S2=[1.1 0.5;0.5 1.2];
[Fi2,L2]=eig(S2);
T2=Fi2*L2^0.5;

X2=zeros(2,N);
for i=1:N
    X2(:,i)=T1*randn(2,1)+M2;
end

figure(1);
plot(X2(1,:),X2(2,:),'ro');
axis equal;
hold on;
plot(X1(1,:),X1(2,:),'b*');
xlabel('X1');
ylabel('X2');
title('Linearni klasifikator metodom zeljenog izlaza');
hold on;

%% Metod zeljenog izlaza a)

Gama=[ones(N,1); 1*ones(N,1)];  
Z=[ones(1,N) -ones(1,N);X1 -X2];  
W=(Z*Z')^(-1)*Z*Gama; %W=pinv(Z')*Gama
v0=W(1);
v=W(2:3);
x=[0,10];
%V(1)x+V(2)Y+VO
y=-v0/v(2)-v(1)/v(2).*x;
figure(1); plot(x,y,'k'); hold on;

%% Metod zeljenog izlaza b)

Gama=[ones(N,1); 3*ones(N,1)];   %Gama zeljeni izlazi,mogu se menjati,veca vaznost nekoj klasi
Z=[ones(1,N) -ones(1,N);X1 -X2];   %definicija Z,1 i -1,X i -X
W=(Z*Z')^(-1)*Z*Gama; %W=pinv(Z')*Gama
%W(1)=v0; [W(2);W(3)]=v;
v0=W(1);
v=W(2:3);
x=[0,10];
%V(1)x+V(2)Y+VO
y=-v0/v(2)-v(1)/v(2).*x;
figure(1); plot(x,y,'g'); hold on;

%% Metod zeljenog izlaza c)

Gama=[3*ones(N,1); 1*ones(N,1)];   %Gama zeljeni izlazi,mogu se menjati,veca vaznost nekoj klasi
Z=[ones(1,N) -ones(1,N);X1 -X2];   %definicija Z,1 i -1,X i -X
W=(Z*Z')^(-1)*Z*Gama;
v0=W(1);
v=W(2:3);
x=[0,10];
y=-v0/v(2)-v(1)/v(2).*x;
figure(1); plot(x,y,'b'); hold off;
legend('Klasa1','Klasa2','Jednaki','Prvi vazniji','Drugi vazniji','Location','SouthEast')

%% Druga iterativna metoda

%Generisanje

clc; clear all; close all;

M1=[2;2];
S1=[1.5 0.5;0.5 1.5];
[Fi1,L1]=eig(S1);
T1=Fi1*L1^0.5;

N=500;

X1=zeros(2,N);
for i=1:N
    X1(:,i)=T1*randn(2,1)+M1;   
end

M2=[9;9];
S2=[1.1 0.5;0.5 1.2];
[Fi2,L2]=eig(S2);
T2=Fi2*L2^0.5;

X2=zeros(2,N);
for i=1:N
    X2(:,i)=T1*randn(2,1)+M2;    
end

figure(1);
plot(X2(1,:),X2(2,:),'ro');
axis equal;
hold on;
plot(X1(1,:),X1(2,:),'b*');
hold on;


Nopt = 500;
v0opt= 0;
sopt= 0;

for s= 0:0.01:1
   V = (s*S1+(1-s)*S2)^(-1)*(M2-M1);
   sigma1 = V'*S1*V;
   sigma2 = V'*S2*V;
   ymin = V'*X1(:,1);
   ymax = V'*X1(:,1);
   for i= 1:N
      if(V'*X1(:,i)<ymin)    %Nalazim minimalno i maksimalno y zbog v0
          ymin = V'*X1(:,i);
      end
      if(V'*X2(:,i)<ymin) 
          ymin = V'*X2(:,i);
      end
      if(V'*X1(:,i)>ymax) 
          ymax = V'*X1(:,i);
      end
      if(V'*X2(:,i)>ymax) 
          ymax = V'*X2(:,i);
      end
   end
   v0 = floor(ymin):ceil(ymax); %Napravim v0 po minimalnom i maksimalnom y
   Nerr = 500;
   v0pom= 0;
   j=1;
   % v0 je oznaka za zapravo -v0
   for v0 = floor(ymin):1:ceil(ymax)
   Npom = 0;
   for i=1:N
      if(V'*X1(:,i)>v0) 
          Npom = Npom+1; 
      end
      if(V'*X2(:,i)<v0) 
          Npom = Npom+1; 
      end
   end
     if(Npom<Nerr)
        Nerr = Npom;
        v0pom = -v0;
     end
   end
   if (Nerr< Nopt)
   Nopt= Nerr;
   v0opt= v0pom;
   sopt= s;
   end
end
Vopt = (sopt*S1+(1-sopt)*S2)^(-1)*(M2-M1);

x= -2:0.1:12;
y= -2:0.1:12;
for i=1:length(x)
   for j=1:length(y)
      X=[x(i); y(j)];
      b(i,j)=Vopt'*X+v0opt;
   end
end
figure(1), hold on
contour(x,y, b', [0 0], 'g');
xlabel('X1');
ylabel('X2');
title('Linearni klasifikator drugom iterativnom metodom');
legend('Klasa1','Klasa2','Klasifikator','Location','NorthWest')
hold off
   
   
%% Kvadratni klasifikator  
   
clc; clear all;
N=500;

%Generisanje prve klase

x10= 5.5;
x20= 6.5;
for i=1:N
alfa = 2*pi*rand;
R= 2*rand;
x1= x10+ R*cos(alfa);
x2= x20+ R*sin(alfa);
X1(:,i)=[x1; x2];
end
figure(1);
plot(X1(1,:), X1(2,:), 'b*');

% Generisanje druge klase

y10=5.2;
y20=7.9;
for i=1:N
alfa =3*pi*rand/4+pi/4;
R= rand+1.5;
y1= y10+ R*cos(alfa);
y2= y20+ R*sin(alfa);
X2(:,i)=[y1; y2];
end
hold on;
plot(X2(1,:), X2(2,:), 'ro');
legend('Klasa 1','Klasa 2','location','SouthEast');   
   
%% Kvadratni a)  

Gama=[1*ones(N,1); 1*ones(N,1)];    %Gama zeljeni izlazi,mogu se menjati
Z=[ones(1,N) -ones(1,N); X1 -X2; X1(1,:).^2 -X2(1,:).^2; 2.*X1(1,:).*X1(2,:) -2.*X2(1,:).*X2(2,:); X1(2,:).^2 -X2(2,:).^2];
W=(Z*Z')^(-1)*Z*Gama;
v0=W(1); v=W(2:3); Q=[W(4) W(5);W(5) W(6)];

syms xp yp
[xp,yp,~,~]=solve(v0+v(1)*xp+v(2)*yp+xp^2*Q(1)+xp*yp*2*Q(2)+yp^2*Q(4),xp,yp,'returnconditions',true);
z=4:0.001:12;
xp=eval(xp); xplot=[xp(1,:) fliplr(xp(2,:))]; %da nastavi crta ceo krug npr
yp=eval(yp); yplot=[yp(1,:) fliplr(yp(2,:))];
xplot1=xplot(imag(xplot)==0);
yplot1=yplot(imag(xplot)==0);%za kompleksne
plot(xplot1,yplot1,'k'); hold on;  
   
%% Kvadratni b)   
   
Gama=[1*ones(N,1); 3*ones(N,1)];    %Gama zeljeni izlazi,mogu se menjati
Z=[ones(1,N) -ones(1,N); X1 -X2; X1(1,:).^2 -X2(1,:).^2; 2.*X1(1,:).*X1(2,:) -2.*X2(1,:).*X2(2,:); X1(2,:).^2 -X2(2,:).^2];
W=(Z*Z')^(-1)*Z*Gama;
v0=W(1); v=W(2:3); Q=[W(4) W(5);W(5) W(6)];

syms xp yp
[xp,yp,~,~]=solve(v0+v(1)*xp+v(2)*yp+xp^2*Q(1)+xp*yp*2*Q(2)+yp^2*Q(4),xp,yp,'returnconditions',true);
z=4:0.001:12;
xp=eval(xp); xplot=[xp(1,:) fliplr(xp(2,:))];
yp=eval(yp); yplot=[yp(1,:) fliplr(yp(2,:))];
xplot1=xplot(imag(xplot)==0);
yplot1=yplot(imag(xplot)==0);%za kompleksne
plot(xplot1,yplot1,'g');  hold on; 
   
   
%% Kvadratni c)

Gama=[3*ones(N,1); 1*ones(N,1)];    %Gama zeljeni izlazi,mogu se menjati
Z=[ones(1,N) -ones(1,N); X1 -X2; X1(1,:).^2 -X2(1,:).^2; 2.*X1(1,:).*X1(2,:) -2.*X2(1,:).*X2(2,:); X1(2,:).^2 -X2(2,:).^2];
W=(Z*Z')^(-1)*Z*Gama;
v0=W(1); v=W(2:3); Q=[W(4) W(5);W(5) W(6)];

syms xp yp
[xp,yp,~,~]=solve(v0+v(1)*xp+v(2)*yp+xp^2*Q(1)+xp*yp*2*Q(2)+yp^2*Q(4),xp,yp,'returnconditions',true);
z=4:0.001:12;
xp=eval(xp); xplot=[xp(1,:) fliplr(xp(2,:))]; %da nastavi crta ceo krug npr
yp=eval(yp); yplot=[yp(1,:) fliplr(yp(2,:))];
xplot1=xplot(imag(xplot)==0);
yplot1=yplot(imag(xplot)==0);%za kompleksne
plot(xplot1,yplot1,'b');
xlabel('X1');
ylabel('X2');
title('Kvadratni klasifikator metodom zeljenog izlaza');
legend('Klasa2','Klasa1','Jednaki','Prvi vazniji','Drugi vazniji','Location','SouthEast')
