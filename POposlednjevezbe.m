
%%
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


Z=[X Y];     %smesti sve u 1
clear X Y;
z=rand(2*N,1);
X=[];
Y=[];
for i=1:2*N
    if z(i)<0.5
        X=[X Z(:,i)];
    else
        Y=[Y Z(:,i)];
    end
end
plot(X(1,:),X(2,:),'ro');    %Prikazi i izmesane sve zajedno
hold on;
plot(Y(1,:),Y(2,:),'bv'); hold on;
legend('Klasa 1','Klasa 2');

M1=mean(X,2);
M2=mean(Y,2);   %Moze i transponovanjem,ali kako vise volimo

lmax=100; %maksimalan broj iteracija
l=1;
a=0;      %ima reklasterizacije
while (l<lmax) && (a==0)
       a=1;
       X1=[];Y1=[];
       for i=1:max(size(X))    %kroz sve odbirke jedne klase
    %d1=pdist([X(:,i)';M1']);  moze formula,moze funkcija pdist,moze vektorski
           d1=((X(:,i)-M1)'*(X(:,i)-M1))^0.5;  
           d2=((X(:,i)-M2)'*(X(:,i)-M2))^0.5;
           if (d1>d2)
               Y1=[Y1 X(:,i)];
               a=0; %doslo je do reklasterizacije
           else
               X1=[X1 X(:,i)];  %ovde ne stavljamo a=1 da ne bismo pregazili info o reklasterizaciji
           end
       end
        for i=1:max(size(Y))    %kroz sve odbirke jedne klase
 %d1=pdist([X(:,i)';M1']);  moze formula,moze funkcija pdist,moze vektorski
        d1=((Y(:,i)-M1)'*(Y(:,i)-M1))^0.5;  
        d2=((Y(:,i)-M2)'*(Y(:,i)-M2))^0.5;
           if (d1>d2)
               Y1=[Y1 Y(:,i)];   %treba da ostane u Y
           else
               X1=[X1 Y(:,i)];
               a=0; %doslo je do reklasterizacije
        end
        end
        X=X1; Y=Y1;
      clear X1 Y1;
    M1=mean(X,2);
    M2=mean(Y,2);
    
    l=l+1;
figure(2);plot(X(1,:),X(2,:),'ro');    %Prikazi i izmesane sve zajedno
hold on;
plot(Y(1,:),Y(2,:),'bv'); hold on;
legend('Klasa 1','Klasa 2');
hold off
l=l+1;
end

%% Iris podaci klasifikacija cvetova

clc; clear all;
%redukcija dimenzija na bazi matrica rasejanja
load 'Iris_podaci.mat'
%U X se nalaze atributi(obelezja,po 4) za svaki od 150 cvetova
%U Y se nalazi informacija o klasi,tj. vrsti cveta

Xt=[X(1:25,:);X(51:75,:);X(101:125,:)];
Yt=[Y(1:25,:);Y(51:75,:);Y(101:125,:)];
Xo=[X(26:50,:);X(76:100,:);X(126:150,:)];
Yo=[Y(26:50,:);Y(76:100,:);Y(126:150,:)];

X1o=Xo(1:25,:);X1t=Xt(1:25,:);     %podela podataka,mogla je i lakse verovatno
X2o=Xo(26:50,:);X2t=Xt(26:50,:);
X3o=Xo(51:75,:);X3t=Xt(51:75,:);

M1=mean(X1o,1)';
S1=cov(X1o);

M2=mean(X2o,1)';
S2=cov(X2o);

M3=mean(X3o,1)';
S3=cov(X3o);

%Matrice rasejanja

Sw=1/3*S1+1/3*S2+1/3*S3;   %unutarklasno rasejanje
M0=1/3*M1+1/3*M2+1/3*M3;
Sb=1/3*(M1-M0)*(M1-M0)' +1/3*(M2-M0)*(M2-M0)'+1/3*(M3-M0)*(M3-M0)';
%Sb medjuklasno rasejanje

[V,D]=eig(Sw^-1*Sb);   %4 sopstvenih vrednosti jer imamo 4 obelezja
%metode za selekciju obelezja,od postojecih izaberu najbolja
%metode za ekstrakciju obelezja,kombinuju vec postojeca obelezja i od njih prave nova

A=V(:,1:2);   %transformaciona matrica
Xo_red=A'*Xo';
Xt_red=A'*Xt';

figure(1);
plot(Xo_red(1,1:25),Xo_red(2,1:25),'ro');hold on;
plot(Xo_red(1,26:50),Xo_red(2,26:50),'bv');hold on;
plot(Xo_red(1,51:75),Xo_red(2,51:75),'m*');hold off;
figure(2);
plot(Xt_red(1,1:25),Xt_red(2,1:25),'ro');hold on;
plot(Xt_red(1,26:50),Xt_red(2,26:50),'bv');hold on;
plot(Xt_red(1,51:75),Xt_red(2,51:75),'m*');hold off;

A=V(:,1:1);   %transformaciona matrica
Xo_red1=A'*Xo';
Xt_red1=A'*Xt';

figure(3);
plot(1:25,Xo_red1(1,1:25),'ro');hold on;
plot(26:50,Xo_red1(1,26:50),'bv');hold on;
plot(51:75,Xo_red1(1,51:75),'m*');hold off;
figure(4);
plot(1:25,Xt_red1(1,1:25),'ro');hold on;
plot(26:50,Xt_red1(1,26:50),'bv');hold on;
plot(51:75,Xt_red1(1,51:75),'m*');hold off;