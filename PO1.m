%%

clc; clear all;


%Klasa 1 je bimodalna Gausova: 
%0.6N(M1,S1) + 0.4N(M2,S2)
%Klasa 2 je unimodalna Gausaova: N(M3,S3)

M1=[4;4]; %Srednja vrednost
S1=[2 -0.5; -0.5 2]; %Kovarijaciona matrica
M2=[4;8];
S2=[0.9 0.7;0.7 0.9]; %12 i 21 jednaki,ali 11 i 22 ne moraju
M3=[-4;6];
S3=[1.5 0.5;0.5 1.5];

%Transformacija bojenja(generisanja odbiraka): 
%T=F*L^0.5

%postoji i eigs i ona vraca sortirane sop. vred. od max kad min,medjutim,
%ne radi uvek u novijim verzijama

[Fi1,L1]=eig(S1);
T1=Fi1*L1^0.5;
[Fi2,L2]=eig(S2);  
T2=Fi2*L2^0.5;
[Fi3,L3]=eig(S3);
T3=Fi3*L3^0.5;

N=1000; %broj elemenata po klasi

%Klasa 2
X2=zeros(2,N);
for i=1:N
    X2(:,i)=T3*randn(2,1)+M3;      %dodajemo srednju vrednost
end

%Moze i prostije u jednom redu
%X2=T3*randn(2,N) + M3*ones(1,N);


%Klasa 1: 0.6 je verovatnoca da je odbirak generisan sa parametrima 
% M1,S1

X1=zeros(2,N);
for i=1:N
    if rand(1)<0.6
        X1(:,i)=T1*randn(2,1)+M1;  %dodamo srednju vrednost
    else
        X1(:,i)=T2*randn(2,1)+M2;  %dodamo srednju vrednost
    end
end

%p=rand(1,N);                      ILI OVAKO BRZE
%pnom=rand(2,N);
%X1=(p<=0.6).*(T1*pnom+M1*ones(1,N))+ (p>0.6).*(T2*pnom+M2*ones(1,N));

figure(1);
plot(X2(1,:),X2(2,:),'ro');
axis equal;
hold all;
plot(X1(1,:),X1(2,:),'b*');
legend('Klasa 2','Klasa 1','Location','SouthWest');
hold on;

%% prikazivanje fgv za K1 i K2
x=-8:0.1:8;
y=-1:0.1:11;
% ispred eksponencijalnog clana 1 sa 2pi na N/2 sa sigma na 1/2,
%a nema je 2,pa se zato izgubilo  1/(2pi)^(N/2)/S^0.5
f1=zeros(length(x),length(y));
f2=zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        X=[x(i);y(j)];
        f11=1/(2*pi)/det(S1)^0.5*exp(-0.5*(X-M1)'*S1^(-1)*(X-M1));
        f12=1/(2*pi)/det(S2)^0.5*exp(-0.5*(X-M2)'*S2^(-1)*(X-M2));
        f1(i,j)=0.6*f11+0.4*f12;
        f2(i,j)=1/(2*pi)/det(S3)^0.5*exp(-0.5*(X-M3)'*S3^(-1)*(X-M3));
    end
end

figure(2);
mesh(x,y,f1');
hold on;
mesh(x,y,f2'); %transponujemo jer se pakuje kolona po kolonu,pa da obrnemo
    
%U unimodalnom slucaju presecamo fgv nekom ravni i dobijamo neku elipsu

%% d2 krive
% za K2 mozemo jednostavnije
%generisemo krugove sa odgovarajucim poluprecnikom
% d=1,2,3 (d2=1,4,9) i primenimo T3
dk=[1 4 9];
for d2=dk
    teta=0:0.01:2*pi;
    Xpom=[d2^0.5*cos(teta);d2^0.5*sin(teta)];
    Y2=T3*Xpom+M3*ones(1,length(teta));  %od belih pravimo obojene,krug u elipsu
    figure(1);
    hold on;
    plot(Y2(1,:),Y2(2,:),'k.');
end

%f=C*exp(-d^2/2)
%fmax=C za d=0;
%fd2=1=fmax*exp(-1/2);...
f2max=max(max(f2));
prag1=f2max*exp(-1/2);
prag2=f2max*exp(-4/2);
prag3=f2max*exp(-9/2);
contour(x,y,f2',[prag1 prag1],'m')  %mora dvaput iz nepoznatog razloga []
contour(x,y,f2',[prag2 prag2],'m')
contour(x,y,f2',[prag3 prag3],'m')  
%sve tacke na elipsi imaju jednako statisticko odstojanje od centra klase

f1max=max(max(f1));
prag1=f1max*exp(-1/2);
prag2=f1max*exp(-4/2);
prag3=f1max*exp(-9/2);   %nadji mi sve tacke koje su dovoljno blizu nivoa koji trazim

Y11=[];
for i=1:length(x)
    for j=1:length(y)
        if abs(f1(i,j)-prag1)<5e-3
        Y11=[Y11 [x(i);y(j)]];
        end
    end
end
figure(1);
plot(Y11(1,:),Y11(2,:),'k.');

Y12=[];
for i=1:length(x)
    for j=1:length(y)
        if abs(f1(i,j)-prag2)<5e-3
        Y12=[Y12 [x(i);y(j)]];
        end
    end
end
figure(1);
plot(Y12(1,:),Y12(2,:),'k.');

Y13=[];
for i=1:length(x)
    for j=1:length(y)
        if abs(f1(i,j)-prag3)<5e-4
        Y13=[Y13 [x(i);y(j)]];
        end
    end
end
figure(1);
plot(Y13(1,:),Y13(2,:),'k.');

contour(x,y,f1',[prag1 prag1],'m')  
contour(x,y,f1',[prag2 prag2],'m')
contour(x,y,f1',[prag3 prag3],'m')

%% Testiranje hipoteza

%p1f1><p2f2
%l=f1/f2>< p2/p1
%h=-lnl=lnf2-lnf1=ln(f2/f1)><ln(p1/p2)
%p1=p2=1/2 sledi h><0

%Bajesov test minimalne verovatnoce greske

x=-8:0.1:8;
y=-1:0.1:11;
f1=zeros(length(x),length(y));
f2=zeros(length(x),length(y));
h=zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        X=[x(i);y(j)];
        f11=1/(2*pi)/det(S1)^0.5*exp(-0.5*(X-M1)'*S1^(-1)*(X-M1));
        f12=1/(2*pi)/det(S2)^0.5*exp(-0.5*(X-M2)'*S2^(-1)*(X-M2));
        f1(i,j)=0.6*f11+0.4*f12;
        f2(i,j)=1/(2*pi)/det(S3)^0.5*exp(-0.5*(X-M3)'*S3^(-1)*(X-M3));
        h(i,j)=log(f2(i,j))-log(f1(i,j));
    end
end
%koje su tacke na h gde ce h biti vece od nule ili manje od nule
figure(3);
plot(X2(1,:),X2(2,:),'ro');
axis equal;
hold all;
plot(X1(1,:),X1(2,:),'b*');
legend('Klasa 2','Klasa 1','Location','SouthWest');
contour(x,y,h',[0 0]);
hold off;
%Moze se desiti kad radimo domaci da nadjemo h i da se pojavi vise linija
%odnodno krivih,ako su Gausovske,dobijaju se krive drugog reda,sto znaci da
%moze doci do pojave hiperbole,parabole, itd. Ako imamo kovarijacione
%matrice iste,dobijamo zaista nesto sto je samo funkcija prvog reda,prava
%ili ravan,zavisi od prostora u kom radimo
%Ako imamo dve klase i jedna je izmedju dve krive,ova prva se vise
%rasipa,njen "rep" je veci,ako uzmemo random tacku negde,fgv1 ce biti vece
%od fgv2
%Sto vise tacaka generisemo,to ce procena biti bolja
%Treba naci gresku,integralimo fgv1 po regiji druge klase i fgv2 po regiji
%prve klase
%Radimo numericku integraciju,ideja je da izracunamo povrsinu koju zauzima
%prva i koju druga klasa
%Oblast 2 je h>0 i integralimo fgv1,za svaku kockicu vidimo vrednost fgv i
%sabiramo

%procena greske
Eps1=0;
Eps2=0;
for i=1:length(x)
    for j=1:length(y)
        if h(i,j)<0  %mozemo uzeti tacku ili tacku i tri oko nje,pa naci srednju vrednost
            Eps2=Eps2+f2(i,j)*0.1*0.1;   
            %Eps2=Eps2+(f2(i,j)+f2(i+1,j)+f2(i+1,j+1)+f2(i,j+1))/4*0.1*0.1;
        else
            Eps1=Eps1+f1(i,j)*0.1*0.1;
        end
    end
end
Eps1  %ispisi mi    1/2000=5.0000e-04
Eps2  %ispisi mi