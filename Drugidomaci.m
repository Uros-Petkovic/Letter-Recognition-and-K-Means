
%% Generisanje klasa i odbiraka i njihovo prikazivanje

clc; clear all; close all;

% Napraviti dve bimodalne Gausove raspodele po ugledu na primer 4.3

%Prva klasa
M11=[1;1];               %Srednja vrednost
S11=[4 1.1; 1.1 2];      %Kovarijaciona matrica
M12=[6;4];
S12=[3 -0.8;-0.8 1.5]; 
P11=0.6; P12=0.4;
%Druga klasa
M21=[7;-4];               %Srednja vrednost
S21=[2 1.1; 1.1 4];       %Kovarijaciona matrica
M22=[6;0];
S22=[3 0.8;0.8 0.5];
P21=0.55; P22=0.45;

%Generisanje odbiraka prve klase transformacijom bojenja

[Fi11,L11]=eig(S11);
T11=Fi11*L11^0.5;
[Fi12,L12]=eig(S12);  
T12=Fi12*L12^0.5;

N=500;
X1=zeros(2,N);
for i=1:N
    if rand(1)<0.6
        X1(:,i)=T11*randn(2,1)+M11;  %dodamo srednju vrednost
    else
        X1(:,i)=T12*randn(2,1)+M12;  %dodamo srednju vrednost
    end
end

%Generisanje odbiraka druge klase transformacijom bojenja

[Fi21,L21]=eig(S21);
T21=Fi21*L21^0.5;
[Fi22,L22]=eig(S22);  
T22=Fi22*L22^0.5;

N=500;
X2=zeros(2,N);
for i=1:N
    if rand(1)<0.55
        X2(:,i)=T21*randn(2,1)+M21;  %dodamo srednju vrednost
    else
        X2(:,i)=T22*randn(2,1)+M22;  %dodamo srednju vrednost
    end
end

figure(1);
plot(X2(1,:),X2(2,:),'ro');
axis equal;
hold on;
plot(X1(1,:),X1(2,:),'b*');
hold on;

%% Generisanje funkcija gustina verovatnoca za obe klase

x=-8:0.1:12;
y=-10:0.1:8;
f1=zeros(length(x),length(y));
f2=zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        X=[x(i);y(j)];
        f11=1/(2*pi)/det(S11)^0.5*exp(-0.5*(X-M11)'*S11^(-1)*(X-M11));
        f12=1/(2*pi)/det(S12)^0.5*exp(-0.5*(X-M12)'*S12^(-1)*(X-M12));
        f1(i,j)=0.6*f11+0.4*f12;
        f21=1/(2*pi)/det(S21)^0.5*exp(-0.5*(X-M21)'*S21^(-1)*(X-M21));
        f22=1/(2*pi)/det(S22)^0.5*exp(-0.5*(X-M22)'*S22^(-1)*(X-M22));
        f2(i,j)=0.55*f21+0.45*f22;
    end
end

figure(2);
mesh(x,y,f1');
hold on;
mesh(x,y,f2');

%% Prikazivanje d2 krivih sa konstantnom fgv

f1max=max(max(f1));
prag1=f1max*exp(-1/2);
prag2=f1max*exp(-4/2);
prag3=f1max*exp(-9/2);   %nadji mi sve tacke koje su dovoljno blizu nivoa koji trazim

Y11=[];
for i=1:length(x)
    for j=1:length(y)
        if abs(f1(i,j)-prag1)<2e-3
        Y11=[Y11 [x(i);y(j)]];
        end
    end
end
figure(1);
plot(Y11(1,:),Y11(2,:),'g.'); hold on;

Y12=[];
for i=1:length(x)
    for j=1:length(y)
        if abs(f1(i,j)-prag2)<7e-4
        Y12=[Y12 [x(i);y(j)]];
        end
    end
end
figure(1);
plot(Y12(1,:),Y12(2,:),'g.'); hold on;
Y13=[];
for i=1:length(x)
    for j=1:length(y)
        if abs(f1(i,j)-prag3)<7e-5
        Y13=[Y13 [x(i);y(j)]];
        end
    end
end
figure(1);
plot(Y13(1,:),Y13(2,:),'g.'); hold on;

contour(x,y,f1',[prag1 prag1],'m')  
contour(x,y,f1',[prag2 prag2],'m')
contour(x,y,f1',[prag3 prag3],'m')


f2max=max(max(f2));
prag1=f2max*exp(-1/2);
prag2=f2max*exp(-4/2);
prag3=f2max*exp(-9/2);  %nadji mi sve tacke koje su dovoljno blizu nivoa koji trazim

Y21=[];
for i=1:length(x)
    for j=1:length(y)
        if abs(f2(i,j)-prag1)<0.5e-2
        Y21=[Y21 [x(i);y(j)]];
        end
    end
end
figure(1);
plot(Y21(1,:),Y21(2,:),'y.'); hold on

Y22=[];
for i=1:length(x)
    for j=1:length(y)
        if abs(f2(i,j)-prag2)<1.5e-3
        Y22=[Y22 [x(i);y(j)]];
        end
    end
end
figure(1);
plot(Y22(1,:),Y22(2,:),'y.'); hold on;

Y23=[];
for i=1:length(x)
    for j=1:length(y)
        if abs(f2(i,j)-prag3)<3e-4
        Y23=[Y23 [x(i);y(j)]];
        end
    end
end
figure(1);
plot(Y23(1,:),Y23(2,:),'y.'); hold on;

contour(x,y,f2',[prag1 prag1],'m')  
contour(x,y,f2',[prag2 prag2],'m')
contour(x,y,f2',[prag3 prag3],'m')

%% Bajesov test minimalne greske

%  p1f1><p2f2
%  l=f1/f2>< p2/p1
%  h=-lnl=lnf2-lnf1=ln(f2/f1)><ln(p1/p2)
%  p1=p2=1/2 sledi h><0

x=-8:0.1:12;
y=-10:0.1:8;
h=zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        X=[x(i);y(j)];
        h(i,j)=log(f2(i,j))-log(f1(i,j));
    end
end
%Koje su tacke na h gde ce h biti vece od nule ili manje od nule
figure(3);
plot(X2(1,:),X2(2,:),'ro');
axis equal;
hold all;
plot(X1(1,:),X1(2,:),'b*');
legend('Klasa 2','Klasa 1','Location','SouthWest');
contour(x,y,h',[0 0]);
hold off;


%Nalazenje greske prvog i drugog reda
Eps11=0;
Eps21=0;
for i=1:length(x)
    for j=1:length(y)
        if h(i,j)<0   %P1=P2=1/2 
            Eps21=Eps21+f2(i,j)*0.1*0.1;   
        else
            Eps11=Eps11+f1(i,j)*0.1*0.1;
        end
    end
end

%% Bajesov test minimalne cene

%  f1/f2 >< (C12-C22)/(C21-C11)*P2/P1
%  h=ln(f2/f1) >< ln(p1(C21-C11)/p2(C12-C22))
%  p1 = p2 = 1/2 sledi h>< ln(C21-C11)/(C22-C12)
%  C11=C22=0;
%  a) C21=C11
%  b) C21=3C11
%  c) C11=3C21


x=-8:0.1:12;
y=-10:0.1:8;
h=zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        X=[x(i);y(j)];
        h(i,j)=log(f2(i,j))-log(f1(i,j));
    end
end

%  a) resenje   C12=C21

P1=0.5; P2=0.5;
C11=0; C22=0;
C21=5; C12=5;   %C12=C21
Temp=real(log((C21-C11)/(C22-C12)*P1/P2));

figure(3);
plot(X2(1,:),X2(2,:),'ro');
axis equal;
hold all;
plot(X1(1,:),X1(2,:),'b*');
legend('Klasa 2','Klasa 1','Location','SouthWest');
contour(x,y,h',[Temp Temp]);
hold on;


%Nalazenje greske prvog i drugog reda
Eps12=0;
Eps22=0;
for i=1:length(x)
    for j=1:length(y)
        if h(i,j)<Temp   %P1=P2=1/2 
            Eps22=Eps22+f2(i,j)*0.1*0.1;   
        else
            Eps12=Eps12+f1(i,j)*0.1*0.1;
        end
    end
end

%  B) resenje   C12=25*C21

P1=0.5; P2=0.5;
C11=0; C22=0;
C21=5; C12=25;   %C12=25*C21
Temp=real(log((C21-C11)/(C22-C12)*P1/P2));

figure(3);
contour(x,y,h',[Temp Temp]);
hold on;


%Nalazenje greske prvog i drugog reda
Eps13=0;
Eps23=0;
for i=1:length(x)
    for j=1:length(y)
        if h(i,j)<Temp   %P1=P2=1/2 
            Eps23=Eps23+f2(i,j)*0.1*0.1;   
        else
            Eps13=Eps13+f1(i,j)*0.1*0.1;
        end
    end
end

%  c) resenje   25*C12=C21

P1=0.5; P2=0.5;
C11=0; C22=0;
C21=25; C12=5;   %25*C12=C21
Temp=real(log((C21-C11)/(C22-C12)*P1/P2));

figure(3);
contour(x,y,h',[Temp Temp],'k');
hold off;


%Nalazenje greske prvog i drugog reda
Eps14=0;
Eps24=0;
for i=1:length(x)
    for j=1:length(y)
        if h(i,j)<Temp   %P1=P2=1/2 
            Eps24=Eps24+f2(i,j)*0.1*0.1;   
        else
            Eps14=Eps14+f1(i,j)*0.1*0.1;
        end
    end
end

%% Wald-ov test

clc; clear all; close all;

%Prva klasa
M11=[1;1];               %Srednja vrednost
S11=[4 1.1; 1.1 2];      %Kovarijaciona matrica
M12=[6;4];
S12=[3 -0.8;-0.8 1.5]; 
P11=0.6; P12=0.4;
%Druga klasa
M21=[7;-4];               %Srednja vrednost
S21=[2 1.1; 1.1 4];       %Kovarijaciona matrica
M22=[6;0];
S22=[3 0.8;0.8 0.5];
P21=0.55; P22=0.45;

[Fi11,L11]=eig(S11);
T11=Fi11*L11^0.5;
[Fi12,L12]=eig(S12);  
T12=Fi12*L12^0.5;

[Fi21,L21]=eig(S21);
T21=Fi21*L21^0.5;
[Fi22,L22]=eig(S22);  
T22=Fi22*L22^0.5;

    
    
e1=10^-8;
e2=10^-8;

A=(1-e1)/e2;
B=e1/(1-e2);
a=-log(A);
b=-log(B);


    
for k=1:200
        odluka=false;
        Smm=0;
        Sm=0;
        i=1;
        while(odluka==false)
            if rand(1)<0.6
               x1=T11*randn(2,1)+M11;  %dodamo srednju vrednost
            else
               x1=T12*randn(2,1)+M12;  %dodamo srednju vrednost
            end
            X=[x1(1,1);x1(2,1)];
            f11=1/(2*pi)/det(S11)^0.5*exp(-0.5*(X-M11)'*S11^(-1)*(X-M11));
            f12=1/(2*pi)/det(S12)^0.5*exp(-0.5*(X-M12)'*S12^(-1)*(X-M12));
            f1=0.6*f11+0.4*f12;
            f21=1/(2*pi)/det(S21)^0.5*exp(-0.5*(X-M21)'*S21^(-1)*(X-M21));
            f22=1/(2*pi)/det(S22)^0.5*exp(-0.5*(X-M22)'*S22^(-1)*(X-M22));
            f2=0.55*f21+0.45*f22;
            Smm=Smm-log(f1/f2);
            Sm(i)=Smm;
            i=i+1;
            if (Smm<a || Smm>b)
                odluka=true;
            end
        end
        figure(1);
        plot(0:length(Sm)-1,Sm,'r'); hold on;
        title('Wald-ov sekvencijalni test');
        xlabel('Broj potrebnih odbiraka');
        ylabel('Sm');
end
    

for k=1:100
        odluka=false;
        Smm=0;
        Sm=0;
        i=1;
        while(odluka==false)
            if rand(1)<0.55
               x1=T21*randn(2,1)+M21;  %dodamo srednju vrednost
            else
               x1=T22*randn(2,1)+M22;  %dodamo srednju vrednost
            end
            X=[x1(1,1);x1(2,1)];
            f11=1/(2*pi)/det(S11)^0.5*exp(-0.5*(X-M11)'*S11^(-1)*(X-M11));
            f12=1/(2*pi)/det(S12)^0.5*exp(-0.5*(X-M12)'*S12^(-1)*(X-M12));
            f1=0.6*f11+0.4*f12;
            f21=1/(2*pi)/det(S21)^0.5*exp(-0.5*(X-M21)'*S21^(-1)*(X-M21));
            f22=1/(2*pi)/det(S22)^0.5*exp(-0.5*(X-M22)'*S22^(-1)*(X-M22));
            f2=0.55*f21+0.45*f22;
            Smm=Smm-log(f1/f2);
            Sm(i)=Smm;
            i=i+1;
            if (Smm<a || Smm>b)
                odluka=true;
            end
        end
        figure(1);
        plot(0:length(Sm)-1,Sm,'b'); hold on;
end

%%   Zavisnosti od gresaka

eks=-20:0.1:0;
e1niz=10.^eks;
e2niz=10.^eks;
    
%e1=10^-5;
e2=10^-5;

%Zavisnost od e1
z=1;
brodb=zeros(1,length(e1niz));
for e1=e1niz
    
   A=(1-e1)/e2;
   B=e1/(1-e2);
   a=-log(A);
   b=-log(B);


   brodbiter=zeros(1,200);
   m=1;
   for k=1:200
   odluka=false;
   Smm=0;
   Sm=0;
   i=1;
   while(odluka==false)
       if rand(1)<0.6
          x1=T11*randn(2,1)+M11;  %dodamo srednju vrednost
       else
          x1=T12*randn(2,1)+M12;  %dodamo srednju vrednost
       end
       X=[x1(1,1);x1(2,1)];
       f11=1/(2*pi)/det(S11)^0.5*exp(-0.5*(X-M11)'*S11^(-1)*(X-M11));
       f12=1/(2*pi)/det(S12)^0.5*exp(-0.5*(X-M12)'*S12^(-1)*(X-M12));
       f1=0.6*f11+0.4*f12;
       f21=1/(2*pi)/det(S21)^0.5*exp(-0.5*(X-M21)'*S21^(-1)*(X-M21));
       f22=1/(2*pi)/det(S22)^0.5*exp(-0.5*(X-M22)'*S22^(-1)*(X-M22));
       f2=0.55*f21+0.45*f22;
       Smm=Smm-log(f1/f2);
       Sm(i)=Smm;
       i=i+1;
       if (Smm<a || Smm>b)
           odluka=true;
       end
   end
   brodbiter(m)=length(Sm);
   m=m+1;
   end
   brodb(z)=mean(brodbiter);
   z=z+1;
end
figure(5);
semilogx(e1niz,brodb,'b');
title('Wald-ov sekvencijalni test');
xlabel('Greska prvog tipa');
ylabel('Broj potrebnih odbiraka');

%Za gresku epsilon 2
z=1;
brodb=zeros(1,length(e2niz));
e1=10^-5;
for e2=e2niz
    
   A=(1-e1)/e2;
   B=e1/(1-e2);
   a=-log(A);
   b=-log(B);


   brodbiter=zeros(1,200);
   m=1;
   for k=1:200
   odluka=false;
   Smm=0;
   Sm=0;
   i=1;
   while(odluka==false)
       if rand(1)<0.6
          x1=T11*randn(2,1)+M11;  %dodamo srednju vrednost
       else
          x1=T12*randn(2,1)+M12;  %dodamo srednju vrednost
       end
       X=[x1(1,1);x1(2,1)];
       f11=1/(2*pi)/det(S11)^0.5*exp(-0.5*(X-M11)'*S11^(-1)*(X-M11));
       f12=1/(2*pi)/det(S12)^0.5*exp(-0.5*(X-M12)'*S12^(-1)*(X-M12));
       f1=0.6*f11+0.4*f12;
       f21=1/(2*pi)/det(S21)^0.5*exp(-0.5*(X-M21)'*S21^(-1)*(X-M21));
       f22=1/(2*pi)/det(S22)^0.5*exp(-0.5*(X-M22)'*S22^(-1)*(X-M22));
       f2=0.55*f21+0.45*f22;
       Smm=Smm-log(f1/f2);
       Sm(i)=Smm;
       i=i+1;
       if (Smm<a || Smm>b)
           odluka=true;
       end
   end
   brodbiter(m)=length(Sm);
   m=m+1;
   end
   brodb(z)=mean(brodbiter);
   z=z+1;
end
figure(6);
semilogx(e2niz,brodb,'b');
title('Wald-ov sekvencijalni test');
xlabel('Greska drugog tipa');
ylabel('Broj potrebnih odbiraka');










