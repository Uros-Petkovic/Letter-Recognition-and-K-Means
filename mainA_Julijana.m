clear all
close all
clc 

N=500;
% N=10000;

M1=[2 2]';
M=[0 0]';
M21=[8 6]';
M22=[7 -6]';

Sig1=[4 -1;-1 2];
[Q, Lam]=eig(Sig1);
Q
Lam
LT1=sqrt(Lam)*Q';

Sig21=[1 -3;-3 10];
[Q, Lam]=eig(Sig21);
Q
Lam
LT21=sqrt(Lam)*Q';

Sig22=[2 -1;-1 2];
[Q, Lam]=eig(Sig22);
Q
Lam
LT22=sqrt(Lam)*Q';

% ODBIRCI IZ PRVE KLASE:

x1=randn(N,2)*LT1+[M1(1)*ones(N,1), M1(2)*ones(N,1)];
p21=0.4;
p22=1-p21;

% ODBIRCI IZ DRUGE KLASE:

x2=zeros(N,2);

for i=1:N
    p=rand(1,1);
    if p<p21
        x2(i,:)=randn(1,2)*LT21+[M21(1), M21(2)];
    else 
        x2(i,:)=randn(1,2)*LT22+[M22(1), M22(2)];
    end
end

% CRTANJE ODBIRAKA DVEJU KLASA

plot(x1(:,1),x1(:,2),'b.');
hold on
plot(x2(:,1),x2(:,2),'ro');
axis([-10 20 -15 20]);
legend('w1','w2'); 
title('Prikaz odbiraka klasa w1 i w2');

% Srednja vrednost
for i=1:N
    M(1)=x1(i,1)+M(1);
    M(2)=x1(i,2)+M(2);
end
M=M/N;

% SIGMA
S=zeros(2,2);
for i=1:N
    S(1,1)=S(1,1)+(x1(i,1)-M(1))*(x1(i,1)-M(1));
    S(2,2)=S(2,2)+(x1(i,2)-M(2))*(x1(i,2)-M(2));
    S(2,1)=S(2,1)+(x1(i,2)-M(2))*(x1(i,1)-M(1));
    S(1,2)=S(1,2)+(x1(i,1)-M(1))*(x1(i,2)-M(2));
end
S=S/N
 

d2array=[1 4 9 16];



% for d2=d2array
% NumPoints=1000;
% teta=linspace(0,2*pi,NumPoints);
% dx1=sqrt(d2)*sin(teta)';
% dy1=sqrt(d2)*cos(teta)';
% Y1=[dx1, dy1]*LT1+ones(NumPoints,1)*M1';
% Y21=[dx1, dy1]*LT21+ones(NumPoints,1)*M21';
% Y22=[dx1, dy1]*LT22+ones(NumPoints,1)*M22';
% plot(Y1(:,1),Y1(:,2),'k-');
% plot(Y21(:,1),Y21(:,2),'k-');
% plot(Y22(:,1),Y22(:,2),'k-');
% end
% 
x=-10:0.5:20;
y=-15:0.5:20;
T=0;
bw=zeros(length(x),length(y),'uint8');
% 
for i=1:length(x)
    for j=1:length(y)
    f1(j,i)=pfunc([x(i);y(j)],Sig1,M1);
    f21(j,i)=pfunc([x(i);y(j)],Sig21,M21);
    f22(j,i)=pfunc([x(i);y(j)],Sig22,M22);
    f2(j,i)=p21*f21(j,i)+p22*f22(j,i);
    h(j,i)=log((f2(j,i))/(f1(j,i)));
    if(h(j,i)>T)
        bw(j,i)=255;
    end
    end
end
figure
mesh(x,y,f1);
hold on 
mesh(x,y,f21);
hold on
mesh(x,y,f22);


figure(1);
hold all
fmax=max(max(f1));
for k=1:4
    hold all
    contour(x,y,f1,[fmax*exp(-k^2/2), fmax*exp(-k^2/2)],'m');
end
fmax=max(max(f2));
for k=1:4
    hold all
    contour(x,y,f2,[fmax*exp(-k^2/2), fmax*exp(-k^2/2)],'g');
end
title('d\^2 krive');
hold off

x1g=0;
x2g=0;

for i=1:length(x1)
    hx=log((pfunc([x1(i,1);x1(i,2)],Sig22,M22)+pfunc([x1(i,1);x1(i,2)],Sig21,M21))/(pfunc([x1(i,1);x1(i,2)],Sig1,M1)));
    if(hx>=T)
        x1g=x1g+1;
        x1(i,:)
    end
end
for i=1:length(x2)
    hx=log((pfunc([x2(i,1);x2(i,2)],Sig22,M22)+pfunc([x2(i,1);x2(i,2)],Sig21,M21))/(pfunc([x2(i,1);x2(i,2)],Sig1,M1)));
    if(hx<T)
        x2g=x2g+1;
        x2(i,:)
    end
end

e1=x1g/N
e2=x2g/N

 %Iscrtavanje klasifikacione linije
 figure(1); hold on;
 contour(x, y, h,[0 0], 'c');
 title('Bayesov klasifikator');

%% ZADATAK POD C

T=0;

%  x1g=0;
% x2g=0;


Nmi=300;

    
%     RACUNANJE PRAGA
brT=1;
for mi=linspace(-50,50,Nmi)
        
%         T(brT)=-log(mi);
        Tmat(brT)=mi;
        
         x1g=0;
         x2g=0;
        
        for i=1:length(x1)
            hx=log((pfunc([x1(i,1);x1(i,2)],Sig22,M22)+pfunc([x1(i,1);x1(i,2)],Sig21,M21))/(pfunc([x1(i,1);x1(i,2)],Sig1,M1)));
            if(hx>=Tmat(brT))
                x1g=x1g+1;
%                 x1(i,:)
            end
        end
        for i=1:length(x2)
            hx=log((pfunc([x2(i,1);x2(i,2)],Sig22,M22)+pfunc([x2(i,1);x2(i,2)],Sig21,M21))/(pfunc([x2(i,1);x2(i,2)],Sig1,M1)));
            if(hx<Tmat(brT))
                x2g=x2g+1;
%                 x2(i,:)
            end
        end

        e1mat(brT)=x1g/N;
        e2mat(brT)=x2g/N;
        
        brT=brT+1;
end

e2mat=fliplr(e2mat);
e1mat=fliplr(e1mat);
Tmat=fliplr(Tmat);

e2niz=linspace(0,0.5,4);

for e2k=1:length(e2niz)

%     e2pam=e2;
%     e1pam=e1;

    [a,b]=min(abs(e2mat-e2niz(e2k)));
    e1=e1mat(b)
    e2=e2mat(b)
    T=Tmat(b);
    Tpam(e2k)=T;
    
%     e2=1;
%     
%     for k=1:length(e2mat)
%         if(abs(e2mat(k)-e2niz(e2k))<=e2)
%             e2=e2mat(k);
%             e1=e1mat(k);
%             T=Tmat(k);
%         end
%     end
    
%      CRTANJE KRIVIH

        x=-10:0.5:20;
        y=-15:0.5:20;
        bw=zeros(length(x),length(y),'uint8');

        for i=1:length(x)
            for j=1:length(y)
            f1(j,i)=pfunc([x(i);y(j)],Sig1,M1);
            f21(j,i)=pfunc([x(i);y(j)],Sig21,M21);
            f22(j,i)=pfunc([x(i);y(j)],Sig22,M22);
            h(j,i)=log((f22(j,i)+f21(j,i))/(f1(j,i)));
            if(h(j,i)>T)
                bw(j,i)=255;
            end
            end
        end


 %Iscrtavanje klasifikacione linije
 figure(1); hold on;
 contour(x, y, h,[Tpam(e2k) Tpam(e2k)], 'o');
 title('Neyman-Perason-ov test za razlicite Epsilon0');
 hold on

end


%% PRVI POD D Waldov sekv. test

    %Generisanje iz prve klase
    
breps=1;
eks=-20:1:-1;
epsniz=10.^eks;
% epsniz=[1e-20 1e-15 1e-10 1e-5 1e-4 1e-3 1e-2 1e-1];
for eps=epsniz
    
e1=eps;
e2=eps;

A=(1-e1)/e2;
B=e1/(1-e2);
a=-log(A);
b=-log(B);


    
    for k=1:400
        odluka=false;
        br(k)=0;
        s=0;
        while(odluka==false)
            x1=randn(1,2)*LT1+[M1(1), M1(2)];
            p1=pfunc([x1(1);x1(2)],Sig1,M1);
            p21=pfunc([x1(1);x1(2)],Sig21,M21);
            p22=pfunc([x1(1);x1(2)],Sig22,M22);
            p2=0.6*p21+0.4*p22;
            s=s-log(p1/p2);
            if (s<a || s>b)
                odluka=true;
            end
            br(k)=br(k)+1;
        end
    end
    brmer(breps)=mean(br);
    breps=breps+1;
end

figure
plot(log10(epsniz),brmer);
title('Wald-ov sekvencijalni test za odbirke prve klase');
xlabel('greska');
ylabel('broj potrebnih odbiraka');

    %Generisanje iz druge klase

breps=1;
eks=-20:1:-1;
epsniz=10.^eks;
% epsniz=[1e-20 1e-15 1e-10 1e-5 1e-4 1e-3 1e-2 1e-1];
for eps=epsniz
    
e1=eps;
e2=eps;

A=(1-e1)/e2;
B=e1/(1-e2);
a=-log(A);
b=-log(B);
    
    for k=1:400
        odluka=false;
        br(k)=0;
        s=0;
        while(odluka==false)
            
            p=rand(1,1);
            if p<p21
                x1=randn(1,2)*LT21+[M21(1), M21(2)];
            else 
                x1=randn(1,2)*LT22+[M22(1), M22(2)];
            end
            
            p1=pfunc([x1(1);x1(2)],Sig1,M1);
            p21=pfunc([x1(1);x1(2)],Sig21,M21);
            p22=pfunc([x1(1);x1(2)],Sig22,M22);
            p2=0.6*p21+0.4*p22;
            s=s-log(p1/p2);
            if (s<a || s>b)
                odluka=true;
            end
            br(k)=br(k)+1;
        end
    end
    brmer(breps)=mean(br);
    breps=breps+1;
end

figure
plot(log10(epsniz),brmer);
title('Wald-ov sekvencijalni test za odbirke druge klase');
xlabel('greska');
ylabel('broj potrebnih odbiraka');