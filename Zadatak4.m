%% Generisanje po 500 odbiraka Gausovki raspodeljenih iz 4 klase

clear all; close all; clc;

M1=[4;4]; S1=[1.5 -0.5; -0.5 1.5];
M2=[-5;-1]; S2=[0.9 0.7; 0.7 0.9];
M3=[-4;6]; S3=[1.5 0.5; 0.5 1.5];
M4=[3;-3]; S4=[1 -0.7; -0.7 1];

[Fi1,L1]=eig(S1); T1=Fi1*L1^0.5;
[Fi2,L2]=eig(S2); T2=Fi2*L2^0.5;
[Fi3,L3]=eig(S3); T3=Fi3*L3^0.5;
[Fi4,L4]=eig(S4); T4=Fi4*L4^0.5;

N=500;

X1=T1*randn(2,N)+M1*ones(1,N);
X2=T2*randn(2,N)+M2*ones(1,N);
X3=T3*randn(2,N)+M3*ones(1,N);
X4=T4*randn(2,N)+M4*ones(1,N);

figure(1);
plot(X1(1,:),X1(2,:),'ro');
axis equal; hold on;
plot(X2(1,:),X2(2,:),'bx'); 
plot(X3(1,:),X3(2,:),'go');
plot(X4(1,:),X4(2,:),'y*'); hold off;
legend('Klasa 1', 'Klasa 2','Klasa 3','Klasa 4');
xlabel('X1');
ylabel('X2');
title('4 dvodimenzionalne klase Gausovski raspodeljene');

%%
Z=[X1 X2 X3 X4];
clear X1 X2 X3 X4
z=rand(4*N,1);
X1=[];
X2=[];
X3=[];
X4=[];
for i=1:4*N
    if z(i)<0.25
        X1=[X1 Z(:,i)];
    elseif z(i)<0.5
        X2=[X2 Z(:,i)];
    elseif z(i)<0.75
        X3=[X3 Z(:,i)];
    else
        X4=[X4 Z(:,i)];
    end
end
plot(X1(1,:),X1(2,:),'bo'); hold on;
plot(X2(1,:),X2(2,:),'rx'); 
plot(X3(1,:),X3(2,:),'yo');
plot(X4(1,:),X4(2,:),'g*'); hold off;
legend('Klasa 1', 'Klasa 2','Klasa 3','Klasa 4');
title('Prikaz inicijalno klasifikovanih odbiraka');
M1=mean(X1,2);
M2=mean(X2,2);
M3=mean(X3,2);
M4=mean(X4,2);
lmax=100; %maksimalan broj iteracija
l=1;
a=0; %ima reklasterizacije

%%
while (l<lmax) && (a==0)
    X11=[]; X21=[]; X31=[]; X41=[];
    a=1;
    for i=1:max(size(X1))
        d1=((X1(:,i)-M1)'*(X1(:,i)-M1))^0.5;
        d2=((X1(:,i)-M2)'*(X1(:,i)-M2))^0.5;
        d3=((X1(:,i)-M3)'*(X1(:,i)-M3))^0.5;
        d4=((X1(:,i)-M4)'*(X1(:,i)-M4))^0.5;
        mind=min([d1,d2,d3,d4]);
        if(mind==d1)
            X11=[X11 X1(:,i)];
        elseif (mind==d2)
            X21=[X21 X1(:,i)];
            a=0;
        elseif (mind==d3)
            X31=[X31 X1(:,i)];
            a=0;
        else
            X41=[X41 X1(:,i)];
            a=0;
        end
    end
    for i=1:max(size(X2))
        d1=((X2(:,i)-M1)'*(X2(:,i)-M1))^0.5;
        d2=((X2(:,i)-M2)'*(X2(:,i)-M2))^0.5;
        d3=((X2(:,i)-M3)'*(X2(:,i)-M3))^0.5;
        d4=((X2(:,i)-M4)'*(X2(:,i)-M4))^0.5;
        mind=min([d1,d2,d3,d4]);
        if(mind==d1)
            X11=[X11 X2(:,i)];
            a=0;
        elseif (mind==d2)
            X21=[X21 X2(:,i)];
        elseif (mind==d3)
            X31=[X31 X2(:,i)];
            a=0;
        else
            X41=[X41 X2(:,i)];
            a=0;
        end
    end
    for i=1:max(size(X3))
        d1=((X3(:,i)-M1)'*(X3(:,i)-M1))^0.5;
        d2=((X3(:,i)-M2)'*(X3(:,i)-M2))^0.5;
        d3=((X3(:,i)-M3)'*(X3(:,i)-M3))^0.5;
        d4=((X3(:,i)-M4)'*(X3(:,i)-M4))^0.5;
        mind=min([d1,d2,d3,d4]);
        if(mind==d1)
            X11=[X11 X3(:,i)];
            a=0;
        elseif (mind==d2)
            X21=[X21 X3(:,i)];
            a=0;
        elseif (mind==d3)
            X31=[X31 X3(:,i)];
        else
            X41=[X41 X3(:,i)];
            a=0;
        end
    end
    for i=1:max(size(X4))
        d1=((X4(:,i)-M1)'*(X4(:,i)-M1))^0.5;
        d2=((X4(:,i)-M2)'*(X4(:,i)-M2))^0.5;
        d3=((X4(:,i)-M3)'*(X4(:,i)-M3))^0.5;
        d4=((X4(:,i)-M4)'*(X4(:,i)-M4))^0.5;
        mind=min([d1,d2,d3,d4]);
        if(mind==d1)
            X11=[X11 X4(:,i)];
            a=0;
        elseif (mind==d2)
            X21=[X21 X4(:,i)];
            a=0;
        elseif (mind==d3)
            X31=[X31 X4(:,i)];
            a=0;
        else
            X41=[X41 X4(:,i)];           
        end
    end    
    X1=X11; X2=X21; X3=X31; X4=X41;
    clear X11 X21 X31 X41;
    M1=mean(X1,2); M2=mean(X2,2);
    M3=mean(X3,2); M4=mean(X4,2);
    l=l+1;    
end
figure(2);
plot(X1(1,:),X1(2,:),'bo'); hold on;
plot(X2(1,:),X2(2,:),'rx'); 
plot(X3(1,:),X3(2,:),'yo');
plot(X4(1,:),X4(2,:),'g*'); hold off;
legend('Klasa 1', 'Klasa 2','Klasa 3','Klasa 4');
title('Prikaz krajnje klasifikovanih odbiraka');

%% u slucaju da ne znamo broj klasa pr. 2

clear X1 X2 X3 X4
z=rand(4*N,1);
X1=[];
X2=[];
for i=1:4*N
    if z(i)<0.5
        X1=[X1 Z(:,i)];
    else
        X2=[X2 Z(:,i)];
    end
end
figure(1);
plot(X1(1,:),X1(2,:),'bo'); hold on;
plot(X2(1,:),X2(2,:),'rx');  hold off;
legend('Klasa 1', 'Klasa 2');
title('Prikaz inicijalno klasifikovanih odbiraka');
M1=mean(X1,2);
M2=mean(X2,2);
lmax=100; %maksimalan broj iteracija
l=1;
a=0; %ima reklasterizacije

%%
while (l<lmax) && (a==0)
    X11=[]; X21=[];
    a=1;
    for i=1:max(size(X1))
        d1=((X1(:,i)-M1)'*(X1(:,i)-M1))^0.5;
        d2=((X1(:,i)-M2)'*(X1(:,i)-M2))^0.5;
        mind=min([d1,d2]);
        if(mind==d1)
            X11=[X11 X1(:,i)];
        else
            X21=[X21 X1(:,i)];
            a=0;
        end
    end
    for i=1:max(size(X2))
        d1=((X2(:,i)-M1)'*(X2(:,i)-M1))^0.5;
        d2=((X2(:,i)-M2)'*(X2(:,i)-M2))^0.5;
        mind=min([d1,d2]);
        if(mind==d1)
            X11=[X11 X2(:,i)];
            a=0;
        else
            X21=[X21 X2(:,i)];        
        end
    end    
    X1=X11; X2=X21;
    clear X11 X21
    M1=mean(X1,2); M2=mean(X2,2);
    l=l+1;
end
figure(3);
plot(X1(1,:),X1(2,:),'bo'); hold on;
plot(X2(1,:),X2(2,:),'rx');  hold off;
legend('Klasa 1', 'Klasa 2');  
title('Prikaz krajnje klasifikovanih odbiraka');

%% pr. 5
clear X1 X2 X3 X4
z=rand(4*N,1);
X1=[];
X2=[];
X3=[];
X4=[];
X5=[];
for i=1:4*N
    if z(i)<1/5
        X1=[X1 Z(:,i)];
    elseif z(i)<2/5
        X2=[X2 Z(:,i)];
    elseif z(i)<3/5
        X3=[X3 Z(:,i)];
    elseif z(i)<4/5
        X4=[X4 Z(:,i)];
    else 
        X5=[X5 Z(:,i)];
    end
end
figure(1);
plot(X1(1,:),X1(2,:),'bo'); hold on;
plot(X2(1,:),X2(2,:),'rx'); 
plot(X3(1,:),X3(2,:),'yo');
plot(X4(1,:),X4(2,:),'g*'); 
plot(X5(1,:),X5(2,:),'ro'); hold off;
legend('Klasa 1', 'Klasa 2','Klasa 3','Klasa 4','Klasa 5');
title('Prikaz inicijalno klasifikovanih odbiraka');
M1=mean(X1,2);
M2=mean(X2,2);
M3=mean(X3,2);
M4=mean(X4,2);
M5=mean(X5,2);
lmax=100; %maksimalan broj iteracija
l=1;
a=0; %ima reklasterizacije

%% 
while (l<lmax) && (a==0)
    X11=[]; X21=[]; X31=[]; X41=[]; X51=[];
    a=1;
    for i=1:max(size(X1))
        d1=((X1(:,i)-M1)'*(X1(:,i)-M1))^0.5;
        d2=((X1(:,i)-M2)'*(X1(:,i)-M2))^0.5;
        d3=((X1(:,i)-M3)'*(X1(:,i)-M3))^0.5;
        d4=((X1(:,i)-M4)'*(X1(:,i)-M4))^0.5;       
        d5=((X1(:,i)-M5)'*(X1(:,i)-M5))^0.5;
        mind=min([d1,d2,d3,d4,d5]);
        if(mind==d1)
            X11=[X11 X1(:,i)];
        elseif (mind==d2)
            X21=[X21 X1(:,i)];
            a=0;
        elseif (mind==d3)
            X31=[X31 X1(:,i)];
            a=0;
        elseif (mind==d4)
            X41=[X41 X1(:,i)];
            a=0;
        else 
            X51=[X51 X1(:,i)];
            a=0;
        end
    end
    for i=1:max(size(X2))
        d1=((X2(:,i)-M1)'*(X2(:,i)-M1))^0.5;
        d2=((X2(:,i)-M2)'*(X2(:,i)-M2))^0.5;
        d3=((X2(:,i)-M3)'*(X2(:,i)-M3))^0.5;
        d4=((X2(:,i)-M4)'*(X2(:,i)-M4))^0.5;
        d5=((X2(:,i)-M5)'*(X2(:,i)-M5))^0.5;
        mind=min([d1,d2,d3,d4,d5]);
        if(mind==d1)
            X11=[X11 X2(:,i)];
            a=0;
        elseif (mind==d2)
            X21=[X21 X2(:,i)];
        elseif (mind==d3)
            X31=[X31 X2(:,i)];
            a=0;
        elseif (mind==d4)
            X41=[X41 X2(:,i)];
            a=0;
        else 
            X51=[X51 X2(:,i)];
            a=0;
        end
    end
    for i=1:max(size(X3))
        d1=((X3(:,i)-M1)'*(X3(:,i)-M1))^0.5;
        d2=((X3(:,i)-M2)'*(X3(:,i)-M2))^0.5;
        d3=((X3(:,i)-M3)'*(X3(:,i)-M3))^0.5;
        d4=((X3(:,i)-M4)'*(X3(:,i)-M4))^0.5;
        d5=((X3(:,i)-M5)'*(X3(:,i)-M5))^0.5;
        mind=min([d1,d2,d3,d4,d5]);
        if(mind==d1)
            X11=[X11 X3(:,i)];
            a=0;
        elseif (mind==d2)
            X21=[X21 X3(:,i)];
            a=0;
        elseif (mind==d3)
            X31=[X31 X3(:,i)];
        elseif (mind==d4)
            X41=[X41 X3(:,i)];
            a=0;
        else
            X51=[X51 X3(:,i)];
            a=0;
        end
    end
    for i=1:max(size(X4))
        d1=((X4(:,i)-M1)'*(X4(:,i)-M1))^0.5;
        d2=((X4(:,i)-M2)'*(X4(:,i)-M2))^0.5;
        d3=((X4(:,i)-M3)'*(X4(:,i)-M3))^0.5;
        d4=((X4(:,i)-M4)'*(X4(:,i)-M4))^0.5;
        d5=((X4(:,i)-M5)'*(X4(:,i)-M5))^0.5;
        mind=min([d1,d2,d3,d4,d5]);
        if(mind==d1)
            X11=[X11 X4(:,i)];
            a=0;
        elseif (mind==d2)
            X21=[X21 X4(:,i)];
            a=0;
        elseif (mind==d3)
            X31=[X31 X4(:,i)];
            a=0;
        elseif (mind==d4)
            X41=[X41 X4(:,i)];  
        else
            X51=[X51 X4(:,i)];
            a=0;
        end
    end    
    for i=1:max(size(X5))
        d1=((X5(:,i)-M1)'*(X5(:,i)-M1))^0.5;
        d2=((X5(:,i)-M2)'*(X5(:,i)-M2))^0.5;
        d3=((X5(:,i)-M3)'*(X5(:,i)-M3))^0.5;
        d4=((X5(:,i)-M4)'*(X5(:,i)-M4))^0.5;
        d5=((X5(:,i)-M5)'*(X5(:,i)-M5))^0.5;
        mind=min([d1,d2,d3,d4,d5]);
        if(mind==d1)
            X11=[X11 X5(:,i)];
            a=0;
        elseif (mind==d2)
            X21=[X21 X5(:,i)];
            a=0;
        elseif (mind==d3)
            X31=[X31 X5(:,i)];
            a=0;
        elseif (mind==d4)
            X41=[X41 X5(:,i)];  
            a=0;
        else
            X51=[X51 X5(:,i)];           
        end
    end    
    X1=X11; X2=X21; X3=X31; X4=X41; X5=X51;
    clear X11 X21 X31 X41 X51;
    M1=mean(X1,2); M2=mean(X2,2);
    M3=mean(X3,2); M4=mean(X4,2);
    M5=mean(X5,2);
    l=l+1;    
end
figure(4);
plot(X1(1,:),X1(2,:),'bo'); hold on;
plot(X2(1,:),X2(2,:),'rx'); 
plot(X3(1,:),X3(2,:),'yo');
plot(X4(1,:),X4(2,:),'g*');
plot(X5(1,:),X5(2,:),'ko');hold off;
legend('Klasa 1', 'Klasa 2','Klasa 3','Klasa 4','Klasa 5');
title('Prikaz krajnje klasifikovanih odbiraka');
    
%% 2. zadatak

clear all; close all; clc;

%PRVA KLASA

N=500;
M1=[6;6];
R1=2;
Tetax=rand(1,N)*2*pi;
Rx=rand(1,N)*R1;
X1=[Rx.*cos(Tetax);Rx.*sin(Tetax)]+M1*ones(1,N);

%DRUGA KLASA

M2=[6;6];
R2=4.5; d=2;
Tetay=rand(1,N)*2*pi;
Ry=rand(1,N)*d+R2;
X2=[Ry.*cos(Tetay);Ry.*sin(Tetay)]+M2*ones(1,N);
figure(5); plot(X1(1,:),X1(2,:),'bo');
hold on;
plot(X2(1,:),X2(2,:),'rv'); hold on;
legend('Klasa 1','Klasa 2');

%%
Z=[X1 X2];
clear X1 X2 
z=rand(2*N,1);
X1=[];
X2=[];
for i=1:2*N
    if z(i)<0.5
        X1=[X1 Z(:,i)];
    else
        X2=[X2 Z(:,i)];
    end
end
figure(6);
plot(X1(1,:),X1(2,:),'bo'); hold on;
plot(X2(1,:),X2(2,:),'rx');  hold off;
legend('Klasa 1', 'Klasa 2');
title('Prikaz inicijalno klasifikovanih odbiraka');

%%
l=1;
lmax=100;
a=0;

while (l<lmax) && (a==0)
   X11=[]; X21=[]; 
   a=1;
   P1=size(X1,2)/(2*N);
   P2=size(X2,2)/(2*N);
   M1=mean(X1,2);
   M2=mean(X2,2);
   S1=cov(X1');
   S2=cov(X2');
   
   for i=1:max(size(X1))
       d1=0.5*(X1(:,i)-M1)'*inv(S1)*(X1(:,i)-M1)+0.5*log(det(S1)/P1);
       d2=0.5*(X1(:,i)-M2)'*inv(S2)*(X1(:,i)-M2)+0.5*log(det(S2)/P2);
       if d1>d2
           X21=[X21 X1(:,i)];
           a=0;
       else
           X11=[X11 X1(:,i)];
       end
   end
   
   for i=1:max(size(X2))
       d1=0.5*(X2(:,i)-M1)'*inv(S1)*(X2(:,i)-M1)+0.5*log(det(S1)/P1);
       d2=0.5*(X2(:,i)-M2)'*inv(S2)*(X2(:,i)-M2)+0.5*log(det(S2)/P2);
       if d1>d2
           X21=[X21 X2(:,i)];
       else
           X11=[X11 X2(:,i)];
           a=0;
       end
   end   
   l=l+1;
   X1=X11; X2=X21;
   clear X11 X21;
   figure(7);
   plot(X1(1,:),X1(2,:),'bo'); hold on;
   plot(X2(1,:),X2(2,:),'rx'); hold off;
   legend('Klasa 1', 'Klasa 2');
   title('Prikaz krajnje klasifikovanih odbiraka');
end

%%
Z=[X1 X2];
clear X1 X2 
z=rand(2*N,1);
X1=[];
X2=[];
X3=[];
for i=1:2*N
    if z(i)<1/3
        X1=[X1 Z(:,i)];
    elseif z(i)<2/3
        X2=[X2 Z(:,i)];
    else
        X3=[X3 Z(:,i)];
    end
end
figure(6);
plot(X1(1,:),X1(2,:),'bo'); hold on;
plot(X2(1,:),X2(2,:),'rx');  
plot(X3(1,:),X3(2,:),'y*'); hold off;
legend('Klasa 1', 'Klasa 2','Klasa 3');
title('Prikaz inicijalno klasifikovanih odbiraka');

%%
l=1;
lmax=100;
a=0;

while (l<lmax) && (a==0)
   X11=[]; X21=[]; X31=[];
   a=1;
   P1=size(X1,2)/(2*N);
   P2=size(X2,2)/(2*N);
   P3=size(X3,2)/(2*N);
   M1=mean(X1,2);
   M2=mean(X2,2);
   M3=mean(X3,2);
   S1=cov(X1');
   S2=cov(X2');
   S3=cov(X3');
   
   for i=1:max(size(X1))
       d1=0.5*(X1(:,i)-M1)'*inv(S1)*(X1(:,i)-M1)+0.5*log(det(S1)/P1);
       d2=0.5*(X1(:,i)-M2)'*inv(S2)*(X1(:,i)-M2)+0.5*log(det(S2)/P2);
       d3=0.5*(X1(:,i)-M3)'*inv(S3)*(X1(:,i)-M3)+0.5*log(det(S3)/P3);
       mind=min([d1,d2,d3]);
       if mind==d2
           X21=[X21 X1(:,i)];
           a=0;
       elseif mind==d1
           X11=[X11 X1(:,i)];
       else
           X31=[X31 X1(:,i)];
           a=0;
       end
   end
   
   for i=1:max(size(X2))
       d1=0.5*(X2(:,i)-M1)'*inv(S1)*(X2(:,i)-M1)+0.5*log(det(S1)/P1);
       d2=0.5*(X2(:,i)-M2)'*inv(S2)*(X2(:,i)-M2)+0.5*log(det(S2)/P2);
       d3=0.5*(X2(:,i)-M3)'*inv(S3)*(X2(:,i)-M3)+0.5*log(det(S3)/P3);
       mind=min([d1,d2,d3]);
       if mind==d2
           X21=[X21 X2(:,i)];           
       elseif mind==d1
           X11=[X11 X2(:,i)];
           a=0;
       else
           X31=[X31 X2(:,i)];
           a=0;
       end
   end 
   for i=1:max(size(X3))
       d1=0.5*(X3(:,i)-M1)'*inv(S1)*(X3(:,i)-M1)+0.5*log(det(S1)/P1);
       d2=0.5*(X3(:,i)-M2)'*inv(S2)*(X3(:,i)-M2)+0.5*log(det(S2)/P2);
       d3=0.5*(X3(:,i)-M3)'*inv(S3)*(X3(:,i)-M3)+0.5*log(det(S3)/P3);
       mind=min([d1,d2,d3]);
       if mind==d2
           X21=[X21 X3(:,i)];
           a=0;
       elseif mind==d1
           X11=[X11 X3(:,i)];
           a=0;
       else
           X31=[X31 X3(:,i)];
       end
   end
   l=l+1;
   X1=X11; X2=X21; X3=X31;
   clear X11 X21 X31;
   figure(7);
   plot(X1(1,:),X1(2,:),'bo'); hold on;
   plot(X2(1,:),X2(2,:),'rx'); 
   plot(X3(1,:),X3(2,:),'y*'); hold off;
   legend('Klasa 1', 'Klasa 2','Klasa 3');
   title('Prikaz krajnje klasifikovanih odbiraka');
end

