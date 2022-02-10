
%% Prepoznavanje oblika 2, obrada slova pisanih na predavanjima

% Ideja je prvo malo iseci sliku tako da ostane na slici samo slovo i da
% nema slobodnog prostora, u stvari isecemo okvir koliko god je to moguce
% U pitanju su slova A,O i U

% U bazi cifara nalaze se slike pod nazivom bazaX001-bazaX120

N=120; % broj slova

%krosvalidacija-kombinovanje sta je test skup,sta trening skup,ali to
%necemo sad
NO=100; %za obucavajuci skup,posto prvih 100 sve obuhvata
        %inace bismo morali malo da izmesamo to krosvalidacijom

for i=1:N                 %radimo odmah sve zajedno
    name=['bazaA' num2str(i,'%03d')];        %03d integer na 3 mesta
    x=imread(name,'bmp');
    X1(:,i)=Obelezja(x);    %treba da uzme sliku i vrati nam broj obelezja i sama obelezja
     
    name=['bazaO' num2str(i,'%03d')];   
    x=imread(name,'bmp');
    X2(:,i)=Obelezja(x);
    
    name=['bazaU' num2str(i,'%03d')];     
    x=imread(name,'bmp');
    X3(:,i)=Obelezja(x);
    
    name=['bazaE' num2str(i,'%03d')];     
    x=imread(name,'bmp');
    X4(:,i)=Obelezja(x);
    
    name=['bazaI' num2str(i,'%03d')];     
    x=imread(name,'bmp');
    X5(:,i)=Obelezja(x);
end

O1=X1(:,1:NO); T1=X1(:,NO+1:end);   %podelimo na obucavajuci i test skup
O2=X2(:,1:NO); T2=X2(:,NO+1:end);       
O3=X3(:,1:NO); T3=X3(:,NO+1:end);      
O4=X4(:,1:NO); T4=X4(:,NO+1:end);
O5=X5(:,1:NO); T5=X5(:,NO+1:end);

plot(O1(1,:),O1(2,:),'ro',O2(1,:),O2(2,:),'bx',O3(1,:),O3(2,:),'y*',O4(1,:),O4(2,:),'g*',O5(1,:),O5(2,:),'k*');
legend('A','O','U','E','I');
%Plotuje obelezja

figure(2);
hist3(X1');  %mogli smo i ovde ispitati pretpostavke o raspodeli,ali necemo da komplikujemo
figure(3);
hist3(X2');
figure(4);
hist3(X3');
figure(5);
hist3(X4');
figure(6);
hist3(X5');

M1=mean(O1')'; S1=cov(O1');                     % ili mean(O1(:,2))
M2=mean(O2')'; S2=cov(O2'); 
M3=mean(O3')'; S3=cov(O3'); 
M4=mean(O4')'; S4=cov(O4'); 
M5=mean(O5')'; S5=cov(O5'); 

%% Testiranje klasifikatora

Mk=zeros(5,5);
for k=1:5
    if k==1
        T=T1;
    elseif k==2
        T=T2;
    elseif k==3
        T=T3;
    elseif k==4
        T=T4;
    else
        T=T5;
    end

    for i=1:(N-NO)
    x=T(:,i);
    f1=1/(2*pi*det(S1)^0.5)*exp(-0.5*(x-M1)'*S1^(-1)*(x-M1));
    f2=1/(2*pi*det(S2)^0.5)*exp(-0.5*(x-M2)'*S2^(-1)*(x-M2));
    f3=1/(2*pi*det(S3)^0.5)*exp(-0.5*(x-M3)'*S3^(-1)*(x-M3));
    f4=1/(2*pi*det(S4)^0.5)*exp(-0.5*(x-M4)'*S4^(-1)*(x-M4));
    f5=1/(2*pi*det(S5)^0.5)*exp(-0.5*(x-M5)'*S5^(-1)*(x-M5));
    m=max([f1,f2,f3,f4,f5]);   %prolazimo kroz K1,sracunamo fgv123 za svaki odbirak,gledamo maks od ova tri
                               %i gledamo jesmo li napravili gresku tako sto cemo
                               %odredjeni clan u konf. matrici povecati za 1
       if m==f1
        Mk(k,1)=Mk(k,1)+1;
        disp(['Oblik' num2str(i) 'iz klase' num2str(k) ' je prepoznat kao slovo A']);
         elseif m==f2
        Mk(k,2)=Mk(k,2)+1;
        disp(['Oblik' num2str(i) 'iz klase' num2str(k) ' je prepoznat kao slovo O']);
       elseif m==f3
        Mk(k,3)=Mk(k,3)+1;
        disp(['Oblik' num2str(i) 'iz klase' num2str(k) ' je prepoznat kao slovo U']);
       elseif m==f4
        Mk(k,4)=Mk(k,4)+1;
        disp(['Oblik' num2str(i) 'iz klase' num2str(k) ' je prepoznat kao slovo E']);
       else
        Mk(k,5)=Mk(k,5)+1;
        disp(['Oblik' num2str(i) 'iz klase' num2str(k) ' je prepoznat kao slovo I']);
       end
    end
end

Mk
greska=(sum(sum(Mk))-trace(Mk))/sum(sum(Mk));
disp(['Ukupna greska klasifikacije je: ' num2str(greska)]);

% ako uzmemo vise obelezja,bice problema sa vizuelizacijom podataka,ali
% mozemo da napravimo npr. drvo prepoznavanja,tri slicna slova sklonimo na
% jednu stranu,druga slicna na drugu i onda radimo podklasifikitore toga
% ako nam je tako lakse i bolje i mozemo prikazati preko matrice konfuzije





