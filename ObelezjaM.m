function P=Obelezja(x)

y=double(x); %uint8 -> double
% y=255-y;     %negativ, da nije njega,obrnuli bismo < i >=
T=0.8;         %prag za binarizaciju zbog toga da nema veze cime smo pisali
y(y<T*max(max(y)))=0;
y(y>=T*max(max(y)))=255;                     %imshow(255*(x<90)); sve sto je
                                      %manje od 90,zakucaj na 255 i pokazi mi
x=uint8(y);

%bez negativa : y(y<=T*max(max(y)))=255;y(y>T*max(max(y)))=0;

%Eliminacija okvira slike,predobrada slike
[nr,nc]=size(x);
poc=1;  %krecemo odozgo
 %mozemo da idemo redom i spustamo se i ako je bela linija,idi dole ili
 %mozemo da saberemo celu vrstu i podelimo sa duzinom i proverimo je li 255
 %i ako je manje vidimo hocemo li stati ili ne,zavisi sta nam treba
 
 
 while (poc<nr) && (sum(x(poc,1:nc))/nc<220)   %isecanje crnog okvira,moze biti greska ako ga imamo samo na kraju
     poc=poc+1;
 end
 
 kraj=nr;                                             %CRNOG
 while (kraj>1) && (sum(x(kraj,1:nc))/nc<220)     %isecanj okvira odozdo
     kraj=kraj-1;
 end
 levo=1;
 while (levo<nc) && (sum(x(1:nr,levo))/nr<220)    %isecanje okvira sleva
     levo=levo+1;
 end
 desno=nc;
 while  (desno>1) && (sum(x(1:nr,desno))/nr<220)   %isecanje okvira zdesna
     desno=desno-1;
 end
 
x=x(poc:kraj,levo:desno);
[nr,nc]=size(x);
poc=1;
levo=1;
desno=nc;
kraj=nr;

                                                 %idemo do slova
 
 
while (poc<nr) && ((sum(x(poc,1:nc))/nc)>250)... %252 robusnost na crne tacke
        || (((sum(x(poc,1:nc))/nc)<=250) && (sum(x(poc+1,1:nc))/nc)>250)
    poc=poc+1; %isecanje belog segmenta
end
while (kraj>1) && ((sum(x(kraj,1:nc))/nc)>250)... %252 robusnost na crne tacke
        || (((sum(x(kraj,1:nc))/nc)<=250) && (sum(x(kraj-1,1:nc))/nc)>250)
    kraj=kraj-1; %isecanje belog segmenta
end
while (levo<nc) && ((sum(x(1:nr,levo))/nr)>250)... %252 robusnost na crne tacke
        || (((sum(x(1:nr,levo))/nr)<=250) && (sum(x(1:nr,levo+1))/nr)>250)
    levo=levo+1; %isecanje belog segmenta
end
while (desno>1) && ((sum(x(1:nr,desno))/nr)>250)... %252 robusnost na crne tacke
        || (((sum(x(1:nr,desno))/nr)<=250) && (sum(x(1:nr,desno-1))/nr)>250)
    desno=desno-1; %isecanje belog segmenta
end
%Kod A posmatramo segment ispod slova,dosta je prazan,moze gore levo i gore
%desno,za O moze centar jer je prazan,kod U dole desno ili gore iznad slova

x=x(poc:kraj,levo:desno); 
[nr,nc]=size(x);
figure(); imshow(x);
%Obelezja
P(1,1)=mean(mean(x(round(1/3*nr):round(2/3*nr),round(1/3*nc):round(2/3*nc))));
P(2,1)=mean(mean(x(round(1/2*nr):end,round(2/3*nc):end)));


end