function P=Obelezja(x)
y=double(x); 
y=y(10:(end-6),1:(end-4));
T=0.8; %prag za binarizaciju
y(y<T*max(max(y)))=0;
y(y>=T*max(max(y)))=255;
x=uint8(y);

%predobrada slike
[nr,nc]=size(x);
poc=1; %krecemo odozgo
while (poc<nr) && ((sum(x(poc,1:nc))/nc)<230) %isecanje okvira
    poc=poc+1;
end
kraj=nr;
while(kraj>1) && (sum(x(kraj,1:nc))/nc<230)
    kraj=kraj-1;
end
levo=1;
while(levo<nc) && (sum(x(1:nr,levo))/nr<230)
    levo=levo+1;
end
desno=nc;
while (desno>1) && (sum(x(1:nr,desno))/nr<230)
    desno=desno-1;
end
x=x(poc:kraj,levo:desno);
[nr,nc]=size(x);
poc=1;
levo=1;
desno=nc;
kraj=nr;
while (poc<nr) && ((sum(x(poc,1:nc))/nc)>251)... %252 robusnost na crne tacke
        || (((sum(x(poc,1:nc))/nc)<=251) && (sum(x(poc+1,1:nc))/nc)>251)
    poc=poc+1; %isecanje belog segmenta
end
while (kraj>1) && ((sum(x(kraj,1:nc))/nc)>251)... %252 robusnost na crne tacke
        || (((sum(x(kraj,1:nc))/nc)<=251) && (sum(x(kraj-1,1:nc))/nc)>251)
    kraj=kraj-1; %isecanje belog segmenta
end
while (levo<nc) && ((sum(x(1:nr,levo))/nr)>251)... %252 robusnost na crne tacke
        || (((sum(x(1:nr,levo))/nr)<=251) && (sum(x(1:nr,levo+1))/nr)>251)
    levo=levo+1; %isecanje belog segmenta
end
while (desno>1) && ((sum(x(1:nr,desno))/nr)>251)... %252 robusnost na crne tacke
        || (((sum(x(1:nr,desno))/nr)<=251) && (sum(x(1:nr,desno-1))/nr)>251)
    desno=desno-1; %isecanje belog segmenta
end

x=x(poc:kraj,levo:desno);
[nr,nc]=size(x);
P(1,1)=mean(mean(x(1:round(1/5*nr),round(1/6*nc):round(5/6*nc))));
P(2,1)=mean(mean(x(round(4/5*nr):end,round(1/6*nc):round(5/6*nc))));
P(3,1)=mean(mean(x(round(1/3*nr):round(2/3*nr),round(1/6*nc):round(5/6*nc))));
P(4,1)=mean(mean(x(round(3/5*nr):end,1:round(2/6*nc))));
end                                                      

