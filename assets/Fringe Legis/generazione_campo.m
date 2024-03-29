%elettrodo di estrazione si trova a 0 mm
clearvars
close all
load('Bfield')
rangex=min(xm):0.001:max(xm);%range da usare per l'interpolazione che va dal minimo al massimo dei punti x con step di 1 mm
rangex=rangex';
B_field=interp1(xm,BT,rangex,'phchp');%interpolazione del campo magnetico
figure
plot(rangex,B_field)
gradB=gradient(B_field);%calcolo il gradiente per la componente radiale
figure
plot(rangex,gradB)
r=0:0.001:0.11;%range di punti lungo r a step di 1 mm (deve essere impostato in base alle dimensioni trasverse della simulazione)
matrice_campo=[];
%calcol ola componente radiale del campo magnetico
for i=1:length(rangex)
    for j=1:length(r)
        Br=-0.5*r(j)*gradB(i);
        matrice_campo=[matrice_campo;rangex(i),r(j),B_field(i),Br];%salvo la matrice per IbSimu con x,r,Bx,Br
    end
end
dlmwrite('matrice_campo.txt', matrice_campo, 'delimiter', '\t')%file di testo della matrice per IbSimu con x,r,Bx,Br gi� ordinato in ordine decrescente
