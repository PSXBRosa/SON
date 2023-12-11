%%% BE Introduction au traitement de la parole
%%%
%%% Manipulation de PSOLA pour 3 signaux artificiels (voir variable 'ccc')
%%%   ccc=1 : signal sinuso�dal de p�riode T0
%%%   ccc=2 : signal p�riodique (somme de sinuso�des) de p�riode T0
%%%   ccc=3 : suite d'impulsions espac�es de T0
%%% S. Rossignol -- 20/01/14 -- 03/12/2021

clear all;
close all;


%%% construction du signal

fe=48000;  %%% fr�quence �chantillonnage
ls1=1330;  %%% longueur du silence en d�but du signal
ls2=1320;  %%% longueur du silence en fin du signal
f0=440;    %%% fr�quence fondamentale du signal vois� (avec 110 Hz, on n'entend pas � cause des HPs des ordis)
nper=100;  %%% nombre de p�riodes dans le signal vois�

lsig=round(fe/f0*nper); %%% longueur du signal utile(=vois�) (en �chantillons)

tfen=round(fe*0.005); %%% taille des fen�tres PSOLA pour les parties non vois�es (en �chantillons)


tfen1=tfen-1;         %%% utilitaire (� ne pas changer)

%%% signal -- il y a 3 signaux possibles
ccc=1;
if ccc==1 || ccc==2
  %%% --------------------------------> ccc=1 : sinus
  xsin=0.5*sin(2*pi*f0*[0:lsig]/fe);

  %%% --------------------------------> ccc=2 : somme de sinus
  if ccc==2
    npartiels=10;
    for ii=2:npartiels
      xsin=xsin+0.5/(ii^1.5)*sin(2*pi*f0*ii*[0:lsig]/fe);
    end;
    xsin=xsin/(max(abs(xsin)))*0.25; %%% on manipule l'amplitude pour qu'il n'y ait pas saturation
                                    %%% et aussi pour que le son ne soit pas aussi d�sagr�able
  end;
elseif ccc==3
  %%% --------------------------------> ccc=3 : suite d'impulsions
  xsin=zeros(1,lsig);
  aa=0;
  ii=0;
  while aa<=lsig
    aa=17+round(ii*fe/f0);
    xsin(aa)=0.125;
    ii=ii+1;
  end;
end;

xx=[zeros(1,ls1) xsin zeros(1,ls2)]; %%% on ajoute les segments de silence en d�but et fin de signal
%sound(xx,fe);


%%% marques PSOLA et voisement correspondant
%%% => bien s�r, id�alement il faudrait placer les marques
%%%    automatiquement et pas � la main comme ici

marques=[];
voise1=[];
nsil=5; %%% nombre de fen�tres PSOLA pour le silence du d�but
        %%% => si "tfen" change, adapter �a !!!!!!!
for ii=1:nsil
  marques = [marques ii*tfen]; %%% silence de d�part
  voise1  = [voise1  0];
end;
if ccc==1
  dd=round(fe/f0/4);  %%% position de la premi�re marque vois�e -- cas o� "npartiels==1";
else
  dd=17;              %%% position de la premi�re marque vois�e -- cas o� "npartiels==10"/impulsions
end;
ll=ls1+dd;
jj=1;
while ll<ls1+lsig
  marques = [marques ll];
  ll = ls1+round(dd+jj*fe/f0);
  jj=jj+1;
end;
ll=ls1+round(fe/f0/4+(jj-2)*fe/f0);
%la=361+round(20*(rand(1,1)-0.5)); %%% position de la derni�re marque vois�e, ajust�e
                                   %%% pour r�duire le SNR/RSB (essai d'optimisation automatique)
                                   %%% => si "tfen" change, adapter �a !!!!!!!
                                   %%% => si "f0" change, adapter �a !!!!!!
%la = 361;
la=109;
marques = [marques ll+la  ll+la+tfen ll+la+2*tfen ll+la+3*tfen ll+la+4*tfen ll+la+5*tfen]; %%% silence de fin
voise2  = [        1      0          0            0            0            0];

voise   = [voise1  ones(1,nper)  voise2]; %%% voisement (accompagnant chacune des marques)


figure(1);
clf;
grid on;
hold on;
for ii=1:length(marques)
  if (voise(ii)==1) %%% marque vois�e en rouge
    plot([marques(ii) marques(ii)],[-1 1],'r');
  else              %%% marque non-vois�e en noir
    plot([marques(ii) marques(ii)],[-1 1],'k');
  end;
  text(marques(ii),-0.9,num2str(ii));
end;
plot(xx);
title('signal original avec les marques PSOLA et leur voisement');
hold off;


%%% examen sur des petites parties du son de la position
%%% des marques et de leur voisement (qui ont �t� fix�s � la main)

doit=0;
if doit==1
  lll=400;  %%% largeur des petites parties du son examin�es
  for ii=0:ceil(length(xx)/lll)
    figure(5);
    clf;
    grid on;
    hold on;
    xlim([ii*lll (ii+1)*lll]);
    ylim([-1.01 1.01]);
    for ii=1:length(marques)
      plot([marques(ii) marques(ii)],[-1 1],'r');
      text(marques(ii),0,num2str(ii));
    end;
    plot(xx);
    title('bout du signal original avec les marques PSOLA et leur voisement');
    hold off;
    pause;    %%% tapez une touche pour examiner la petite partie du son suivante
  end;
end;


%%% essai de reconstruction du signal de d�part
%%% => il y a un SNR non infini parce que les marques
%%%    et les fen�tres ne sont pas parfaitement comme il faut
%%%    (notamment � cause des troncatures de fen�tres)

sigi=zeros(length(xx)+tfen+1,1); %%% tailles l�g�rement allong�es pour �tre
x1=[xx zeros(1,tfen)];         %%% s�r que �a se passe bien pour la derni�re trame
                               %%% => on coupe � la fin (voir "diff") le suppl�ment �ventuel
for ii=1:length(marques)
  if (voise(ii)==0)
    sig=x1(marques(ii)-tfen1:marques(ii)+tfen1).*hanning(2*tfen-1)';
    sigi(marques(ii)-tfen1:marques(ii)+tfen1)=sigi(marques(ii)-tfen1:marques(ii)+tfen1)+sig';
  else
    %%% dans le cas vois�, il faut adapter la taille de la fen�tre au pitch
    %% 1. on force une fen�tre sym�trique
    %tfena=round ( (marques(ii+1) - marques(ii-1))/2 );
    %sig=x1(marques(ii)-tfena:marques(ii)+tfena).*hanning(2*tfena+1)';

    %% 2. on permet une fen�tre asym�trique
    ll1=2*(marques(ii)-marques(ii-1));
    fen1=hanning(ll1); %%% premi�re partie de la fen�tre
    fen1=fen1(1:ll1/2);

    ll2=2*(marques(ii+1)-marques(ii));
    fen2=hanning(ll2); %%% deuxi�me partie de la fen�tre
    fen2=fen2(ll2/2+1:end);

    sig=x1(marques(ii)-ll1/2:marques(ii)+ll2/2).*[fen1' 1 fen2'];

    sigi(marques(ii)-ll1/2:marques(ii)+ll2/2)=sigi(marques(ii)-ll1/2:marques(ii)+ll2/2)+sig';
  end;
end;

diff = xx-sigi(1:length(xx))';
snr = 10.*log10(sum(xx.^2)/sum(diff.^2));
fprintf(1,'bruit apres reconstitution du son : snr=%f dB\n\n\n', snr);

figure(2);
clf;
plot(diff,'r');
hold on;
plot(xx,'b');
title('signal de depart et difference entre lui et le signal reconstruit par PSOLA');
hold off;
drawnow;


%%% modification du pitch -> son plus aigu

doit=1;

if doit==1
  sigi=zeros(2*length(xx),1); %%% signal reconstruit
  x1=[xx zeros(1,tfen)];      %%% 'xx' : signal original ; on ajoute un silence � la fin pour assurer
                              %%% que �a se passe bien pour la derni�re fen�tre
  tlim=0; %%% taille du signal reconstruit
  dec=0;  %%% d�calage de la position temporelle courante dans le signal reconstruit (d� au changement
          %%% de pitch �ventuellement requis)
  for ii=1:length(marques)
    if (voise(ii)==0)
      sig=x1(marques(ii)-tfen1:marques(ii)+tfen1).*hanning(2*tfen-1)';
      sigi(marques(ii)-tfen1+dec:marques(ii)+tfen1+dec)=sigi(marques(ii)-tfen1+dec:marques(ii)+tfen1+dec)+sig';
      tlim=marques(ii)+tfen1+dec;
    else
      %%% dans le cas vois�, il faut adapter la taille de la fen�tre au pitch
      %% 1. on force une fen�tre sym�trique
      %tfena=round ( (marques(ii+1) - marques(ii-1))/2 );
      %dec=dec-round(tfena*0.1); %%% on contrôle le changement de pitch
      %sig=x1(marques(ii)-tfena:marques(ii)+tfena).*hanning(2*tfena+1)';
      %sigi(marques(ii)-tfena+dec:marques(ii)+tfena+dec)=sigi(marques(ii)-tfena+dec:marques(ii)+tfena+dec)+sig';

      %% 2. on permet une fen�tre asym�trique
      ll1=2*(marques(ii)-marques(ii-1));
      fen1=hanning(ll1); %%% premi�re partie de la fen�tre
      fen1=fen1(1:ll1/2);

      ll2=2*(marques(ii+1)-marques(ii));
      fen2=hanning(ll2); %%% deuxi�me partie de la fen�tre
      fen2=fen2(ll2/2+1:end);

      chgt=0.1;                             %%% chgt*100 donne en % le changement de pitch
      chch=num2str(round(chgt*100));
      dec=dec+round((ll1/2+ll2/2)/2*-chgt); %%% on contr�le le changement de pitch
                                            %%% => "dec" indique le d�calage entre les num�ros d'�chantillon
                                            %%%   du fichier original et les num�ros d'�chantillon
                                            %%%   du fichier d'arriv�e

      sig=x1(marques(ii)-ll1/2:marques(ii)+ll2/2).*[fen1' 1 fen2'];

      sigi(marques(ii)-ll1/2+dec:marques(ii)+ll2/2+dec)=sigi(marques(ii)-ll1/2+dec:marques(ii)+ll2/2+dec)+sig';
      tlim=marques(ii)+ll2/2+dec;
    end;
  end;
  sigi=sigi(1:tlim); %%% on coupe la fin �ventuelle

  fprintf(1,'son original\n');
  sound(xx,fe);
  pause(length(xx)/fe+1); %%% pour attendre que le son soit fini d'�tre jou�
  audiowrite('fic1.wav',xx,fe);

  fprintf(1,'son au pitch augmente d environ %s pc\n',chch);
  sound(sigi,fe);
  audiowrite('fic2.wav',sigi,fe);


  %%% spectres
  tfft=2^16;              %%% taille FFT
  fff=[0:tfft-1]/tfft*fe; %%% vecteur des fr�quences

  ppt=round(length(xx)/2); %%% position de d�part (le centre du signal)
  tana=0.2; %%% taille de la fen�tre d'analyse (en seconde)
            %%% => grande ici parce que les signaux sont stationnaires (car artificiels)

  figure(3);
  clf;
  grid on;
  hold on;
  plot(xx(ppt:ppt+round(0.04*fe)));       %%% et pas "tana", pour mieux voir car "f0" est grand
  plot(sigi(ppt:ppt+round(0.04*fe)),'r'); %%% et pas "tana", pour mieux voir car "f0" est grand
  title('bouts du signal original et du signal modifie');
  hold off;

  spe1=abs(fft(xx(ppt:ppt+round(tana*fe)),tfft));
  spe2=abs(fft(sigi(ppt:ppt+round(tana*fe)),tfft));

  T0=1/f0;
  f00=1/(T0*(1-chgt));

  figure(4);
  clf;
  grid on;
  hold on;
  plot(fff,spe1);
  plot(fff,spe2,'r');
  plot([f0 f0],[0 max(spe1)],'k');
  plot([f00 f00],[0 max(spe1)],'k');
  xlim([f0-50 f0*1.1+50]);
  title('spectre des 2 signaux -- basses frequences -- avec position des pitchs en noir');
  legend('signal original','signal modifie','pitch original','pitch (theorique) modifie');
  xlabel('frequence en Hz');
  hold off;
end;

