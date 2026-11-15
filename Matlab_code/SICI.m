
So_start=1370; % upphafsgildi Sólfastans. Stuðlar ebm eru stilltir þannig
               % að þetta gildi valdi ekki ísi hulinni jörð.

% Kanna efrihluta "íslausu" greinarinnar.
SoUpp=linspace(So_start,So_start*1.30,40); % S gildi til að nota við
                                          % könnun efri hluta
                                          % "íslausu" greinarinnar
Tupp=zeros(size(SoUpp)); % vigur til að geyma meðalhitagildi
					  
ebm % búa til T til að starta með
 
disp('Kanna efrihluta íslausu greinarinnar')

for i=1:length(SoUpp)
  Talb=IterTalb(T,SoUpp(i));
  Tupp(i)=mean(Talb);
end

% Kanna neðrihluta íslausu greinarinnar og minnka sólfastann þangað
% til jörðin verður ísihulinn (veit að ég þarf að fara niður fyrir
% 1277). Ef ég fer 15% niður verður það yfrið nóg.

SoDown=linspace(So_start,So_start*(1-0.15),71);
Tdown=zeros(size(SoDown));

disp('Kanna neðri hluta íslausu greinarinnar niður á isihuldu greinina')

for i=1:length(SoDown)
  Talb=IterTalb(T,SoDown(i));
  Tdown(i)=mean(Talb);
end



% Kanna ísilögðu greinina upp frá 1250W
% til jörðin verður ísihulinn (veit að ég þarf að fara niður fyrir
% 1277). Ef ég fer 15% niður verður það yfrið nóg.

So_cold=1250;
Tcold=IterTalb(T,So_cold);

SoUpp2=linspace(So_cold,So_start*1.30,150);
Tupp2=zeros(size(SoUpp2));

disp('Kanna ísihuldu greinina upp á þá íslausu')
for i=1:length(SoUpp2)
  Talb=IterTalb(Tcold,SoUpp2(i));
  Tupp2(i)=mean(Talb);
end

figure
plot(SoUpp,Tupp)
hold
plot(SoDown,Tdown)
plot(SoUpp2,Tupp2)
