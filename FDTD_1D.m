%% Implementation of FDTD method

close all; clear all;
%% Parameters of the propagation

% Param�tres physiques
L = 40e-6; % Dimension de propagation
l = 800e-9; % Longueur d'onde
c = 3e8; % Vitesse de la lumi�re dans le vide
eps0 = 8.854e-12; % Permittivit� �lectrique du vide
espR = 1; % Permittivit� relative du milieu
mu0 = 4*pi*10^(-7); % Perm�abilit� magn�tique du vide
muR = 1; % Perm�abilit� relative du milieu
eta0 = sqrt(mu0/eps0); % Imp�dance du vide
tau = 10e-15;

% Param�tres num�riques
dz = l/20; % Dimension unitaire de la grille
dt = dz/c; % Dimension temporelle de la simulation
N = 1001; % R�solution spatiale de la simulation
Ni = 2000;
z = linspace(0,L,N);
t = 0:dt:Ni*dt;

% Propri�t�s d'espace
EpsR = ones(1,N); % Matrice de permittivit� �lectrique
MuR = ones(1,N); % Matrice de perm�abilit� magn�tique
sigma = zeros(1,N); % Conductivit�

% Profil pour les variables d'espaces
profile = 5;
for i=1:N
   w1 = 0.5*1.e-6;
   w2 = 1.5*1.e-6;
   if (profile==1) %dielectric window
      if (abs(z(i)-L/2)<0.5e-6) EpsR(i)=1; end
   end
   if (profile==2)%dielectric window with smooth transition
   	if (abs(z(i)-L/2)<1.5*1.e-6) EpsR(i)=1+3*(1+cos(pi*(abs(z(i)-L/2)-w1)/(w2-w1)))/2; end
      if (abs(z(i)-L/2)<0.5*1.e-6) EpsR(i)=4; end
   end
   if (profile==3) %dielectric discontinuity
      if (z(i)>L/2) EpsR(i) = 9; end
   end
  	if (profile==4) %dielectric disontinuity with 1/4-wave matching layer
      if (z(i)>(L/2-0.1443)) EpsR(i) = 3; end
      if (z(i)>L/2) EpsR(i) = 9; end
   end
   if (profile==5) %conducting half space
      if (z(i)>L/2) sigma(i) = 5e3; end
   end
   if (profile==6) %sinusoidal dielectric
      EpsR(i) = 1+sin(2*pi*z(i)/L)^2;
   end
   if (profile==7) %sinusoidal dimagnetic
      MuR(i) = 1+sin(2*pi*z(i)/L)^2;
   end

end
eta = eta0.*sqrt(MuR./EpsR); % Imp�dance du milieu
%% Field definition

% Initialisation des champs
E0 = 2; % Amplitude du champ (prise en compte des solutions de Maxwell)
E = zeros(size(z,2),size(t,2));
H = zeros(size(z,2),size(t,2));

% Terme de source
Es = @(ti,t0) E0*exp(-((ti-t0)/(tau)).^2).*sin(2*pi.*ti*c/l);

for ti = 2:numel(t)-1
    H(1, :) = H(2, :); % Conditions aux limites (Absorption sur les bords � gauche)
    for zi = 2:numel(z)-1
        % Mise � jour des param�tres d'espace
        ae = (dt)/((dz*eps0*EpsR(zi))*(1+dt*sigma(zi)/(2*eps0*EpsR(zi))));
        am = (dt)/(dz*mu0*muR);
        as = (1 - dt*sigma(zi)/(2*eps0*EpsR(zi)))/(1 + dt*sigma(zi)/(2*eps0*EpsR(zi)));

        % Mise � jour des champs
        if zi == 2
            H(zi,ti) = H(zi,ti-1) - am*(E(zi+1,ti-1) - E(zi,ti-1));
            E(zi,ti) = as*E(zi,ti-1) - ae*(H(zi,ti) - H(zi-1,ti)) + Es(t(ti),3*tau);
        else
            H(zi,ti) = H(zi,ti-1) - am*(E(zi+1,ti-1) - E(zi,ti-1));
            E(zi,ti) = as*E(zi,ti-1) - ae*(H(zi,ti) - H(zi-1,ti));
        end
    end
    E(numel(z), :) = E(numel(z)-1, :); % Conditions aux limites (Absorption sur les bords � droite)
end
%% Plotting of the fields
i = 1;
figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])% Pour mettre la figure en plein �cran
set(gcf,'doublebuffer','on'); %Pour des figures plus adoucies
for ti = 1:numel(t)
    plot(E(:,ti)); hold on;
    plot(eta*H(:,ti)); hold off;
    ylim([-4 4])
    axis off
    print(sprintf('%d.png',i),'-dpng');
    i = i + 1;
end
close all;
