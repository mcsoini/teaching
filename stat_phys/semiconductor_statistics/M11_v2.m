clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%       Distribution functions and densities of states:   %%%%%%%%%%%%
%%%%%%        Semiconductor statistics                   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1

E=-0.5:0.001:1.5    %eV

k=8.62*10^-5    %eV/K
m=5.69*10^-16   %eV*s^2/cm^2
h=4.14*10^-15   %eV*s

T=300   %K

EV=0    %eV
EC=1.12 %eV
mu=EC/2 %eV

%Fermi-Dirac distribution:
FD=1./(exp((E-mu)./(T*k))+1);
%Conduction band density of states
g(1,:)=pi*(8*m)^1.5/(2*h^3)*real(sqrt(E-EC));
%Valence band density of states
g(2,:)=pi*(8*m)^1.5/(2*h^3)*real(sqrt(-(E-EV)));

%%%%%%%%%%%%%%%%%%%
figure(1)
[AX,H1,H2] = plotyy(E,FD,E,g)
title({'#1: Densities of states and FD distribution'}, ...
    'fontsize', 14, 'fontweight','b')
set(get(AX(1),'Ylabel'),'String','Fermi-Dirac distribution','fontsize', 16,'fontweight','b')
set(get(AX(2),'Ylabel'),'String','Densities of states [1/(eV*cm^3)]','fontsize', 16,'fontweight','b')
legend('FD distribution','CB density of states','VB Density of states', 'location','best')
xlabel('Energy [eV]','fontsize', 16,'fontweight','b')
yLimit = get(gca,'YLim');
line([EC EC], yLimit, 'color', 'black', 'lineStyle', '-')
line([EV EV], yLimit, 'color', 'black', 'lineStyle', '-')
line([EC/2 EC/2], yLimit, 'color', 'black', 'lineStyle', '--')
%%%%%%%%%%%%%%%%%%%
%% 3

%Density of electrons:
nperE=g(1,:).*FD

figure(2)
plot(E,nperE)
xlim([1,1.5])
title({'#3: Electron density'}, ...
    'fontsize', 14, 'fontweight','b')
xlabel('Energy [eV]','fontsize', 16,'fontweight','b')
ylabel('Electron density [1/(eV*cm^3)]','fontsize', 16,'fontweight','b')
set(gca,'fontsize',16,'fontweight','normal')
yLimit = get(gca,'YLim');
line([EC EC], yLimit, 'color', 'black', 'lineStyle', '-')

%% 4

%Electron concentration
trapz(E,nperE)

%Result: 9.8506e+09 1/cm^3

%% 5

T=100:100:1500
[TM,EM]=meshgrid(T,E);

g=zeros(size(TM));
FD=zeros(size(TM));

FD=1./(exp((EM-mu)./(TM*k))+1);
gT=pi*(8*m)^1.5/(2*h^3)*real(sqrt(EM-EC));
%%%%%%%%%%%%%%%%%%%
figure(3)
plot(E,FD(:,1),E,FD(:,end))
title({'#5.1: FD distribution for different T'}, ...
    'fontsize', 14, 'fontweight','b')
yLimit = get(gca,'YLim');
line([EC EC], yLimit, 'color', 'black', 'lineStyle', '-')
line([EV EV], yLimit, 'color', 'black', 'lineStyle', '-')
xlim([-0.1,1.2])
legend('T=100K','T=1500K', 'location','best')
xlabel('Energy [eV]','fontsize', 16,'fontweight','b')
ylabel('Fermi-Dirac distribution','fontsize', 16,'fontweight','b')
set(gca,'fontsize',16,'fontweight','normal')
%%%%%%%%%%%%%%%%%%%


nfunceT=gT.*FD;
%%%%%%%%%%%%%%%%%%%
figure(4)
semilogy(T,trapz(E,nfunceT),'-or')
title({'#5.2: Temperature dependent electron density'}, ...
    'fontsize', 14, 'fontweight','b')
set(gca,'fontsize',16,'fontweight','normal')
xlabel('Temperature [K]','fontsize', 16,'fontweight','b')
ylabel('Electron density [1/cm^3]','fontsize', 16,'fontweight','b')
%%%%%%%%%%%%%%%%%%%


%% 6
mu=-0:0.01:1.5;
ED=0.044; %Phosphorus
ND=10^18;

T0=300

%calculate n:
[EM,muM]=meshgrid(E,mu);

noccT=1./(exp((EM-muM)./(T0*k))+1);

gT=pi*(8*m)^1.5/(2*h^3)*real(sqrt(EM-EC));
nfunceT=gT.*noccT;
n=trapz(E,nfunceT,2)

%calculate p:
noccT=1./(exp(-(EM-muM)./(T0*k))+1);
gT=pi*(8*m)^1.5./(2*h^3)*real(sqrt(EV-EM));
nfunceT=gT.*noccT;
p=trapz(E,nfunceT,2)

NDp=ND*(1+2*exp(-(EC-ED-mu)./(k*T0))).^-1

%%%%%%%%%%%%%%%%%%%
figure(5)
semilogy(mu,p.'+NDp, mu,n)

yLimit = get(gca,'YLim');
title({'#6: Charge concentrations as a function of \mu'}, ...
    'fontsize', 14, 'fontweight','b')
line([EC EC],yLimit, 'color', 'r', 'lineStyle', '-')
line([EV EV],yLimit, 'color', 'r', 'lineStyle', '-')
line([EC/2 EC/2],yLimit, 'color', 'r', 'lineStyle', '--')
line([EC-ED EC-ED],yLimit, 'color', 'r', 'lineStyle', '--')
legend('Positive charges','Negative charges', 'location','best')
xlabel('Chemical potential [eV]','fontsize', 16,'fontweight','b')
ylabel('Charge concentrations [1/cm^3]','fontsize', 16,'fontweight','b')
set(gca,'fontsize',16,'fontweight','normal')
%%%%%%%%%%%%%%%%%%%

[~,nmin]=min(abs(p.'+NDp-n.'))

mu(nmin)
n(nmin)

%% 7

Tvec=50:50:1000

for i=1:length(Tvec)
    T=Tvec(i)
    
    mu=0:0.001:EC;
    ED=0.044; %Phosphorus
    ND=10^16;

    %calculate n:
    [EM,muM]=meshgrid(E,mu);
    noccT=1./(exp((EM-muM)./(T*k))+1);
    gT=pi*(8*m)^1.5/(2*h^3)*real(sqrt(EM-EC));
    nfunceT=gT.*noccT;
    n=trapz(E,nfunceT,2);

    %calculate p:
    noccT=1./(exp(-(EM-muM)./(T*k))+1);
    gT=pi*(8*m)^1.5/(2*h^3)*real(sqrt(EV-EM));
    nfunceT=gT.*noccT;
    p=trapz(E,nfunceT,2);

    NDp=ND*(1+2*exp(-((EC-ED)-mu)/(k*T))).^-1;

    %semilogy(mu,p.'+NDp, mu,n)
    
    [~,nmin]=min(abs(p.'+NDp-n.'))

    mufin(i)=mu(nmin)
    nfin(i)=n(nmin)
    
end

%%%%%%%%%%%%%%%%%%%
figure(6)
plot(Tvec,mufin,'o-')
title({'#7.1: Chemical potential as a function of temperature'}, ...
    'fontsize', 14, 'fontweight','b')
yLimit = get(gca,'XLim');
line(yLimit,[EC EC], 'color', 'r', 'lineStyle', '-')
line(yLimit,[EV EV], 'color', 'r', 'lineStyle', '-')
line(yLimit,[EC/2 EC/2], 'color', 'r', 'lineStyle', '--')
line(yLimit,[EC-ED EC-ED], 'color', 'r', 'lineStyle', '--')
xlabel('Temperature [K]','fontsize', 16,'fontweight','b')
ylabel('Chemical potential [eV]','fontsize', 16,'fontweight','b')
set(gca,'fontsize',16,'fontweight','normal')
legend('Chemical potential', 'location','best')

figure(7)
semilogy(1./Tvec,nfin,'o-')
title({'#7.2: Electron density as a function of T'}, ...
    'fontsize', 14, 'fontweight','b')
xlabel('Inverse Temperature [1/K]','fontsize', 16,'fontweight','b')
ylabel('Electron density [1/cm^3]','fontsize', 16,'fontweight','b')
set(gca,'fontsize',16,'fontweight','normal')
%%%%%%%%%%%%%%%%%%%




