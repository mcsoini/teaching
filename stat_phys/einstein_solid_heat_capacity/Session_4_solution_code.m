%% Matlab Session 4 - March 7, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1

clear all
clc
N=25;
q=0:100;
kb=1.38*10^-23; %J/K
sok=gammaln(q+N)-gammaln(q+1)-...
    ones(1,length(q)).*gammaln(N)

figure(1)
plot(q,sok,'b')
xlabel ('Number of energy units','fontname','arial','fontsize',16,'fontweight','bold')
ylabel ('Entropy S/k_b','fontname','arial','fontsize',16,'fontweight','bold')
%legend('plot1','location','best')
Q1_title1 = {'- 1 -';'Entropy vs Energy units'};
title(Q1_title1,'fontsize',16,'fontweight','bold')
set(gca,'fontname','arial','fontsize',16,'fontweight','normal')
hold on
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2

%Energy in units of epsilon
U=q;
T=diff(U)./diff(sok)

%Interpolate the vector U
Uint=(diff(U)*0.5+U(1:(end-1)))

figure(2)
plot(Uint,T,'b.-')
xlabel ('Energy U/\epsilon','fontname','arial','fontsize',16,'fontweight','bold')
ylabel ('Temperature T/[\epsilon/k_b]','fontname','arial','fontsize',16,'fontweight','bold')
%legend('plot1','location','best')
Q1_title1 = {'- 2 -';'Temperature vs Energy'};
title(Q1_title1,'fontsize',16,'fontweight','bold')
set(gca,'fontname','arial','fontsize',16,'fontweight','normal')
hold on
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3
%Interpolate the vector sok
sokint=(diff(sok)*0.5+sok(1:(end-1)))

figure(3)
plot(T,sokint,'.-')
xlabel ('Temperature T/[\epsilon/k_b]','fontname','arial','fontsize',16,'fontweight','bold')
ylabel ('Entropy S/k_b','fontname','arial','fontsize',16,'fontweight','bold')
%legend('plot1','location','best')
Q1_title1 = {'- 3 -';'Entropy vs Temperature'};
title(Q1_title1,'fontsize',16,'fontweight','bold')
set(gca,'fontname','arial','fontsize',16,'fontweight','normal')


%%
%%Problem 4: Heat capacity dU/dT
C=diff(Uint)./diff(T)
Tint=(diff(T)*0.5+T(1:(end-1)))
Uint2=(diff(Uint)*0.5+Uint(1:(end-1)))

figure(4)
plot(Tint,C,'b.-')
xlabel ('Temperature T/[\epsilon/k_b]','fontname','arial','fontsize',16,'fontweight','bold')
ylabel ('Heat capacity C_V/k_b','fontname','arial','fontsize',16,'fontweight','bold')
legend('N=25','N=50','location','best')
Q1_title1 = {'- 4/2 -';'Heat capacity vs Temperature'};
title(Q1_title1,'fontsize',16,'fontweight','bold')
set(gca,'fontname','arial','fontsize',16,'fontweight','normal')
hold on

figure(5)
plot(Uint2,C,'b.-')
xlabel ('Energy U/\epsilon','fontname','arial','fontsize',16,'fontweight','bold')
ylabel ('Heat capacity C_V/k_b','fontname','arial','fontsize',16,'fontweight','bold')
legend('N=25','N=50','location','best')
Q1_title1 = {'- 4/1 -';'Heat capacity vs Energy'};
title(Q1_title1,'fontsize',16,'fontweight','bold')
set(gca,'fontname','arial','fontsize',16,'fontweight','normal')
hold on

%%
%%Problem 5
clear all
clc

A_Pb=importdata('lead.dat');
C_Pb=A_Pb(:,2);
T_Pb=A_Pb(:,1);

A_Al=importdata('aluminum.dat');
C_Al=A_Al(:,2);
T_Al=A_Al(:,1);

figure(6)
plot(T_Al,C_Al,'b*',T_Pb,C_Pb,'ro')
xlabel ('Temperature T/[K]','fontname','arial','fontsize',16,'fontweight','bold')
ylabel ('Heat capacity C_V/[J/K]','fontname','arial','fontsize',16,'fontweight','bold')
legend('Aluminum','Lead','location','best')
Q1_title1 = {'- 5 -';'Experimental data - heat capacity'};
title(Q1_title1,'fontsize',16,'fontweight','bold')
set(gca,'fontname','arial','fontsize',16,'fontweight','normal')

%%
%%Problem 6:
T_exp=T_Pb;
C_exp=C_Pb;
%T_exp=T_Al;
%C_exp=C_Al;

Na=6.022*10^23;
kb=1.38*10^-23; %J/K
ev=1.6*10^-19;
C_func=@(eps) sum((C_exp-3*Na*kb.*(eps./(kb.*T_exp)).^2.*exp(eps./(kb.*T_exp))./(exp(eps./(kb.*T_exp))-1).^2).^2)

epsf=fminsearch(C_func,0.01*1.6*10^-19)/ev

%%
%%Problem 7
eps_Pb=0.0055*ev
eps_Al=0.0248*ev;
T=0:500
C_fit_Pb=3*Na*kb.*(eps_Pb./(kb.*T)).^2.*exp(eps_Pb./(kb.*T))./(exp(eps_Pb./(kb.*T))-1).^2;
C_fit_Al=3*Na*kb.*(eps_Al./(kb.*T)).^2.*exp(eps_Al./(kb.*T))./(exp(eps_Al./(kb.*T))-1).^2;
plot(T_Pb,C_Pb,'or',T_Al,C_Al,'ob',T,C_fit_Pb,'-r',T,C_fit_Al,'-b')
xlabel ('Temperature T/[K]','fontname','arial','fontsize',16,'fontweight','bold')
ylabel ('Heat capacity C_V/[J/K]','fontname','arial','fontsize',16,'fontweight','bold')
legend('Experiment Pb','Experiment Al','Fit Pb','Fit Al','location','best')
Q1_title1 = {'- 7 -';'Fitted heat capacities'};
title(Q1_title1,'fontsize',16,'fontweight','bold')
set(gca,'fontname','arial','fontsize',16,'fontweight','normal')







