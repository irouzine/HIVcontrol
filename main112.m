% dynamics of the HIV model with 2 steady states 
% to predict progressors and controllers
%  includes helper cells, one type of CTL, 
% latency, indirect killing, ART
% 
% simplifed version of main.m with one type of E


global b p k c I0H I0E H0 dT dI dH dE Iin w tARTstart tARTend effic r0 pL0 E0L rmax Em

%% Most important parameters
R0=8;           % basic reproduction number    
I0E=1e-3; %0.03      % avidity th. of ED in I, [T0]
I0H=3e-3 ;   % avidity th. of H in I, [T0]
w= 10;            % waste, indirect killing factor
time=200;       % total time in days
col='r';
yesphase=1;
%% ART parameters
tARTstart=6; % days
tARTend=12; 
effic=    0.97;

%% Latency parameters
r0=0.001;       % reactivation without E
pL0=0.4;        % max latency prob without E
E0L=0.01;%1     % value of E that halv p0
rmax=0.2*dI;    % maximal reactivation rate

%% Time rate parameters  
c=2;            % max division rate of E and H, 1/day
dT=0.3;         % death rate of Target cells,  1/day
dI=1;           % death rate of Infected, 1/day
dE=0.15;         % death rate of Effectors, 1/day
dH=0.15;         % death rate of Helper cells, 1/day

%% Scaling parameters
T0=1;           % normal target cell #, a.u., scales T and I
b=T0*dT;        % linear source of Target cells, cell/day
%R0= p*b/dT/dI
p=R0*dT*dI/b;   % infection coefficient, 1/cell/day
H0=0.01;%1;           % avidity of EH in H, a.u., scales H 
E0=0.01;%1;           % value E which doubles death of I, a.u., scales E 
k= dI/E0;       % killing coefficient
%Em=100;         % ceiling for E in %
Em=1;

%% Initial values
% Hin=1e-3;       % naive cells
% Ein=1e-3;       % naive cells 
 Hin=1e-5;       % naive cells
 Ein=1e-5;       % naive cells 
Iin=1e-11*T0;   % := one cell ; for E and H, 100*Iin 
% I=Iin/2 is the low cutoff where I:=0
Lin=Iin; % one latent cell

y0=[T0 Iin Lin Hin Ein]';
    
%% running dynamics
% Set accuracy
RelTol=0.001;  
Tol=0.001;
options=odeset('RelTol',RelTol,'AbsTol',...
    Tol*[Iin Iin Lin Hin Ein],'NonNegative',[1,2,3,4,5],...
    'InitialStep',1e-3,'Jacobian',@Jacob112);
% options=odeset('RelTol',RelTol,...
%      'NonNegative',[1,2,3,4,5],...
%     'InitialStep',1e-3,'Jacobian',@Jacob112);
[t,y]=ode15s(@odefun112,[0,time],y0,options);
 
T=y(:,1);I=y(:,2);L=y(:,3);H=y(:,4);E=y(:,5);

%% predicted points in I
% low state (dH/dt=0)
I1=0;
I1=-I0H*log(1-(p*w*I1+dH)/c);
I1=-I0H*log(1-(p*w*I1+dH)/c); 
% high state (dE/dt=0 at H=0)
I2=-I0E*log(1-dE/c);      
Istar=(c-dH)/p/w;             
Istar=(c*(1-exp(-Istar/I0H))-dH)/p/w; % upper border (dH/dt=0)

%% plots
% subplot(2,1,1)
% loglog(I,H,'r',L,H,'c')
% ylabel('H/H0');xlabel('I/T0 red and L/T0 cyan')
% title(sprintf('R0=%g,c=%g,I0E=%g,I0H=%g,dT=%g,dI=%g,dH=%g,dE=%g,w=%g, effic=%g\n RelTol=%g,Tol=%g,Iin=%g,r0=%g,pL0=%g,E0L=%g',R0,c,I0E,I0H,dT,dI,dH,dE,w,effic,RelTol,Tol,Iin,r0,pL0,E0L))
% axis([Iin 10 Hin/100 10^5])
% hold on; 
% loglog([I1 I1],[Hin/100 10^3],'--k',[I2 I2],[Hin/100 10^3],'ok',[Istar Istar],[Hin/100 10^5],'-.m'); hold off

subplot(2,1,1)
semilogy(t,L,'c',t,H,'b',t,E,'g',t,T,'k',t,I,'r')
hold on
ylabel('Cell numbers, a.u.'); xlabel('Time, days')
axis([0 time 10^(-8) 10^1])
semilogy([0 time],[I1 I1],'--k',[0 time],[I2 I2],'ok',[0 time], [Istar Istar],'-.m'); hold off
title(sprintf('R0=%g c=%g I0E=%g I0H=%g dT=%g dI=%g dH=%g dE=%g w=%g eff =%g\n RelTol=%g Tol=%g Iin=%g r0=%g pL0=%g E0L=%g rmax=%g Em=%g\n T black I red L cyan H blue E green',...
    R0,c,I0E,I0H,dT,dI,dH,dE,w,effic,RelTol,Tol,Iin,r0,pL0,E0L,rmax,Em))
grid

%phase_diagram
if yesphase
phase_diagram2
end

% adding data
yesdata=0;
if yesdata==1
data_avidity

% rescale and scatter in Y axis
A=1e-5; 
HIC=HIC/A; HAART= HAART/A; VIR= VIR/A;
yHIC=8e-2*rand(1,length(HIC));
yHAART=8e-2*rand(1,length(HAART));
yVIR=8e-2*rand(1,length(VIR));
hold on
loglog(HIC,yHIC,'or',VIR,yVIR,'*m',HAART,yHAART,'ob')
hold off
end