function rhs=odefun112(t,x)
global b p k c I0H I0E H0 dT dI dH dE Iin w tARTstart tARTend effic r0 pL0 E0L rmax Em

% version of odefun11, with a single phenotype of E 


T=x(1); I=x(2); L=x(3); H=x(4); E=x(5); 

% control parameters
alphaH=1-exp(-I/I0H);
alphaH=alphaH*(1-H/Em)^(1/3); % limit on H
 sigma=1-exp(-alphaH*H/H0-I/I0E); % added 2 signals
% sigma=1-exp(-H/H0-I/I0E); % added 2 signals
sigma=sigma*(1-E/Em)^(1/3); % limit on E
r=r0+rmax*E/(E0L+E);
pL=pL0*E0L/(E0L+E);

% cutoff of infection below 1 cell; H and E are measured in % max
p1=p*(I>Iin/2); 
% alphaH=alphaH*(H>100*Iin/2); 
% sigma=sigma*(E>100*Iin/2);
alphaH=alphaH*(H>Iin/2); 
sigma=sigma*(E>Iin/2);


% ART
if t > tARTstart && t < tARTend
    p1=p1*(1-effic);
end
 
rhs=zeros(5,1);
rhs(1)=b - p1*I*T - dT*T; % derivative T
rhs(2)=(1-pL)*p1*I*T - k*I*E - dI*I+r*L; %I
rhs(3)=pL*p1*I*T - r*L; %L
rhs(4)=c*alphaH*H - p1*w*I*H - dH*H; %H
rhs(5)=c*sigma*E - dE*E; %E
end