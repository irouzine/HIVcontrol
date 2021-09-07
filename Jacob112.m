function J=Jacob112(t,x)
global p k c I0H I0E H0 dT dI dH dE Iin w tARTstart tARTend effic r0 pL0 E0L rmax Em

T=x(1); I=x(2); L=x(3); H=x(4); E=x(5); 

% control functions and their derivatives
alphaH=1-exp(-I/I0H);
alphaH=alphaH*(1-H/Em)^(1/3); % limit on H
sigma=1-exp(-alphaH*H/H0-I/I0E);
%  sigma=1-exp(-H/H0-I/I0E);
sigma=sigma*(1-E/Em)^(1/3); % limit on E

alphaHder=exp(-I/I0H)/I0H;
alphaHderH=-alphaH*(1/3)/Em/(1-H/Em);
%alphaEder=exp(-I/I0E)/I0E;
 %sigmaderH=(1-sigma)/alphaH/H0; %??? must be *alphaH ??
  sigmaderH=(1-sigma)*alphaH/H0; 
%sigmaderH=(1-sigma)/H0;
 sigmaderI=(1-sigma)*(1/I0E+alphaHder/H0);
%sigmaderI=(1-sigma)*(1/I0E);
sigmaderE=-sigma*(1/3)/Em/(1-E/Em);

r=r0+rmax*E/(E0L+E);
pL=pL0*E0L/(E0L+E);
r_der=rmax*E0L/(E0L+E)^2;
pL_der=-pL0*E0L/(E0L+E)^2;

 % cutoff of infection and proliferation below 1 cell, H and E measured in
 % % max, 
p1=p*(I>Iin/2); 
% alphaH=alphaH*(H>100*Iin/2); 
% sigma=sigma*(E>100*Iin/2);
 alphaH=alphaH*(H>Iin/2); 
 sigma=sigma*(E>Iin/2);

% ART
if t > tARTstart && t < tARTend
    p1=p1*(1-effic);
end

J=zeros(5,5);
%bb - p1*w*I*T - dT*T; %T
J(1,:)=[-p1*I-dT,-p1*T,0,0,0];
%(1-pL)*p1*I*T - k*I*E - dI*I+r*L %I
J(2,:)=[(1-pL)*p1*I, (1-pL)*p1*T-k*E-dI,r,0,-k*I-pL_der*p1*I*T];
%pL*p1*I*T-r*L; %L
J(3,:)=[pL*p1*I,pL*p1*T,-r,0,pL_der*p1*I*T-r_der*L];
%c*alphaH*H - p1*w*I*H - dH*H  %H
J(4,:)=[0,-p1*w*H+c*H*alphaHder,0,c*alphaH-p1*w*I-dH+c*H*alphaHderH,0];
%c*sigma*E - dE*E %EH
J(5,:)=[0,c*E*sigmaderI,0,c*E*sigmaderH,c*sigma+c*sigmaderE*E-dE];

end