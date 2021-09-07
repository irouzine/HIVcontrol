% plotting analytic and numeric phase diagram w vs H: 
% bistable state, only high state, only low state

global b p k c I0H I0E H0 dT dI dH dE Iin w tARTstart tARTend effic r0 pL0 E0L


%% analytic diagram for fixed w

x=-6:0.01:0; IE=10.^x; 

%% existence border of helper state
IH_high = c*T0/w/R0/dI*(1-sqrt(2*dH/c));

%% I1=I2 border, stability of helper state
I2=-IE*log(1-dE/c); 
IH_sta=-I2.*log(1-(p*w*I2+dH)/c).^(-1);
ii=find(1-(p*w*I2+dH)/c > 0);
IH_sta=IH_sta(ii);IE=IE(ii);
[~,imax]=max(IH_sta); % critical point
IH_sta=IH_sta(1:imax);IE=IE(1:imax);

%% w < w_low, where H are not depleted
[Imax,mm]=max(I);                               % virus peak 
[~,m2]=min(abs(Imax/10-I(mm+1:end))); m2=m2+mm; % integration limit right of peak
intI=sum((t(2:m2)-t(1:m2-1)).*I(2:m2));         % virus integral


w_low=T0/(intI*R0*dI)*(15.0+(c-dH)/dI/R0*log(Imax*IH_sta.^(-1)));
%log(Hin/Iin/T0*1e-3)=15??

if all(w_low < w)
    label='post-treatment control';
elseif all(w_low > w)
    label='elite control';
else 
    label='elite or post-treatment';
end
    
%% plotting
if ~(tARTstart < tARTend && effic > 0)
subplot(2,1,2)

IEmin=1e-7;IEmax=1;
col='m';
loglog(IH_sta,IE, col)
hold on
%loglog([IH_sta(1) IH_sta(end)], IE(end)*[1 1],col) % horizontal
loglog(IH_high*[1 1],[IEmin IEmax],col)  % existence limit

% labeling
im=round(length(IE)/2);
text(IH_sta(1),IE(im),label)
text(IH_sta(1),2*IEmin,'transient control')
text(IH_high,2*IEmin,sprintf('w=%g',w))
xlabel('Helper cell threshold, I_{0H}/T_0'); 
ylabel('Effector cell threshold, I_{0E}/T_0')
title(sprintf('R0=%g,w=%g,c/dI=%g,c/dE=%g,(c-dH)/dI=%g,T0/(intI R0 dI)=%g',R0,w,c/dI,c/dE,(c-dH)/dI,T0/(intI*R0*dI)))
axi=axis; axi(1)=IH_sta(1); axi(2)=10; axi(3)=IEmin; axi(4)=IEmax; axis(axi);
hold off

%im=round(length(iii)/2);
%text(IH(im),w_low(im),10,sprintf('%g',R0))
%text(IH(im),w_low(im),10,sprintf('%g',c))
%text(max(IHc),10,sprintf('%g',dH))
%text(IH(im),w_low(im),sprintf('%g',dI))
%text(IH(im),w_low(im),sprintf('%g',dE))
%text(IH(im),w_low(im),sprintf('%g',dT))
%text(IH(im),w_low(im),sprintf('%g',Hin))

end
