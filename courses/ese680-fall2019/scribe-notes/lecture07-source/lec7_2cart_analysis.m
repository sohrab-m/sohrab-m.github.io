%Triple lead compenstaor
s = zpk('s'); %Laplace s
C = 100*ss((s+1)/(.001*s+1))^3;
%Uncertain parameters
k = ureal('k',1,'percent',20);
m1 = ureal('m1',1,'percent',20);
m2 = ureal('m2',1,'percent',20);
%Uncertain cart models
G1 = 1/s^2/m1;
G2 = 1/s^2/m2;
%Uncertain model of closed-loop systems
% Spring-less inner block F(s)
F = [0;G1]*[1 -1]+[1;-1]*[0,G2];
P = lft(F,k);
% Uncertain open-loop model is
L = P*C;
T = feedback(L,1)
%Nominal plant
Pnom = zpk(P.nominal)
%Nominal closed-loop stability
Tnom = zpk(T.nominal);
maxrealpole = max(real(pole(Tnom)))

%Robust stability
% Show report and compute sensitivity
opt = robOptions('Display','on','Sensitivity','on');
[StabilityMargin,wcu] = robstab(T,opt);

%Worst case performance
[PeakGain,wcu] = wcgain(T);
PeakGain
%Substitute worst case param into T to compter wc closed loop TF, Twc
Twc = usubs(T,wcu); 

Trand = usample(T,4);         % 4 random samples of uncertain model T
clf
subplot(311), bodemag(Trand,'b',Twc,'r',{10 1000});  % plot Bode response
subplot(312), step(Trand,'b',Twc,'r',0.2);           % plot step response
figure()
wcsigma(T, {10 1000});           % plot step response