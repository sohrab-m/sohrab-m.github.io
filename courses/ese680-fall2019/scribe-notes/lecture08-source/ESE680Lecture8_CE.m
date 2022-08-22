%Kendall J. Queen
%ESE680 Lecture 8 Scribe
%Two-mass-one-spring system Computational Example

clear;

s = tf('s'); %laplace s
 
a = 0.1;
A0 = [0 0 1 0
      0 0 0 1
      -1.25 1.25 0 0
      1.25 -1.25 0 0];
A1 = [0;0;-sqrt(2)/2; sqrt(2)/2];
A2 = [(3*sqrt(2))/4; (-3*sqrt(2))/4; 0; 0];

B1 = [0;0;1;0];
B2 = [0; 0; 0; -1];

Ac = [0 -0.7195 1 0 
      0, -2.9732 0 1
      -2.5133 4.8548 -1.7287 -0.9616
      1.0063 -5.4097 -0.0081 0.0304];
Bc = [0.72; 2.973; -3.37; 4.419];
Cc = [-1.506 0.494 -1.738 -0.932];

%[num,denom] = ss2tf(A0,[1;1;1;1], [1,1,1,1], 0);
Gp = inv([s s s-1 s; s s s s-1; s+1.25 s-1.25 s s; s-1.25 s+1.25 s s]);

%[num2,denom2] = ss2tf(Ac,Bc,Cc, 0);
z = [s s+0.7195 s-1 s];
y = [s s+2.9732 s s-1];
x = [s+2.5133 s-4.8548 s+1.7287 s+0.9616];
w = [s-1.0063, s+5.4097, s+0.0081, s-0.0304];
pre1 = inv([z;y;x;w]);

pre2 = pre1*Bc;
Gc = Cc * pre2;


abst_init_iqc;  %Initialize IQC-environment

%Define basic signals
w1=signal;
w2=signal;
w3=signal;
u=signal;

%Relate basic signals to the remaining signals
x=Gp*(A1*w3+B1*(u+w1)+B2*w2);
v=transpose(A2)*x;
x2=x(2); %Also x2 = C*x
uc=Gc*x2;



%Define IQC
M=variable; %creates LMI toolbox variable 
M>0;
v'*(a^2*M)*v-w3'*M*w3>0; %descibes the IQC using toolbox
u == uc;
%run solver to find values for g and M
g = iqc_gain_tbx([w1;w2],x2)

iqc_value
value_iqc(M)


