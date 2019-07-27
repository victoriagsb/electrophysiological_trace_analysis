function [Rs, Ri, C]=stepanalysis(step)
s=size(step,1);
%average trace
avg_step=zeros(1,s);
for j=1:s
avg_step (j) = mean(step(j,:));
end
%Calculting series resistance and input resistance in MOhm
%peak current
amp_step=min(avg_step);
bl=mean(avg_step(1:1000));
I=abs(amp_step-bl);
amp_r=mean(avg_step(1340:1355));
Ir=abs(amp_r-bl);
%Ohm's law
Rs=10/I*1000;
R=10/Ir*1000;
Ri=R-Rs;
%Calculating capacitance in pF
%integrate are under curve
Cstep=avg_step(1156:1356);
Cstep=(Cstep-amp_r).*10^(-12);
Cstep=Cstep(Cstep<=0);
Ctime=[0:0.00005:(0.00005*(length(Cstep)-1))];
Q = trapz(Ctime,Cstep);
%Q=CV in picoFarad
C=abs(Q)/(0.010*10^(-12));
end