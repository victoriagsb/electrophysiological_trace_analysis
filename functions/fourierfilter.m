function [rectransform1]=fourierfilter(recording, L, Fs,cutofffreq)
%ephys scripts
%Fs=Sampling frequency; 
%cutofffreq=filtering cutoff frequency  
%L=Length of signal
T = 1/Fs;             %   Sampling period       
t = (0:L-1)*T;

Y = fft(recording(:,1));

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

cutoff = find(f==cutofffreq);
cutoff2=length(P2)-(cutoff-1);
rectransform1=zeros(size(recording));
nsweeps=length(recording(1,:));
for p=1:nsweeps
ftrans = fft(recording(:,p));
ftrans(cutoff:cutoff2)=0;
rectransform1(:,p)=real(ifft(ftrans));
end

time=[0:T/1000:T/1000*(L-1)];
for n=1:nsweeps
figure(n)
plot(time,recording(:,n), 'k')
hold on
plot(time,rectransform1(:,n),'r')
waitfor n
end
end