function [avg, amp_avg, amp_singletrial, onsetamp_avg,onsetval_avg, FWHM_avg, FWHM_onset_avg, FWHM_singletrial, FWHM_onset_singletrial, RT_avg, RTon_avg, decay_avg, decay8020_avg]=analysis_ephys_2exp(recording, stim1, onsetval)
s=size(recording,1);
nsweeps=size(recording,2);
%average trace
avg=zeros(1,s);
for j=1:s
avg (j) = mean(recording(j,:));
end
%amplitude
[amp_singletrial_raw, ~]=max(recording(stim1:(stim1+200),:));
amp_singletrial=amp_singletrial_raw-mean(recording((stim1-400):stim1,:));
[amp_avg_raw, peakpos_avg]=max(avg(stim1:(stim1+200)));
peakpos_avg=peakpos_avg+stim1-1;
amp_avg=amp_avg_raw-mean(avg((stim1-400):stim1));
%FWHM from baseline to peak
peakprofile_avg=avg(stim1:(stim1+300));
time=[0:0.05:(0.05*(length(peakprofile_avg)-1))];
HM_avg=amp_avg/2+mean(avg((stim1-400):stim1));
aboveHM_avg=peakprofile_avg>HM_avg;
FWHM_intpos1_avg = find(aboveHM_avg, 1, 'first');
FWHM_intpos2_avg = find(aboveHM_avg, 1, 'last');
FWHM_avg=time(FWHM_intpos2_avg)-time(FWHM_intpos1_avg);
FWHM_time1_avg = time(FWHM_intpos1_avg)+stim1*0.05-0.05;
FWHM_time2_avg = time(FWHM_intpos2_avg)+stim1*0.05-0.05;
FWHM_singletrial=zeros(1,nsweeps);
for p=1:nsweeps
peakprofile=recording(stim1:(stim1+300),p);
HM=amp_singletrial(p)/2+mean(recording((stim1-400):stim1,p));
aboveHM=peakprofile>HM;
FWHM_intpos1 = find(aboveHM, 1, 'first');
FWHM_intpos2 = find(aboveHM, 1, 'last');
FWHM_singletrial(p)=time(FWHM_intpos2)-time(FWHM_intpos1);
end
%find onset on average
onprofile=avg((peakpos_avg-50):(peakpos_avg-2));
ontime=0:0.05:0.05*(length(onprofile)-1);
%derive 3 times and find local maxima on 3rd derivative
d1=diff(onprofile)./diff(ontime);
d2=diff(d1)./diff(ontime(1:48));
d3=diff(d2)./diff(ontime(1:47));
[~,locs] = findpeaks(d3);
%plot rising phase and located maxima
fig=figure;
h(1)=subplot(2,1,1);
plot(ontime(4:end),onprofile(4:end),'b')
xlim([ontime(4) ontime(end)])
h(2)=subplot(2,1,2);
p1=plot(ontime(1:46),d3,'k');
hold on
p2=plot(ontime(locs),d3(locs),'r*','MarkerSize',20);
xlim([ontime(1) ontime(46)])
%select points using brush tool
hbrush=brush;
set(hbrush,'Enable','on');
hstart = uicontrol('Position',[250 10 80 20],'String','Ok',...
              'Callback','uiresume(gcbf)');
uiwait(gcf); 
%retrieve position of the selected point in the vector containing all maxima (logical vector), then find position value within rising phase profile, then correct for time difference in original profile  
brushed_locs = get(p2, 'BrushData');
onsetpos=locs(brushed_locs==1);
onset_avg=onsetpos+1+(peakpos_avg-51);
onsetval_avg=avg(onset_avg);
onsetamp_avg=onsetval_avg-mean(avg((stim1-400):stim1));
close(fig)
%FWHM from onset to peak
HM_onset_avg=(amp_avg_raw-onsetval_avg)/2+onsetval_avg;
aboveHM_onset_avg=peakprofile_avg>HM_onset_avg;
FWHM_intpos1_onset_avg = find(aboveHM_onset_avg, 1, 'first');
FWHM_intpos2_onset_avg = find(aboveHM_onset_avg, 1, 'last');
FWHM_onset_avg=time(FWHM_intpos2_onset_avg)-time(FWHM_intpos1_onset_avg);
FWHM_onset_time1_avg = time(FWHM_intpos1_onset_avg)+stim1*0.05-0.05;
FWHM_onset_time2_avg = time(FWHM_intpos2_onset_avg)+stim1*0.05-0.05;
FWHM_onset_singletrial=zeros(1,nsweeps);
for p=1:nsweeps
peakprofile=recording(stim1:(stim1+200),p);
HM=(amp_singletrial_raw(p)-onsetval(p))/2+onsetval(p);
aboveHM=peakprofile>HM;
FWHM_intpos1 = find(aboveHM, 1, 'first');
FWHM_intpos2 = find(aboveHM, 1, 'last');
FWHM_onset_singletrial(p)=time(FWHM_intpos2)-time(FWHM_intpos1);
end
%Calculate 20-80% rise time on the average
riseprofile_avg=peakprofile_avg(1:(peakpos_avg-stim1+1));
[~, t20] = min(abs(riseprofile_avg-(0.2*amp_avg+mean(avg((stim1-400):stim1)))));
[~, t80] = min(abs(riseprofile_avg-(0.8*amp_avg+mean(avg((stim1-400):stim1)))));
i20=peakprofile_avg(t20);
i80=peakprofile_avg(t80);
RT_avg=time(t80)-time(t20);
t20_avg = time(t20)+stim1*0.05-0.05;
t80_avg = time(t80)+stim1*0.05-0.05;
%Calculate 20-80% rise time from onset on the average
[~, t20on] = min(abs(riseprofile_avg-(0.2*(amp_avg_raw-onsetval_avg)+onsetval_avg)));
[~, t80on] = min(abs(riseprofile_avg-(0.8*(amp_avg_raw-onsetval_avg)+onsetval_avg)));
i20on=peakprofile_avg(t20on);
i80on=peakprofile_avg(t80on);
RTon_avg=time(t80on)-time(t20on);
t20on_avg = time(t20on)+stim1*0.05-0.05;
t80on_avg = time(t80on)+stim1*0.05-0.05;
%Calculate 80-20% decay time on the average
decayprofile_avg=peakprofile_avg((peakpos_avg-stim1+1):end);
[~, t20d] = min(abs(decayprofile_avg-(0.2*amp_avg+mean(avg((stim1-400):stim1)))));
[~, t80d] = min(abs(decayprofile_avg-(0.8*amp_avg+mean(avg((stim1-400):stim1)))));
i20d=decayprofile_avg(t20d);
i80d=decayprofile_avg(t80d);
decay8020_avg=time(t20d)-time(t80d);
t20d_avg = time(t20d)+peakpos_avg*0.05-0.05;
t80d_avg = time(t80d)+peakpos_avg*0.05-0.05;
%fit double exponential to average decay profile and calculate decay tau
decayprofile=avg((peakpos_avg+5):(stim1+400));
%decayprofile=peakprofile_avg((peakpos_avg-stim1+6):end);
decayprofile=transpose(decayprofile);
decaytime=[0:0.05:(0.05*(length(decayprofile)-1))];
decaytime=transpose(decaytime);
decayfit = fit(decaytime, decayprofile,'exp2');
decayfig=decayfit(decaytime);
decay_coeff = coeffvalues(decayfit);
decay_avg=zeros(1,2);
decay_avg(1)=-1/decay_coeff(2);
decay_avg(2)=-1/decay_coeff(4);
%plot it all to check
time=[0:0.05:(0.05*(length(avg)-1))];
figure
plot(time,recording,':k')
hold on
plot(time,avg,'r', 'LineWidth',1.5)
hold on
plot(time(peakpos_avg),avg(peakpos_avg),'p','MarkerSize',10, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm')
hold on
line([FWHM_time1_avg FWHM_time2_avg], [HM_avg HM_avg], 'Color', 'b', 'LineWidth', 2);
hold on
line([FWHM_onset_time1_avg FWHM_onset_time2_avg], [HM_onset_avg HM_onset_avg], 'Color', 'b', 'LineWidth', 2);
hold on
line([t20_avg t80_avg], [i20 i80], 'Color', 'c', 'LineWidth', 2);
hold on
line([t20on_avg t80on_avg], [i20on i80on], 'Color', 'y', 'LineWidth', 2);
hold on
line([t20d_avg t80d_avg], [i20d i80d], 'Color', 'y', 'LineWidth', 2);
hold on
plot(time((peakpos_avg+5):(stim1+400)),decayfig,'c', 'LineWidth',2)
xlim([time(stim1) time(stim1+300)])
%save everything
%folder_name = uigetdir;
%oldFolder = cd(folder_name);
%save stuff
%csvwrite('avg.csv', avg);
%csvwrite('amp_singletrial.csv',amp_singletrial);
%csvwrite('FWHM_singletrial.csv',FWHM_singletrial);
%csvwrite('FWHM_onset_singletrial.csv',FWHM_onset_singletrial);
%cd(oldFolder);
end