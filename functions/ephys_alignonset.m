function [alignedtopeak_rec, shift_peak, peakpos, alignedtoonset_rec, shift_onset, onset, onsetval, onsetvalamp]=ephys_alignonset(recording, stim1, folder_name)
nsweeps=size(recording,2);
%define time vector
time=[0:0.05:(0.05*(length(recording)-1))];
%find peak values
[~, peakpos]=max(recording(stim1:(stim1+200),:));
peakpos=peakpos+(stim1-1);
%run loop to find peaks in 3rd derivative of the rising phase of the AP and
%choose the peak corresponding to the onset.
onset=zeros(1,nsweeps);
onsetval=zeros(1,nsweeps);
onsetvalamp=zeros(1,nsweeps);
disp('For each sweep, click on peak to choose onset timepoint. Press ok when you have selected the peak.');
for p=1:nsweeps
onprofile=transpose(recording((peakpos(p)-50):(peakpos(p)-2),p));
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
onset(p)=onsetpos+1+(peakpos(p)-51);
onsetval(p)=recording(onset(p),p);
onsetvalamp(p)=onsetval(p)-mean(recording((stim1-300):stim1,p));
close(fig)
end
%plot recording with detected peaks and onset points to check.
for p=1:nsweeps
hfig = figure(p);
plot(time,recording(:,p),'k')
hold on
plot(time(peakpos(p)),recording(peakpos(p),p),'r*','MarkerSize',20)
hold on
plot(time(onset(p)),recording(onset(p),p),'g*','MarkerSize',20)
xlim([time(stim1) (time(stim1)+10)])
waitfor(hfig)
end
%align to everything by eliminating a certain number of points before the
%alignment position and adding the same number of 0 in the end as points
%deleted

%align to peak
alignedtopeak_rec=recording;
shift_peak=zeros(1,nsweeps);
    for p=1:nsweeps
    temp_peak1=recording(:,p);
    shift_peak(p)=peakpos(p)- min(peakpos);
    temp_peak1(1:shift_peak(p))=[];
    temp_peak1(length(temp_peak1):(length(temp_peak1)+shift_peak(p)))=0;
    alignedtopeak_rec(:,p)=temp_peak1;
    end
    
    %align to onset
alignedtoonset_rec=recording;
shift_onset=zeros(1,nsweeps);
    for p=1:nsweeps
    temp_peak1=recording(:,p);
    shift_onset(p)=onset(p)-min(onset);
    temp_peak1(1:shift_onset(p))=[];
    temp_peak1(length(temp_peak1):(length(temp_peak1)+shift_onset(p)))=0;
    alignedtoonset_rec(:,p)=temp_peak1;
    end
    figure
subplot(1,3,1)
plot(recording(stim1:(stim1+500),:))
subplot(1,3,2)
plot(alignedtopeak_rec(stim1:(stim1+500),:))
subplot(1,3,3)
plot(alignedtoonset_rec(stim1:(stim1+500),:))
oldFolder = cd(folder_name);
%save stuff
csvwrite('peakpos_soma.csv',peakpos);
csvwrite('onset_soma.csv',onset);
csvwrite('onsetval_soma.csv',onsetval);
csvwrite('onsetvalamp_soma.csv',onsetvalamp);
csvwrite('shift_peak_soma.csv',shift_peak);
csvwrite('shift_onset_soma.csv',shift_onset);
csvwrite('alignedtoonset_recsoma.csv',alignedtoonset_rec);
csvwrite('alignedtopeak_recsoma.csv',alignedtopeak_rec);
csvwrite('recsoma.csv',recording);
cd(oldFolder);
end