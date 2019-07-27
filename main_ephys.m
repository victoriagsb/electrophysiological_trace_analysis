%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%CODE FOR ANALYSIS OF ELECTROPHYSIOLOGICAL RECORDINGS
%
%AUTHOR: Victoria Gonzalez Sabater
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%FILE:main_ephys.m
%
%CODE DESCRIPTION:Code for analysis of electrohpysiological traces, 
%including a membrane test and repeated time locked stimulation.
%
%get resitance and capacitante values from  -10mV step recording
voltage_step=yourmembranetest;
[Rs, Ri, C]=stepanalysis(voltage_step);

%upload your ephys trace and indicate stimulation timepoint
recording=yourephystrace;
stim1=25885;

%Optional - filter recording with a low-pass 3Khz filter

[rec]=fourierfilter(recording,200000, 20000, 3000);

%analyse recording

folder_name = uigetdir;

%Align traces to peak and onset of the action potential
[alignedtopeak_rec, shift_peak, peakpos, alignedtoonset_rec, shift_onset, onset, onsetval, onsetvalamp]=ephys_alignonset(rec, stim1, folder_name);
%Extract average trace and action potential parameters from not aligned,
%aligned to peak and aligned to onset
[avg1, amp_avg1, amp_singletrial1, onsetamp_avg1,onsetval_avg1, FWHM_avg1, FWHM_onset_avg1, FWHM_singletrial1, FWHM_onset_singletrial1, RT_avg1, RTon_avg1, decay_avg1, decay8020_avg1]=analysis_ephys_2exp(rec, stim1, onsetval);
[avg_peak, amp_avg_peak, amp_singletrial_peak, onsetamp_avg_peak,onsetval_avg_peak, FWHM_avg_peak, FWHM_onset_avg_peak, FWHM_singletrial_peak, FWHM_onset_singletrial_peak, RT_avg2, RTon_avg_peak, decay_avg_peak, decay8020_avg_peak]=analysis_ephys_2exp(alignedtopeak_rec, (stim1-50), onsetval);
[avg_onset, amp_avg_onset, amp_singletrial_onset, onsetamp_avg_onset,onsetval_avg_onset, FWHM_avg_onset, FWHM_onset_avg_onset, FWHM_singletrial_onset, FWHM_onset_singletrial_onset, RT_avg_onset, RTon_avg_onset, decay_avg_onset, decay8020_avg_onset]=analysis_ephys_2exp(alignedtoonset_rec, (stim1-50), onsetval);
avg_all=zeros(size(avg1,2),3);
avg_all(:,1)=avg1;
avg_all(:,2)=avg_peak;
avg_all(:,3)=avg_onset;
results=zeros(3,9);
results(1,:)=[amp_avg1 FWHM_avg1 FWHM_onset_avg1 RT_avg1 RTon_avg1 min(abs(decay_avg1)) decay8020_avg1 onsetamp_avg1 onsetval_avg1];
results(2,:)=[amp_avg_peak FWHM_avg_peak FWHM_onset_avg_peak RT_avg2 RTon_avg_peak min(abs(decay_avg_peak)) decay8020_avg_peak onsetamp_avg_peak onsetval_avg_peak];
results(3,:)=[amp_avg_onset FWHM_avg_onset FWHM_onset_avg_onset RT_avg_onset RTon_avg_onset min(abs(decay_avg_onset)) decay8020_avg_onset onsetamp_avg_onset onsetval_avg_onset];

oldFolder = cd(folder_name);
%save stuff
csvwrite('amp_singletrial.csv',amp_singletrial1);
csvwrite('FWHM_singletrial.csv',FWHM_singletrial1);
csvwrite('FWHM_onset_singletrial.csv',FWHM_onset_singletrial1);
csvwrite('amp_singletrial_peak.csv',amp_singletrial_peak);
csvwrite('FWHM_singletrial_peak.csv',FWHM_singletrial_peak);
csvwrite('FWHM_onset_singletrial_peak.csv',FWHM_onset_singletrial_peak);
csvwrite('amp_singletrial_onset.csv',amp_singletrial_onset);
csvwrite('FWHM_singletrial_onset.csv',FWHM_singletrial_onset);
csvwrite('FWHM_onset_singletrial_onset.csv',FWHM_onset_singletrial_onset);
csvwrite('avg_all.csv',avg_all);
csvwrite('rec.csv',rec);
csvwrite('results.csv',results);
cd(oldFolder); 
