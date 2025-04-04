%{
Script to reproduce Figure 2 of de Cheveigné 2023, Is EEG better left
alone?, submitted to Scientific Reports.
%}

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Check for required data and toolboxes.
%}
datapath='/Users/adc/data/ds003061/';
iSub=1; iRun=1;
eegname=['sub-',num2str(iSub, '%03d'), '/eeg/', ...
    'sub-',num2str(iSub, '%03d'), '_task-P300_run-',num2str(iRun), '_eeg.set'];
if 2~=exist([datapath,eegname])
    error('Download data from https://doi.org/10.18112/openneuro.ds003061.v1.1.2 and adjust path in this script.')
end
if 2~=exist('nt_version')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and adjust Matlab path');
end
if 2~=exist('ft_read_header')
    error('Download FieldTrip from https://www.fieldtriptoolbox.org and ajust Matlab path to include fileio/, or use other data-reading tool');
end
if 2~=exist('ft_plot_topo')
    error('Download FieldTrip from https://www.fieldtriptoolbox.org and ajust Matlab path to include plotting/');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Panel (a) of Figure 2 illustrates slow drift in EEG and the effect of
trial-based detrending.
%}
figure(1); clf;
set(gcf, 'position', [616 115 600 200])

% plot raw data for 2 subjects to illustrate various kinds of 'slow drift'
subplot 131;
iSub=1; iRun=1;
[x,sr]=getdata(datapath,iSub,iRun);
plot((1:size(x,1))/sr, x/1000); 
ylim([-4 4])
set(gca,'fontsize',14);
set(gca,'xtick', [0 300 600 900]); xlabel('time (s)');
set(gca, 'ytick', [-4 -2 0 2 4]);; ylabel('mV'); title('Subject 1')
xlim([0 size(x,1)/sr])

subplot 132;
iSub=12; iRun=2;
[x,sr]=getdata(datapath,iSub,iRun);
plot((1:size(x,1))/sr, x/1000); 
ylim([-4 4]);
set(gca,'fontsize',14);
set(gca,'xtick', [0 300 600 900]); xlabel('time (s)'); title('Subject 12');
set(gca, 'ytick', []); 
xlim([0 size(x,1)/sr])
plot_tweak([-0.04 0.09 0 -.14])

% plot trial-averaged waveform of one channel, before and after detrending
iSub=4; iRun=3; % choose subject for best visual clarity
[x,sr,info]=getdata(datapath,iSub,iRun);
pre=2; % s
post=2; %s
[standard,oddball,button,standard_nodetrend]=gettrials(x,info,pre,post,sr);
t=linspace(-1,2,size(standard,1));

% plot without detrend
subplot 233;
ch=11; % this channel looks good
a=mean(standard_nodetrend(:,ch,:),3);
plot(t,a, 'k:'); hold on
a(t<0)=nan; a(t>1)=nan; % clip to [0 1]
plot(t,a, 'k', 'linewidth', 1);
set(gca, 'fontsize',14)
set(gca, 'xgrid', 'on', 'xtick', [ 0 1]);
ylim([-8 8]);
xlim([-1 2])
ylabel('\mu V', 'interpreter', 'tex');
set(gca, 'xtick', [0 1], 'xticklabel', [], 'xgrid', 'on')
title('raw')
plot_tweak([0 .05, 0 -.1])

% plot with detrend
subplot 236
b=mean(standard(:,ch,:),3); 
plot(t,b, 'r:'); hold on
b(t<0)=nan; b(t>1)=nan;
plot(t,b, 'r', 'linewidth', 2);
set(gca, 'fontsize',14)
set(gca, 'xgrid', 'on', 'xtick', [ 0 1]);
ylim([-8  8]);
xlim([-1 2])
set(gca, 'xtick', [0 1], 'xgrid', 'on')
ylabel('\mu V', 'interpreter', 'tex');
set(gca, 'xtick', [0 1], 'xgrid', 'on')
xlabel('time (s)');
title('detrended')
plot_tweak([0 .1, 0 -.1])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Panel (b) of Figure 2 illustrates power line artifact removal
%}
iSub=9; iRun=2; % choose subject for best visual clarity
[x,sr,info]=getdata(datapath,iSub,iRun);
pre=2; % s
post=2; %s
[standard,oddball,button]=gettrials(x,info,pre,post,sr);

p=[];
p.fig2=figure(2); clf;
set(gcf, 'position', [616 115 600 200])
p.nfft=256; % default was 512, widen to accommodate line frequency variations
fline=50; % Hz, line frequency
nremove=3; % number of components to remove
nt_zapline(standard(:,:,1:30),fline/sr, nremove, p);

drawnow; % or matlab gets confused

subplot 121;
set(gca,'fontsize', 14); 
set(gca, 'xticklabel', [0 50 100 150]); xlabel('frequency (Hz)')
plot_tweak([0 .1 0 -.1])

subplot 122;
set(gca,'fontsize', 14); 
set(gca, 'xticklabel', [0 50 100 150]); xlabel('frequency (Hz)')
plot_tweak([0 .1 0 -.1])
c=get(gca,'children'); set(c(2), 'linewidth', 0.5);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Panel (c) of Figure 2 shows the histogram of reaction times and ERPs for
oddball and button-press, across all subjects and trials, and for one
subject/trial.
%}
if ~isdir('./tmp')
    % extract reaction times of all subjects and trials
    mkdir tmp
    a=[];
    for iSub=1:13
        rtime=[];
        for iRun=1:3
            [x,sr,info]=getdata(datapath,iSub,iRun);
            if iSub==12 && iRun==1; continue; end  % corrupt, skip            
            sample=info.sample;
            trial_type=info.trial_type;
            response_time=info.response_time;
            sample=sample(3:end-3); % remove events too close to beginning or end 
            trial_type=trial_type(3:end-3,:); 
            response_time=response_time(3:end-3,:);
            ntrials=numel(sample);
            idxResp=[]; 
            lastwasresp=0;
            for iTrial=1:ntrials
                if lastwasresp
                    lastwasresp=0;
                    continue; % skip stims that immediately follow a button press
                else
                    if strcmp('response',trial_type(iTrial,:)); 
                        idxResp=[idxResp,iTrial];
                        lastwasresp=1;
                    end
                end
            end
            rtime=[rtime; str2num(response_time(idxResp-1,:))];
        end
        a=[a; hist(rtime, 0:10:1000)];
    end
    a=a';
    a=a(1:100,:);
    aa=reshape(a,[4,25,13]); 
    aa=squeeze(sum(aa)); 
    save ./tmp/RTs aa
else 
    load ./tmp/RTs
end
figure(3); clf;
set(gcf, 'position', [616 115 600 350])
subplot 243
plot((0:40:960)/1000, aa, 'k', 'color', [1 1 1]*0.8); 
hold on
plot ((0:40:960)/1000, aa(:,3),'k', 'linewidth', 2);
set(gca,'fontsize', 14);
set(gca,'xtick',0:0.5:1, 'xgrid', 'on')
xlabel('time (s)'); 
h=ylabel('count'); set(h,'position', [-0.25 50 -1])
title('RT');

%{
Panel (c) also shows ERPs calculated based on the first component(s) of a 
DSS analysis designed to emphasize components that reliably repeat across
trials, time-aligned to both the oddball stimulus and the button press. 

The cross-correlation matrix between time-averaged components and
time-averaged electrodes is calculated, and each component is projected
back to the electrode for which the correlation is largest in absolute
value. This is plotted as an ERP with microvolt units. 

The cross-correlation matrix is also used to plot the topography associated
with the first (best) component.

Of principal interest is the similarity of the stimulus ERP to the button
press ERP, suggesting that the stimulus ERP is dominated by
response-related activity.
%}

% concatenate trial matrices
oddball=[]; button=[];
iSub=1; 
for iRun=1:3
    [x,sr,info,hdr]=getdata(datapath, iSub, iRun);
    pre=2; post=2;
    [a,b,c]=gettrials(x,info,pre,post,sr);
    oddball=cat(3,oddball,b);
    button=cat(3,button,c);
end

% button press
xx=button;
xx=nt_demean2(xx,1:(1.5*sr));   % baseline correction: subtract pre -0.5s mean
xx=nt_rereference(xx);          % to get nice topos
[todss,pwr0,pwr1]=nt_dss1(xx);  % DSS to find most reliable components
zz=nt_mmat(xx,todss);
zz=nt_normcol(zz);
[iBest,y,topo]=nt_back(mean(xx,3),mean(zz(:,1:3,:),3)); % 
t=linspace(-pre, post, size(xx,1));
subplot 223
plot(t,y(:,2:3,:), 'color', [1 1 1]*0.8); % components 2 & 3 in gray
hold on
h=plot(t, y(:,1,:), 'b', 'linewidth', 2); % component 1 in blue
set(gca,'xgrid','on');
set(gca,'fontsize', 14); 
xlabel('time (s)'); ylabel('\mu V', 'interpreter', 'tex');
xlim([-1 1])
ylim([-40 20]);
title('button press ERP');
legend(h,hdr.label{iBest(1)}, 'location', 'southeast'); legend boxoff

% topography of normalized correlation
axes('position', [.1 .11 .2 .2]); 
axis off; axis equal
xlim([-1 1]); ylim([-1 1])
[~,~,topo]=nt_back(nt_normcol(mean(xx,3)),mean(zz(:,1,:),3)); % --> correlation
nt_topoplot('biosemi64', topo); 
set(gca, 'clim',[-1 1])
colorbar('location', 'south')
hh=title('correlation');
set(hh, 'position',[0 .6 0], 'fontsize',11)

% oddball
xx=oddball;
xx=nt_demean2(xx,1:(1.5*sr));   % baseline correction: subtract pre -0.5s mean
xx=nt_rereference(xx);          % to get nice topos
[todss,pwr0,pwr1]=nt_dss1(xx);  % DSS to find most reliable components
zz=nt_mmat(xx,todss);
zz=nt_normcol(zz);
[iBest,y,topo]=nt_back(mean(xx,3),mean(zz(:,1:3,:),3)); % 
t=linspace(-pre, post, size(xx,1));
subplot 221
plot(t,y(:,2:3,:), 'color', [1 1 1]*0.8); % components 2 & 3 in gray
hold on
h=plot(t, y(:,1,:), 'b', 'linewidth', 2); % component 1 in blue
set(gca,'xgrid','on');
set(gca,'fontsize', 14); 
set(gca,'xticklabel',[]);
ylabel('\mu V', 'interpreter', 'tex');
xlim([-1 1])
ylim([-40 20]);
title('oddball ERP');
legend(h,hdr.label{iBest(1)}, 'location', 'southeast'); legend boxoff

% cross-correlation between trial-averaged components and electrodes ->
% topography
axes('position', [.1 .58 .2 .2]); axis off; axis equal
axis off; axis equal
xlim([-1 1]); ylim([-1 1])
[~,~,topo]=nt_back(nt_normcol(mean(xx,3)),mean(zz(:,1,:),3)); % --> correlation
nt_topoplot('biosemi64', topo); 
set(gca, 'clim',[-1 1])
colorbar('location', 'south')
hh=title('correlation');
set(hh, 'position',[0 .6 0], 'fontsize',11)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Panel (d) of Figure 2 shows the histogram of reaction times and ERPs for
oddball and button-press, across all subjects and trials, and for one
subject/trial.
%}

figure(4); clf
set(gcf, 'position', [316 115 600 300])

% concatenate trial matrices across runs
standard=[]; oddball=[]; button=[];
iSub=1; 
for iRun=1:3
    [x,sr,info,hdr]=getdata(datapath, iSub, iRun);
    x=x(:,1:64);
    pre=2; post=2;
    [a,b,c]=gettrials(x,info,pre,post,sr);
    standard=cat(3,standard,a);
    oddball=cat(3,oddball,b);
    button=cat(3,button,c);
end

t=linspace(-pre,post,size(xx,1));

% DSS to create denoising matrix to remove motor-related components
xx=button;
xx=nt_demean2(xx,t<-0.5); % baseline correction: subtract pre -0.5s mean
[todss,pwr0,pwr1]=nt_dss1(xx);
fromdss=pinv(todss);
NREMOVE=3;
D=todss(:,NREMOVE+1:end)*fromdss(NREMOVE+1:end,:); % denoising matrix
FOCUS=find(t>-0.5 & t<1); % focus on a shorter epoch for sensory responses

% sensory response to standard (no denoise)
xx=standard;
xx=xx(FOCUS,:,:);
xx=nt_demean2(xx, find(t(FOCUS)<0)); % baseline correction
[todss,pwr0,pwr1]=nt_dss1(xx);
zz=nt_mmat(xx,todss);
[iBest1,y_standard_nodenoise,topo]=nt_back(mean(xx,3),mean(zz(:,1,:),3)); % 
    
% sensory response to standard (w denoise)
xx=standard;
xx=nt_mmat(xx,D); % apply denoising matrix
xx=xx(FOCUS,:,:);
xx=nt_demean2(xx, find(t(FOCUS)<0 )); % baseline correction
[todss_standard,pwr0,pwr1]=nt_dss1(xx);
zz=nt_mmat(xx,todss_standard);
[iBest2,y_standard,topo]=nt_back(mean(xx,3),mean(zz(:,1,:),3)); % 
xx_standard=xx;
    
% sensory response to oddball
xx=oddball;
xx=nt_mmat(xx,D); % apply denoising matrix
xx=nt_demean2(xx);
xx=nt_detrend(xx,1); 
xx=xx(FOCUS,:,:);
xx=nt_demean2(xx, find(t(FOCUS)<0)); % baseline correction: subtract pre 0 mean
[todss_oddball,pwr0,pwr1]=nt_dss1(xx);
xx_oddball=xx;    
    
% apply both DSS matrices to both conditions, concatenate for CCA
z11=nt_mmat(xx_standard,todss_standard);
z12=nt_mmat(xx_standard,todss_oddball);
z22=nt_mmat(xx_oddball,todss_oddball);
z21=nt_mmat(xx_oddball,todss_standard);
z1=cat(3,z11,z21);
z2=cat(3,z12,z22);
    
% select best components for each, use CCA to align them
NKEEP1=3;
NKEEP2=3;
[A,B,R]=nt_cca(mean(z1(:,1:NKEEP1,:),3), mean(z2(:,1:NKEEP2,:),3));
y1=nt_mmat(z11(:,1:NKEEP1,:),A); 
y2=nt_mmat(z22(:,1:NKEEP2,:),B);

[iBestA,yA,topoA]=nt_back(mean(xx_standard,3),mean(y1(:,1,:),3)); % standard
[iBestB,yB,topoB]=nt_back(mean(xx_oddball,3),mean(y2(:,1,:),3)); % oddball

% % backprojection matrices
% from_cc1=nt_xcov(nt_normcol(mean(y1,3)), mean(xx_standard,3)) / size(xx_standard,1); 
% from_cc2=nt_xcov(nt_normcol(mean(y2,3)), mean(xx_oddball,3)) / size(xx_standard,1); 
% [~,idxBest]=max(abs(from_cc1(1,:))+abs(from_cc2(1,:)));
% %idxBest=47;
% yy1=nt_mmat(nt_normcol(mean(y1(:,1,:),3)),from_cc1(1,idxBest));
% yy2=nt_mmat(nt_normcol(mean(y2(:,1,:),3)),from_cc2(1,idxBest));

% standard: compare with and without denoising (should be similar)
subplot 222
plot([y_standard_nodenoise,y_standard])
plot(t(FOCUS),y_standard_nodenoise, 'k:') 
hold on; 
plot(t(FOCUS),y_standard, 'k') 
xlim([-0.2 0.7]); ylim([-5 5])
set(gca,'fontsize',14)
set(gca,'xgrid','on');
set(gca,'xticklabel',[], 'xtick',[-0.2 0 0.2 0.4 0.6]);
ylabel('\mu V', 'interpreter', 'tex')
legend('raw', 'denoised', 'location', 'southeast'); legend boxoff
title('standard, DSS 1', 'fontweight','normal', 'fontsize',12)
xlim([-0.2 0.7]); 

subplot 121
plot(t(FOCUS),yA, 'k') ; hold on
plot(t(FOCUS),yB, 'r') 
xlim([-0.2 0.7]); ylim([-15 15])
set(gca,'fontsize',14)
set(gca,'xgrid','on');
set(gca,'ytick',[-10 0 10]);
xlabel('time (s)'); ylabel('\mu V', 'interpreter', 'tex')
legend('standard (P10)', 'oddball (P10)', 'location', 'northwest'); legend boxoff
plot_tweak([0 0.05 0 -0.05])

subplot 224
plot(t(FOCUS),-mean(z22(:,1,:),3), 'r'); hold on 
plot(t(FOCUS),mean(z22(:,2:4,:),3)) 
set(gca,'fontsize',14)
set(gca,'xgrid','on');
set(gca,'ytick',[], 'xtick',[-0.2 0 0.2 0.4 0.6]);
ylim([-2 1.5])
xlabel('time (s)'); ylabel('a.u')
title('oddball, DSS 1-4', 'fontweight','normal', 'fontsize',12)
plot_tweak([0 0.05 0 0])
legend({'1','2','3','4'}', 'orientation', 'horizontal', 'location','south'); legend boxoff
xlim([-0.2 0.7]); 

% plot topographies
[iBestA,~,topoA]=nt_back(nt_normcol(mean(xx_standard,3)),nt_normcol(mean(yA(:,1,:),3))); % standard
[iBestB,~,topoB]=nt_back(nt_normcol(mean(xx_oddball,3)),nt_normcol(mean(yB(:,1,:),3))); % oddball
axes('position', [.08 .17 .25 .25]); 
axis off; axis equal
xlim([-1 1]); ylim([-1 1])
nt_topoplot('biosemi64', topoA); ; set(gca, 'clim',[-1 1])
h=title('standard', 'fontweight','normal'); set(h,'position', [0 0.7 0])
colorbar('location', 'south');
axes('position', [.27 .17 .25 .25]); 
axis off; axis equal
xlim([-1 1]); ylim([-1 1])
nt_topoplot('biosemi64', topoB); set(gca, 'clim',[-1 1])
h=title('oddball', 'fontweight','normal'); set(h,'position', [0 0.7 0])
colorbar('location', 'south');

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Read data from file.
%}
function [x,sr,info,hdr]=getdata(datapath,iSub,iRun)
%
%  x: data (time X channels)
%  sr: sampling rate
%  info: trial information
%  hdr: file header
%
%  datapath: path of directory containing EEG
%  iSub: subject number (1 to 13)
%  iRun: run number (1 to 3)

    eegname=['sub-',num2str(iSub, '%03d'), '/eeg/', ...
        'sub-',num2str(iSub, '%03d'), '_task-P300_run-',num2str(iRun), '_eeg.set'];
    hdr=ft_read_header([datapath, eegname]);
    sr=hdr.Fs;
    x=ft_read_data([datapath, eegname])';
    x=x(:,1:64);
    x=nt_demean(x);
    infoname=['sub-',num2str(iSub, '%03d'), '/eeg/', ...
        'sub-',num2str(iSub, '%03d'), '_task-P300_run-',num2str(iRun), '_events.tsv'];

    info=tdfread([datapath, infoname]);
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
Cut data into trials.
%}
function [standard,oddball,button,standard_nodetrend]=gettrials(x,info,pre,post,sr)
%
%  standard, oddball, button: trial matrices 
%
%  x: raw data
%  info: metadata
%  pre, post: (s) peri-event epoch
%  sr: sampling rate

% events
onset=info.onset;
sample=info.sample;
trial_type=info.trial_type;
value=info.value;
ntrials=size(sample,1);
ntrials=ntrials-5; % seems to be missing data at end
idxStim=[]; idxResp=[]; idxStan=[]; idxOdd=[];
value=value(:,1:20);
for iTrial=1:ntrials
    if strcmp('stimulus',trial_type(iTrial,:)); idxStim=[idxStim,iTrial]; end
    if strcmp('response',trial_type(iTrial,:)); idxResp=[idxResp,iTrial]; end
    if strcmp('standard            ',value(iTrial,:)); idxStan=[idxStan,iTrial]; end
    if strcmp('oddball_with_reponse',value(iTrial,:)); idxOdd=[idxOdd,iTrial]; end
    % ignore 'oddball' and 'standard_with_response'
end
% standards preceded by a response have a big pre-onset artifact --> remove
a=idxStan;
for iTrial=1:numel(idxStan)
    if any(idxResp==idxStan(iTrial)-1 | idxResp==idxStan(iTrial)+1  | idxResp==idxStan(iTrial)-2  | idxResp==idxStan(iTrial)+2 )
        a(iTrial)=0;
    end
end
idxStan(find(a==0))=[];

% make trial matrix
epochsize=(pre+post)*sr; 
advance=pre*sr;
xx=zeros(epochsize,64,ntrials);
for iTrial=2:ntrials
    start=round(onset(iTrial)*sr) - advance;
    xx(:,:,iTrial)=x(start+(0:epochsize-1),:);
end
standard=xx(:,:,idxStan); %
oddball=xx(:,:,idxOdd); %
button=xx(:,:,idxResp); %

standard_nodetrend=nt_demean2(standard); % need this for (a)
order=1; % linear detrend
standard=nt_detrend( nt_demean2(standard),order);
oddball=nt_detrend(nt_demean2(oddball),order);
button=nt_detrend(nt_demean2(button),order);

end % function               
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
Twean position of plot.
%}
function plot_tweak(tweak, ax)
%plot_tweak(tweak,ax) - tweak position and size of plot
%
%  tweak: vector of values to adjust plot: [x_shift, y_shift, x_expand, y_expand]
%  ax: handle to axes
if nargin<2; ax=gca; end
axes(ax)
tweak(end+1:4)=0;
tweak=tweak(:)';
set(ax,'position',get(ax,'position')+tweak);
end % function

