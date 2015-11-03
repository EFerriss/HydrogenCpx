% Elizabeth Ferriss, 2013
% Diffusion of water in clinopyroxene project
clear all; close all;

%% database of thicknesses in microns
% K5
th_K5a=3450;                            % crystal thickness in microns of K5 // a* 
th_K5b=[1614 1609 1612 1602 1607];      % crystal thickness in microns of K5 // b
th_K5c=[1756 1759 1748 1763 1759];      % crystal thickness in microns of K5 // c
% K4
th_K4a=7000;                            % crystal thickness in microns of K4 // a* 
th_K4b=[2185 2190 2188 2185 2188];      % crystal thickness in microns of K4 // b
th_K4c=[1546 1551 1536 1548 1548];      % crystal thickness in microns of K4 // c
% K3
th_K3a=[1998 1988 2002 1998 1997];      % crystal thickness in microns of K3 // a*
th_K3b=[1481 1484 1479 1482 1482];      % crystal thickness in microns of K3 // b
th_K3c=[1803 1801 1800 1796 1804];      % crystal thickness in microns of K3 // c
% J1
% ORIENTATION NOTE
% only original labels are used throughout program to avoid confusion
% BUT, now I think a is c; c is b; and b is a
th_J1a=[2465 2460 2467 2470 2452];      % crystal thickness in microns of J1 // a*
th_J1b=[4367 4367 4366 4365 4364];      % crystal thickness in microns of J1 // b
th_J1c=[3218 3236 3190 3232 3231];      % crystal thickness in microns of J1 // c
% PMR P1
th_PMR=[885 867 873 878 879];

%% user input section

% data
r=1;   sig(r,:,:)=load('P_0_unpol.CSV'); thick(r)=mean2(th_PMR); desc(r)={'initial'}; hr(r)=0;
r=r+1; sig(r,:,:)=load('P_1_unpol.CSV'); thick(r)=mean2(th_PMR); desc(r)={'15 m'}; hr(r)=0.25;
r=r+1; sig(r,:,:)=load('P_2_unpol.CSV'); thick(r)=mean2(th_PMR); desc(r)={'30 m'}; hr(r)=0.5;
       sigmix=load('P_2_unpol2.CSV'); 
       sig(r,:,2)=(sig(r,:,2)+sigmix(:,2)')./2;
r=r+1; sig(r,:,:)=load('P_3_unpol.CSV'); thick(r)=mean2(th_PMR); desc(r)={'45 m'}; hr(r)=0.75;
       sigmix1=load('P_3_unpol_sp2.CSV'); 
       sigmix2=load('P_3_unpol_sp3.CSV'); 
       sig(r,:,2)=(sig(r,:,2)+sigmix1(:,2)'+sigmix2(:,2)')./3;
r=r+1; sig(r,:,:)=load('P_4_unpol.CSV'); thick(r)=mean2(th_PMR); desc(r)={'1 hr'};   hr(r)=1;
r=r+1; sig(r,:,:)=load('P_5_unpol.CSV'); thick(r)=mean2(th_PMR); desc(r)={'2 hr'};   hr(r)=2;
       sigmix1=load('P_5_unpol_sp2.CSV'); 
       sigmix2=load('P_5_unpol_2.CSV'); 
       sig(r,:,2)=(sig(r,:,2)+sigmix1(:,2)'+sigmix2(:,2)')./3;
r=r+1; sig(r,:,:)=load('P_6_unpol.CSV'); thick(r)=mean2(th_PMR); desc(r)={'3 hr'};   hr(r)=3;
       sigmix1=load('P_6_unpol_sp2.CSV'); 
       sigmix2=load('P_6_unpol_2.CSV'); 
       sig(r,:,2)=(sig(r,:,2)+sigmix1(:,2)'+sigmix2(:,2)')./3;

ndata=r; %Total number of spectra to be evaluated

% wave number range to plot 
% in bulk water movie
low=3000;
high=4000;
% in peakfit movie
lowFIT=3200;
highFIT=3700;

% label positions
% MakePeakMovie
tagx=3620; % where the labels show up along x direction relative to spectra
xshift=200; % wavenumber shift of labels
yshift=0; % shift labels up or down
tagy=0.5;  % where the labels show up in y relative to the uppermost label
% peak shift movie
plaby=4;

tit='800°C, 1-3 unpolarized spectra each'; % title

% name for autosaving (.txt, .jpg., .avi)
name='PMR_unpol';

% required baseline info
bwni=3150;  % initial wavenumber
bwnf=3700;  % final wavenumber
bwnm=3500;  % middle wavenumber
bam=-0.2;   % how much to deviate from linear baseline for quadratic fitting
nfit=2;     % polynomial order for the baseline fit

% axis limits
% in bulk water movie
top=4; % y axis max;
bottom=0; % y axis min;
% in peakfit movie
topFIT=4; % y axis max;
bottomFIT=0; % y axis min;
% in bulk D figure
topD=3e5; % sqrt(time)/length max
perTopD=0.32; % for x location of labels

% description location in peakfit images
dxpos=3250; 
dypos=topFIT-(0.1*topFIT); 

% number of peaks to look for in peak "deconvolution"
nfitparam=3;
peakshape=2; % 1=Guassian; 2=Lorentzian 
MFE=7; % MaxFunEvals in peakfit.m

%% Make movie, final figure, bulk H area info
[peakarea,remain]=MakePeakMovie(ndata,sig,thick,desc,hr,low,high,tagx,xshift,tagy,tit,bwni,bwnf,bwnm,bam,nfit,name,bottom,top,yshift);

%% Peak fitting 
[DR,Derror]=FitFTIRpeaks(ndata,sig,thick,desc,hr,lowFIT,highFIT,topFIT,bottomFIT,tagx,xshift,tagy,tit,bwni,bwnf,bwnm,bam,nfit,name,nfitparam,peakshape,MFE,dxpos,dypos,plaby);

% issue: peak order in results (DR) is not necessarily consistent!
% solution: always sort resulting bands in order of wavenumber

for k=1:ndata
    DRsort(k,:,:)=sortrows(squeeze(DR(k,:,:)),2);
end

%% fit bulk D
fignum=r+2; % label of figure
titD='PMR-53 bulk water dehydration';
FitD_1D(thick(1),name,hr,remain,topD,perTopD,fignum,titD);

%% fit D for each peak

for k=1:nfitparam
    namepeak=strcat(name,'_peak',num2str(k));
    titD=['PMR-53 peak at ~' num2str(mean2(DRsort(:,k,2)),4) ' cm^{-1}']
    remainpeak=(100*DRsort(:,k,5))./DRsort(1,k,5);
    FitD_1D(thick(1),namepeak,hr,remainpeak',topD,perTopD,k+r+2,titD);
end
