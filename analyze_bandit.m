function output=analyze_bandit(cell,plots,region,give_keyboard)
% os: 0=windows 1=linux 2=mac

if ~exist('plots')
    plots=[];
end

if ~exist('give_keyboard')
    give_keyboard=0;
end

if ~exist('region')
    region=5;
end

try
open_bandit

event5=Event005_1;
event6=Event006_1;


switch unit,
    case 1,
        sig=sig001a;
    case 2,
        sig=sig001b;
    case 3,
        sig=sig001c;
    case 4,
        sig=sig001d;
end

plx_starttime=event5(2);
plx_juicetime=event6(2);
mat_starttime=data(1).trial_start_time;
matlab_offset=mat_starttime-plx_starttime;

%load up vars
i=1; %i counts successful trials
j=1; %j counts broken fixation trials
rwd_last_trial(1)=150;
for x=1:length(data),
    if ~strcmp(data(x).result,'broke fixation')&&(data(x).location~=0)
        trial_start_time(i)=data(x).trial_start_time;
        trial_start_time(i)=trial_start_time(i)-matlab_offset;
        plx_open(i)=trial_start_time(i)-1;
        plx_close(i)=trial_start_time(i)+5;
        spks=sig(find(sig>plx_open(i)&sig<plx_close(i)));
        spks=spks-plx_open(i);
        if ~isempty(spks)
            rasterline=histc(spks,[0:0.01:5]);
            raster(i,:)=rasterline;
        else
            raster(i,:)=zeros(1,numel([0:0.01:5]));
        end
    
        % my code, set up vars
        rwd_all(i,:)=1000*data(x).payouts;
        location(i)=data(x).location;
        if location(i)~=0
            rwd(i)=rwd_all(i,location(i));
            if i>1
                rwd_last_trial(i)=rwd(i-1);
            end
            chose_which(i,location(i))=1;
        else
            rwd(i)=0;
            if i>1
                rwd_last_trial(i)=rwd(i-1);
            end
        end
        rt(i)=data(x).rt;
        juice_ind=find(strcmp(data(x).ev,'juice_rec')==1);
        if ~isempty(juice_ind)
            juice_start(i)=data(x).evt(juice_ind)-rwd(i); %offset to *start* of juice delivery
            juice_end(i)=data(x).evt(juice_ind);
        end
        go_ind=find(strcmp(data(x).ev,'go_cue')==1);
        if ~isempty(go_ind)
            go_time(i)=data(x).evt(go_ind);
        end
        iti_ind=find(strcmp(data(x).ev,'iti_begin')==1);
        if ~isempty(iti_ind)
            iti_time(i)=data(x).evt(iti_ind);
        end
        over_ind=find(strcmp(data(x).ev,'trial_over')==1);
        if ~isempty(over_ind)
            over_time(i)=data(x).evt(over_ind);
        end
                
        i=i+1;
    else
        trial_start_time_bf(j)=data(x).trial_start_time;
        trial_start_time_bf(j)=trial_start_time_bf(j)-matlab_offset;
        plx_open_bf(j)=trial_start_time_bf(j)-1;
        plx_close_bf(j)=trial_start_time_bf(j)+4;
        spks=sig(find(sig>plx_open_bf(j)&sig<plx_close_bf(j)));
        spks=spks-plx_open_bf(j);
        if ~isempty(spks)
            rasterline=histc(spks,[0:0.01:5]);
            raster_bf(j,:)=rasterline;
        end
        
        j=j+1;
    end
end

%define explore/exploit relative to index-guided behavior
[b_advice,b_flips,Q_b]=bayes_advice(rwd_all,location); %greedy behavior
exploit=(b_advice==location);

%define empirically instead:
%define explore vs exploit (need to refine; here, exploit is four
%choices of the same target in a row, defined crudely)
% crit1=zeros(1,numel(rwd));
% crit2=zeros(1,numel(rwd));
% crit3=zeros(1,numel(rwd));
% crit4=zeros(1,numel(rwd));
% 
% for i=1:numel(rwd)
%     if i>3 % previous three
%         crit1(i)=(location(i)==location(i-1))&&(location(i)==location(i-2))&&(location(i)==location(i-3));
%     end
%     if (i>2)&&(i<numel(rwd)) % previous two and one after
%         crit2(i)=(location(i)==location(i-2))&&(location(i)==location(i-1))&&(location(i)==location(i+1));
%     end
%     if (i>1)&&(i<(numel(rwd)-1)) % previous one and two after
%         crit3(i)=(location(i)==location(i-1))&&(location(i)==location(i+1))&&(location(i)==location(i+2));
%     end
% end
% exploit=crit1|crit2;
        
explore=~exploit;

%define alternate explore/exploit according to heuristic
[h_advice,h_flips,Q_h]=heur_advice(rwd_all,location,130); 
for i=1:size(h_advice,1)
    h2_advice(i)=find(h_advice(i,:)==1); % reduce to vector
end
h_exploit=(location==h2_advice);
h_explore=~h_exploit;

%define explore/exploit according to Kalman filter
warnstate=warning('off','all'); %turn off warnings to suppress copious fmincon output
%search options
options=optimset('MaxFunEvals',1200);
options=optimset('Display','off');


fk=@(beta)(-1)*sum(LL_kalman(rwd_all,location,beta));
lb=[0 40 0 0];
ub=[2 260 0.999 260];
[kalman_beta,fval]=fmincon(fk,[0.1 100 0.01 100],[],[],[],[],lb,ub,[],options);

[k_advice,k_flips,Q_k]=kalman_advice(rwd_all,location,kalman_beta); %greedy behavior
k_exploit=(k_advice==location);
k_explore=~k_exploit;

%define explore/exploit according to bias-improved Kalman filter

fk2=@(beta)(-1)*sum(LL_kalman2(rwd_all,location,beta));
lb=[0 40 0 0 0.001 0.001 0.001 0.001];
ub=[2 260 0.999 260 1000 1000 1000 1000];
[kalman2_beta,fval]=fmincon(fk2,[0.1 100 0.01 100 1 1 1 1],[],[],[],[],lb,ub,[],options);

[k2_advice,k2_flips,Q_k2]=kalman_advice(rwd_all,location,kalman2_beta); %greedy behavior
k2_exploit=(k2_advice==location);
k2_explore=~k2_exploit;

% %yet another Kalman filter
% %in an emergency, pull relevant fits from aic.mat, use those to generate
% %vars
% load aic2.mat
% load neuron_list
% 
% %nlist is a list of units to feed to this function
% %blist is a list of indices in nlist corresponding to distinct behavioral
% %sessions
% nlist_ind=find(cell==nlist);
% blist_ind=find(nlist_ind>=blist,1,'last'); %find highest blist index <= this one
% 
% %get the right model from list
% modlist=aic(1).label;
% rightmodel= strcmp(modlist,'kalman3');
% kalman3_beta=aic(blist_ind).beta{rightmodel}; 
% [foo,kalman3_vars]=LL_kalman3(rwd_all,location,kalman3_beta); %retrieve variables

fk3=@(beta)(-1)*sum(LL_kalman3(rwd_all,location,beta));
lb=[0 40 0 0 0.1 0.1 0.01 0.01 0.01];
ub=[2 260 0.999 260 30 30 100 100 100];
[kalman3_beta,fval]=fmincon(fk3,[0.1 100 0.01 100 5 1 1 1 1],[],[],[],[],lb,ub,[],options);



%analyze in terms of Bayesian jump model
fj=@(beta)(-1)*sum(LL_nwhg(rwd_all,location,beta));
lb=[0 1e-5 0.1 1 0 0 0 0];
ub=[1 0.9999 1000 1000 1 1000 0.99 1];
b0=[0.1 0.05 150 10 0.5 5 0.05 0.9];
multmin=1;
numiter=10;

if multmin
    for ind=1:numiter
        bthis=lb+rand(1,numel(ub)).*(ub-lb); %seed random starting point
        [nwhg_beta(ind,:),fval(ind)]=fmincon(fj,bthis,[],[],[],[],lb,ub,[],options);
    end
    [dummy,best_ind]=min(fval); %find best fit
    nwhg_beta=nwhg_beta(best_ind,:); %get best betas
    [foo,nwhg_vars]=LL_nwhg(rwd_all,location,nwhg_beta); %retrieve variables
else     
    [nwhg_beta,fval]=fmincon(fj,b0,[],[],[],[],lb,ub);
end

warning(warnstate); %restore warnings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get expected reward
for i=1:length(rwd)
    b_rwd(i)=Q_b(i,location(i));
    h_rwd(i)=Q_h(i,location(i));
    k_rwd(i)=Q_k(i,location(i));
    k2_rwd(i)=Q_k2(i,location(i));
end

%define strategy switch trials
switchup(1)=0;
stratswitch(1)=0;
for i=2:numel(rwd)
    switchup(i)=(location(i)~=location(i-1));
    stratswitch(i)=(exploit(i)~=exploit(i-1));
end

% control trials for target switch vs explore/exploit
lore_switch=explore&switchup;
lore_noswitch=explore&~switchup;
hlore_switch=h_explore&switchup;
hlore_noswitch=h_explore&~switchup;
klore_switch=k_explore&switchup;
klore_noswitch=k_explore&~switchup;
switch_lore=switchup&explore;
switch_loit=switchup&exploit;
switch_hlore=switchup&h_explore;
switch_hloit=switchup&h_exploit;
switch_klore=switchup&k_explore;
switch_kloit=switchup&k_exploit;
loit_noswitch=exploit&~switchup;
hloit_noswitch=h_exploit&~switchup;
kloit_noswitch=k_exploit&~switchup;

%define time axis
t_axis=-1000:10:4000;
t_rwdalign=-2000:10:2000;

%for rwd-aligned psth's
raster_rwd_align=nan(size(raster,1),401);

for i=1:size(juice_start,2)
    before_juice=find(t_axis<=juice_start(i));
    align_bins(i)=before_juice(end);
    if align_bins(i)<201
        pre=1:align_bins(i)-1; %2s pre interval
    else
        pre=(align_bins(i)-200):(align_bins(i)-1);
    end
    if align_bins(i)>301
        post=(align_bins(i)+1):501; %2s post interval
    else
        post=(align_bins(i)+1):(align_bins(i)+200);
    end
    
    t_rwd_axes=[pre align_bins(i) post];
    offset=201-align_bins(i);
    raster_rwd_align(i,t_rwd_axes+offset)=raster(i,t_rwd_axes);
    if size(raster_rwd_align,2)>401
        keyboard
    end
end

raster_rwd_align2=nan(size(raster,1),401);

for i=1:size(juice_start,2)
    before_juice2=find(t_axis<=juice_end(i));
    align_bins2(i)=before_juice2(end);
    if align_bins2(i)<201
        pre=1:align_bins2(i)-1; %2s pre interval
    else
        pre=(align_bins2(i)-200):(align_bins2(i)-1);
    end
    if align_bins2(i)>301
        post=(align_bins2(i)+1):501; %2s post interval
    else
        post=(align_bins2(i)+1):(align_bins2(i)+200);
    end
    
    t_rwd_axes2=[pre align_bins2(i) post];
    offset=201-align_bins2(i);
    raster_rwd_align2(i,t_rwd_axes2+offset)=raster(i,t_rwd_axes2);
    if size(raster_rwd_align2,2)>401
        keyboard
    end
end

%additional data manipulations
juice_start=juice_start(find(juice_start<2000)); %clear out the worst outliers
juice_end=juice_end(find(juice_end<2300));
mean_juice_start=mean(juice_start);
std_juice_start=std(juice_start);
mean_juice_end=mean(juice_end);
std_juice_end=std(juice_end);
go_time=go_time(find(go_time<2000));
mean_go_time=mean(go_time);
iti_time=iti_time(find(iti_time<3000));
mean_iti_time=mean(iti_time); %mean time of *iti_begin*
over_time=over_time(find(iti_time<5000));
mean_over_time=mean(over_time);

%make data selections

start_bin=find(t_axis==0);
juice_start_bin=find(t_axis<=(mean_juice_start-std_juice_start));
juice_start_bin=juice_start_bin(end);
juice_end_bin=find(t_axis<=(mean_juice_end+std_juice_end));
juice_end_bin=juice_end_bin(end);
iti_begin_bin=find(t_axis<=mean_iti_time);
iti_begin_bin=iti_begin_bin(end);
end_bin=find(t_axis<=mean_over_time);
end_bin=end_bin(end);

R_prior=1:start_bin;
R_decision=start_bin:juice_start_bin;
R_juice_to_iti=juice_end_bin:iti_begin_bin;
R_iti=iti_begin_bin:end_bin;
R_pre=union(R_prior,R_decision);
R_post=union(R_juice_to_iti,R_iti);
R_tot=1:end_bin;
R_rwd=(juice_start_bin-10):(juice_end_bin+10);

roi{1}=R_prior;
roi{2}=R_decision;
roi{3}=R_juice_to_iti;
roi{4}=R_iti;
roi{5}=R_pre;
roi{6}=R_post;
roi{7}=R_tot;
roi{8}=R_rwd;

R=roi{region};

%write some things out to file
output.basic.h_explore=h_explore;
output.basic.k_explore=k_explore;
output.basic.k2_explore=k2_explore;
output.basic.explore=explore;
output.basic.switch=switchup;
output.basic.raster=raster;
output.basic.raster_rwd_align=raster_rwd_align;
output.basic.raster_bf=raster_bf;
output.basic.rwd=rwd;
output.basic.rwd_all=rwd_all;
output.basic.b_rwd=b_rwd;
output.basic.h_rwd=h_rwd;
output.basic.k_rwd=k_rwd;
output.basic.k2_rwd=k2_rwd;
output.basic.rt=rt;
output.basic.rwd_last_trial=rwd_last_trial;
output.basic.location=location;
output.basic.roi=roi;
output.basic.rasttimes=[mean_go_time, mean_juice_start, mean_juice_end, mean_iti_time, mean_over_time, std_juice_start, std_juice_end];
output.basic.kalman3=kalman3_vars;
output.basic.nwhg=nwhg_vars;

%smaller epochs (100ms bins)
sides=1:10:end_bin;
for i=1:(numel(sides)-1)
    epoch{i}=sides(i):(sides(i+1)-1);
end


%behavioral analysis vars
k=length(data(1,1).payouts);
numtrials=length(rwd);

max_reward=260;
min_reward=40;
mean_reward=150;
central_tendency=0.05;
step_std=10;
rwd_axis=min_reward:5:max_reward;

%load colors
colors{1}=[255 0 0];
colors{2}=[255 0 255];
colors{3}=[255 255 0];
colors{4}=[128 0 255];
colors{5}=[0 0 255];
colors{6}=[255 128 128]; 
colors{7}=[0 255 128];
colors{8}=[255 128 0];
gray=[0.5 0.5 0.5];

plot_colors=[];
for i=1:k
    plot_colors=[plot_colors;colors{i}];
end
plot_colors=[plot_colors; [255 255 255]];
plot_colors=plot_colors./255; %translate to Matlab color intensities

if give_keyboard
    keyboard
end
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(plots)
    if ~isempty(find(9==plots))
        numpanes=ceil((numel(plots)-1)/4);
    else
        numpanes=ceil(numel(plots)/4);
    end
    currpane=1;
    currplot=1;
    
    for i=1:numpanes
        fig_handle(i)=figure;
    end
end 


%%%%%%%%%% plot 1 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing: explore and exploit
q1=find(exploit==1);
q2=find(exploit==0);
c1=100*mean(raster(q1,:));
c2=100*mean(raster(q2,:));
q1_next=q1-1; %list of trials *before* exploit trials
q1_next=q1_next(q1_next~=0); % removes 0 
q2_next=q2-1; 
q2_next=q2_next(q2_next~=0);

%t-tests
for i=1:size(roi,2)
%     [h,pval1(1,i)]=ttest2(c1(roi{i}),c2(roi{i}));
    sample=100*mean(raster(:,roi{i}),2);
    [h,pval1(1,i)]=ttest2(sample(q1),sample(q2));
    [h,pval1(2,i),ksstat]=kstest2(c1(roi{i}),c2(roi{i}));
    cv(1,i)=std(c1(roi{i}))/mean(c1(roi{i}));
    cv(2,i)=std(c2(roi{i}))/mean(c2(roi{i}));
end

% does firing this trial distinguish between lore and loit *next* trial
sample=100*mean(raster(:,roi{6}),2);
[h,ppred]=ttest2(sample(q1_next),sample(q2_next));

for i=1:size(epoch,2) % epoch p-vals (much smaller time windows)
    sample=100*mean(raster(:,epoch{i}),2);
    [h,pval_ep(1,i)]=ttest2(sample(q1),sample(q2));
end

output.plot1.tpvals=pval1(1,:);
output.plot1.ppred=ppred;
output.plot1.ksvals=pval1(2,:);
output.plot1.tcv=cv(1,:);
output.plot1.kscv=cv(2,:);
output.plot1.loitpsth=c1;
output.plot1.lorepsth=c2;
output.plot1.ep_pvals=pval_ep(1,:);

if ~isempty(find(1==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    hold on
    plot(t_axis,gsmooth(c1,5),'r')
    plot(t_axis,gsmooth(c2,5),'k')
    plot([t_axis(1) t_axis(end)],[mean(c1(R)) mean(c1(R))],'r')
    plot([t_axis(1) t_axis(end)],[mean(c2(R)) mean(c2(R))],'k')
    vert=ylim;
    vert=vert(2);
    plot([0 0],[0 vert],'linestyle','--','color','k')
    text(0,-1.5,'Begin','HorizontalAlignment','center')
    plot([mean_go_time mean_go_time],[0 vert],'linestyle','--','color','g');
    text(mean_go_time,-1.5,'go','color','g','HorizontalAlignment','center')
    plot([mean_juice_start mean_juice_start],[0 vert],'r')
    text(mean_juice_start,-1.5,'juice','color','r','HorizontalAlignment','center')
    plot([mean_juice_start-std_juice_start mean_juice_start-std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_juice_start+std_juice_start mean_juice_start+std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_iti_time mean_iti_time],[0 vert],'linestyle','--','color',gray)
    text(mean_iti_time,-1.5,'Target off','color','k','HorizontalAlignment','center')
    plot([mean_over_time mean_over_time],[0 vert],'linestyle','--','color','k')
    text((mean_iti_time+mean_over_time)/2,-1.5,'ITI','color','k','HorizontalAlignment','center')
    text(mean_over_time,-1.5,'End','color','k','HorizontalAlignment','center')
    legend({'exploit','explore'},'Location','SouthWest')
    title(sprintf('p=%f %% explore=%f',pval1(1,region),mean(explore)))
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 14 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing: explore and exploit (heuristic definition)
clear q1 q2 c1 c2
q1=find(h_exploit==1);
q2=find(h_exploit==0);
c1=100*mean(raster(q1,:));
c2=100*mean(raster(q2,:));
q1_next=q1-1; %list of trials *before* exploit trials
q1_next=q1_next(q1_next~=0); % removes 0 
q2_next=q2-1; 
q2_next=q2_next(q2_next~=0);


%t-tests
for i=1:size(roi,2)
    sample=100*mean(raster(:,roi{i}),2);
    [h,pval1(4,i)]=ttest2(sample(q1),sample(q2));
%     [h,pval1(2,i),ksstat]=kstest2(c1(roi{i}),c2(roi{i}));
%     cv(1,i)=std(c1(roi{i}))/mean(c1(roi{i}));
%     cv(2,i)=std(c2(roi{i}))/mean(c2(roi{i}));
end

% does firing this trial distinguish between lore and loit *next* trial
sample=100*mean(raster(:,roi{6}),2);
[h,ppred]=ttest2(sample(q1_next),sample(q2_next));

for i=1:size(epoch,2)
    sample=100*mean(raster(:,epoch{i}),2);
    [h,pval_ep(2,i)]=ttest2(sample(q1),sample(q2));
end

output.plot14.tpvals=pval1(4,:);
output.plot14.ppred=ppred;
output.plot14.loitpsth=c1;
output.plot14.lorepsth=c2;
output.plot14.ep_pvals=pval_ep(2,:);

if ~isempty(find(14==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    hold on
    plot(t_axis,gsmooth(c1,5),'r')
    plot(t_axis,gsmooth(c2,5),'k')
    plot([t_axis(1) t_axis(end)],[mean(c1(R)) mean(c1(R))],'r')
    plot([t_axis(1) t_axis(end)],[mean(c2(R)) mean(c2(R))],'k')
    vert=ylim;
    vert=vert(2);
    plot([0 0],[0 vert],'linestyle','--','color','k')
    text(0,-1.5,'Begin','HorizontalAlignment','center')
    plot([mean_go_time mean_go_time],[0 vert],'linestyle','--','color','g');
    text(mean_go_time,-1.5,'go','color','g','HorizontalAlignment','center')
    plot([mean_juice_start mean_juice_start],[0 vert],'r')
    text(mean_juice_start,-1.5,'juice','color','r','HorizontalAlignment','center')
    plot([mean_juice_start-std_juice_start mean_juice_start-std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_juice_start+std_juice_start mean_juice_start+std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_iti_time mean_iti_time],[0 vert],'linestyle','--','color',gray)
    text(mean_iti_time,-1.5,'Target off','color','k','HorizontalAlignment','center')
    plot([mean_over_time mean_over_time],[0 vert],'linestyle','--','color','k')
    text((mean_iti_time+mean_over_time)/2,-1.5,'ITI','color','k','HorizontalAlignment','center')
    text(mean_over_time,-1.5,'End','color','k','HorizontalAlignment','center')
    legend({'h_exploit','h_explore'},'Location','SouthWest')
    title(sprintf('p=%f %% explore=%f',pval1(4,region),mean(h_explore)))
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end


%%%%%%%%%% plot 20 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing: explore and exploit (kalman definition)
clear q1 q2 c1 c2
q1=find(k_exploit==1);
q2=find(k_exploit==0);
c1=100*mean(raster(q1,:));
c2=100*mean(raster(q2,:));
q1_next=q1-1; %list of trials *before* exploit trials
q1_next=q1_next(q1_next~=0); % removes 0 
q2_next=q2-1; 
q2_next=q2_next(q2_next~=0);


%t-tests
for i=1:size(roi,2)
    sample=100*mean(raster(:,roi{i}),2);
    [h,pval1(4,i)]=ttest2(sample(q1),sample(q2));
%     [h,pval1(2,i),ksstat]=kstest2(c1(roi{i}),c2(roi{i}));
%     cv(1,i)=std(c1(roi{i}))/mean(c1(roi{i}));
%     cv(2,i)=std(c2(roi{i}))/mean(c2(roi{i}));
end

% does firing this trial distinguish between lore and loit *next* trial
sample=100*mean(raster(:,roi{6}),2);
[h,ppred]=ttest2(sample(q1_next),sample(q2_next));

for i=1:size(epoch,2)
    sample=100*mean(raster(:,epoch{i}),2);
    [h,pval_ep(2,i)]=ttest2(sample(q1),sample(q2));
end

output.plot20.tpvals=pval1(4,:);
output.plot20.ppred=ppred;
output.plot20.loitpsth=c1;
output.plot20.lorepsth=c2;
output.plot20.ep_pvals=pval_ep(2,:);

if ~isempty(find(20==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    hold on
    plot(t_axis,gsmooth(c1,5),'r')
    plot(t_axis,gsmooth(c2,5),'k')
    plot([t_axis(1) t_axis(end)],[mean(c1(R)) mean(c1(R))],'r')
    plot([t_axis(1) t_axis(end)],[mean(c2(R)) mean(c2(R))],'k')
    vert=ylim;
    vert=vert(2);
    plot([0 0],[0 vert],'linestyle','--','color','k')
    text(0,-1.5,'Begin','HorizontalAlignment','center')
    plot([mean_go_time mean_go_time],[0 vert],'linestyle','--','color','g');
    text(mean_go_time,-1.5,'go','color','g','HorizontalAlignment','center')
    plot([mean_juice_start mean_juice_start],[0 vert],'r')
    text(mean_juice_start,-1.5,'juice','color','r','HorizontalAlignment','center')
    plot([mean_juice_start-std_juice_start mean_juice_start-std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_juice_start+std_juice_start mean_juice_start+std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_iti_time mean_iti_time],[0 vert],'linestyle','--','color',gray)
    text(mean_iti_time,-1.5,'Target off','color','k','HorizontalAlignment','center')
    plot([mean_over_time mean_over_time],[0 vert],'linestyle','--','color','k')
    text((mean_iti_time+mean_over_time)/2,-1.5,'ITI','color','k','HorizontalAlignment','center')
    text(mean_over_time,-1.5,'End','color','k','HorizontalAlignment','center')
    legend({'h_exploit','h_explore'},'Location','SouthWest')
    title(sprintf('p=%f %% explore=%f',pval1(4,region),mean(h_explore)))
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 16 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing: explore and exploit (aligned to reward onset)
clear q1 q2 c1 c2 d1 d2 d3 d4
q1=find(exploit==1);
q2=find(exploit==0);
q3=find(h_exploit==1);
q4=find(h_exploit==0);
q5=find(k_exploit==1);
q6=find(k_exploit==0);
c1=100*nanmean(raster_rwd_align(q1,:));  % aligned to reward onset
c2=100*nanmean(raster_rwd_align(q2,:));
c3=100*nanmean(raster_rwd_align2(q1,:)); % aligned to reward offset
c4=100*nanmean(raster_rwd_align2(q2,:));
hc1=100*nanmean(raster_rwd_align(q3,:));  % aligned to reward onset
hc2=100*nanmean(raster_rwd_align(q4,:));
hc3=100*nanmean(raster_rwd_align2(q3,:)); % aligned to reward offset
hc4=100*nanmean(raster_rwd_align2(q4,:));
kc1=100*nanmean(raster_rwd_align(q5,:));  % aligned to reward onset
kc2=100*nanmean(raster_rwd_align(q6,:));
kc3=100*nanmean(raster_rwd_align2(q5,:)); % aligned to reward offset
kc4=100*nanmean(raster_rwd_align2(q6,:));

r_pre=1:200;
r_post=201:401;
if region==5
    R_align=r_pre;
elseif region==6
    R_align=r_post;
else
    R_align=200;
end
d1=100*nanmean(raster_rwd_align(q1,R_align),2);
d2=100*nanmean(raster_rwd_align(q2,R_align),2);
d3=100*nanmean(raster_rwd_align2(q1,R_align),2);
d4=100*nanmean(raster_rwd_align2(q2,R_align),2);
hd1=100*nanmean(raster_rwd_align(q3,R_align),2);
hd2=100*nanmean(raster_rwd_align(q4,R_align),2);
hd3=100*nanmean(raster_rwd_align2(q3,R_align),2);
hd4=100*nanmean(raster_rwd_align2(q4,R_align),2);
kd1=100*nanmean(raster_rwd_align(q5,R_align),2);
kd2=100*nanmean(raster_rwd_align(q6,R_align),2);
kd3=100*nanmean(raster_rwd_align2(q5,R_align),2);
kd4=100*nanmean(raster_rwd_align2(q6,R_align),2);

%t-tests
% pre-juice
[h,pval_align(1)]=ttest2(mean(raster_rwd_align(q1,1:200),2),mean(raster_rwd_align(q2,1:200),2));
[h,pval_align(3)]=ttest2(mean(raster_rwd_align(q3,1:200),2),mean(raster_rwd_align(q4,1:200),2));
[h,pval_align(5)]=ttest2(mean(raster_rwd_align(q5,1:200),2),mean(raster_rwd_align(q6,1:200),2));
% post-juice
[h,pval_align(2)]=ttest2(mean(raster_rwd_align(q1,202:401),2),mean(raster_rwd_align(q2,202:401),2));
[h,pval_align(4)]=ttest2(mean(raster_rwd_align(q3,202:401),2),mean(raster_rwd_align(q4,202:401),2));
[h,pval_align(6)]=ttest2(mean(raster_rwd_align(q6,202:401),2),mean(raster_rwd_align(q6,202:401),2));

output.plot16.tpvals=pval_align;
output.plot16.loitpsth=[c1;hc1;kc1];   % aligned to onset
output.plot16.lorepsth=[c2;hc2;kc2];
output.plot16.loitpsth2=[c3;hc3;kc3];  % aligned to offset
output.plot16.lorepsth2=[c4;hc4;kc4];
output.plot16.loitrates={d1,hd1,kd1};
output.plot16.lorerates={d2,hd2,kd2};
output.plot16.loitrates2={d3,hd3,kd3};
output.plot16.lorerates2={d4,hd4,kd4};

if ~isempty(find(16==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    hold on
    plot(t_rwdalign,gsmooth(c1,5),'r')
    plot(t_rwdalign,gsmooth(c2,5),'k')
    plot([t_rwdalign(1) t_rwdalign(end)],[mean(c1) mean(c1)],'r')
    plot([t_rwdalign(1) t_rwdalign(end)],[mean(c2) mean(c2)],'k')
    vert=ylim;
    vert=vert(2);
    
    plot([0 0],[0 vert],'r')
    text(0,-1.5,'Reward','color','r','HorizontalAlignment','center')
    
    legend({'exploit','explore'},'Location','SouthWest')
%     title(sprintf('p=%f %% explore=%f',pval1(5,region),mean(h_explore)))
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 10 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing: strategy switches
q1=find(switchup==1);
q2=find(switchup==0);
c1=100*mean(raster(q1,:));
c2=100*mean(raster(q2,:));

%t-tests
for i=1:size(roi,2)
    sample=100*mean(raster(:,roi{i}),2);
    [h,pval1(3,i)]=ttest2(sample(q1),sample(q2));
    %[h,pval1(3,i)]=ttest2(c1(roi{i}),c2(roi{i}));
%     [h,pval10(2,i),ksstat]=kstest2(c1(roi{i}),c2(roi{i}));
%     cv(1,i)=std(c1(roi{i}))/mean(c1(roi{i}));
%     cv(2,i)=std(c2(roi{i}))/mean(c2(roi{i}));
end

output.plot10.tpvals=pval1(3,:);

if ~isempty(find(10==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    hold on
    plot(t_axis,gsmooth(c1,1),'r')
    plot(t_axis,gsmooth(c2,1),'k')
    plot([t_axis(1) t_axis(end)],[mean(c1(R)) mean(c1(R))],'r')
    plot([t_axis(1) t_axis(end)],[mean(c2(R)) mean(c2(R))],'k')
    vert=ylim;
    vert=vert(2);
    plot([0 0],[0 vert],'linestyle','--','color','k')
    text(0,-1.5,'Begin','HorizontalAlignment','center')
    plot([mean_go_time mean_go_time],[0 vert],'linestyle','--','color','g');
    text(mean_go_time,-1.5,'go','color','g','HorizontalAlignment','center')
    plot([mean_juice_start mean_juice_start],[0 vert],'r')
    text(mean_juice_start,-1.5,'juice','color','r','HorizontalAlignment','center')
    plot([mean_juice_start-std_juice_start mean_juice_start-std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_juice_start+std_juice_start mean_juice_start+std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_iti_time mean_iti_time],[0 vert],'linestyle','--','color',gray)
    text(mean_iti_time,-1.5,'Target off','color','k','HorizontalAlignment','center')
    plot([mean_over_time mean_over_time],[0 vert],'linestyle','--','color','k')
    text((mean_iti_time+mean_over_time)/2,-1.5,'ITI','color','k','HorizontalAlignment','center')
    text(mean_over_time,-1.5,'End','color','k','HorizontalAlignment','center')
    legend({'switch','perseverate'},'Location','SouthWest')
    title(sprintf('p=%f %% switchup=%f',pval1(3,region),mean(switchup)))
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 17 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing: explore disentangled from switching
clear q1 q2 c1 c2
q1=find(lore_switch==1);
q2=find(lore_noswitch==1);
q3=find(hlore_switch==1);
q4=find(hlore_noswitch==1);
q5=find(switch_lore==1);
q6=find(switch_loit==1);
q7=find(switch_hlore==1);
q8=find(switch_hloit==1);
q9=find(loit_noswitch==1);
q10=find(hloit_noswitch==1);
q11=find(klore_switch==1);
q12=find(klore_noswitch==1);
q13=find(switch_klore==1);
q14=find(switch_kloit==1);
q15=find(kloit_noswitch==1);

c1=100*nanmean(raster(q1,:),1);
c2=100*nanmean(raster(q2,:),1);
c3=100*nanmean(raster(q3,:),1);
c4=100*nanmean(raster(q4,:),1);
c5=100*nanmean(raster(q11,:),1);
c6=100*nanmean(raster(q12,:),1);

d1=100*nanmean(raster_rwd_align(q1,:),1);
d2=100*nanmean(raster_rwd_align(q2,:),1);
d3=100*nanmean(raster_rwd_align(q3,:),1);
d4=100*nanmean(raster_rwd_align(q4,:),1);
d5=100*nanmean(raster_rwd_align(q11,:),1);
d6=100*nanmean(raster_rwd_align(q12,:),1);

e1=100*nanmean(raster(q5,:),1);
e2=100*nanmean(raster(q6,:),1);
e3=100*nanmean(raster(q7,:),1);
e4=100*nanmean(raster(q8,:),1);
e5=100*nanmean(raster(q13,:),1);
e6=100*nanmean(raster(q14,:),1);

f1=100*nanmean(raster_rwd_align(q5,:),1);
f2=100*nanmean(raster_rwd_align(q6,:),1);
f3=100*nanmean(raster_rwd_align(q7,:),1);
f4=100*nanmean(raster_rwd_align(q8,:),1);
f5=100*nanmean(raster_rwd_align(q13,:),1);
f6=100*nanmean(raster_rwd_align(q14,:),1);

g1=100*nanmean(raster(q9,:),1);
g2=100*nanmean(raster(q10,:),1);
g3=100*nanmean(raster_rwd_align(q9,:),1);
g4=100*nanmean(raster_rwd_align(q10,:),1);
g5=100*nanmean(raster(q15,:),1);
g6=100*nanmean(raster_rwd_align(q15,:),1);

%t-tests
for i=1:size(roi,2)
    [h,pval1(6,i)]=ttest2(c1(roi{i}),c2(roi{i}));
    [h,pval1(7,i)]=ttest2(c3(roi{i}),c4(roi{i}));
    [h,pval1(10,i)]=ttest2(c5(roi{i}),c6(roi{i}));
    
    %take means over rois and do t-tests over trials
    m1=100*nanmean(raster(q1,roi{i}),2);
    m2=100*nanmean(raster(q2,roi{i}),2);
    m3=100*nanmean(raster(q3,roi{i}),2);
    m4=100*nanmean(raster(q4,roi{i}),2);
    m5=100*nanmean(raster(q11,roi{i}),2);
    m6=100*nanmean(raster(q12,roi{i}),2);
    [h,pval1(8,i)]=ttest2(m1,m2);
    [h,pval1(9,i)]=ttest2(m3,m4);
    [h,pval1(11,i)]=ttest2(m5,m6);
end

% mean firing rates for each trial in epochs of interest
aa1=100*nanmean(raster(:,R),2);
aa2=100*nanmean(raster_rwd_align(:,R_align),2);
cc1=100*nanmean(raster(q1,R),2);
cc2=100*nanmean(raster(q2,R),2);
cc3=100*nanmean(raster(q3,R),2);
cc4=100*nanmean(raster(q4,R),2);
cc5=100*nanmean(raster(q11,R),2);
cc6=100*nanmean(raster(q12,R),2);
dd1=100*nanmean(raster_rwd_align(q1,R_align),2);
dd2=100*nanmean(raster_rwd_align(q2,R_align),2);
dd3=100*nanmean(raster_rwd_align(q3,R_align),2);
dd4=100*nanmean(raster_rwd_align(q4,R_align),2);
dd5=100*nanmean(raster_rwd_align(q11,R_align),2);
dd6=100*nanmean(raster_rwd_align(q12,R_align),2);
ee1=100*nanmean(raster(q5,R),2);
ee2=100*nanmean(raster(q6,R),2);
ee3=100*nanmean(raster(q7,R),2);
ee4=100*nanmean(raster(q8,R),2);
ee5=100*nanmean(raster(q13,R),2);
ee6=100*nanmean(raster(q14,R),2);
ff1=100*nanmean(raster_rwd_align(q5,R_align),2);
ff2=100*nanmean(raster_rwd_align(q6,R_align),2);
ff3=100*nanmean(raster_rwd_align(q7,R_align),2);
ff4=100*nanmean(raster_rwd_align(q8,R_align),2);
ff5=100*nanmean(raster_rwd_align(q13,R_align),2);
ff6=100*nanmean(raster_rwd_align(q14,R_align),2);
gg1=100*nanmean(raster(q9,R),2);
gg2=100*nanmean(raster(q10,R),2);
gg3=100*nanmean(raster_rwd_align(q9,R_align),2);
gg4=100*nanmean(raster_rwd_align(q10,R_align),2);
gg5=100*nanmean(raster(q15,R),2);
gg6=100*nanmean(raster_rwd_align(q15,R_align),2);

output.plot17.tpvals=pval1(6,:);
output.plot17.tpvals_alt=pval1(8,:);
output.plot17.thpvals=pval1(7,:);
output.plot17.thpvals_alt=pval1(9,:);
output.plot17.tkpvals=pval1(10,:);
output.plot17.tkpvals_alt=pval1(11,:);
output.plot17.loreswitch_psths=[c1;c2;c3;c4;c5;c6];
output.plot17.loreswitch_psths_align=[d1;d2;d3;d4;d5;d6];
output.plot17.switchlore_psths=[e1;e2;e3;e4;e5;e6];
output.plot17.switchlore_psths_align=[f1;f2;f3;f4;f5;f6];
output.plot17.loitnoswitch_psths=[g1;g2;g5];
output.plot17.loitnoswitch_psths_align=[g3;g4;g6];
output.plot17.all_means=aa1;
output.plot17.all_means_align=aa2;
output.plot17.loreswitch_means.cc1=cc1;
output.plot17.loreswitch_means.cc2=cc2;
output.plot17.loreswitch_means.cc3=cc3;
output.plot17.loreswitch_means.cc4=cc4;
output.plot17.loreswitch_means.cc5=cc5;
output.plot17.loreswitch_means.cc6=cc6;
output.plot17.loreswitch_means_align.dd1=dd1; 
output.plot17.loreswitch_means_align.dd2=dd2; 
output.plot17.loreswitch_means_align.dd3=dd3; 
output.plot17.loreswitch_means_align.dd4=dd4; 
output.plot17.loreswitch_means_align.dd5=dd5; 
output.plot17.loreswitch_means_align.dd6=dd6; 
output.plot17.switchlore_means.ee1=ee1;
output.plot17.switchlore_means.ee2=ee2;
output.plot17.switchlore_means.ee3=ee3;
output.plot17.switchlore_means.ee4=ee4;
output.plot17.switchlore_means.ee5=ee5;
output.plot17.switchlore_means.ee6=ee6;
output.plot17.switchlore_means_align.ff1=ff1;
output.plot17.switchlore_means_align.ff2=ff2;
output.plot17.switchlore_means_align.ff3=ff3;
output.plot17.switchlore_means_align.ff4=ff4;
output.plot17.switchlore_means_align.ff5=ff5;
output.plot17.switchlore_means_align.ff6=ff6;
output.plot17.loitnoswitch_means.gg1=gg1;
output.plot17.loitnoswitch_means.gg2=gg2;
output.plot17.loitnoswitch_means.gg3=gg3;
output.plot17.loitnoswitch_means.gg4=gg4;
output.plot17.loitnoswitch_means.gg5=gg5;
output.plot17.loitnoswitch_means.gg6=gg6;

if ~isempty(find(17==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    hold on
    plot(t_axis,gsmooth(c1,5),'r')
    plot(t_axis,gsmooth(c2,5),'k')
    plot(t_axis,gsmooth(c3,5),'g')
    plot(t_axis,gsmooth(c4,5),'b')
    plot([t_axis(1) t_axis(end)],[mean(c1(R)) mean(c1(R))],'r')
    plot([t_axis(1) t_axis(end)],[mean(c2(R)) mean(c2(R))],'k')
    plot([t_axis(1) t_axis(end)],[mean(c3(R)) mean(c3(R))],'g')
    plot([t_axis(1) t_axis(end)],[mean(c4(R)) mean(c4(R))],'b')
    vert=ylim;
    vert=vert(2);
    plot([0 0],[0 vert],'linestyle','--','color','k')
    text(0,-1.5,'Begin','HorizontalAlignment','center')
    plot([mean_go_time mean_go_time],[0 vert],'linestyle','--','color','g');
    text(mean_go_time,-1.5,'go','color','g','HorizontalAlignment','center')
    plot([mean_juice_start mean_juice_start],[0 vert],'r')
    text(mean_juice_start,-1.5,'juice','color','r','HorizontalAlignment','center')
    plot([mean_juice_start-std_juice_start mean_juice_start-std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_juice_start+std_juice_start mean_juice_start+std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_iti_time mean_iti_time],[0 vert],'linestyle','--','color',gray)
    text(mean_iti_time,-1.5,'Target off','color','k','HorizontalAlignment','center')
    plot([mean_over_time mean_over_time],[0 vert],'linestyle','--','color','k')
    text((mean_iti_time+mean_over_time)/2,-1.5,'ITI','color','k','HorizontalAlignment','center')
    text(mean_over_time,-1.5,'End','color','k','HorizontalAlignment','center')
    legend({'lore switch','lore noswitch','hlore switch','hlore noswitch'},'Location','SouthWest')
    title(sprintf('p=%f p_h=%f',pval1(6,region),pval1(7,region)))
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 2 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing: completed and broke fixation trials
c1=100*mean(raster);
c2=100*mean(raster_bf);

%t-tests
for i=1:size(roi,2)
    [h,pval2(i)]=ttest2(mean(raster,2),mean(raster_bf,2));
%     [h,pval2(i)]=ttest2(c1(roi{i}),c2(roi{i}));
end

output.plot2.tpvals=pval2;

if ~isempty(find(2==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    hold on
    plot(t_axis, gsmooth(c1,1),'r')
    plot(t_axis, gsmooth(c2,1),'k')
    plot([t_axis(1) t_axis(end)],[mean(c1(R)) mean(c1(R))],'r')
    plot([t_axis(1) t_axis(end)],[mean(c2(R)) mean(c2(R))],'k')
    vert=ylim;
    vert=vert(2);
    plot([0 0],[0 vert],'linestyle','--','color','k')
%     plot([t_axis(R(end)) t_axis(R(end))],[0 vert],'linestyle','--','color','k')
    plot([mean_go_time mean_go_time],[0 vert],'linestyle','--','color','g');
    plot([mean_juice_start mean_juice_start],[0 vert],'r')
    plot([mean_juice_start-std_juice_start mean_juice_start-std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_juice_start+std_juice_start mean_juice_start+std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_iti_time mean_iti_time],[0 vert],'linestyle','--','color',gray)
    plot([mean_over_time mean_over_time],[0 vert],'linestyle','--','color','k')
    legend({'saccade complete','broke fixation'})
    title(pval2(region))
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 3 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing vs target chosen
q1=find(chose_which(:,1)==1);
q2=find(chose_which(:,2)==1);
q3=find(chose_which(:,3)==1);
q4=find(chose_which(:,4)==1);
c1=100*mean(raster(q1,:),1);
c2=100*mean(raster(q2,:),1);
c3=100*mean(raster(q3,:),1);
c4=100*mean(raster(q4,:),1);
fr_bytarg=100*mean(raster,2);

targ_tune=[mean(c1(R)) mean(c2(R)) mean(c3(R)) mean(c4(R))];
[pval3,table3,stats3]=anova1(fr_bytarg,location,'off');

output.plot3.tuning=targ_tune;
output.plot3.pval=pval3;
output.plot3.table=table3;
output.plot3.stats=stats3;

if ~isempty(find(3==plots))
    figure(fig_handle(currpane))
    targ_axes=subplot(2,2,currplot);
    
    hold on
    set(targ_axes,'ColorOrder',plot_colors);
    plot(t_axis,gsmooth(c1,1),'color',plot_colors(1,:))
    plot(t_axis,gsmooth(c2,1),'color',plot_colors(2,:))
    plot(t_axis,gsmooth(c3,1),'color',plot_colors(3,:))
    plot(t_axis,gsmooth(c4,1),'color',plot_colors(4,:))
   
    plot([t_axis(1) t_axis(end)],[mean(c1(R)) mean(c1(R))],'color',plot_colors(1,:))
    plot([t_axis(1) t_axis(end)],[mean(c2(R)) mean(c2(R))],'color',plot_colors(2,:))
    plot([t_axis(1) t_axis(end)],[mean(c3(R)) mean(c3(R))],'color',plot_colors(3,:))
    plot([t_axis(1) t_axis(end)],[mean(c4(R)) mean(c4(R))],'color',plot_colors(4,:))
    vert=ylim;
    vert=vert(2);
    plot([0 0],[0 vert],'linestyle','--','color','k')
%     plot([t_axis(R(end)) t_axis(R(end))],[0 vert],'linestyle','--','color','k')
    plot([mean_go_time mean_go_time],[0 vert],'linestyle','--','color','g');
    plot([mean_juice_start mean_juice_start],[0 vert],'r')
    plot([mean_juice_start-std_juice_start mean_juice_start-std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_juice_start+std_juice_start mean_juice_start+std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_iti_time mean_iti_time],[0 vert],'linestyle','--','color',gray)
    plot([mean_over_time mean_over_time],[0 vert],'linestyle','--','color','k')
    legend({'T1','T2','T3','T4'},'Location','SouthWest')
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end



%%%%%%%%%% plot 4 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing vs reward on last trial (or this trial; region-dependent)
clear c;
clear sem;
clear raw_c;
clear raw_x;
raw_c=[];
raw_x=[];
last_regions=[1,2,5];
this_regions=[3,4,6];
clear q1 q2 q3 q4 c1 c2 c3 c4
q1=(h_exploit==1);
q2=(h_exploit==0);
q3=(exploit==1);
q4=(exploit==0);
q5=(k_exploit==1);
q6=(k_exploit==0);
q1_next=circshift(q1,[1 -1]); %logical list of trials *before* exploit trials
q2_next=circshift(q2,[1 -1]);
q3_next=circshift(q3,[1 -1]);
q4_next=circshift(q4,[1 -1]);
q5_next=circshift(q5,[1 -1]);
q6_next=circshift(q6,[1 -1]);

llpv=nan(numel(rwd_axis),size(roi,2));
hllpv=nan(numel(rwd_axis),size(roi,2));
kllpv=nan(numel(rwd_axis),size(roi,2));
llpv_pred=nan(numel(rwd_axis),size(roi,2));
hllpv_pred=nan(numel(rwd_axis),size(roi,2));
kllpv_pred=nan(numel(rwd_axis),size(roi,2));

for i=1:length(rwd_axis);
    parsed(i,:)=(rwd==rwd_axis(i));
    parsed_prev(i,:)=(rwd_last_trial==rwd_axis(i));
end

for i=1:length(rwd_axis)
    if ~isempty(find(region==last_regions,1))
        query{i}=find(rwd_last_trial==rwd_axis(i));
    else
        query{i}=find(rwd==rwd_axis(i));
    end
    if ~isempty(query{i})
        d=100*mean(raster(query{i},R),2);
        c(i)=mean(d,1); % averaged over trials
        sem(i)=std(d)/sqrt(length(d)); % std over trials
        %for regression:
        raw_c=[raw_c;d]; %observed firing rates
        raw_x=[raw_x;rwd_axis(i)*ones(length(d),1)]; %matrix of inputs to linear model
    else
        %        c(i)=0;
        %        sem(i)=0;
        c(i)=NaN;
        sem(i)=NaN;
    end

    % now, do t-tests of explore/exploit firing rates in each epoch, controlled for
    % reward
    for j=1:size(roi,2);
        if ~isempty(find(j==last_regions,1))
            sel=parsed_prev;
        else
            sel=parsed;
        end
        if ~isempty(find(sel(i,:)&q3,1))&&~isempty(find(sel(i,:)&q4,1))
            sample=100*mean(raster(:,roi{j}),2);
            [h,llpv(i,j)]=ttest2(sample(sel(i,:)&q3),sample(sel(i,:)&q4));
        end
        if ~isempty(find(sel(i,:)&q1,1))&&~isempty(find(sel(i,:)&q2,1))
            sample=100*mean(raster(:,roi{j}),2);
            [h,hllpv(i,j)]=ttest2(sample(sel(i,:)&q1),sample(sel(i,:)&q2));
        end
        if ~isempty(find(sel(i,:)&q5,1))&&~isempty(find(sel(i,:)&q6,1))
            sample=100*mean(raster(:,roi{j}),2);
            [h,kllpv(i,j)]=ttest2(sample(sel(i,:)&q5),sample(sel(i,:)&q6));
        end
        %corrected significance of strategy selected *next* trial 
        % (only valid for regions 3,4, and 6 (at best))
        if ~isempty(find(sel(i,:)&q3_next,1))&&~isempty(find(sel(i,:)&q4_next,1))
            sample=100*mean(raster(:,roi{j}),2);
            [h,llpv_pred(i,j)]=ttest2(sample(sel(i,:)&q3_next),sample(sel(i,:)&q4_next));
        end
        if ~isempty(find(sel(i,:)&q1_next,1))&&~isempty(find(sel(i,:)&q2_next,1))
            sample=100*mean(raster(:,roi{j}),2);
            [h,hllpv_pred(i,j)]=ttest2(sample(sel(i,:)&q1_next),sample(sel(i,:)&q2_next));
        end
        if ~isempty(find(sel(i,:)&q5_next,1))&&~isempty(find(sel(i,:)&q6_next,1))
            sample=100*mean(raster(:,roi{j}),2);
            [h,kllpv_pred(i,j)]=ttest2(sample(sel(i,:)&q5_next),sample(sel(i,:)&q6_next));
        end
    end

end

% another set of t-tests, now for chunks of reward axis
rwd_subset=1:5:45;
parsed2=zeros(numel(rwd_subset)-1,numel(rwd));
parsed2_prev=zeros(numel(rwd_subset)-1,numel(rwd));
llpv2=nan(numel(rwd_subset)-1,size(roi,2));
hllpv2=nan(numel(rwd_subset)-1,size(roi,2));
kllpv2=nan(numel(rwd_subset)-1,size(roi,2));
for i=1:numel(rwd)
    k=find(rwd(i)>rwd_axis(rwd_subset),1,'last');
    parsed2(k,i)=1;
    k=find(rwd_last_trial(i)>rwd_axis(rwd_subset),1,'last');
    parsed2_prev(k,i)=1;
end
%now do t-tests
for i=1:(numel(rwd_subset)-1)
    for j=1:size(roi,2);
        if ~isempty(find(j==last_regions,1))
            sel=parsed2_prev;
        else
            sel=parsed2;
        end
        if ~isempty(find(sel(i,:)&q3,1))&&~isempty(find(sel(i,:)&q4,1))

            sample=100*mean(raster(:,roi{j}),2);
            [h,llpv2(i,j)]=ttest2(sample(sel(i,:)&q3),sample(sel(i,:)&q4));
        end
        if ~isempty(find(sel(i,:)&q1,1))&&~isempty(find(sel(i,:)&q2,1))
            sample=100*mean(raster(:,roi{j}),2);
            [h,hllpv2(i,j)]=ttest2(sample(sel(i,:)&q1),sample(sel(i,:)&q2));
        end
        if ~isempty(find(sel(i,:)&q5,1))&&~isempty(find(sel(i,:)&q6,1))
            sample=100*mean(raster(:,roi{j}),2);
            [h,kllpv2(i,j)]=ttest2(sample(sel(i,:)&q5),sample(sel(i,:)&q6));
        end
    end
end

stats=regstats(raw_c,raw_x,'quadratic',{'beta','tstat','fstat','rsquare'});
b=stats.beta;

p4=stats.tstat.pval(2);
p4q=stats.tstat.pval(3);
%[h,p4_dumb]=ttest(raw_c,mean(raw_c)); %probability tuning curve is just a flat line

tuning=c;
pisflat1=~(p4<0.05)||~(p4q<0.05); %either linear or quadratic term is significantly nonzero
pisflat2=~(stats.fstat.pval<0.05);

output.plot4.tuning=c;
output.plot4.stats=stats;
output.plot4.pisflat1=pisflat1;
output.plot4.pisflat2=pisflat2;

output.plot4.parsed=parsed;
output.plot4.parsed_prev=parsed_prev;

output.plot4.llpv=llpv;
output.plot4.hllpv=hllpv;
output.plot4.kllpv=kllpv;
output.plot4.llpv_pred=llpv_pred;
output.plot4.hllpv_pred=hllpv_pred;
output.plot4.kllpv_pred=kllpv_pred;
output.plot4.llpv_grouped=llpv2;
output.plot4.hllpv_grouped=hllpv2;
output.plot4.kllpv_grouped=kllpv2;

if ~isempty(find(4==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    
    errorbar(rwd_axis,c,sem,'d')
    if ~isempty(find(region==last_regions))
        xlabel('Reward Last Trial')
    else
        xlabel('Reward This Trial')
    end
    hold on
    %plot(raw_x(:,2)',raw_c','rd')
    plot(rwd_axis,x2fx(rwd_axis','quadratic')*b,'r')
    title(sprintf('slope=%f p=%f',b(2),p4))
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 5 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing vs rwd change on last trial
clear c;
clear sem;
clear raw_c;
clear raw_x;
raw_c=[];
raw_x=[];
rwd_change_axis=-30:5:30;
rwd_change(1)=0;
rwd_change(2)=0;
if ~isempty(find(region==last_regions))
    for j=3:length(rwd)
        rwd_change(j)=rwd_last_trial(j)-rwd_last_trial(j-1);
    end
else
    for j=2:length(rwd)
        rwd_change(j)=rwd(j)-rwd_last_trial(j);
    end
end
for i=1:length(rwd_change_axis)
   query{i}=find((rwd_change==rwd_change_axis(i)));
   if ~isempty(query{i})
       d=100*mean(raster(query{i},R),2);
       c(i)=mean(d,1); % averaged over trials
       sem(i)=std(d)/sqrt(length(d)); % std over trials
       %for regression:
       raw_c=[raw_c;d]; %observed firing rates
       raw_x=[raw_x;rwd_change_axis(i)*ones(length(d),1)]; %matrix of inputs to linear model 
   else
       c(i)=NaN;
       sem(i)=NaN;
   end
end
stats=regstats(raw_c,raw_x,'linear',{'beta','tstat'});
b=stats.beta;
p5=stats.tstat.pval(2);

ouput(5).stats=stats;

if ~isempty(find(5==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    
    errorbar(rwd_change_axis,c,sem,'d')
    if ~isempty(find(region==last_regions))
        xlabel('Change in Reward (last trial)')
    else
        xlabel('Change in Reward (this trial)')
    end
    hold on
    %plot(raw_x(:,2)',raw_c','rd')
    plot(rwd_change_axis,[ones(length(rwd_change_axis),1),rwd_change_axis']*b,'r')
    title(sprintf('slope=%f p=%f',b(2),p5))
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 6 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing vs rwd change on last trial (when exploiting) 
clear c;
clear sem;
clear raw_c;
clear raw_x;
raw_c=[];
raw_x=[];
rwd_change_axis=-30:5:30;
rwd_change(1)=0;
rwd_change(2)=0;
if region==1
    for j=3:length(rwd)
        rwd_change(j)=rwd_last_trial(j)-rwd_last_trial(j-1);
    end
else
    for j=2:length(rwd)
        rwd_change(j)=rwd(j)-rwd_last_trial(j);
    end
end
for i=1:length(rwd_change_axis)
   query{i}=find((rwd_change==rwd_change_axis(i))&exploit);
   if ~isempty(query{i})
       d=100*mean(raster(query{i},R),2);
       c(i)=mean(d,1); % averaged over trials
       sem(i)=std(d)/sqrt(length(d)); % std over trials
       %for regression:
       raw_c=[raw_c;d]; %observed firing rates
       raw_x=[raw_x;rwd_change_axis(i)*ones(length(d),1)]; %matrix of inputs to linear model 
   else
       c(i)=NaN;
       sem(i)=NaN;
   end
end
stats=regstats(raw_c,raw_x,'linear',{'beta','tstat'});
b=stats.beta;
p6=stats.tstat.pval(2);

output.plot6.stats=stats;

if ~isempty(find(6==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    
    errorbar(rwd_change_axis,c,sem,'d')
    xlabel('Change in Reward (while exploiting)')
    hold on
    %plot(raw_x(:,2)',raw_c','rd')
    plot(rwd_change_axis,[ones(length(rwd_change_axis),1),rwd_change_axis']*b,'r')
    title(sprintf('slope=%f p=%f',b(2),p6))
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 7 %%%%%%%%%%%%%%%%%%%%%%%%%%
% raster plot
if ~isempty(find(7==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    
    imagesc(raster)
    colormap gray
    hold on
    vert=size(raster,1);
    plot([200 200],[0 vert],'linestyle','--','color','r')
    plot([300 300],[0 vert],'linestyle','--','color','r')
    % line(200*ones(1,vert+1),0:vert,'linestyle','--','color','r')
    % line(300*ones(1,vert+1),0:vert,'linestyle','--','color','r')
    ylabel('trials')
    xlabel('bins')
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 8 %%%%%%%%%%%%%%%%%%%%%%%%%%
% firing: first explore trial vs all others
last_trial_explore(1)=0;
for i=2:length(explore)
    last_trial_explore(i)=explore(i-1);
end
q1=find(exploit==1);
q2=find((explore==1)&(last_trial_explore==0));
c1=100*mean(raster(q1,:));
c2=100*mean(raster(q2,:));

%t-tests
for i=1:size(roi,2)
    sample=100*mean(raster(:,roi{i}),2);
    [h,pval7(i)]=ttest2(sample(q1),sample(q2));
%     [h,pval7(i)]=ttest2(c1(roi{i}),c2(roi{i}));
end

output.plot8.tpvals=pval7;

if ~isempty(find(8==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    
    hold on
    plot(t_axis,gsmooth(c1,1),'r')
    plot(t_axis,gsmooth(c2,1),'k')
    plot([t_axis(1) t_axis(end)],[mean(c1(R)) mean(c1(R))],'r')
    plot([t_axis(1) t_axis(end)],[mean(c2(R)) mean(c2(R))],'k')
    vert=ylim;
    vert=vert(2);
    plot([0 0],[0 vert],'linestyle','--','color','k')
%     plot([t_axis(R(end)) t_axis(R(end))],[0 vert],'linestyle','--','color','k')
    plot([mean_go_time mean_go_time],[0 vert],'linestyle','--','color','g');
    plot([mean_juice_start mean_juice_start],[0 vert],'r')
    plot([mean_juice_start-std_juice_start mean_juice_start-std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_juice_start+std_juice_start mean_juice_start+std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_iti_time mean_iti_time],[0 vert],'linestyle','--','color',gray)
    plot([mean_over_time mean_over_time],[0 vert],'linestyle','--','color','k')
    legend({'all others','first explore trial'},'Location','SouthWest') 
    title(pval7(region))
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 9 %%%%%%%%%%%%%%%%%%%%%%%%%%
% behavior vs firing rate
if ~isempty(find(9==plots))
    firingvsbehavior_plot(raster,location,b_advice,rwd,1:numel(rwd),exploit,R,tuning,targ_tune,plot_colors,rt);
end

%%%%%%%%%% plot 11 %%%%%%%%%%%%%%%%%%%%%%%%%
% derive filter going back n trials
if ~isempty(find(11==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    
    filter=findfilter(raster,location,rwd,1:numtrials,R,15);
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 12 %%%%%%%%%%%%%%%%%%%%%%%%%
% probability of switching as a function of firing rate

fr=100*mean(raster(:,R),2);
numbins=20;
[counts centers]=hist(fr,numbins);
spacing=0.5*(centers(2)-centers(1));
edges=[-inf centers+spacing inf];
edges(numel(edges)-1)=[]; % remove extra edge, or we have one too many bins

for i=1:(numel(edges)-1)
    trials=find((fr>=edges(i))&(fr<edges(i+1)));
    if ~isempty(find(trials==1))
        trials(find(trials==1))=[]; %exclude first trial
    end
    if ~isempty(find(trials==numel(rwd)))
        trials(find(trials==numel(rwd)))=[]; %exclude last trial
    end
    av_loit(i)=nanmean(exploit(trials+1));
    sem_loit(i)=nanstd(exploit(trials+1))/max(sqrt(numel(trials)),1);
    av_switch(i)=nanmean(switchup(trials+1));
    sem_switch(i)=nanstd(switchup(trials+1))/max(sqrt(numel(trials)),1);
    av_strat(i)=nanmean(stratswitch(trials+1));
    sem_strat(i)=nanstd(stratswitch(trials+1))/max(sqrt(numel(trials)),1);
    av_hloit(i)=nanmean(h_exploit(trials+1));
    sem_hloit(i)=nanstd(h_exploit(trials+1))/max(sqrt(numel(trials)),1);
    av_kloit(i)=nanmean(k_exploit(trials+1));
    sem_kloit(i)=nanstd(k_exploit(trials+1))/max(sqrt(numel(trials)),1);
end

R_corr=corrcoef([centers' av_loit' av_switch'],'rows','complete');


sh_exploit=exploit;
sh_switchup=switchup;
sh_stratswitch=stratswitch;
sh_hexploit=h_exploit;
sh_kexploit=k_exploit;
sh_fr=fr;

if ~isempty(find(region==last_regions)) % regressing against *next* trial, so shift binaries
    sh_exploit(1)=[];
    sh_switchup(1)=[];
    sh_stratswitch(1)=[];
    sh_hexploit(1)=[];
    sh_kexploit(1)=[];
    sh_fr(end)=[];
end
    
stats1=regstats(sh_exploit,sh_fr,'linear',{'beta','tstat','fstat','rsquare'});
b1=stats1.beta;
loit_pv=stats1.tstat.pval(2);
stats2=regstats(sh_switchup,sh_fr,'linear',{'beta','tstat','fstat','rsquare'});
b2=stats2.beta;
switch_pv=stats2.tstat.pval(2);
stats3=regstats(sh_stratswitch,sh_fr,'linear',{'beta','tstat','fstat','rsquare'});
b3=stats3.beta;
strat_pv=stats3.tstat.pval(2);
stats4=regstats(sh_hexploit,sh_fr,'linear',{'beta','tstat','fstat','rsquare'});
b4=stats4.beta;
hloit_pv=stats4.tstat.pval(2);
stats5=regstats(sh_kexploit,sh_fr,'linear',{'beta','tstat','fstat','rsquare'});
b5=stats4.beta;
kloit_pv=stats5.tstat.pval(2);

pv_pred(1,:)=[loit_pv,switch_pv,strat_pv,hloit_pv,kloit_pv];

output.plot12.corrs=R_corr;
output.plot12.pvals=pv_pred(1,:);
output.plot12.centers=centers;
output.plot12.probs=[av_loit;av_switch;av_strat;av_hloit;av_kloit];
output.plot12.sem=[sem_loit;sem_switch;sem_strat;sem_hloit;sem_kloit];
output.plot12.slopes=[b1(2),b2(2),b3(2),b4(2),b5(2)];

if ~isempty(find(12==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)

    hold on
    errorbar(centers,av_loit,sem_loit,'bs')
    errorbar(centers,av_switch,sem_switch,'rd')
    errorbar(centers,av_strat,sem_strat,'kh')
    errorbar(centers,av_hloit,sem_hloit,'gx')
    errorbar(centers,av_kloit,sem_kloit,'y^')
    plot(fr,x2fx(fr,'linear')*b1,'b')
    plot(fr,x2fx(fr,'linear')*b2,'r')
    plot(fr,x2fx(fr,'linear')*b3,'k')
    plot(fr,x2fx(fr,'linear')*b4,'g')
    plot(fr,x2fx(fr,'linear')*b5,'y')
    xlabel('firing rate')
    ylabel('probability next trial');
    legend('exploitation','target switch','strategy switch','heuristic','kalman','location','southoutside')
    title(sprintf('p-values: exploit=%f target switch=%f \n strat switch=%f heuristic=%f heuristic=%f',loit_pv,switch_pv,strat_pv,hloit_pv,kloit_pv))

    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 13 %%%%%%%%%%%%%%%%%%%%%%%%%
% probability of switching as a function of firing rate (corrected for
% dependency of firing on rwd and targ_tune

fr=100*mean(raster(:,R),2);

if output.plot4.pisflat2<0.05 %if rwd tuning cannot be considered flat
    isflat_rwd=0;
else
    isflat_rwd=1;
end
if output.plot3.pval<0.05 %if target tuning has main effect of location by F-test
    isflat_targ=0;
else
    isflat_targ=1;
end

if isflat_targ && isflat_rwd %tuning for neither
    reg_code=0;
elseif isflat_targ && ~isflat_rwd %tuning for rwd only
    reg_code=1;
elseif ~isflat_targ && isflat_rwd %tuning for target only
    reg_code=2;
elseif ~isflat_targ && ~isflat_rwd %tuning for both
    reg_code=3;
end

switch reg_code
    case 1,
        for i=1:numel(rwd)
            rwd_ind=find(rwd(i)==rwd_axis);
            inputs(i,:)=tuning(rwd_ind);
        end
    case 2,
        for i=1:numel(rwd)
            rwd_ind=find(rwd(i)==rwd_axis);
            inputs(i,:)=targ_tune(location(i));
        end
    case 3,
        for i=1:numel(rwd)
            rwd_ind=find(rwd(i)==rwd_axis);
            inputs(i,:)=[tuning(rwd_ind) targ_tune(location(i))];
        end
end

if reg_code ~= 0
    stats0=regstats(fr,inputs,'linear',{'beta','tstat','fstat','rsquare'});
    b0=stats0.beta;

    fr_rev=fr-x2fx(inputs,'linear')*b0;
else
    fr_rev=fr;
    stats0=0;
end

reg_out=(fr-fr_rev)/100; % remember, we have to multiply rasters by 100

for i=1:size(raster,1)
    raster_rev(i,:)=raster(i,:)-reg_out(i);
end

%%% calculate effect of target tuning alone
if ((reg_code==2)||(reg_code==3))
    for i=1:numel(rwd)
        rwd_ind=find(rwd(i)==rwd_axis);
        inputs2(i,:)=targ_tune(location(i));
    end
    stats00=regstats(fr,inputs2,'linear',{'beta','tstat','fstat','rsquare'});
    b00=stats00.beta;
    fr_rev2=fr-x2fx(inputs2,'linear')*b00;
else
    fr_rev2=fr;
    stats00=0;
end

reg_out2=(fr-fr_rev2)/100; % remember, we have to multiply rasters by 100

for i=1:size(raster,1)
    raster_rev2(i,:)=raster(i,:)-reg_out2(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numbins=20;
[counts centers]=hist(fr_rev,numbins);
spacing=0.5*(centers(2)-centers(1));
edges=[-inf centers+spacing inf];
edges(numel(edges)-1)=[]; % remove extra edge, or we have one too many bins

for i=1:(numel(edges)-1)
    trials=find((fr_rev>=edges(i))&(fr_rev<edges(i+1)));
    if ~isempty(find(trials==1))
        trials(find(trials==1))=[]; %exclude first trial
    end
    if ~isempty(find(trials==numel(rwd)))
        trials(find(trials==numel(rwd)))=[]; %exclude last trial
    end
    av_loit(i)=nanmean(exploit(trials+1));
    sem_loit(i)=nanstd(exploit(trials+1))/max(sqrt(numel(trials)),1);
    av_switch(i)=nanmean(switchup(trials+1));
    sem_switch(i)=nanstd(switchup(trials+1))/max(sqrt(numel(trials)),1);
    av_strat(i)=nanmean(stratswitch(trials+1));
    sem_strat(i)=nanstd(stratswitch(trials+1))/max(sqrt(numel(trials)),1);
    av_hloit(i)=nanmean(h_exploit(trials+1));
    sem_hloit(i)=nanstd(h_exploit(trials+1))/max(sqrt(numel(trials)),1);
    av_kloit(i)=nanmean(k_exploit(trials+1));
    sem_kloit(i)=nanstd(k_exploit(trials+1))/max(sqrt(numel(trials)),1);
end

R_corr=corrcoef([centers' av_loit' av_switch'],'rows','complete');

sh_fr_rev=fr_rev;

if ~isempty(find(region==last_regions)) % regressing against *next* trial, so shift binaries
    sh_fr_rev(end)=[];
end

stats1=regstats(sh_exploit,sh_fr_rev,'linear',{'beta','tstat','fstat','rsquare'});
b1=stats1.beta;
loit_pv=stats1.tstat.pval(2);
stats2=regstats(sh_switchup,sh_fr_rev,'linear',{'beta','tstat','fstat','rsquare'});
b2=stats2.beta;
switch_pv=stats2.tstat.pval(2);
stats3=regstats(sh_stratswitch,sh_fr_rev,'linear',{'beta','tstat','fstat','rsquare'});
b3=stats3.beta;
strat_pv=stats3.tstat.pval(2);
stats4=regstats(sh_hexploit,sh_fr_rev,'linear',{'beta','tstat','fstat','rsquare'});
b4=stats4.beta;
hloit_pv=stats4.tstat.pval(2);
stats5=regstats(sh_kexploit,sh_fr_rev,'linear',{'beta','tstat','fstat','rsquare'});
b5=stats5.beta;
hloit_pv=stats5.tstat.pval(2);

pv_pred(2,:)=[loit_pv,switch_pv,strat_pv,hloit_pv,kloit_pv];

output.plot13.reg_code=reg_code;
output.plot13.stats0=stats0;
output.plot13.stats00=stats00;
output.plot13.corrs=R_corr;
output.plot13.pvals=pv_pred(2,:);
output.plot13.centers=centers;
output.plot13.probs=[av_loit;av_switch;av_strat;av_hloit;av_kloit];
output.plot13.sem=[sem_loit;sem_switch;sem_strat;sem_hloit;sem_kloit];
output.plot13.slopes=[b1(2),b2(2),b3(2),b4(2),b5(2)];


if ~isempty(find(13==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)

    hold on
    errorbar(centers,av_loit,sem_loit,'bs')
    errorbar(centers,av_switch,sem_switch,'rd')
    errorbar(centers,av_strat,sem_strat,'kh')
    errorbar(centers,av_hloit,sem_hloit,'gx')
    errorbar(centers,av_hloit,sem_hloit,'y^')
    plot(fr_rev,x2fx(fr_rev,'linear')*b1,'b')
    plot(fr_rev,x2fx(fr_rev,'linear')*b2,'r')
    plot(fr_rev,x2fx(fr_rev,'linear')*b3,'k')
    plot(fr_rev,x2fx(fr_rev,'linear')*b4,'g')
    plot(fr_rev,x2fx(fr_rev,'linear')*b5,'y')
    xlabel('firing rate (corrected)')
    ylabel('probability next trial');
    legend('exploitation','target switch','strategy switch','heuristic','kalman','location','southoutside')
    title(sprintf('p-values: exploit=%f target switch=%f \n strat switch=%f heuristic=%f kalman=%f',loit_pv,switch_pv,strat_pv,hloit_pv,kloit_pv))

    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 15 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing: explore and exploit (bayesian and heuristic, adjusted for reward,targ_tune)
clear q1 q2 q3 q4 c1 c2 c3 c4
q1=find(h_exploit==1);
q2=find(h_exploit==0);
q3=find(exploit==1);
q4=find(exploit==0);
q5=find(k_exploit==1);
q6=find(k_exploit==0);
%use raster_rev from plot13
c1=100*mean(raster_rev(q1,:));
c2=100*mean(raster_rev(q2,:));
c3=100*mean(raster_rev(q3,:));
c4=100*mean(raster_rev(q4,:));
c5=100*mean(raster_rev(q5,:));
c6=100*mean(raster_rev(q6,:));

%t-tests
%only good for R, since raster_rev uses the regression in 13, done over R
sample=mean(raster_rev(:,R),2);
[h,pval15a]=ttest2(sample(q1),sample(q2));
[h,pval15b]=ttest2(sample(q3),sample(q4));
[h,pval15e]=ttest2(sample(q5),sample(q6));

sample=mean(raster_rev2(:,R),2); % corrected for target tuning alone
[h,pval15c]=ttest2(sample(q1),sample(q2));
[h,pval15d]=ttest2(sample(q3),sample(q4));
[h,pval15f]=ttest2(sample(q5),sample(q6));

output.plot15.tpval_h=pval15a;
output.plot15.tpval=pval15b;
output.plot15.tpval_k=pval15e;
output.plot15.ftpval_h=pval15c;
output.plot15.ftpval=pval15d;
output.plot15.ftpval_k=pval15f;
frate=mean(raster(:,R),2);
for i=1:numel(rwd)
    if rwd(i) <= rwd_axis(15)
        rwd_level(i)=1;
    elseif rwd(i) <= rwd_axis(30)
        rwd_level(i)=2;
    else
        rwd_level(i)=3;
    end
end

p_ind=anovan(frate,{h_exploit rwd_level location},'display','off'); %independent factors
p_int=anovan(frate,{h_exploit rwd_level location},'model','interaction','display','off'); %interaction model

output.plot15.anova.p_ind=p_ind;
output.plot15.anova.p_int=p_int;

if ~isempty(find(15==plots))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    hold on
    plot(t_axis,gsmooth(c1,1),'r')
    plot(t_axis,gsmooth(c2,1),'k')
    plot([t_axis(1) t_axis(end)],[mean(c1(R)) mean(c1(R))],'r')
    plot([t_axis(1) t_axis(end)],[mean(c2(R)) mean(c2(R))],'k')
    vert=ylim;
    vert=vert(2);
    plot([0 0],[0 vert],'linestyle','--','color','k')
    text(0,-1.5,'Begin','HorizontalAlignment','center')
    plot([mean_go_time mean_go_time],[0 vert],'linestyle','--','color','g');
    text(mean_go_time,-1.5,'go','color','g','HorizontalAlignment','center')
    plot([mean_juice_start mean_juice_start],[0 vert],'r')
    text(mean_juice_start,-1.5,'juice','color','r','HorizontalAlignment','center')
    plot([mean_juice_start-std_juice_start mean_juice_start-std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_juice_start+std_juice_start mean_juice_start+std_juice_start],[0 vert],'linestyle','--','color','r')
    plot([mean_iti_time mean_iti_time],[0 vert],'linestyle','--','color',gray)
    text(mean_iti_time,-1.5,'Target off','color','k','HorizontalAlignment','center')
    plot([mean_over_time mean_over_time],[0 vert],'linestyle','--','color','k')
    text((mean_iti_time+mean_over_time)/2,-1.5,'ITI','color','k','HorizontalAlignment','center')
    text(mean_over_time,-1.5,'End','color','k','HorizontalAlignment','center')
    legend({'h_exploit','h_explore'},'Location','SouthWest')
    title(sprintf('p=%f %% explore=%f',pval15a,mean(h_explore)))
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

%%%%%%%%%% plot 19 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing vs reward on last trial (or this trial; region-dependent)
%CORRECTED FOR TARGET TUNING
clear c;
clear sem;
clear raw_c;
clear raw_x;
raw_c=[];
raw_x=[];
last_regions=[1,2,5];
this_regions=[3,4,6];
clear q1 q2 q3 q4 q5 q6 c1 c2 c3 c4 c5 c6
q1=(h_exploit==1);
q2=(h_exploit==0);
q3=(exploit==1);
q4=(exploit==0);
q5=(k_exploit==1);
q6=(k_exploit==0);
q1_next=circshift(q1,[1 -1]); %logical list of trials *before* exploit trials
q2_next=circshift(q2,[1 -1]);
q3_next=circshift(q3,[1 -1]);
q4_next=circshift(q4,[1 -1]);
q5_next=circshift(q5,[1 -1]);
q6_next=circshift(q6,[1 -1]);

llpv=nan(numel(rwd_axis),size(roi,2));
hllpv=nan(numel(rwd_axis),size(roi,2));
kllpv=nan(numel(rwd_axis),size(roi,2));
llpv_pred=nan(numel(rwd_axis),size(roi,2));
hllpv_pred=nan(numel(rwd_axis),size(roi,2));
kllpv_pred=nan(numel(rwd_axis),size(roi,2));

for i=1:length(rwd_axis);
    parsed(i,:)=(rwd==rwd_axis(i));
    parsed_prev(i,:)=(rwd_last_trial==rwd_axis(i));
end

for i=1:length(rwd_axis)
    % do t-tests of explore/exploit firing rates in each epoch, controlled for
    % reward
    for j=1:size(roi,2);
        if ~isempty(find(j==last_regions,1))
            sel=parsed_prev;
        else
            sel=parsed;
        end
        if ~isempty(find(sel(i,:)&q3,1))&&~isempty(find(sel(i,:)&q4,1))
            sample=100*mean(raster_rev2(:,roi{j}),2);
            [h,llpv(i,j)]=ttest2(sample(sel(i,:)&q3),sample(sel(i,:)&q4));
        end
        if ~isempty(find(sel(i,:)&q1,1))&&~isempty(find(sel(i,:)&q2,1))
            sample=100*mean(raster_rev2(:,roi{j}),2);
            [h,hllpv(i,j)]=ttest2(sample(sel(i,:)&q1),sample(sel(i,:)&q2));
        end
        if ~isempty(find(sel(i,:)&q5,1))&&~isempty(find(sel(i,:)&q6,1))
            sample=100*mean(raster_rev2(:,roi{j}),2);
            [h,kllpv(i,j)]=ttest2(sample(sel(i,:)&q5),sample(sel(i,:)&q6));
        end
        %corrected significance of strategy selected *next* trial 
        % (only valid for regions 3,4, and 6 (at best))
        if ~isempty(find(sel(i,:)&q3_next,1))&&~isempty(find(sel(i,:)&q4_next,1))
            sample=100*mean(raster_rev2(:,roi{j}),2);
            [h,llpv_pred(i,j)]=ttest2(sample(sel(i,:)&q3_next),sample(sel(i,:)&q4_next));
        end
        if ~isempty(find(sel(i,:)&q1_next,1))&&~isempty(find(sel(i,:)&q2_next,1))
            sample=100*mean(raster_rev2(:,roi{j}),2);
            [h,hllpv_pred(i,j)]=ttest2(sample(sel(i,:)&q1_next),sample(sel(i,:)&q2_next));
        end
        if ~isempty(find(sel(i,:)&q5_next,1))&&~isempty(find(sel(i,:)&q6_next,1))
            sample=100*mean(raster_rev2(:,roi{j}),2);
            [h,kllpv_pred(i,j)]=ttest2(sample(sel(i,:)&q5_next),sample(sel(i,:)&q6_next));
        end
    end

end

% another set of t-tests, now for chunks of reward axis
rwd_subset=1:5:45;
parsed2=zeros(numel(rwd_subset)-1,numel(rwd));
parsed2_prev=zeros(numel(rwd_subset)-1,numel(rwd));
llpv2=nan(numel(rwd_subset)-1,size(roi,2));
hllpv2=nan(numel(rwd_subset)-1,size(roi,2));
kllpv2=nan(numel(rwd_subset)-1,size(roi,2));
for i=1:numel(rwd)
    k=find(rwd(i)>rwd_axis(rwd_subset),1,'last');
    parsed2(k,i)=1;
    k=find(rwd_last_trial(i)>rwd_axis(rwd_subset),1,'last');
    parsed2_prev(k,i)=1;
end
%now do t-tests
for i=1:(numel(rwd_subset)-1)
    for j=1:size(roi,2);
        if ~isempty(find(j==last_regions,1))
            sel=parsed2_prev;
        else
            sel=parsed2;
        end
        if ~isempty(find(sel(i,:)&q3,1))&&~isempty(find(sel(i,:)&q4,1))

            sample=100*mean(raster_rev2(:,roi{j}),2);
            [h,llpv2(i,j)]=ttest2(sample(sel(i,:)&q3),sample(sel(i,:)&q4));
        end
        if ~isempty(find(sel(i,:)&q1,1))&&~isempty(find(sel(i,:)&q2,1))
            sample=100*mean(raster_rev2(:,roi{j}),2);
            [h,hllpv2(i,j)]=ttest2(sample(sel(i,:)&q1),sample(sel(i,:)&q2));
        end
        if ~isempty(find(sel(i,:)&q5,1))&&~isempty(find(sel(i,:)&q6,1))
            sample=100*mean(raster_rev2(:,roi{j}),2);
            [h,kllpv2(i,j)]=ttest2(sample(sel(i,:)&q5),sample(sel(i,:)&q6));
        end
    end
end

output.plot19.llpv=llpv;
output.plot19.hllpv=hllpv;
output.plot19.kllpv=kllpv;
output.plot19.llpv_pred=llpv_pred;
output.plot19.hllpv_pred=hllpv_pred;
output.plot19.kllpv_pred=kllpv_pred;
output.plot19.llpv_grouped=llpv2;
output.plot19.hllpv_grouped=hllpv2;
output.plot19.kllpv_grouped=kllpv2;

%%%%%%%%%% plot 18 %%%%%%%%%%%%%%%%%%%%%%%%%%
%firing: explore and exploit (bayesian and heuristic): correlation with
%reaction time
clear q1 q2 q3 q4 c1 c2 c3 c4
q1=find(h_exploit==1);
q2=find(h_exploit==0);
q3=find(exploit==1);
q4=find(exploit==0);
q5=find(k_exploit==1);
q6=find(k_exploit==0);

%t-tests, kstests (distributions probably not normal) for significant rt differences in 
%explore vs exploit
[h,rt_pv]=ttest2(rt(q3),rt(q4));
[h,rt_pvh]=ttest2(rt(q1),rt(q2));
[h,rt_pvk]=ttest2(rt(q5),rt(q6));
[h,rt_pks]=kstest2(rt(q3),rt(q4));
[h,rt_pksh]=kstest2(rt(q1),rt(q2));
[h,rt_pksk]=kstest2(rt(q5),rt(q6));

%regression: firing vs rt (not clear which is independent or dependent, but
%all that matters is significant correlation)
sample=100*mean(raster(:,R),2);
stats=regstats(sample,rt,'linear',{'beta','tstat','fstat','rsquare'});
[rval,corrp]=corrcoef(sample,rt);

output.plot18.tpval=rt_pv;
output.plot18.tpval_h=rt_pvh;
output.plot18.tpval_k=rt_pvk;
output.plot18.kspval=rt_pks;
output.plot18.kspval_h=rt_pksh;
output.plot18.kspval_k=rt_pksk;
output.plot18.loit=mean(rt(q3));
output.plot18.lore=mean(rt(q4));
output.plot18.hloit=mean(rt(q1));
output.plot18.hlore=mean(rt(q2));
output.plot18.kloit=mean(rt(q5));
output.plot18.klore=mean(rt(q6));
output.plot18.fvsrt=stats;
output.plot18.corr=rval;
output.plot18.pcorr=corrp;

if ~isempty(find(18==plots, 1))
    figure(fig_handle(currpane))
    subplot(2,2,currplot)
    hold on
    % plot stuff goes here %%%%%%%%%
    
    currplot=currplot+1;
    if currplot>4
        currpane=currpane+1;
        currplot=mod(currplot,4);
    end
end

if give_keyboard
    keyboard
end
%}
catch
    q=lasterror;
    keyboard
end