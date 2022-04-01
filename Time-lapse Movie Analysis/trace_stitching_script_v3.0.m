%this script does mitosis calling and stitching of traces for silencing and
%calls silencing events




%part I identify mitotic events to break stitching up into fragments
%involves generating gradient traces to detect change in slope for H2B

f=smoothrows(cit_good,2);
f3=smoothrows(mch_good,2);
c1=[14 189 109]./255; %green
c2=[187 23 167]./255; %maroon

%  subplot(2,2,1)
%  plot(f(222,:)),title('raw trace (I)'),grid on,hold on
 
 for j=1:size(f,1) %loop over number of sc traces 
     
    for n=1:size(f,2)-3  %loop over timepoints
     
      if abs(f(j,n+1)/f(j,n))<0.4  %if signal drop after division << half signal,fix the point using neighbors +/- 1hr window
         f(j,n+1)=f(j,n+3);
         f(j,n+2)=f(j,n+3);
      end
     end
 
 end
 
 
 filtered_sc=[];
 
 for j=1:size(f,1)
     
    f2=smooth((f(j,:)));    
     
 
        spikes=find(gradient(f2(1:end-4)));
       
            if max(diff(spikes)<4) || numel(spikes)>4
 
                filtered_sc(j,:)=f2;
 
                else
                f2(find(gradient(f2(1:end-4))>1e4):find(gradient(f2(1:end-4))>1e4)+3)=f(find(gradient(f2)>1e4)+3);
 
                filtered_sc(j,:)=f2;
            end

 end



 
 %call division events in the traces. variable called splits
splits={};

for q=1:size(filtered_sc,1)
    
  % p=findchangepts(filtered_sc(q,:),'Statistic','linear','MaxNumChanges',5);
  

     [peaks,p]=findpeaks(sgolayfilt(-1*gradient(filtered_sc(q,:)),1,3),'MinPeakDistance',30,'MinPeakHeight',0.008); %0.01
     
     
    if (filtered_sc(q,p(find(diff(p)<20))))<filtered_sc(q,p(find(diff(p)<20)+1));
        p(find(diff(p)<20)+1)=[];
    else
    p(find(diff(p)<20))=[];
    end

                splits{q}=p;


end
  
 
 
 %%

%generate stitched traces from the data for citrine
yf={};
stitched_lineage=[];
yf_fit={};
yf_fit_sil=[];
fs={};
ts={};
c5=[200 200 200]./255;


slope_threshold=4e-3;
grad_threshold=0.7;


for f=1:200 %1465:1469%f=779%size(filtered_sc,1)  %f is number of tracked cells to loop over 
    
    %size(filtered_sc,1)%pick cell index


close all
subplot(3,2,1)
plot(cit_good(f,:),'-','color',c1),title('Raw Citrine'),grid on,
if ~isempty(splits{f})
vline(splits{f},'k'),axis tight
end
windowsize=4; %window around division where to look for max/min points for stitching
test=cit_good(f,:);

divisions=splits{f};
trace_stitched=test;
frames_stitched=1:length(test);

if min(divisions)<= windowsize 
    
    %|| max(divisions)>= length(test)-windowsize;
    yf_fit{f}=spline_dat;
    sil_fun=NaN;
    yf_fit_sil(f)=sil_fun;
    f=f+1;
        
  
    
    
else
    
    
%stitch across divisions, from the last to the first
for n=1:length(divisions)
    
    %account for cells that divide super early or late (frame<windowsize)
     divisions(divisions>length(test)-windowsize)=length(test)-windowsize;
      divisions(divisions<windowsize)=windowsize;

    
    i=length(divisions)-n+1; %go backwards with stitching, otherwise indixes get messed up
    [maxf(i),maxind(i)]=max(filtered_sc(f,divisions(i)-windowsize:divisions(i))); %find max point right before division
    [minf(i),minind(i)]=min(filtered_sc(f,divisions(i):divisions(i)+windowsize)); %find min point right after division
    delta(i)=maxf(i)-minf(i); %change in fluorescence


%shift the entire trace after this division up, throw away points between
%max and min, variable dependeing on division, but smaller than windowsize
    trace_stitched=[trace_stitched(1:divisions(i)-windowsize+maxind(i)-1),trace_stitched(divisions(i)+minind(i):end)+delta(i)];
%keep track of frames we discarded
    frames_stitched=[frames_stitched(1:divisions(i)-windowsize+maxind(i)-1),frames_stitched(divisions(i)+minind(i):end)];
end


end

fs{f}=frames_stitched;
ts{f}=trace_stitched;

yf{f}=trace_stitched;
stitched_lineage{f}=lineage_data1(f);
subplot(3,2,3)
plot(frames_stitched./3,trace_stitched, '-','color',c1,'linewidth',1),axis tight,hold on,

% polyfitting method for plots

% r=polyfit(frames_stitched,trace_stitched,2);
% x1 = linspace(0,360);
% y1 = polyval(r,x1);

%spline fitting for plots
ydata=yf{f};
ydata=ydata';
xdata=(1:1:size(ydata,1))';
f1 = fit(xdata,ydata,'smoothingspline','SmoothingParam',0.0001);

spline_dat=medfilt1(feval(f1,xdata),3);

%if spline_dat(end)>2.2||spline_dat(72)<0.7 % set to 2.2 and 1.1 for krab these are traces that were not stitched properly and need to be ommitted 
if spline_dat(end)>2.2||spline_dat(72)<0.7 % set to 2.2 and 1.1 for krab these are traces that were not stitched properly and need to be ommitted 
 
    
    yf_fit{f}=spline_dat;
    sil_fun=NaN;
    yf_fit_sil(f)=sil_fun;
    
 plot(frames_stitched./3,spline_dat,'k','linewidth',1),hold on,
  yyaxis right,plot(frames_stitched./3,sgolayfilt(gradient(spline_dat),3,35),'r') %for splinefit
  ylabel('Spline Gradient')
    
    
else
    
    
%plot(x1,y1,'k','linewidth',1),axis tight,hold on   %for polyfit
plot(frames_stitched./3,spline_dat,'k','linewidth',1),hold on,
%sil_fun=min(find(gradient(y1)<0.2*max(gradient(y1)))); %for polyfit
sil_fun=find(gradient(spline_dat)<grad_threshold*max(gradient(spline_dat)));
sil_fun=min(sil_fun(sil_fun>72));
%plot(x1(sil_fun),y1(sil_fun),'ko','MarkerFaceColor', 'k') for polyfit
 yyaxis right,plot(frames_stitched./3,sgolayfilt(gradient(spline_dat),3,35),'b'),hold on %for splinefit
ylabel('Spline Gradient')

y1=spline_dat;

sil_fun
yf_fit{f}=y1;
if (isempty(sil_fun)||sil_fun<1)
    sil_fun=NaN;
    
    
else
    
    sil_fun=floor(sil_fun/(size(spline_dat,1))*size(filtered_sc,2)); %correction based on missing frames of sititching (pushes silencing time a bit back)

end



if size(find(splits{f}>sil_fun),2)>1 %can check cells only if they have another division after silencing (most cells do)
slope_check=median(gradient(spline_dat(splits{f}(min(find(splits{f}>sil_fun))):end)));

if slope_check > slope_threshold
    sil_fun=NaN;
end
end




yf_fit_sil(f)=sil_fun;
end

if ~isnan(sil_fun)
    
vline(frames_stitched(sil_fun)./3),title('Stitched Gradient')


subplot(3,2,5)

plot(time,cit_good(f,:),'-','color',c5,'linewidth',1),title('filtered trace (II)'),vline(splits{f}./3),axis tight
hold on, yyaxis right,
 plot(frames_stitched./3,spline_dat,'-o','color',c1,'MarkerSize',3),axis tight,hold on,
  plot(frames_stitched(sil_fun)./3,spline_dat(sil_fun),'ko','MarkerFaceColor', 'k','MarkerSize',8),title('Stitched Citrine')
  legend('Raw','Stitched')
  
  
  
sil_data_tracked(f).cit_stitched=horzcat((frames_stitched')./3,spline_dat); %record frames/stitched trace
sil_data_tracked(f).cit_sil=(frames_stitched(sil_fun))./3; %record silencing time
sil_data_tracked(f).cit_sil_frame=(sil_fun); %store silencing frame in movie
sil_data_tracked(f).cit_raw=cit_good(f,:); %store silencing time
sil_data_tracked(f).cit_raw_stitched=trace_stitched; %store raw stitched trace
sil_data_tracked(f).cit_grad=sgolayfilt(gradient(spline_dat),3,35); %store gradient of stitiched citrine
sil_data_tracked(f).cit_splits=(splits{f})./3; %store mitosis points


else %for cells that did not silence
    

    
subplot(3,2,5)

plot(time,cit_good(f,:),'-','color',c5,'linewidth',1),title('filtered trace (II)'),vline(splits{f}./3),axis tight
hold on, yyaxis right,
 plot(frames_stitched./3,spline_dat,'-o','color',c1,'MarkerSize',3),axis tight,hold on,
 % plot(frames_stitched(sil_fun)./3,spline_dat(sil_fun),'ko','MarkerFaceColor', 'k','MarkerSize',8),title('Stitched Citrine')
  legend('Raw','Stitched')
  
  
  
sil_data_tracked(f).cit_stitched=horzcat((frames_stitched')./3,spline_dat); %record frames/stitched trace
sil_data_tracked(f).cit_sil=NaN; %record silencing time
sil_data_tracked(f).cit_sil_frame=(sil_fun); %store silencing frame in movie
sil_data_tracked(f).cit_raw=cit_good(f,:); %store silencing time
sil_data_tracked(f).cit_raw_stitched=trace_stitched; %store raw stitched trace
sil_data_tracked(f).cit_grad=sgolayfilt(gradient(spline_dat),3,35); %store gradient of stitiched citrine
sil_data_tracked(f).cit_splits=(splits{f})./3; %store mitosis points



end
end
yf_lin=yf;
yf=yf(~cellfun('isempty',yf));
stitched_lineage=stitched_lineage(~cellfun('isempty',stitched_lineage));
%%

%generate stitched traces from the data for mcherry

%modify parameters 
%clear sil_data_tracked %delete this when not debugginf
slope_threshold=4e-3;
grad_threshold=0.7;



c4=[232 196 245]./255
mcf={};
stitched_lineage=[];
mcf_fit={};
mcf_fit_sil=[];


for f=1:200%1356%1:1000 %size(filtered_sc,1)%pick cell index

close all
subplot(3,2,1)
plot(mch_good(f,:),'-','color',c2),title('Raw mCherry'),grid on,

if ~isempty(splits{f})

vline(splits{f},'k'),axis tight
end
windowsize=4; %window around division where to look for max/min points for stitching
test=mch_good(f,:);

divisions=splits{f};
trace_stitched=test;
frames_stitched=1:length(test);

if isempty (splits{f}) 
    f=f+1;
elseif min(divisions)<= windowsize
%omit traces with very eary division even
    f=f+1;
else
    
    
%stitch across divisions, from the last to the first
for n=1:length(divisions)
    
    
      %account for cells that divide super early or late (frame<windowsize)
     divisions(divisions>length(test)-windowsize)=length(test)-windowsize;
      divisions(divisions<windowsize)=windowsize;
      
      
    i=length(divisions)-n+1; %go backwards with stitching, otherwise indixes get messed up
    [maxf(i),maxind(i)]=max(f3(f,divisions(i)-windowsize:divisions(i))); %find max point right before division
    [minf(i),minind(i)]=min(f3(f,divisions(i):divisions(i)+windowsize)); %find min point right after division
    delta(i)=maxf(i)-minf(i); %change in fluorescence


%shift the entire trace after this division up, throw away points between
%max and min, variable dependeing on division, but smaller than windowsize
    trace_stitched=[trace_stitched(1:divisions(i)-windowsize+maxind(i)-1),trace_stitched(divisions(i)+minind(i):end)+delta(i)];
%keep track of frames we discarded
    frames_stitched=[frames_stitched(1:divisions(i)-windowsize+maxind(i)-1),frames_stitched(divisions(i)+minind(i):end)];
end


end

fs2{f}=frames_stitched;
ts2{f}=trace_stitched;

frames_stitched=frames_stitched./3;
mcf{f}=trace_stitched;
subplot(3,2,3)
plot(frames_stitched,trace_stitched, '-','color',c2,'linewidth',1),axis tight


mcf{f}=trace_stitched;
stitched_lineage{f}=lineage_data1(f);
subplot(3,2,3)
plot(frames_stitched,trace_stitched, '-','color',c2,'linewidth',1),axis tight
hold on,

% polyfit
%  r=polyfit(frames_stitched,trace_stitched,3);
%  x1 = linspace(0,360);
% y1 = polyval(r,x1);
% 

%spline fitting for plots
ydata=medfilt1(mcf{f},3);  %medfilt???????
ydata=ydata';
xdata=(1:1:size(ydata,1))';
f1 = fit(xdata,ydata,'smoothingspline','SmoothingParam',0.0001);
spline_dat=feval(f1,xdata);
%spline_dat=sgolayfilt(ydata,2,19);
%sil_fun=min(find(gradient(spline_dat)<0.45*max(gradient(spline_dat))));
1


 plot(frames_stitched,spline_dat,'k','linewidth',1),axis tight,hold on,


d=find(gradient(spline_dat(1:end))==max(gradient(spline_dat(1:end))));    
peak=max(gradient(spline_dat(72:end))); %delete
h=find(gradient(spline_dat(72:end))==peak);%delete



index=find((gradient(spline_dat(72:end)))==max(gradient(spline_dat(72:end)))); %delete
%check downstream of spike peak and 




 
 
 if isempty(index)
     sil_fun=min(plats)+1;

 else
%sil_fun=max(plats(index));

sil_fun=min(find(gradient(spline_dat(h+72:end))<0.7*peak))+h+72 ;%delete
if ~isempty(sil_fun)
    
if size(find(splits{f}>sil_fun),2)>1 %can check cells only if they have another division after silencing (most cells do)
slope_check=median(gradient(spline_dat(splits{f}(min(find(splits{f}>sil_fun))):end)));

if slope_check > slope_threshold
    sil_fun=NaN;
end
end
end

 end
 
 %check if in addion to gradient that protein level also flattens by next
 %mitosis event
 
 
 
 
 
 
% 
%plot(x1(sil_fun),y1(sil_fun),'ko','MarkerFaceColor', 'k')


y1=spline_dat;


mcf_fit{f}=y1;
if isempty(sil_fun)
    sil_fun=NaN;

% else   %removed this correction since we are re-includng lost time 
%     sil_fun=floor(sil_fun/(size(spline_dat,1))*size(filtered_sc,2)); %correction based on missing frames of sititching (pushes silencing time a bit back)
%    
%     if sil_fun>size(frames_stitched,2)
%         sil_fun=size(frames_stitched,2);
%     else
%     end
    %sil_fun=floor(sil_fun/(size(spline_dat,1))*size(filtered_sc,2))
end

mcf_fit_sil(f)=sil_fun;

if isnan(sil_fun)
    yyaxis right, plot(frames_stitched,gradient(spline_dat),'b'),axis tight,ylabel('gradient spline') %for splinefit
else
    yyaxis right, plot(frames_stitched,smooth(gradient(spline_dat),29),'b'),hold on,ylabel('gradient spline') %vline(frames_stitched(sil_fun)),axis tight %for splinefit
end


if ~isnan(sil_fun) && max(frames_stitched)>sil_fun;

subplot(3,2,5)
plot(time,mch_good(f,:),'-','color',c5,'linewidth',1),title('Stitched mCherry'),vline(splits{f}./3),axis tight
hold on, yyaxis right,
 plot(frames_stitched,spline_dat,'-o','color',c2,'MarkerSize',3),axis tight,hold on,
  plot(frames_stitched(sil_fun),spline_dat(sil_fun),'ko','MarkerFaceColor', 'k','MarkerSize',8)
  legend('Raw','Stitched')

sil_data_tracked(f).mch_stitched=horzcat(frames_stitched',spline_dat);%store stitched trace
sil_data_tracked(f).mch_sil=(frames_stitched(sil_fun)); %store silencing time
sil_data_tracked(f).mch_sil_frame=(sil_fun); %store silencing frame in movie
sil_data_tracked(f).mch_raw=mch_good(f,:); %store silencing time  
sil_data_tracked(f).mch_raw_stitched=trace_stitched; %store raw stitched trace
sil_data_tracked(f).mch_grad=sgolayfilt(gradient(spline_dat),3,35); %store gradient of stitiched citrine
sil_data_tracked(f).mch_splits=(splits{f})./3; %store mitosis points
  
else 
   subplot(3,2,5)
 
plot(time,mch_good(f,:),'-','color',c5,'linewidth',1),title('Stitched mCherry'),vline(splits{f}./3),axis tight
hold on, yyaxis right,
 plot(frames_stitched,spline_dat,'-o','color',c2,'MarkerSize',3),axis tight,hold on,
  legend('Raw','Stitched')

sil_data_tracked(f).mch_stitched=horzcat(frames_stitched',spline_dat);%store stitched trace
sil_data_tracked(f).mch_sil=NaN; %store silencing time
sil_data_tracked(f).mch_sil_frame=(sil_fun); %store silencing frame in movie
sil_data_tracked(f).mch_raw=mch_good(f,:); %store silencing time  
sil_data_tracked(f).mch_raw_stitched=trace_stitched; %store raw stitched trace
sil_data_tracked(f).mch_grad=sgolayfilt(gradient(spline_dat),3,35); %store gradient of stitiched citrine
sil_data_tracked(f).mch_splits=(splits{f})./3; %store mitosis points 
    

%     
%  
% sil_data_tracked(f).mch_stitch=horzcat(frames_stitched',spline_dat);%store stitched trace
% sil_data_tracked(f).mch_sil=(frames_stitched(sil_fun))./3; %store silencing time
% sil_data_tracked(f).mch_raw=mch_good(f,:); %store silencing time  
     
    
end
f



end

mcf_lin=mcf;
mcf=mcf(~cellfun('isempty',mcf));

%%

%check cumulative silencing with the parameters (ie. no dox should have
%very little silencin, and plus dox should be >80% at endpoint

%%

% generate single cell stitched plots (of all traces) 

%192 175 37
subplot(2,2,1)
yfp_ss=[];
stitched_time=min(cellfun('size',yf,2)); %note this should be shortest trace generated after shaving off window size frames after stitiching
for d=1:size(yf,2)  %yfp traces

 if isempty(yf{d})
     d=d+1;
 else
yfp_ss(d,:)=(yf{d}(1:stitched_time));
     patchline(((1:size(yf{d},2))./3),sgolayfilt(yf{d}(1:end),1,25),'edgecolor', [255 169 0]./255,'linewidth',1,'edgealpha',0.1),grid on;
%plot(yf{d}(1:end),'color', c1),hold on,axis tight,grid on
 end
end
mean_yfp_sil=sgolayfilt(mean(yfp_ss),1,25);
hold on,plot((1:stitched_time)./3,mean_yfp_sil,'-k','linewidth',4)
% yfp_pt_sil=min(findchangepts(mean_yfp_sil,'MaxNumChanges',2,'Statistic','std'));
yfp_pt_sil=find(gradient(mean_yfp_sil)<0.3*(max(gradient(mean_yfp_sil))),1,'first');
hold on,
plot((yfp_pt_sil/3),(mean_yfp_sil(floor(yfp_pt_sil))),'o', 'MarkerSize',10,'MarkerFaceColor', 'k','MarkerEdgeColor','k'),axis tight%%

subplot(2,2,3)
mcf_ss=[];

for d=1:size(mcf,2) %mcherry traces

 if isempty(mcf{d})
     d=d+1;
 else
     
     mcf_ss(d,:)=mcf{d}(1:stitched_time);
     patchline(((1:size(mcf{d},2))./3),sgolayfilt(mcf{d}(1:end),1,25),'edgecolor', [  250, 29, 29]./255,'linewidth',1,'edgealpha',0.1),grid on,axis tight;
 end
end

 mean_mcf_sil=sgolayfilt(median(mcf_ss),1,25);
 smooth_mcf_sil=sgolayfilt(median(mcf_ss),1,25);
 
hold on,plot((1:stitched_time)./3,mean_mcf_sil,'-k','linewidth',4)


 %mcf_pt_sil=max(findchangepts(med_mcf_sil,'MaxNumChanges',2,'Statistic','std'));
 mcf_pt_sil=find(gradient(smooth_mcf_sil(150:end))<0.6*(max(gradient(smooth_mcf_sil(150:end)))),1,'first');
 mcf_pt_sil=150+mcf_pt_sil;
 
 hold on,
plot((mcf_pt_sil/3),(mean_mcf_sil(floor(mcf_pt_sil))),'o', 'MarkerSize',10,'MarkerFaceColor', 'k','MarkerEdgeColor','k')




%%

%generate traces for unique cells (random picked from single mother
%lineage)



unique_cells=unique(cell2mat(stitched_lineage))'; %unique parent trace IDs
rand_origin_cells={};

for i=1:size(unique_cells,1)
pickcell=(unique_cells(i)==cell2mat(stitched_lineage)')>0;
index=find(pickcell>0);
pos = randi(size(index,1));
rand_origin_cells{i}=index(pos);

% for j=1:size(index)
% plot(sgolayfilt(mcf_ss(index(j),:),1,19)),hold on

end

 rand_origin_cells=cell2mat(rand_origin_cells);
 mcf_ss_origin=mcf_ss(rand_origin_cells,:);
 yfp_ss_origin=yfp_ss(rand_origin_cells,:);
 mean_mcf_sil_origin=sgolayfilt(mean(mcf_ss_origin),2,11);
  mean_yfp_sil_origin=sgolayfilt(mean(yfp_ss_origin),2,11);

 
%  mcf_pt_sil_origin=max(findchangepts(mean_mcf_sil_origin,'MaxNumChanges',2,'Statistic','std'));
 mcf_pt_sil_origin=find(gradient(mean_mcf_sil(150:end))<0.3*(max(gradient(mean_mcf_sil(150:end)))),1,'first');
 mcf_pt_sil_origin=150+mcf_pt_sil_origin;

 
 yfp_pt_sil_origin=find(gradient(mean_yfp_sil)<0.4*(max(gradient(mean_yfp_sil))),1,'first');

subplot(2,2,1)
for k=1:size(rand_origin_cells,2)
        
      tstep=(1:size(mcf_ss(rand_origin_cells(k),:),2))./3;
      patchline(tstep,sgolayfilt(mcf_ss_origin(k,:),1,11),'edgecolor', [  250, 29, 29]./255,'linewidth',1,'edgealpha',0.2),hold on;
      plot(tstep,sgolayfilt(mean(mcf_ss_origin),1,11),'-k','linewidth',4),hold on,
      plot((mcf_pt_sil_origin/3),(mean_mcf_sil_origin(floor(mcf_pt_sil_origin))),'o', 'MarkerSize',10,'MarkerFaceColor', 'k','MarkerEdgeColor','k'),axis tight%%


end

 

subplot(2,2,3)
 for k=1:size(rand_origin_cells,2)
     
     tstep=(1:size(yfp_ss(rand_origin_cells(k),:),2))./3;
      patchline(tstep,sgolayfilt(yfp_ss_origin(k,:),1,31),'edgecolor', [255 169 0]./255,'linewidth',1,'edgealpha',0.2),hold on;
      plot(tstep,sgolayfilt(mean_yfp_sil_origin,1,31),'-k','linewidth',4),hold on,
      plot((yfp_pt_sil_origin/3),(mean_yfp_sil_origin(floor(yfp_pt_sil_origin))),'o', 'MarkerSize',10,'MarkerFaceColor', 'k','MarkerEdgeColor','k'),axis tight%%

% 
 end
 
 
 
%%
[E,index] = sortrows(yfp_ss,309);
mcf_sort=mcf_ss(index,:);
subplot(2,2,2),imagesc(smoothrows(sortrows(yfp_ss,309))),colorbar,caxis([0,2])
subplot(2,2,4),imagesc(smoothrows(sortrows(mcf_ss,309))),colorbar,caxis([0,3])

%%
%plot single cell citrine and mcherry cumulative traces
timestep=(1:309)./3;

for k=1:49
subplot(7,7,k)
plot(timestep,sgolayfilt(yfp_ss(k,:),2,19),'linewidth',2),hold on,yyaxis right, plot(timestep,sgolayfilt(mcf_ss(k,:),2,19),'linewidth',2),grid on,axis tight
end
%%



selectcells = [1:309];

citrine_ns=yfp_ss(selectcells,:);
cherry_ns=mcf_ss(selectcells,:);


timestep=(1:309)./3;

for k=1:25
subplot(5,5,k)
plot(timestep,sgolayfilt(citrine_ns(k,:),1,11),'linewidth',1,'color', [ 255, 202, 24 ]./255),hold on,yyaxis right, plot(timestep,sgolayfilt(cherry_ns(k,:),1,11),'linewidth',1,'color', [ 254, 81, 81 ]./255),grid on,xlim([0 105])
end


%%








