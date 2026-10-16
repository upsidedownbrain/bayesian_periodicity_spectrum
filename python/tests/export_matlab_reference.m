referenceDir=fileparts(mfilename('fullpath'));
addpath(fullfile(referenceDir,'..','..','matlab'),'-begin');
rng(314); Fs=40; signal=randn(1201,1)+sin(2*pi*7*(0:1200)'/Fs);
cases=cell(1,5);
for k=1:5
    method='STFT'; window=40; overlapLength=39; FFTLength=40; VoicesPerOctave=16; window_BF=100;
    if k==2, window=41;overlapLength=34;FFTLength=63;end
    if k==3, overlapLength=30;end
    if k>=4, method='CWT';end
    if k==5, VoicesPerOctave=12;end
    if strcmp(method,'STFT')
        [~,f]=stft(signal,Fs,'Window',hann(window,'periodic'),'OverlapLength',overlapLength,'FFTLength',FFTLength);
        f=f(f>0&f<Fs/2);
    else
        [~,f]=cwt(signal,'amor',Fs,'VoicesPerOctave',VoicesPerOctave);f=f(f>=1);
    end
    base=struct('mean_psi_non_osc_fitted',linspace(.1,.5,numel(f))');
    result=periodicity_analysis(signal,Fs,window_BF,method,base,'window',window,'overlapLength',overlapLength,'FFTLength',FFTLength,'VoicesPerOctave',VoicesPerOctave,'CI',true);
    result_uncorrected=periodicity_analysis(signal,Fs,window_BF,method,base,'window',window,'overlapLength',overlapLength,'FFTLength',FFTLength,'VoicesPerOctave',VoicesPerOctave,'CI',false,'BF_corr',false);
    cases{k}=struct('method',method,'window',window,'overlapLength',overlapLength,'FFTLength',FFTLength,'VoicesPerOctave',VoicesPerOctave,'window_BF',window_BF,'base',base,'result',result,'uncorrected',result_uncorrected);
end
Tgrid=[-1e20 -100 -12 -2 0 2 12 100 1e20]; Ngrid=[100 30 20 10 20 10 20 30 100]; tails={'right','left','both'};
direct=cell(1,3);
for k=1:3
    [b,p,t,n,l,d]=ttest_bf('T',Tgrid,'N',Ngrid,'tail',tails{k});
    direct{k}=struct('b',b,'p',p,'t',t,'n',n,'l',l,'d',d);
end
X=randn(12,4)+.4;Y=randn(12,4);X(2,1)=NaN;Y(3,2)=NaN;
data=cell(1,3);
for k=1:3
    if k==1,[b,p,t,n,l]=ttest_bf(X,'rho1',[.5 -.5 .99 0],'tail','right');
    elseif k==2,[b,p,t,n,l]=ttest_bf(X,Y,'rho1',[.5 -.5 .99 0],'tail','left');
    else,[b,p,t,n,l]=ttest_bf(X,Y(1:8,:),'paired',false,'tail','both');end
    data{k}=struct('b',b,'p',p,'t',t,'n',n,'l',l);
end
CI_X=randn(9,4);[mu,ci,df,rho]=credible_interval_individual_spectrum(CI_X,1:4);
CI_ref=struct('mu',mu,'ci',ci,'df',df,'rho',rho);
oddsX=[-20 1 3;2 -3 4;0 7 -5];oddsF=[1;3;9];
smooth=posterior_odds_smoothing(oddsX,oddsF,'KernelSigma',2);
temperature=posterior_odds_temperature(oddsX,'Temp',.7);
comparison=spectra_comparison(X.',Y.','right');
noise=randn(601,10);simulation_ref=cell(1,2);
for k=1:2
    method='STFT';if k==2,method='CWT';end
    raw=[];
    for j=1:10
        sn=noise(:,j);sn=(sn-mean(sn))/std(sn);
        if k==1
            [z,f]=stft(sn,20,'Window',hann(20,'periodic'),'OverlapLength',17,'FFTLength',20);keep=f>0&f<10;z=z(keep,:);f=f(keep);hop=3;
        else
            [z,f]=cwt(sn,'amor',20,'VoicesPerOctave',16);keep=f>=1;z=z(keep,:);f=f(keep);hop=1;
        end
        omega=2*pi*f/20;dphi=diff(angle(z),1,2);fi=((omega+(mod(dphi-omega*hop+pi,2*pi)-pi)/hop)*20/(2*pi)).';
        count=floor(size(fi,1)*hop/50);edges=ceil((0:count)*50/hop);psi=[];
        for q=1:count,psi(:,q)=1./std(fi(edges(q)+1:edges(q+1),:)).';end
        raw(:,j)=mean(psi,2);
    end
    mean_psi=mean(raw,2);
    if k==1,fit_mean=mean(mean_psi(2:end-1));fitted=fit_mean*ones(size(mean_psi));
    else,model=fit(log(f),log(mean_psi),'poly6');fit_mean=coeffvalues(model);fitted=exp(feval(model,log(f)));end
    simulation_ref{k}=struct('freqs',f,'mean_psi_non_osc',mean_psi,'fit_mean',fit_mean,'mean_psi_non_osc_fitted',fitted);
end
save(fullfile(referenceDir,'python_reference.mat'),'signal','Fs','cases','Tgrid','Ngrid','direct','X','Y','data','CI_X','CI_ref','oddsX','oddsF','smooth','temperature','comparison','noise','simulation_ref','-v7');
fprintf('MATLAB reference exported: five pipelines, three tails, ESS, missing data, CI and postprocessing.\n');
