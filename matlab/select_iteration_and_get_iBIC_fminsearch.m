clear

%dtset = {'gonogo_data_valentina'};

dtset = {'gonogo_data_jana','gonogo_data_janagroup1','gonogo_data_janagroup2'}
dt = '190605_EM_';

fm{1} = 'llba'; % alpha and beta
fm{2} = 'llbax'; % alpha, beta and noise
fm{3} = 'llbaxb'; % alpha, beta, noise and go bias
fm{4} = 'llbaepxb'; % alpha, beta, noise go bias and pavlovian bias
fm{5} = 'll2baxb';  % alpha, 2 betas, noise, constant go bias
fm{6} = 'll2baepxb'; % alpha, 2 betas (sensitive to reward/punishment), noise, constant go bias, pavlovian
fm{7} = 'll2baepcxb'; % alpha, 2 betas, pavlovian (constant value after encountering win/lose), noise, go bias
% fm{8} = 'll2baxbkwinspos'; % alpha, beta, noise, constant go bias, kappa for wins
% fm{9} = 'll2baxbkpos'; % alpha, 2 betas, noise, constant go bias, kappa for win+lose
% fm{10} = 'll2baxbkwinseppos'; % alpha, 2 betas, noise, constant go bias, kappa for win and pavlovian
% fm{11} = 'll2baxbkeppos'; % alpha, 2 betas, noise, constant go bias, kappa for win+lose and pavlovian
% fm{12} = 'll2baepcxbkpos'; % alpha, 2 betas, pavlovian (constant value after encountering win/lose), noise, go bias
% fm{13} = 'll2baepcxbkwinspos'; % alph

% fm{1} = 'll2baxbkwins'; % alpha, beta, noise, constant go bias, kappa for wins
% fm{2} = 'll2baxbk'; % alpha, 2 betas, noise, constant go bias, kappa for win+lose
% fm{3} = 'll2baxbkwinsep'; % alpha, 2 betas, noise, constant go bias, kappa for win and pavlovian
% fm{4}= 'll2baxbkep'; % alpha, 2 betas, noise, constant go bias, kappa for win+lose and pavlovian
% %     %model 11-13 were prompted by a reviewer cmoment to betts et al
% fm{5}= 'll2baepcxb'; % alpha, 2 betas, pavlovian (constant value after encountering win/lose), noise, go bias
% fm{6}= 'll2baepcxbk'; % alpha, 2 betas, pavlovian (constant value after encountering win/lose), noise, go bias
% fm{7}= 'll2baepcxbkwins';

% Npar=[2 3 4 5 5 6 6 6 7 7 6 7 7 6 6 7 7 7 7 5 5];
% Npar=[2 3 4 5 5 6 6 6 6 7 7 7 7];
%Npar=[6 6 6 7 7 7 7];
Npar=[2 3 4 5 5 6 6];

docomp = 1;
Nsample = 2000;

for setind = 2:length(dtset)
    load(dtset{setind})
    for i = 1:length(R)
        r=R{i};
        r=r./abs(r);
        r(isnan(r))=0;
        R{i}=r;
    end
    dtsetnm = strsplit(dtset{setind}, '_');
    dtsetnm = dtsetnm{3};
    
    
    for model=1:size(fm,2)
        pl=[];
        for ite=1:10
            load ([dtsetnm filesep dt '-' fm{model} '-ite' num2str(ite)]);
            pl(ite)=sum(PL(1:end));
        end
        
        [x y]=min(pl);
        SelectIte(model)=y;
        load ([dtsetnm filesep dt '-' fm{model} '-ite' num2str(y)]);
        Np = Npar(model);
        Nsj = length(A);
        ld = [dt '-' fm{model} '-ite' num2str(y)];
        
        if docomp
            oo = ones(1,Nsample);
            muo = mu*oo; nuo = nu*oo;
            
            rand('seed',sum(100*clock));
            LLi=zeros(Nsample,1);
            ddl = zeros(Np,Nsj);
            
            for sj=1:length(A);
                fprintf('subject %i \r',sj)
                
                a=A{sj};
                r=R{sj};
                s=S{sj};
                
                est = diag(sqrt(nu))*randn(Np,Nsample)+mu*ones(1,Nsample);
                for k=1:Nsample;
                    LLi(k)=feval(fm{model},est(:,k),a,r,s,Z,0);
                    if ~mod(k,50);fprintf('model %i subject %i sample %i\r',model,sj,k);end
                end
                iL(sj) = log(sum(exp(-LLi))/Nsample);
            end
            
            
            % compute integrated BIC
            for k=1:length(A); Nch(k)=length(A{k});end; Nch = Nch(:);
            bici  = -2*sum(iL)   + 2*Np*log(sum(Nch)) % integrated BIC
            
            eval(['save ' dtsetnm filesep ld '_meanerror bici iL mu nu E ld fm model']);
            clear iL bici mu nu E ld
        end
        
    end
    clear A
end