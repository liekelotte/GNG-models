%%%%%%%%%%%%%%%%%
%%% Written by Quentin J. M. Huys, UCL, 2011, modified by Lieke de
%%% Boer
%%% References:
%%% Guitart-Masip M, Quentin JM, Fuentemilla LL, Dayan P, Duzel E, Dolan RJ (2012)
%%% Go and no-go learning in reward and punishment: Interaction between affect and effect NeuroImage doi:10.1016/j.neuroimage.2012.04.024
%%% de Boer, L., J Axelsson, R Chowdhury, K Riklund, RJ Dolan, L Nyberg, L
%%% B�ckman, M Guitart-Masip (2019).Dorsal striatal dopamine D1 receptor availability predicts an
%%% instrumental bias in action learning. PNAS doi: 10.1073/pnas.1816704116

clear;

%% to start parfor
% c=parcluster;
% c.NumWorkers=16 ; % 22 is the maxnumber of workers you're allowed
% parpool(16) ; % same thing

dt = '190605_EM_';

dtset = {'gonogo_data'}; % put in the name of your datafile(s) (more than one is allowed). should contain A, R and S, and 
                        % (ideally) a subject list. 


for setind = 1:length(dtset)
    
    load(dtset{setind}) %make all rewards scale to 1, 0 and -1
    for i = 1:length(R)
        r=R{i};
        r=r./abs(r);
        r(isnan(r))=0;
        R{i}=r;
    end
    dtsetnm = strsplit(dtset{setind}, '_');
    dtsetnm = dtsetnm{3};
    if ~exist(dtsetnm, 'file')
        mkdir(dtsetnm)
    end
    
    dosave = 1;
    docomp = 0;
    docheck = 0;
    Nsample = 2000;
    
    ff{1} = 'llba'; % alpha and beta
    ff{2} = 'llbax'; % alpha, beta and noise
    ff{3} = 'llbaxb'; % alpha, beta, noise and go bias
    ff{4} = 'llbaepxb'; % alpha, beta, noise go bias and pavlovian bias
    ff{5} = 'll2baxb';  % alpha, 2 betas, noise, constant go bias
    ff{6} = 'll2baepxb'; % alpha, 2 betas (sensitive to reward/punishment), noise, constant go bias, pavlovian
    ff{7} = 'll2baepcxb'; % alpha, 2 betas, pavlovian (constant value after encountering win/lose), noise, go bias
    ff{8} = 'll2baxbkwins'; % alpha, beta, noise, constant go bias, kappa for wins
    ff{9} = 'll2baxbk'; % alpha, 2 betas, noise, constant go bias, kappa for win+lose
    ff{10} = 'll2baxbkwinsep'; % alpha, 2 betas, noise, constant go bias, kappa for win and pavlovian
    ff{11} =  'll2baxbkep'; % alpha, 2 betas, noise, constant go bias, kappa for win+lose and pavlovian
    ff{12} = 'll2baxbkwinsepc'; % alpha, 2 betas, noise, constant go bias, kappa for win and pavlovian (constant)
    ff{13} =  'll2baxbkepc'; % alpha, 2 betas, noise, constant go bias, kappa for win+lose and pavlovian (constant)
    
    Npar=[2 3 4 5 5 6 6 6 6 7 7 7 7]; % specify the number of parameters for each model
    
    options=optimset('display','off','DerivativeCheck','on');
    warning('off','optim:fminunc:SwitchingMethod')
    for whichinf=1:size(ff,2)
        for ite=1:10
            
            exx=[]; E=[]; V=[]; PL=[]; mu=[]; nu=[]; par=[]; lt=[]; et=[]; LLi=[]; iL=[]; bici=[];
            
            Np = Npar(whichinf);
            
            ld = [dt '-' ff{whichinf} '-ite' num2str(ite)];
            
            Nsj=length(A);
            Z.mu=zeros(Np,1);
            Z.nui=eye(Np);
            init=.1*randn(Np,1);
            
            E=zeros(Np,Nsj);
            V=zeros(Np,Nsj);
            PL=zeros(1,Nsj);
            exx=zeros(1,Nsj);
            LL=zeros(1,Nsj);
            
            init
            
            emit=0;
            while 1;emit=emit+1;sj=0;
                
                % E step......................................................
                
                parfor sj=1:Nsj
                    fprintf('%2d\n',sj);
                    if isempty(A{sj})==0;
                        a=A{sj};
                        r=R{sj};
                        s=S{sj};
                        
                        
                        init=.1*randn(Np,1);
                        if emit>1; init=E(:,sj); end
                        warning('off','optim:fminunc:SwitchingMethod')
                        [est,fval,ex]=fminsearch(@(x) feval(ff{whichinf},x, a, r, s, Z,1),init,options);
                        
                        mf=fval;mE=est;mx=ex;
                        for i=1:0
                            [est,fval,nex]=fminsearch(@(x) ...
                                feval(ff{whichinf},x, a, r, s, Z,1),init+.1*randn(Np,1),options);
                            if (fval < mf)
                                mf=fval;
                                mE=est;
                                mx=nex;
                            end
                        end
                        
                        fval=mf;
                        est=mE;
                        ex=mx;
                        
                        % 		    if ex<0 ; tmp=tmp+1; fprintf('didn''t converge %i times exit status %i\r',tmp,ex); end
                        
                        hess=NumHessian(@(x) feval(ff{whichinf},x, a, r, s, Z,1),est);
                        
                        exx(sj)=ex;
                        E(:,sj)=est;			% Subjets' parameter estimate
                        V(:,sj) = max(diag(inv(hess)),1e-5);	% inverse of Hessian = variance
                        PL(sj) = fval;
                        LL(sj)=feval(ff{whichinf},est, a, r, s, Z,0);
                        
                        fprintf('Emit=%i subject %i model %i iteration %i exit status=%i\r',emit,sj,whichinf,ite,exx(sj))
                    end
                end
                
                % M step using factorized posterior .................................
                mu = mean(E,2);
                %    nu = sqrt(sum(E.^2 + V,2)/Nsj - mu .^2); %change back to the line below lieke
                nu = sqrt(sum(E.^2,2)/Nsj + trimmean(V,10,2) - mu .^2); % fix to eliminate extreme values
                Z.mu = mu; Z.nui = inv(diag(nu));
                
                whichinf
                [mu nu]
                
                par(emit,:) = [sum(LL) sum(PL) mu' nu(:)'];						% save stuff
                et(emit,:,:) = E;
                lt(emit,:,:) = [LL; PL];
                
                if dosave;
                    eval(['save ' dtsetnm filesep ld ' mu nu E V LL PL par et lt Z exx emit ld ff whichinf;']);
                end
                % check convergence ...........................................
                if emit>1;if sum(abs(par(emit,2:end)-par(emit-1,2:end)))<1e-2;fprintf('\n *** converged *** \n');break;end;end
            end
            
        end
    end
end
% quit




