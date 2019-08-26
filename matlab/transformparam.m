% script that extracts BICs for each model that's compared and then
% transforms and saves the parameters for that
clear
%correctdir='/home/ALDRECENTRUM/lieke.deboer/Experiments/DAD2/kappa/adams/';
%correctdir='\\PSYG8\lieke.deboer\Experiments\DAD2\kappa\lieke\';

%cd(correctdir)
%fbest = dir([correctdir, '*meanerror.mat']);
fbest = dir(['*meanerror.mat']);
fbdate = cellstr(char(fbest.date));
fbest = cellstr(char(fbest.name));
% load /home/ALDRECENTRUM/lieke.deboer/Experiments/DAD2/modelling_gonogo_2017/matgp/llb2aepxb_allgood.mat
%fbeh= dir([correctdir, '*gonogo_data*']);
fbeh= dir(['*gonogo_data*']);

fbeh = char(fbeh.name);
load(fbeh)

for i = 1:length(fbest)
    d = load(fbest{i});
    BICs(i) = d.bici;
    LL(i) = sum(d.iL);
    b=strsplit(fbest{i}, '-');
    nam{i}=char(b(2));
    decpt = sum(cellfun(@length, A));
    llalt = decpt * log(0.5);
    llmd = sum(d.iL);
    rsqr(i) = (llalt-llmd)/llalt;
end

[BICsort,order] = sort(BICs);
modellist = nam(order)'; % list of the models, lowest BIC first, highest last
Rsqr = rsqr(order)';
clear b d BICs nam
modLL = LL(order)';
datemodelled = fbdate(order);
BIComp = [BICsort' ([NaN; diff(BICsort)'])];  % the actual BICs

deltBIC = [NaN; cumsum(BIComp(2:end,2))];
BICstats = [BIComp deltBIC modLL Rsqr];

format shortG
BICstats %#ok<NOPTS> % print the model stats in the console window
modellist %#ok<NOPTS>

save modcomp modellist BICstats
format short

if ~exist('list', 'var')
    list= 'this variable is empty because the original file did not have a list of subjects in it';
end

for i = 1:length(fbest)
    load(fbest{i});
    b=strsplit(fbest{i}, '-');
    nam=char(b(2));
    
    if strcmp('llba', nam)
        beta  = exp(E(1,:));            % sensitivity to reward
        alfa  = 1./(1+exp(-E(2,:)));
        save(nam, 'beta','alfa','list')
        clear(nam, 'beta','alfa')
        
    elseif strcmp('llbax', nam)
        beta  = exp(E(1,:));            % sensitivity to reward
        alfa  = 1./(1+exp(-E(2,:)));
        g     = 1./(1+exp(-E(3,:)));
        save(nam, 'beta','alfa', 'g', 'list')
        clear(nam, 'beta','alfa', 'g')
        
        
    elseif strcmp('llbaxb', nam)
        beta  = exp(E(1,:));            % sensitivity to reward
        alfa  = 1./(1+exp(-E(2,:)));
        g     = 1./(1+exp(-E(3,:)));
        bias  = E(4,:);
        save(nam, 'beta','alfa', 'g', 'bias', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias')
        
    elseif strcmp('llbaepxb', nam)
        beta  = exp(E(1,:));            % sensitivity to reward
        alfa  = 1./(1+exp(-E(2,:)));
        pav   = exp(E(3,:));
        g     = 1./(1+exp(-E(4,:)));
        bias  = E(5,:);
        save(nam, 'beta','alfa', 'g', 'bias', 'pav', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'pav')
        
    elseif strcmp('ll2baxb', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));
        bias      = E(5,:);
        save(nam, 'beta','alfa', 'g', 'bias', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias')
        
    elseif strcmp('ll2baepxb', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        pav       = exp(E(4,:));
        g         = 1./(1+exp(-E(5,:)));
        bias      = E(6,:);
        save(nam, 'beta','alfa', 'g', 'bias', 'pav', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'pav')
        
    elseif strcmp('ll2baxbk', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        alfanogo  = 1./(1+exp(-(E(3,:)-E(6,:))));
        alfago    = 1./(1+exp(-(E(3,:)+E(6,:)))); % instrumental bias on alfa
        kappa     = E(6,:);
        save(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'alfago', 'alfanogo', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'alfago', 'alfanogo')
        
        
    elseif strcmp('ll2baxbkep', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        alfanogo  = 1./(1+exp(-(E(3,:)-E(6,:))));
        alfago    = 1./(1+exp(-(E(3,:)+E(6,:)))); % instrumental bias on alfa
        kappa     = E(6,:);
        pav       = exp(E(7,:));
        save(nam, 'beta','alfa', 'g', 'bias', 'pav', 'kappa', 'alfago', 'alfanogo', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'pav', 'kappa', 'alfago', 'alfanogo')
        
        
    elseif strcmp('ll2baxbkwins', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        alfago    = 1./(1+exp(-(E(3,:)+E(6,:)))); % instrumental bias on alfa
        kappa     = E(6,:);
        save(nam, 'beta','alfa', 'g', 'bias','kappa', 'alfago', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias','kappa', 'alfago')
        
        
    elseif strcmp('ll2baxbkwinsep', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        alfago    = 1./(1+exp(-(E(3,:)+E(6,:)))); % instrumental bias on alfa
        kappa     = E(6,:);
        pav       = exp(E(7,:));
        save(nam, 'beta','alfa', 'g', 'bias', 'pav', 'kappa', 'alfago', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'pav', 'kappa', 'alfago')
        
    elseif strcmp('ll2baxbbehact_gobonus', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        gobonus   = E(6,:);
        
        save(nam, 'beta','alfa', 'g', 'bias', 'gobonus', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias',  'gobonus')
        
    elseif strcmp('ll2baxbbehact_goandnogo', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        gobonus     = E(6,:);
        nogobonus     = E(7,:);
        
        save(nam, 'beta','alfa', 'g', 'bias', 'gobonus', 'nogobonus','list')
        clear(nam, 'beta','alfa', 'g', 'bias',  'gobonus', 'nogobonus')
        
    elseif strcmp('ll2baepcxbkwins', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        alfago    = 1./(1+exp(-(E(3,:)+E(6,:)))); % instrumental bias on alfa
        kappa     = E(6,:);
        pav       = exp(E(7,:));
        save(nam, 'beta','alfa', 'g', 'bias', 'pav', 'kappa', 'alfago', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'pav', 'kappa', 'alfago')
        
    elseif strcmp('ll2baepcxb', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        pav       = exp(E(4,:));         % irreduceable noise
        g	      = 1./(1+exp(-E(5,:)));                   % go bias
        bias      = E(6,:);
        save(nam, 'beta','alfa', 'g', 'bias', 'pav', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'pav')
        
        
    elseif strcmp('ll2baepcxbk', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        alfanogo  = 1./(1+exp(-(E(3,:)-E(6,:))));
        alfago    = 1./(1+exp(-(E(3,:)+E(6,:)))); % instrumental bias on alfa
        kappa     = E(6,:);
        pav       = exp(E(7,:));
        save(nam, 'beta','alfa', 'g', 'bias', 'pav', 'kappa', 'alfago', 'alfanogo', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'pav', 'kappa', 'alfago', 'alfanogo')
        
    elseif strcmp('ll2baxbkwinspos', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        alfago    = 1./(1+exp(-(E(3,:)+E(6,:)))); % instrumental bias on alfa
        kappa     = exp(E(6,:));
        save(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'alfago', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'alfago')
        
    elseif strcmp('ll2baxbkpos', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        alfago    = 1./(1+exp(-(E(3,:)+E(6,:)))); % instrumental bias on alfa
        alfanogo    = 1./(1+exp(-(E(3,:)-E(6,:)))); % instrumental bias on alfa
        kappa     = exp(E(6,:));
        save(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'alfago', 'alfanogo', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'alfago','alfanogo')
        
    elseif strcmp('ll2baxbkwinseppos', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        alfago    = 1./(1+exp(-(E(3,:)+E(6,:)))); % instrumental bias on alfa
        kappa     = exp(E(6,:));
        pav       = exp(E(7,:));
        save(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'pav', 'alfago', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'pav', 'alfago')
        
    elseif strcmp('ll2baxbkeppos', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        alfago    = 1./(1+exp(-(E(3,:)+E(6,:)))); % instrumental bias on alfa
        alfanogo  = 1./(1+exp(-(E(3,:)-E(6,:)))); % instrumental bias on alfa
        kappa     = exp(E(6,:));
        pav       = exp(E(7,:));
        save(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'alfago', 'alfanogo','pav', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'alfago','alfanogo','pav')
    
     elseif strcmp('ll2baepcxbkpos', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        alfago    = 1./(1+exp(-(E(3,:)+E(6,:)))); % instrumental bias on alfa
        alfanogo  = 1./(1+exp(-(E(3,:)-E(6,:)))); % instrumental bias on alfa
        kappa     = exp(E(6,:));
        pav       = exp(E(7,:));
        save(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'alfago', 'alfanogo','pav', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'alfago','alfanogo','pav')
        
        
    elseif strcmp('ll2baepcxbkwinspos', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        g         = 1./(1+exp(-E(4,:)));         % irreduceable noise
        bias	  = E(5,:);                   % go bias
        alfago    = 1./(1+exp(-(E(3,:)+E(6,:)))); % instrumental bias on alfa
        kappa     = exp(E(6,:));
        pav       = exp(E(7,:));
        save(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'pav', 'alfago', 'list')
        clear(nam, 'beta','alfa', 'g', 'bias', 'kappa', 'pav', 'alfago')
        
    elseif strcmp('ll2baepcb', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        bias	  = E(5,:);                   % go bias
        pav       = exp(E(4,:));
        save(nam, 'beta','alfa', 'bias',  'pav','list')
        clear(nam, 'beta','alfa', 'bias', 'pav')
       
    elseif strcmp('ll2baepb', nam)
        beta 	  = exp(E(1:2,:));            % sensitivity to reward
        alfa 	  = 1./(1+exp(-E(3,:)));     % learning rate
        bias	  = E(5,:);                   % go bias
        pav       = exp(E(4,:));
        save(nam, 'beta','alfa', 'bias',  'pav','list')
        clear(nam, 'beta','alfa', 'bias', 'pav')
        
    else
        warning(['asked for a transformation without defining model for: ' nam ])
    end
    
end
