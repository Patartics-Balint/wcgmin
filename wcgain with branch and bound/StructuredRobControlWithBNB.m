%% Structured Robust Control
%  Algorithm to approximate the proposed approach in da Silva, et al.
close all
clear
clc

%% Create a Plant
excase = 3;
switch excase
    case 1
        % Random Example
        rng(0);
        M = rss(7,7,6);
        del1 = ureal('del1',7,'Range',[6 8]);
        del2 = ureal('del2',3,'Range',[2 4]);
        Del = blkdiag(del1*eye(2),del2);
        P = lft(Del,M);
        ny = 1;
        nu = 1;
        
        nxK = 4;
    case 2
        % Flex Aircraft Example from da Silva, et al
        daSilva_Apkaraian_example;
        ny = n_meas;
        nu = n_con;
        nxK = 6;
    case 3
        % Examples from da Silva, et al
        load example3;
        Nblk = size(blk,1);
        Delta = [];
        Nr = 0;
        Nd = 0;
        block_structure = [];
        for i=1:Nblk
            blki = blk(i,:);
            if blki(1)<0
                % UREAL
                Nr = Nr+1;
                I = eye( abs(blki(1,1)) );
                Deltai = ureal(['delta_real_' int2str(Nr)],0)*I;
                block_structure(i, :) = [blki(1, 1), 0];
            else
                % ULTIDYN
                Nd = Nd +1;
                Deltai = ultidyn(['delta_dyn_' int2str(Nd)],blki);
                block_structure(i, :) = [blki(1, 1), blki(1, 1)];
            end
            Delta = blkdiag(Delta,Deltai);
        end
        P0 = P;
        P = lft(Delta,P);
        [nyP, nuP] = size(P);
        nz = nyP - ny; 
        nw = nuP - nu;
        nxK = 4;
end

%% Create a controller structure
K = tunableSS('K', nxK, nu, ny);

%% Initialize the active scenarios
% Start with AS corrresponding to nominal parameter values
AS = P.Uncertainty;
uname = fieldnames(AS);
for i=1:numel(uname)
    ui = AS.(uname{i});
    AS.(uname{i}) = ui.Nominal;
end

%% Iterative HINSTRUCT
tol = 1e-2;
maxcnt = 20;
block_structure = [block_structure; [nw, nz]];

opt = hinfstructOptions;
opt.Display = 'off';

cnt = 0;
GAM = 0;
WCGUB = inf;
fprintf('\n');
while (cnt<maxcnt) && ((abs(GAM-WCGUB)>tol*GAM) || ~isfinite(GAM))
    % Update iteration count
    cnt = cnt+1;
    
    % Create array of active models
    Parr = usubs(P,AS);
    
    %% HINFSTRUCT on Array of Models
    [Kd,GAM,INFO] = hinfstruct(Parr,K,opt);
    Kd = ss(Kd);
    
    %% Robustness analysis
    CLunc = lft(P,Kd);
    %  CLunc = ufrd(CLunc,logspace(-2,2,100));  % Does this call BNB code?
    %wcopt = wcgainOptions('MaxOverFrequency','off');
    %[WCG,WCU,INFO] = wcgain(CLunc,wcopt);
    %WCGLB = max(WCG.LowerBound.ResponseData(:));
    %WCGUB = max(WCG.UpperBound.ResponseData(:));
%     WCF = WCG.CriticalFrequency;
%     if isinf(WCGUB)
%       [SM, WCU] = robstab(CLunc);  
%       WCF = SM.CriticalFrequency;
%     end

    [WCG, WCU] = wcgainWithBNB(CLunc);
%     [~, WCU] = robstab(CLunc);
    WCGLB = WCG.LowerBound;
    WCGUB = WCG.UpperBound;
    WCF = WCG.CriticalFrequency;
    if isinf(WCGUB)
      [SM, WCU] = robstab(CLunc);  
      WCF = SM.CriticalFrequency;
    end
    
    
    %% Add WCU to Active Scenarios
    AS = [AS WCU];
    
    %% Analysis
    %plot_iter_step(P0, Kd, cnt, block_structure, WCU);
    
%     figure(1);
%     semilogx(INFO.Frequency,INFO.Bounds,'b');
%     ax = ylim;
%     hold on; semilogx(wcf*[1 1],ax,'r'); hold off;
%     title(sprintf('WCF = %4.3f',wcf));
%     ylim([0 10]);
%     drawnow;
    
    %% Display Iteration Info
    fprintf('cnt = %d \t GAM = %4.3f \t WCGLB = %4.3f \t WCGUB = %4.3f \t WCF = %4.3f',...
        cnt,GAM,WCGLB,WCGUB,WCF);
    fprintf('\n');
    freqresp(WCU.delta_dyn_1, WCF)
    %pause;
end
% hold off
fprintf('\n');
% plot_iter_step(P0, Kd, cnt, block_structure, WCU);

function plot_iter_step(P0, Kd, cnt, block_structure, WCU)
    n_signals = 4;
    cmap = flip(copper(n_signals));
    freq = logspace(-4, 5, 1000);
    M = lft(P0, Kd);
    mubounds = mussv(frd(M, freq), block_structure);
    if rem(cnt-1, n_signals) == 0
        figure();        
        semilogx(mubounds(1).Frequency([1, end]), [1, 1]);
        hold on
        title([num2str(cnt), ' - ', num2str(cnt+n_signals-1)]);
        xlim(freq([1, end]));
        grid('on');
        xlabel('frequency [rad/s]');
        ylabel('\mu');        
    end
    plot_mu_bound = @(mubounds, idx, cnt, n_signals, cmap_mu)(semilogx(mubounds(1).Frequency,...
        squeeze(mubounds(idx).ResponseData), '-', 'Color', cmap_mu(rem(cnt-1, n_signals)+1, :),...
        'LineWidth', 2));
    plot_mu_bound(mubounds, 1, cnt, n_signals, cmap);
    plot_mu_bound(mubounds, 2, cnt, n_signals, cmap);
%     resp_s = sigma(M, freq);
%     semilogx(freq, max(resp_s), '--', 'Color', cmap(rem(cnt-1, n_signals)+1, :), 'LineWidth', 2);
%     drawnow;       
    resp_u = sigma(ss(struct2array(WCU)), freq);
    semilogx(freq, resp_u, '-', 'Color', cmap(rem(cnt-1, n_signals)+1, :));
%     ylim([0, 1.1*max(squeeze(mubounds(1).ResponseData))]);
%     legend('1', '\mu_{LB}', '\mu_{UB}', '\sigma(M)', '\sigma(\Delta_{WC})');
    legend('1', '\mu_{LB}', '\mu_{UB}', '\sigma(\Delta_{WC})');
    drawnow;
end