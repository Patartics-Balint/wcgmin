function [WCG, WCU] = wcgainWithBNB(P, varargin)

%% Get names of uncertainties
uname = fieldnames(P.Uncertainty);
Nu = numel(uname);
Nr = 0;
for i=1:Nu
    ui = P.Uncertainty.(uname{i});
    if isa(ui,'ureal')
        Nr=Nr+1;
    end
end

%% Initialize Branch and Bound Parameters
reltol = 0.02;
abstol = 1e-3;
maxLB = 1e7;
count = 0;
LB = 0;
UB = inf;
WCU = [];
CritFreq = NaN;

if Nr==0
    % Don't run branch and bound if all uncertainty is complex
    maxcount = 1;
else
    maxcount = 10;
    %maxcount = 1;
end 

%% Run Branch and Bound
%fprintf('\n')
while (UB-LB)>reltol*LB+abstol && count<maxcount  && LB<maxLB
    count=count+1;
    if count==1
        % Initialize active cubes (for UREAL only/No splitting on ULTIDYN)
        % Note: Active Cubes are stored as USS
        AC = {P};
    else
        % Update Active Cubes: Select cube and side to split         
        
        % Select cube with large UB
        cube2split = AC{UBmaxidx};
        
        % Select largest sidelength
        Range = zeros(Nu,2);
        for i=1:Nu
            ui = cube2split.Uncertainty.(uname{i});
            if isa(ui,'ureal')
                Range(i,:) = ui.Range;
            end
        end
        sidelength = Range(:,2)-Range(:,1);        
        [~,SLmaxidx] = max(sidelength);      
        Range = Range(SLmaxidx,:);
        m = (Range(2)+Range(1))/2;

        Range1 = [Range(1) m];
        del1 = ureal( uname{SLmaxidx}, mean(Range1),'Range',Range1);
        AC1 = usubs(cube2split, uname{SLmaxidx}, del1);
        
        
        Range2 = [m Range(2)];
        del2 = ureal( uname{SLmaxidx}, mean(Range2),'Range',Range2);
        AC2 = usubs(cube2split, uname{SLmaxidx}, del2);
        
        % Add new active cubes to list
        AC(UBmaxidx) = [];
        AC(end+1) = {AC1};
        AC(end+1) = {AC2};        
    end
    
    % Run WCGAIN on all cubes
    Ncube = numel(AC);
    UBcube = zeros(Ncube,1);
    for i=1:Ncube
        % Create plant corresponding to cube i
        Pi = AC{i};
        
        % Call worst-case gain on cube i
        [WCGi,WCUi,INFOi] = wcgain(Pi,varargin{:});
        if WCGi.LowerBound>LB
            LB = WCGi.LowerBound;
            CritFreq = WCGi.CriticalFrequency;
            WCU = WCUi;
        end
        UBcube(i) = WCGi.UpperBound;
    
        % Use improved lower bound
        % (Evaluate only on bad frequencies identified by WCGAIN?)
        if isfinite(LB)
            badFreq = INFOi.Frequency;
            badFreq( isinf(badFreq) ) = [];   % Remove inf frequency
            if ~isempty(badFreq)
                [WCGLB2,WCU2,CritFreq2] = wcgainlb(Pi,badFreq);
                %[WCGLB2,WCU2,CritFreq2] = wcgainlb(Pi,INFOi.Frequency);
                if WCGLB2>LB
                    LB = WCGLB2;
                    CritFreq = CritFreq2;
                    WCU = WCU2;
                end
            end
        end
    
    end
    
    % Compute global upper bound
    % Alternative logic here, e.g. if multiple cubes have max we
    % are simply choosing first one?
    [UB,UBmaxidx] = max( UBcube );    
    
%     fprintf('\n Branch = %d LB = %4.3f UB = %4.3f',count,LB,UB);
end
% fprintf('\n')

%% Store final outputs
WCG.LowerBound = LB;
WCG.UpperBound = UB;
WCG.CriticalFrequency = CritFreq; 