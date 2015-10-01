% Simuluation of 2IFC discrimination task

% This is a simulation of the two-interval forced choice orientation 
% discrimination task. This simulation will look at whether the slope of 
% the psychometric function is perfectly predicted by the criterion 
% threshold. It will also look at which staircase (QUEST or PSI) produce
% more accurate estimates of the parameters.

% This code is a simulation of Experiment 1 of WN's PhD
% This code is a simulation using the PSI staircase to estimate the
% parameters of the psychometric function.

% WN started writing this 10 June 2015

% ----------------------------------------------------------------------- %

% clear all;
% IDinVWM_Simulation;

% Set up simulation parameters

simulation.nStaircases = 1;
simulation.nTrialsPerStaircase = 40;
simulation.nTrials = simulation.nTrialsPerStaircase*simulation.nStaircases;

simulation.nIntervals = 2;     % 2 pairs of Gabors shown
simulation.nItems = 2;         % 2 Gabors shown on each interval

% Set up data structures

gabor.baseOrientation = NaN(nSDValues,simulation.nTrials,simulation.nIntervals);
allOrientations = NaN(nSDValues,simulation.nTrials,simulation.nIntervals,simulation.nItems);
allOrientationDifferences = NaN(nSDValues,simulation.nStaircases, simulation.nTrialsPerStaircase);
allDifferent = NaN(nSDValues,simulation.nStaircases,simulation.nTrials);
allEstimates = NaN(nSDValues,simulation.nStaircases,simulation.nTrials,simulation.nIntervals,simulation.nItems);
allResponses = NaN(nSDValues,simulation.nStaircases,simulation.nTrials);
allCorrect = NaN(nSDValues,simulation.nStaircases,simulation.nTrials);
allThresholds = NaN(nSDValues,simulation.nStaircases);


for thisSD = 1:nSDValues

    % Define probability density functions

    mu = 0;     % Mean of the pdf (centred at orientation of stimulus)
    %sd = 10;    % Standard deviation of the pdf

    % set up Weibull PF parameters

    PFparams = [alphaValues(thisSD) betaValues(thisSD) 0.5 lapseValues(thisSD)];

    % Set up staircase parameters

    staircase.alphaRange = linspace(PAL_Weibull(PFparams,.5,'Inverse'),...
        PAL_Weibull(PFparams,.9999,'Inverse'),simulation.nTrialsPerStaircase);              % Range of threshold for orientation differences to be considered in the posterior distribution 
    staircase.betaRange = linspace(log10(15),log10(200),simulation.nTrialsPerStaircase);    % Range for the slope of the Weibull function      
    staircase.gammaRange = .5;                                  % Guess rate for 2IFC
    staircase.lambdaRange = linspace(-.05,.05,simulation.nTrialsPerStaircase);                                     

    staircase.stimRange = linspace(1,90,100);                  % Stimulus range for orientation differences
    staircase.pThreshold = .82; % Optimal for 2-AFC

    for thisStaircase = 1:simulation.nStaircases

        allPSIData(thisSD,thisStaircase) = PAL_AMPM_setupPM('PF', @PAL_Weibull, 'stimRange', staircase.stimRange,...
            'numTrials', simulation.nTrialsPerStaircase, 'gammaEQlambda', 0, ...
            'priorAlphaRange', single(staircase.alphaRange), 'priorBetaRange', single(staircase.betaRange),...
            'priorGammaRange', single(staircase.gammaRange), 'priorLambdaRange', single(staircase.lambdaRange),...
            'marginalize', [-1 -2 3 4]);

        allDifferent(thisSD,thisStaircase,:) = mod(randperm(simulation.nTrialsPerStaircase),2)+1;    % Set up which interval will be different on each trial
        
    end

    % Simulation

    for thisStaircase = 1:simulation.nStaircases

        for thisTrial = 1:simulation.nTrialsPerStaircase

            whichIsDifferent = allDifferent(thisStaircase,thisTrial);     % Determine which interval is different on this trial
            allThresholds(thisSD,thisStaircase,thisTrial) = PAL_Weibull(PFparams,staircase.pThreshold, 'Inverse');
            
           % Get orientation difference

            thisOrientationDifference = allPSIData(thisStaircase).xCurrent;
            allOrientationDifferences(thisSD,thisStaircase,thisTrial) = thisOrientationDifference;

            % Generate response to this trial from Weibull pdf

            pSuccess = PAL_Weibull(PFparams,thisOrientationDifference); % Probability of correct response
            isCorrect = random('bino',1,pSuccess); % Determines if correct response on trial
            allCorrect(thisSD,thisStaircase,thisTrial) = isCorrect;

            if isCorrect == 0

                response(thisTrial) = 3-whichIsDifferent;       % If wrong, generate wrong response

            elseif isCorrect == 1

                response(thisTrial) = whichIsDifferent;         % If right, generate right response

            end

            allResponses(thisStaircase,thisTrial)=response(thisTrial);

            % Update staircase
            allPSIData(thisSD,thisStaircase) = PAL_AMPM_updatePM(allPSIData(thisSD,thisStaircase), isCorrect);
            
            
        end

        nCorrect(thisStaircase) = sum(allResponses(thisStaircase,:)==allDifferent(thisStaircase,:));
        pCorr(thisStaircase) = nCorrect(thisStaircase)/simulation.nTrialsPerStaircase;
        allPSIEstimate(thisSD,thisStaircase,:) = allPSIData(thisSD,thisStaircase).threshold;
        
    end

%     for thisStaircase = 1:simulation.nStaircases
%     
%     % Plot staircase results
%     
%     figure('Color', 'White', 'Name', 'Orientation Differences');
%     hold on;
%     plot(1:simulation.nTrialsPerStaircase,allOrientationDifferences(thisStaircase,:),'k.');
%     
%     figure('Color','White','Name','Staircase Threshold');
%     plot(allPSIData(thisStaircase).threshold,'k-');
%            
%     end
    
end







