% Simuluation of 2IFC discrimination task QUEST

% This is a simulation of the two-interval forced choice orientation 
% discrimination task. This simulation will look at whether the slope of 
% the psychometric function is perfectly predicted by the criterion 
% threshold. It will also look at which staircase (QUEST or PSI) produce
% more accurate estimates of the parameters.

% This code is a simulation of Experiment 1 of WN's PhD
% This code is a simulation using the QUEST staircase.

% WN started writing this 10 June 2015
% This code was changed to generate the responses from the parameters of
% the Weibull psychometric function fit from the simulation. The
% probability is the derivative of the cdf.

% ----------------------------------------------------------------------- %

clear all;
IDinVWM_Simulation;

simulation.nTrialsPerStaircase = 100;
allRsquared = NaN(1,simulation.nTrialsPerStaircase);

simulation.nStaircases = 100;
simulation.nTrials = simulation.nTrialsPerStaircase*simulation.nStaircases;

simulation.nIntervals = 2;                      % 2 pairs of Gabors shown
simulation.nItems = 2;                          % 2 Gabors shown on each interval

% Set up data matrices

gabor.baseOrientation = NaN(nSDValues,simulation.nTrials,simulation.nIntervals);
allOrientations = NaN(nSDValues,simulation.nTrials,simulation.nIntervals,simulation.nItems);
allOrientationDifferences = NaN(nSDValues,simulation.nStaircases, simulation.nTrialsPerStaircase);
allDifferent = NaN(nSDValues,simulation.nStaircases,simulation.nTrialsPerStaircase);
allEstimates = NaN(nSDValues,simulation.nStaircases,simulation.nTrials,simulation.nIntervals,simulation.nItems);
allResponses = NaN(nSDValues,simulation.nStaircases,simulation.nTrials);
allCorrect = NaN(nSDValues,simulation.nStaircases,simulation.nTrials);
allThresholds = NaN(nSDValues,simulation.nStaircases);

for thisSD = 1:nSDValues

    % Define probability density functions

    mu = 0;     % Mean of the pdf (centred at orientation of stimulus)
    %sd = 10;    % Standard deviation of the pdf

    % Set up simulation parameters


    % Set up Weibull PF parameters

    PFparams = [alphaValues(thisSD) betaValues(thisSD) 0.5 lapseValues(thisSD)];

    % Initialise staircase parameters
    staircase.pThreshold = .82; % Optimal for 2-AFC
    staircase.tGuess = log10(PAL_Weibull(PFparams,staircase.pThreshold, 'Inverse')); % Gives the threshold from the Weibull PF fit
    staircase.tGuessSD = 0.5; % Not sure if this is appropriate
    staircase.beta = 3.5; % Standard
    staircase.delta = lapseValues(sd); % Lapse rate
    staircase.gamma = .5; % Guess rate, 2-AFC

    for thisStaircase = 1:simulation.nStaircases

        allQUESTData(thisSD,thisStaircase) = QuestCreate(staircase.tGuess,staircase.tGuessSD,staircase.pThreshold,staircase.beta,staircase.delta,staircase.gamma);
        allDifferent(thisSD,thisStaircase,:) = mod(randperm(simulation.nTrialsPerStaircase),2)+1;    % Set up which interval will be different on each trial
        
    end

    % Simulation

    for thisStaircase = 1:simulation.nStaircases

        for thisTrial = 1:simulation.nTrialsPerStaircase

        whichIsDifferent = allDifferent(thisSD,thisStaircase,thisTrial);     % Determine which interval is different on this trial
        allThresholds(thisSD,thisStaircase,thisTrial) = staircase.tGuess;
        
        % Get Orientation Difference

        thisLogOrientationDifference = QuestMean(allQUESTData(thisSD,thisStaircase));
        thisOrientationDifference = 10^thisLogOrientationDifference;

        if thisLogOrientationDifference > log10(90)

            thisOrientationDifference = 90;
            thisLogOrientationDifference = log10(90);

        end

        allOrientationDifferences(thisSD,thisStaircase,thisTrial) = thisOrientationDifference;
        allOrientationDifferencesPlot(thisStaircase,thisTrial) = thisOrientationDifference;
        
        % Generate response to this trial from Weibull pdf

        pSuccess = PAL_Weibull(PFparams,thisOrientationDifference); % Probability of correct response

        if pSuccess > 1 % Sometimes pSuccess returns as > 1 due to Weibull fit + lapseRate

        pSuccess = 1;

        end

        isCorrect = random('bino',1,pSuccess); % Determines if correct response on trial
        allCorrect(thisSD,thisStaircase,thisTrial) = isCorrect;

        if isCorrect == 0

            response(thisSD,thisStaircase,thisTrial) = 3-whichIsDifferent;       % If wrong, generate wrong response

        elseif isCorrect == 1

            response(thisSD,thisStaircase,thisTrial) = whichIsDifferent;         % If right, generate right response

        end

        allResponses(thisSD,thisStaircase,thisTrial)=response(thisSD,thisStaircase,thisTrial);

        % Update staircase
        allQUESTData(thisSD,thisStaircase) = QuestUpdate(allQUESTData(thisSD,thisStaircase), thisLogOrientationDifference, isCorrect);
        allQUESTMean(thisSD,thisStaircase, thisTrial) = QuestMean(allQUESTData(thisSD,thisStaircase));
        
        end

        % Plot staircase results

        % Plot of orientation differences tested

%         plotName = ['Orientation Differences tested for Staircase: ' num2str(thisStaircase) ' at SD: ' num2str(thisSD)];
%         figure('Color', 'White', 'Name', plotName);
%         hold on;
%         plot(1:simulation.nTrialsPerStaircase,allOrientationDifferencesPlot(thisStaircase,:),'k.');

    end

end

for thisTrial = 1:simulation.nTrialsPerStaircase
    
    % Get R^2 value

    x = reshape(10.^(allQUESTMean(:,:,thisTrial)),1,nSDValues*simulation.nStaircases);
    y = reshape(10.^(allThresholds(:,:,thisTrial)),1,nSDValues*simulation.nStaircases);
    model = fitlm(x,y);
    allRSquared(1,thisTrial) = model.Rsquared.Ordinary;
    
end

% 
% % Plot of Quest Mean (threshold) for each staircase
% 
% plotNameOverall = ['Quest Mean for each Staircase. Number of Trials = ' num2str(simulation.nTrialsPerStaircase)];
% figure('Color', 'White', 'Name', plotNameOverall);
% plot(sdValues, 10.^allQUESTMean, 'k.');
% hold on;
% plot(sdValues, mean(10.^allQUESTMean'),'r-');



