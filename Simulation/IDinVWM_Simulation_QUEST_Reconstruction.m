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

% Define probability density functions

mu = 0;     % Mean of the pdf (centred at orientation of stimulus)
sd = 10;    % Standard deviation of the pdf

% Set up simulation parameters

simulation.nTrialsPerStaircase = 100;
simulation.nStaircases = 1;
simulation.nTrials = simulation.nTrialsPerStaircase*simulation.nStaircases;

simulation.nIntervals = 2;                      % 2 pairs of Gabors shown
simulation.nItems = 2;                          % 2 Gabors shown on each interval

% Set up data matrices

gabor.baseOrientation = NaN(simulation.nTrials,simulation.nIntervals);
allOrientations = NaN(simulation.nTrials,simulation.nIntervals,simulation.nItems);
allOrientationDifferences = NaN(simulation.nStaircases, simulation.nTrialsPerStaircase);
allDifferent = NaN(simulation.nStaircases,simulation.nTrialsPerStaircase);
allEstimates = NaN(simulation.nTrials,simulation.nIntervals,simulation.nItems);
allResponses = NaN(1,simulation.nTrials);

% Initialise staircase parameters
staircase.tGuess = 20;
staircase.tGuessSd = 1.0;
staircase.pThreshold = .82; % Optimal for 2-AFC
staircase.beta = 3.5; % Standard
staircase.delta = .01; % Lapse rate
staircase.gamma = .5; % Guess rate, 2-AFC

for thisStaircase = 1:simulation.nStaircases

    allData(thisStaircase) = QuestCreate(staircase.tGuess,staircase.tGuessSd,staircase.pThreshold,staircase.beta,staircase.delta,staircase.gamma);
    allDifferent(thisStaircase,:) = mod(randperm(simulation.nTrialsPerStaircase),2)+1;    % Set up which interval will be different on each trial
    
end

% Simulation

% Randomise base orientations

for thisInterval = 1:simulation.nIntervals

    gabor.baseOrientation(:,thisInterval) = round(rand(1,simulation.nTrials)*180);

end

% To generate responses from normal distributions, use this:

for thisStaircase = 1:simulation.nStaircases
    
    for thisTrial = 1:simulation.nTrialsPerStaircase

    whichIsDifferent = allDifferent(thisStaircase,thisTrial);     % Determine which interval is different on this trial

    % Get Orientation Difference

    thisOrientationDifference = QuestMean(allData(thisStaircase));
    allOrientationDifferences(thisStaircase,thisTrial) = thisOrientationDifference;

    if thisOrientationDifference > 180
        thisOrientationDifference = 180;
        thisLogOrientationDifference = log10(180);
    end

    for thisInterval = 1:simulation.nIntervals

        % Generate orientations of both items

        for thisItem = 1:simulation.nItems-1

            orientation(thisItem) = gabor.baseOrientation(thisTrial,thisInterval);
            allOrientations(thisTrial,thisInterval,thisItem) = orientation(thisItem);

        end

        if whichIsDifferent==thisInterval

            % This interval is "different". Generate a different second orientation

            orientation(simulation.nItems) = orientation(thisItem)+thisOrientationDifference;
            allOrientations(thisTrial,thisInterval,simulation.nItems) = orientation(simulation.nItems);

        else

            % This interval is "same'. Generate the same orientation

            orientation(simulation.nItems) = orientation(thisItem);
            allOrientations(thisTrial,thisInterval,simulation.nItems) = orientation(simulation.nItems);

        end

        for thisItem = 1:simulation.nItems

            thisEstimate(thisItem) = (randn*sd)+(orientation(thisItem));
            allEstimates(thisTrial,thisInterval,thisItem) = thisEstimate(thisItem);

        end

        % Find difference between these estimates

        estimateDiff(thisTrial,thisInterval) = abs(thisEstimate(1)-thisEstimate(2));

    end

    % Generate response to this trial

    if estimateDiff(thisTrial,1) > estimateDiff(thisTrial,2)

        response(thisTrial) = 1;        % First interval had larger difference from estimates

    elseif estimateDiff(thisTrial,1) < estimateDiff(thisTrial,2)

        response(thisTrial) = 2;        % Second interval had larger difference from estimates

    end

    responseCorrect = response(thisTrial)==allDifferent(thisTrial);
    allResponses(thisStaircase,thisTrial)=response(thisTrial);

    % Update staircase
    allData(thisStaircase) = QuestUpdate(allData(thisStaircase), thisOrientationDifference, responseCorrect);

    end
    
end

for thisStaircase = 1:simulation.nStaircases
    
    % Plot staircase results
    
    figure('Color', 'White', 'Name', 'QuestMean');
    hold on;
    plot(1:simulation.nTrialsPerStaircase,allOrientationDifferences(thisStaircase,:),'k.');

end


