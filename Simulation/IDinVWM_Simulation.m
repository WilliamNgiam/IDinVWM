% Simuluation of 2IFC discrimination task

% This is a simulation of the two-interval forced choice orientation 
% discrimination task. This simulation will look at whether the slope of 
% the psychometric function is perfectly predicted by the criterion 
% threshold. It will also look at which staircase (QUEST or PSI) produce
% more accurate estimates of the parameters.

% This code is a simulation of Experiment 1 of WN's PhD
% This code runs the experiment using "method of constants", but the
% constants aren't randomised in each block since this is a simulation.

% WN started writing this 10 June 2015

% ----------------------------------------------------------------------- %

clear all;

% Define probability density functions

mu = 0;     % Mean of the pdf (centred at orientation of stimulus)
% sd = 12;     % Standard deviation of the pdf
maxSD = 3;

sdValues = .3:.1:maxSD;
nSDValues = length(sdValues);
parameterValues = NaN(nSDValues,4);
LLValues = NaN(nSDValues,1);
exitFlagValues = NaN(nSDValues,1);


for thisSD = 1:nSDValues
    
    sd = sdValues(thisSD);
    
    % Set up simulation parameters

    simulation.nBlocks = 5;
    simulation.nTrialsPerBlock = 1000;
    simulation.nTrialsPerDifference = simulation.nBlocks * simulation.nTrialsPerBlock;

    simulation.nIntervals = 2;                      % 2 pairs of Gabors shown
    simulation.nItems = 2;                          % 2 Gabors shown on each interval
    simulation.nOrientationDifferences = 20;        % Number of orientation differences to test

%     orientationsToTest = linspace(0,90,simulation.nOrientationDifferences);
    
    maxDifference = 7*sd;
    orientationsToTest = linspace(0,maxDifference,simulation.nOrientationDifferences);
    
    if maxDifference > 90
        
        maxDifference = 90;
        orientationsToTest = linspace(0,maxDifference,simulation.nOrientationDifferences);

    end

    % Set up data structures

    gabor.baseOrientation = NaN(simulation.nOrientationDifferences,simulation.nBlocks,simulation.nTrialsPerBlock,simulation.nIntervals);
    allOrientations = NaN(simulation.nOrientationDifferences,simulation.nBlocks,simulation.nTrialsPerBlock,simulation.nIntervals,simulation.nItems);
    allDifferent = NaN(simulation.nOrientationDifferences,simulation.nBlocks,simulation.nTrialsPerBlock);
    allEstimates = NaN(simulation.nOrientationDifferences,simulation.nBlocks,simulation.nTrialsPerBlock,simulation.nIntervals,simulation.nItems);
    allResponses = NaN(simulation.nOrientationDifferences,simulation.nBlocks,simulation.nTrialsPerBlock);

    % Simulation

    for thisDifference = 1:simulation.nOrientationDifferences

        orientationDifference = orientationsToTest(thisDifference);

        for thisBlock = 1:simulation.nBlocks

            allDifferent(thisDifference,thisBlock,:) = mod(randperm(simulation.nTrialsPerBlock),2)+1;    % Set up which interval will be different on each trial

            % Randomise base orientations

            for thisInterval = 1:simulation.nIntervals

            gabor.baseOrientation(thisDifference,thisBlock,:,thisInterval) = round(rand(1,simulation.nTrialsPerBlock)*180);

            end

            response = NaN(1,simulation.nTrialsPerBlock);

            for thisTrial = 1:simulation.nTrialsPerBlock

                whichIsDifferent = allDifferent(thisDifference,thisBlock,thisTrial);     % Determine which interval is different on this trial

                for thisInterval = 1:simulation.nIntervals

                    % Generate orientations of both items

                    for thisItem = 1:simulation.nItems-1

                        orientation(thisDifference,thisBlock,thisTrial,thisInterval,thisItem) = gabor.baseOrientation(thisDifference,thisBlock,thisTrial,thisInterval);
                        allOrientations(thisDifference,thisBlock,thisTrial,thisInterval,thisItem) = orientation(thisDifference,thisBlock,thisTrial,thisInterval,thisItem);

                    end

                    if whichIsDifferent==thisInterval

                        % This interval is "different". Generate a different second orientation

                        orientation(thisDifference,thisBlock,thisTrial,thisInterval,simulation.nItems) = orientation(thisDifference,thisBlock,thisTrial,thisInterval,thisItem)+orientationDifference;
                        allOrientations(thisDifference,thisBlock,thisTrial,thisInterval,simulation.nItems) = orientation(thisDifference,thisBlock,thisTrial,thisInterval,simulation.nItems);

                    else

                        % This interval is "same'. Generate the same orientation

                        orientation(thisDifference,thisBlock,thisTrial,thisInterval,simulation.nItems) = orientation(thisDifference,thisBlock,thisTrial,thisInterval,thisItem);
                        allOrientations(thisDifference,thisBlock,thisTrial,thisInterval,simulation.nItems) = orientation(thisDifference,thisBlock,thisTrial,thisInterval,simulation.nItems);

                    end

                    for thisItem = 1:simulation.nItems

                        % Generate pdf for both orientations

    %                    pdf(thisItem) = makedist('Normal', mu+orientation(thisDifference,thisBlock,thisTrial,thisInterval,thisItem), sd);

                        % Generate estimate from pdf

    %                     thisEstimate(thisItem) = random(pdf(thisItem));
    %                     allEstimates(thisDifference,thisBlock,thisTrial,thisInterval,thisItem) = thisEstimate(thisItem);

                        % Randomly estimate from a linearly transformed normal distribution

                        thisEstimate(thisItem) = (randn*sd)+(orientation(thisDifference,thisBlock,thisTrial,thisInterval,thisItem));
                        allEstimates(thisDifference,thisBlock,thisTrial,thisInterval,thisItem) = thisEstimate(thisItem);

                    end

                    % Find difference between these estimates

                    estimateDiff(thisTrial,thisInterval) = abs(thisEstimate(1)-thisEstimate(2));

                end

                % Generate response to this trial

                if estimateDiff(thisTrial,1) > estimateDiff(thisTrial,2) == 1

                    response(thisTrial) = 1;        % First interval had larger difference from estimates

                elseif estimateDiff(thisTrial,1) > estimateDiff(thisTrial,2) == 0

                    response(thisTrial) = 2;        % Second interval had larger difference from estimates

                end

                allResponses(thisDifference,thisBlock,thisTrial)=response(thisTrial);

            end

        end

    end

    nCorrect = NaN(1,simulation.nOrientationDifferences);
    pCorr = NaN(1,simulation.nOrientationDifferences);

    for thisDifference = 1:simulation.nOrientationDifferences

        nCorrect(thisDifference) = sum(allResponses(thisDifference,:)==allDifferent(thisDifference,:));
        pCorr(thisDifference) = nCorrect(thisDifference)/(simulation.nBlocks*simulation.nTrialsPerBlock);

    end

    % figureName = ['Standard deviation = ' num2str(sd)];
    % figure('Color','White','Name',figureName);
    % plot(orientationsToTest,pCorr,'k+');

    % Psychometric model fitting

    options = PAL_minimize('options');
    PF = @PAL_CumulativeNormal; % PAL_Weibull
    StimLevels = orientationsToTest;
    NumPos = nCorrect;
    OutOfNum = simulation.nTrialsPerDifference.*ones(size(orientationsToTest));
    searchGrid.alpha = [0:1:90];    %structure defining grid to
    searchGrid.beta = 10.^[-1:.01:2]; %search for initial values
    searchGrid.gamma = [.5];
    searchGrid.lambda = [0];
    paramsFree = [1 1 0 1];
    [paramsValues, LL, exitFlag, output] = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum, searchGrid, paramsFree, PF, 'searchOptions', options); 

    % Draw model fit plot
    StimLevelsFineGrain=[min(StimLevels):max(StimLevels)./1000:max(StimLevels)];
    ProportionCorrectModel = PF(paramsValues,StimLevelsFineGrain);

    figure('name','Maximum Likelihood Psychometric Function Fitting');
    plot(StimLevels,pCorr,'k+','markersize',20);
    set(gca, 'fontsize',10);
    set(gca, 'Xtick',StimLevels);
    axis([min(StimLevels) max(StimLevels) 0.4 1]);
    hold on;
    plot(StimLevelsFineGrain,ProportionCorrectModel,'g-','linewidth',4);

    % Save parameter values

    parameterValues(thisSD,:) = paramsValues;
    LLValues(thisSD) = LL;
    exitFlagValues(thisSD) = exitFlag;

end

for thisSD = 1:nSDValues
    
    if exitFlagValues(thisSD) == 1
        
    alphaValues(thisSD) = parameterValues(thisSD,1);
    betaValues(thisSD) = parameterValues(thisSD,2);
    lapseValues(thisSD) = parameterValues(thisSD,4);
    
    end
    
end

figure('name', 'Alpha Values for PF Fit');
plot(sdValues,alphaValues);
set(gca,'fontsize',10);
axis([0 maxSD+1.5 min(alphaValues)-.5 max(alphaValues)+.5]);
xlabel('Standard Deviation of Normal Distribution (degrees)');
ylabel('Alpha of Weibull PF fit');

figure('name', 'Beta Values for PF Fit');
plot(sdValues, betaValues);
set(gca,'fontsize',10);
axis([0 maxSD+1 min(betaValues) max(betaValues)]);
xlabel('Standard Deviation of Normal Distribution (degrees)');
ylabel('Beta of Weibull PF fit');