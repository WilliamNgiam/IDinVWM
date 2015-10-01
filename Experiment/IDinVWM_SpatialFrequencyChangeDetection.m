% Individual differences in perceptual ability and visual working memory

% This study aims to examine whether an individual's ability to
% discriminate between different orientation and spatial frequency of Gabor
% patches predicts the individual's visual working memory (VWM) ability on 
% a change-detection task.

% Run_IDinVWM.m should be run before this
% This code is for the spatial frequency discrimination task in Experiment 
% 1 of WN's PhD

% WN started writing this June 2015

% -------------------------------------------------------------------------

% Participant parameters received from Run_IDinVWM.m

% Set up experiment parameters

experiment.nCalibrationTrials = 10;
experiment.nPracticeTrials = 10;

experiment.maxDifference = 1;   % Maximum difference in spatial frequency (cycles per degree) between Gabor patches
experiment.minDifference = .1;  % Minimum difference in spatial frequency (cycles per degree) between Gabor patches

% Set up equipment parameters

equipment.viewDist = 500;                           % Viewing distance in mm
equipment.ppm = 2.7;                                % Pixels per mm (HP P1120, 1024 x 768, 120 Hz) % Remeasure
equipment.gammaVals = 1.0./[2.6434 2.2312 2.171];   % Gamma values for CRT in GT519 % Recalibrate

% Set up colour parameters

colour.blackVal = 0;
colour.greyVal = 0.5;
colour.whiteVal = 1;

colour.fixVal = 1;
colour.textVal = 0;

% Set up audio parameters

audio.toneLength = .1;
audio.toneFreq = 880;
audio.toneAmplitude = 1.0;
audio.sampleRate = 48000;        % Default audio sample rate
defaultTone = audio.toneAmplitude*sin(linspace(0,2*pi*audio.toneFreq*audio.toneLength,audio.sampleRate*audio.toneLength));

% Set up stimulus parameters

stimulus.size_dva = 2.5;            % Item size in degrees of visual angle
stimulus.eccentricity_dva = 4;      % Eccentricity in degrees of visual angle
stimulus.nArrayItems = 2;           % Number of items shown in array
stimulus.nItems = 2;                % Number of items in discrimination task
stimulus.maxDiff = 90;              % Maximum difference (degrees)

stimulus.fixationOn = 1;            % Fixation on (1) or off (0)
stimulus.fixationSize_dva = .25;    % Fixation size in degrees of visual angle

stimulus.orientation_deg = 90;      % Orientation of Gabor patches (in degrees from cardinal position)
stimulus.sc = 1/2.5;                % Spatial constant of Gabor patches (in degrees)
stimulus.contrast = 1;              % Contrast of Gabor patches
stimulus.aspectratio = 1;           % Aspect ratio of Gabor patches
stimulus.spatialFrequencyMax = 10;  % Max spatial frequency of Gabor patches shown (cycles per degree)
stimulus.spatialFrequencyMin = 1;   % Min spatial frequency to Gabor patches shown(cycles per degree)

% Set up temporal parameters

timing.fixationDuration = 0.5;      % Duration fixation point is displayed
timing.displayDuration = 2;         % Duration Gabors are displayed for each interval
timing.intervalDuration = 1;        % Duration between intervals
timing.blankDuration = 2;           % Duration of inter-trial interval (blank)

% Set up staircase parameters

staircase.pThreshold = .82; % Optimal for 2-AFC
staircase.tGuess = 0; %log10(PAL_Weibull(PFparams,staircase.pThreshold, 'Inverse')); % I think determine this from average of pilot
staircase.tGuessSD = 1.5; % Not sure if this is appropriate
staircase.beta = 2; % From simulation
staircase.delta = .01;%lapseValues(sd); % Lapse rate
staircase.gamma = .5; % Guess rate, 2-AFC 
staircase.grain = .01; % Step size
staircase.range = 4; % Range is tguess +(-range/2:grain:range/2)

% Initialise PsychSound

InitializePsychSound;
ptbAudioPort = PsychPortAudio('Open', [], [], 0, audio.sampleRate, 1);
PsychPortAudio('FillBuffer', ptbAudioPort, defaultTone);

% Set up PsychToolbox Pipeline

screenID = max(Screen('Screens'));
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'EnablePseudoGrayOutput');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
Screen('Preference','SkipSyncTests',2);
[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, colour.greyVal);
PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
[screenWidth, screenHeight] = RectSize(winRect);
screenCentreX = round(screenWidth/2);
screenCentreY = round(screenHeight/2);
flipInterval = Screen('GetFlipInterval', ptbWindow);

% Calculate equipment parameters

equipment.mpd = (equipment.viewDist)*tan(deg2rad(2*stimulus.eccentricity_dva))/stimulus.eccentricity_dva; % Calculate mm per degree of visual angle to the ecccentricity of the stimuli
equipment.ppd = equipment.ppm*equipment.mpd;

% Calculate spatial parameters

stimulus.size_pix = round(stimulus.size_dva*equipment.ppd);                     % Item size in pixels
stimulus.eccentricity_pix = round(stimulus.eccentricity_dva*equipment.ppd);     % Eccentricity of stimulus in pixels
stimulus.fixationSize_pix = stimulus.fixationSize_dva*equipment.ppd;           % Fixation cross size in pixels

stimulus.sc_pix = round(stimulus.sc*equipment.ppd);                            % Spatial constant of Gabor patches in pixels

% Instruction Text

taskExplanationText = ['For this block, you will be shown a pair of stimuli at two different intervals.\n' ...
    'In one pair, the spatial frequency of the stimuli will be different to each other.\n' ...
    'While the other pair will have the same spatial frequency.\n' ...
    'You need to respond with which pair you think was different.\n' ...
    'Press any key for the next set of instructions'];
    
responseInstructionText = ['If you think the first pair was different, respond with a left arrow keypress.\n' ...
    'If you think the second pair was different, respond with a right arrow keypress.\n' ...
    'Press a response key to continue.'];

startCalibrationText = ['You will begin with some calibration trials.\n' ...
    'Press a response key to start the calibration trials.'];

startStaircaseText = ['You will now begin block ' num2str(thisStaircase) ' of ' num2str(experiment.nStaircases) ' blocks.\n' ...
    'Press a response key to start the block.'];

startBlockText = ['Press any key to start the block'];

% ----------------------------------------------------------------------- %

%                           Begin Experiment

% ----------------------------------------------------------------------- %

%                     Spatial Frequency Discrimination                    %

% Practice trials?

allDifferences = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase); % For saving the orientation difference used on each trial

% Create staircases

for thisStaircase = 1:experiment.nStaircases
    
    allData(thisStaircase) = QuestCreate(staircase.tGuess,staircase.tGuessSD,staircase.pThreshold,staircase.beta,staircase.delta,staircase.gamma,staircase.grain,staircase.range);
    
end

% Estimate staircase parameters by testing across range of differences

% Set up "calibration" parameters

differencesToTest = linspace(experiment.maxDifference, experiment.minDifference, experiment.nCalibrationTrials); % In cycles per degree
differencesToTest_cpp = differencesToTest/equipment.ppd; % In cycles per pixels

calibration.allChanges = mod(randperm(experiment.nCalibrationTrials),stimulus.nArrayItmes)+1;   % Determines which Gabor will change in the array
calibration.allChangeDirections = mod(randperm(experiment.nCalibrationTrials),2)+1;             % Determines which direction the change will occur (increase spatial frequency or decrease spatial frequency)

% Need to also ensure that there is more than 1 cycle within the aperture

calibration.allSpatialFrequencies_cpd = (rand(stimulus.nArrayItems,experiment.nCalibrationTrials)*9+1); % Spatial frequency in cycles per degree (Using the Robson (1966) paper CSF to determine appropriate range - About 1 cpd to 10 cpd)
calibration.allSpatialFrequencies_cpp = calibration.allSpatialFrequencies/equipment.ppd; % First spatial frequency in cycles per pixel

% Randomise phase on each trial to prevent it cueing the spatial frequency

calibration.FirstPhase = rand(1,experiment.nCalibrationTrials)*180; % Randomise phase for when the Gabor is drawn (degrees)
calibration.SecondPhase = rand(1,experiment.nCalibrationTrials)*180;

% Set up location rects of Gabors for every trial

calibration.itemBaseTheta_rad = repmat(linspace(0,2*pi-(2*pi/stimulus.nArrayItems), stimulus.nArrayItems)',1,experiment.nCalibrationTrials);
calibration.arrayJitter_rad = ((2*pi)/stimulus.nArrayItems)*repmat(rand(1,experiment.nCalibrationTrials), stimulus.nArrayItems, 1);
calibration.allItemTheta_rad = calibration.itemBaseTheta_rad+calibration.arrayJitter_rad;

% Instruction Text

DrawFormattedText(ptbWindow,taskExplanationText, 'center', 'center', colour.textVal);
Screen('Flip',ptbWindow);
waitResponse = 1;

while waitResponse

    KbWait(-1,2);
    waitResponse = 0;

end

% Response Instruction Text

DrawFormattedText(ptbWindow,responseInstructionText, 'center', 'center', colour.textVal);
Screen('Flip',ptbWindow);
waitResponse = 1;

while waitResponse

    [beginTime, keyCode] = KbWait(-1,2);
    pressedKey = find(keyCode);

    if (pressedKey==79)||(pressedKey==80) % 79 for right arrow key and 80 for left arrow key
        waitResponse = 0;
    end

end

% Wait for keypress response to start calibration

DrawFormattedText(ptbWindow,startCalibrationText, 'center', 'center', colour.textVal);
Screen('Flip',ptbWindow);
waitResponse = 1;

while waitResponse

    [startCalibrationTime, keyCode] = KbWait(-1,2);
    pressedKey = find(keyCode);

    if (pressedKey==79)||(pressedKey==80) % 79 for right arrow key and 80 for left arrow key
        waitResponse = 0;
    end

end
    
responseTime = startCalibrationTime;

% Prevent Vernier discrimination task (orientation is always 90 degrees so prevent
% Gabors being located at 0 degrees and 180 degrees location

for thisCalibrationTrial = 1:experiment.nCalibrationTrials

    stimulus.orientation_rad = deg2rad(stimulus.orientation_deg);     % Minimum difference between orientation and position in radians
    minimumDifference_rad = pi/4;
    
    if calibration.allItemTheta_rad(1,thisCalibrationTrial) < minimumDifference_rad % Between 0 and 45 degrees

        calibration.allItemTheta_rad(:,thisCalibrationTrial) = calibration.allItemTheta_rad(:,thisCalibrationTrial)+minimumDifference_rad; 

    elseif calibration.allItemTheta_rad(1,thisCalibrationTrial) > (2*pi-minimumDifference_rad) % Between 315 and 360
        
        calibration.allItemTheta_rad(:,thisCalibrationTrial) = calibration.allItemTheta_rad(:,thisCalibrationTrial)-minimumDifference_rad; 
        
    elseif calibration.allItemTheta_rad(1,thisCalibrationTrial) > (pi-minimumDifference_rad) && calibration.allItemTheta_rad(1,thisCalibrationTrial) < (pi+minimumDifference_rad) % Between 135 and 225
        
        diceRoll = randi(2);
        
        if diceRoll == 1
            
            calibration.allItemTheta_rad(:,thisCalibrationTrial) = calibration.allItemTheta_rad(:,thisCalibrationTrial)+(minimumDifference_rad*2);
        
        elseif diceRoll == 2
        
            calibration.allItemTheta_rad(:,thisCalibrationTrial) = calibration.allItemTheta_rad(:,thisCalibrationTrial)-(minimumDifference_rad*2);

        end

    end
    
end

% Create Gabor textures

[gaborid, gaborrect] = CreateProceduralGabor(ptbWindow,stimulus.size_pix,stimulus.size_pix, 0, [0.5 0.5 0.5 0], 1, .5);
       
allCalibrationResponses = NaN(1,experiment.nCalibrationTrials);
allCalibrationCorrect = NaN(1,experiment.nCalibrationTrials);
allCalibrationChange = NaN(1,experiment.nCalibrationTrials);

for thisCalibrationTrial = 1:experiment.nCalibrationTrials
    
    % Create rects
    
    fixRect = [0 0 stimulus.fixationSize_pix stimulus.fixationSize_pix];
    fixRect = CenterRectOnPoint(fixRect, screenCentreX, screenCentreY);

    itemRect = [0 0 stimulus.size_pix stimulus.size_pix];

    theseItemTheta = calibration.allItemTheta_rad(:,thisCalibrationTrial);
    [theseItemX, theseItemY] = pol2cart(theseItemTheta,stimulus.eccentricity_pix*ones(stimulus.nArrayItems,1));

    itemRects = NaN(4,stimulus.nArrayItems);

    for thisItem = 1:stimulus.nArrayItems
    
        itemRects(:,thisItem) = CenterRectOnPoint(itemRect,theseItemX(thisItem)+screenCentreX,theseItemY(thisItem)+screenCentreY)';
        
    end

    % Determine what the size of the difference in spatial frequency will be
    
    thisSpatialFrequencyDifference = differencesToTest_cpp(thisCalibrationTrial);
    thisLogSpatialFrequencyDifference = log10(thisSpatialFrequencyDifference);
    
    % Display fixation cross

    if stimulus.fixationOn

        Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

    end
    
    % Display fixation with tone
    
    startTime = Screen('Flip', ptbWindow, responseTime);
    
    % Play tone to signify start of trial

     PsychPortAudio('Start', ptbAudioPort, 1, startTime, 0);
            
    % Draw first set of Gabors
    
    if stimulus.fixationOn

        Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

    end
        
    for thisGabor = 1:stimulus.nItems

        Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), stimulus.orientation_deg, [], [], [], [], kPsychDontDoRotation, [calibration.FirstPhase(thisCalibrationTrial), calibration.FirstSpatialFrequency_cpp(thisCalibrationTrial), stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

    end
        
    % Display first intervals of Gabor
    
    firstIntervalTime = Screen('Flip', ptbWindow, startTime + timing.blankDuration);
    
    % Display blank (with fixation) for interval

    if stimulus.fixationOn

        Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

    end
    
    blankTime = Screen('Flip', ptbWindow, firstIntervalTime + timing.displayDuration);
    
    % Draw second set of Gabors
    
    if stimulus.fixationOn

        Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

    end

    for thisGabor = 1:stimulus.nArrayItems

        if thisGabor ~= calibration.allChanges(thisCalibrationTrial)
        
        Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), stimulus.orientation_deg, [], [], [], [], kPsychDontDoRotation, [calibration.SecondPhase(thisCalibrationTrial), calibration.SecondSpatialFrequency_cpp(thisCalibrationTrial), stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0])

        elseif thisGabor == calibration.allChanges(thisCalibrationTrial)
               
            if calibration.allChangeDirections(thisCalibrationTrial) == 1

                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), stimulus.orientation_deg, [], [], [], [], kPsychDontDoRotation, [calibration.SecondPhase(thisCalibrationTrial), calibration.SecondSpatialFrequency_cpp(thisCalibrationTrial)+thisSpatialFrequencyDifference, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

            elseif calibration.allChangeDirections(thisCalibrationTrial) == 2

                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), stimulus.orientation_deg, [], [], [], [], kPsychDontDoRotation, [calibration.SecondPhase(thisCalibrationTrial), calibration.SecondSpatialFrequency_cpp(thisCalibrationTrial)-thisSpatialFrequencyDifference, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

            end
            
        end
        
    end

    secondIntervalTime = Screen('Flip', ptbWindow, blankTime + timing.intervalDuration);
    
    % Flip to blank
    
    finishTime = Screen('Flip', ptbWindow, secondIntervalTime + timing.displayDuration);
   
    % Get click response, participant indicates which Gabor patch they
    % think has changed
    
    % Get response and store

    waitResponse = 1;
    ShowCursor;
    SetMouse(screenCentreX,screenCentreY,ptbWindow); % Move mouse cursor to centre of screen
    
    CheckResponse = zeros(1,stimulus.nArrayItems);

    while ~any(CheckResponse)

        [~,xClickResponse,yClickResponse] = GetClicks(ptbWindow, 0); % Get x and y location of first mouse click on test array
        clickSecs = GetSecs;

        for thisGabor = 1:stimulus.nArrayItems;

            CheckResponse(thisGabor) = IsInRect(xClickResponse, yClickResponse, itemRects(:,thisGabor)); % Test if mouse click is in aperture of each successive item
                    % CheckResponse returns 0 for no, 1 for yes

        end    

    end

    ClickedTarget = find(CheckResponse);
    calibration.allResponses(thisCalibrationTrial) = ClickedTarget;
    
    HideCursor;

    if ClickedTarget == allChanges(thisCalibrationTrial)
        
        ifCorrect = 1;
        
    elseif ClickedTarget ~= allChanges(thisCalibrationTrial)
        
        ifCorrect = 0;
        
    end
    
    allCalibrationCorrect(thisCalibrationTrial) = ifCorrect;
    
    % Update all staircases
    
    for thisStaircase = 1:experiment.nStaircases

        allData(thisStaircase) = QuestUpdate(allData(thisStaircase), thisLogOrientationDifference, ifCorrect);

    end

end

allFirstOrientations = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);
allSecondOrientations = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);
allDifferences = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);
allResponses = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);
allCorrect = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);

for thisStaircase = 1:experiment.nStaircases

    % Set up staircase parameters
    
    whichEstimateInterval = mod(randperm(experiment.nTrialsPerStaircase),2)+1;          % Determines which pair will have different Gabors
    staircase.allChangeDirections = mod(randperm(experiment.nTrialsPerStaircase,2)+1);   % Determines whether the change will be an increase or decrease in spatial frequency

    staircase.FirstSpatialFrequency_cpd = round(rand(1,experiment.nTrialsPerStaircase)*20)*equipment.ppd;
    staircase.SecondSpatialFrequency_cpd = round(rand(1,experiment.nTrialsPerStaircase)*20)*equipment.ppd;

    % Randomise phase on each trial to prevent it cueing the spatial
    % frequency
    
    staircase.FirstPhase = round(rand(1,experiment.nTrialsPerStaircase)*180); % Phase in degrees
    staircase.SecondPhase = round(rand(1,experiment.nTrialsPerStaircase)*180);

    % Set up location rects of Gabors for every trial

    staircase.itemBaseTheta_rad = repmat(linspace(0,2*pi-(2*pi/stimulus.nArrayItems), stimulus.nArrayItems)',1,experiment.nTrialsPerStaircase);
    staircase.arrayJitter_rad = ((2*pi)/stimulus.nArrayItems)*repmat(rand(1,experiment.nTrialsPerStaircase), stimulus.nArrayItems, 1);
    staircase.allItemTheta_rad = staircase.itemBaseTheta_rad+staircase.arrayJitter_rad;
    
    % Prevent Vernier discrimination task
    
    stimulus.orientation_rad = deg2rad(stimulus.orientation_deg);     % Minimum difference between orientation and position in radians
    minimumDifference_rad = pi/4;
    
    if staircase.allItemTheta_rad(thisStaircaseTrial) < minimumDifference_rad % Between 0 and 45 degrees

        staircase.allItemTheta_rad(thisstaircaseTrial) = staircase.allItemTheta_rad(thisStaircaseTrial)+minimumDifference_rad; 

    elseif staircase.allItemTheta_rad(thisStaircaseTrial) > (2*pi-minimumDifference_rad) % Between 315 and 360
        
        staircase.allItemTheta_rad(thisStaircaseTrial) = staircase.allItemTheta_rad(thisStaircaseTrial)-minimumDifference_rad; 
        
    elseif staircase.allItemTheta_rad(thisStaircaseTrial) > (pi-minimumDifference_rad) && staircase.allItemTheta_rad(thisStaircaseTrial) < (pi+minimumDifference_rad) % Between 135 and 225
        
        diceRoll = randi(2);
        
        if diceRoll == 1
            
            staircase.allItemTheta_rad(thisStaircaseTrial) = staircase.allItemTheta_rad+(minimumDifference_rad*2);
        
        elseif diceRoll == 2
        
            staircase.allItemTheta_rad(thisStaircaseTrial) = staircase.allitemTheta_rad-(minimumDifference_rad*2);

        end

    end
    
    % Create Gabor textures

    [gaborid, gaborrect] = CreateProceduralGabor(ptbWindow,stimulus.size_pix,stimulus.size_pix, 0, [0.5 0.5 0.5 0], 1, .5);

    % Wait for keypress response to start staircase

    DrawFormattedText(ptbWindow,startStaircaseText, 'center', 'center', colour.textVal);
    Screen('Flip',ptbWindow);
    waitResponse = 1;
    
    while waitResponse

    [startTime, keyCode] = KbWait(-1,2);
    pressedKey = find(keyCode);

        if (pressedKey==79)||(pressedKey==80) % 79 for right arrow key and 80 for left arrow key
            waitResponse = 0;
        end

    end
    
    responseTime = startTime;
    
    % Set up trial parameters

    whichInterval = mod(randperm(experiment.nTrialsPerStaircase),2)+1;    % Determines which pair will have different Gabors

    for thisStaircaseTrial = 1:experiment.nTrialsPerStaircase

        % Create rects

        fixRect = [0 0 stimulus.fixationSize_pix stimulus.fixationSize_pix];
        fixRect = CenterRectOnPoint(fixRect, screenCentreX, screenCentreY);

        itemRect = [0 0 stimulus.size_pix stimulus.size_pix];

        theseItemTheta = staircase.allItemTheta_rad(:,thisStaircaseTrial);
        [theseItemX, theseItemY] = pol2cart(theseItemTheta,stimulus.eccentricity_pix*ones(stimulus.nArrayItems,1));

        itemRects = NaN(4,stimulus.nArrayItems);

        for thisItem = 1:stimulus.nArrayItems

            itemRects(:,thisItem) = CenterRectOnPoint(itemRect,theseItemX(thisItem)+screenCentreX,theseItemY(thisItem)+screenCentreY)';

        end

        % Determine which interval will have the different pair

        whichOne = whichEstimateInterval(thisStaircaseTrial);

        % Determine what the size of the difference in orientation will be

        thisLogSpatialFrequencyDifference = QuestMean(allData(thisStaircase));
        thisSpatialFrequencyDifference = 10.^thisLogOrientationDifference;
        allDifferences(thisStaircase,thisStaircaseTrial) = thisSpatialFrequencyDifference;

        % Display fixation cross

        if stimulus.fixationOn

            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

        end

        % Display fixation with tone

        startTrialTime = Screen('Flip', ptbWindow, responseTime);

        % Play tone to signify start of trial

         PsychPortAudio('Start', ptbAudioPort, 1, startTrialTime, 0);

        % Draw first set of Gabors

        if stimulus.fixationOn

            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

        end

        if whichOne == 1 % First interval is different

            for thisGabor = 1:stimulus.nItems-1

                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), stimulus.orientation_deg, [], [], [], [], kPsychDontDoRotation, [180, staircase.FirstSpatialFrequency_cpd, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

            end

            if staircase.allChangeDirections(thisStaircaseTrial) == 1
                
                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), stimulus.orientation_deg, [], [], [], [], kPsychDontDoRotation, [180, staircase.FirstSpatialFrequency_cpd + thisSpatialFrequencyDifference, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

            elseif staircase.allChangeDirections(thisStaircaseTrial) == 2
                
                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), stimulus.orientation_deg, [], [], [], [], kPsychDontDoRotation, [180, staircase.FirstSpatialFrequency_cpd - thisSpatialFrequencyDifference, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);
                
            end
            
        else

            for thisGabor = 1:stimulus.nItems

                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), stimulus.orientation_deg, [], [], [], [], kPsychDontDoRotation, [180, staircase.FirstSpatialFrequency_cpd, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

            end

        end

        % Display first intervals of Gabor

        firstIntervalTime = Screen('Flip', ptbWindow, startTime + timing.blankDuration);

        % Display blank (with fixation) for interval

        if stimulus.fixationOn

            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

        end

        blankTime = Screen('Flip', ptbWindow, firstIntervalTime + timing.displayDuration);

        % Draw second set of Gabors

        if stimulus.fixationOn

            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

        end

        if whichOne == 2 % Second interval is different

            for thisGabor = 1:stimulus.nItems-1

                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), stimulus.orientation_deg, [], [], [], [], kPsychDontDoRotation, [180, staircase.SecondSpatialFrequency_cpd, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0])
            end
            
            if staircase.allChangeDirections == 1
            
                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), stimulus.orientation_deg, [], [], [], [], kPsychDontDoRotation, [180, staircase.SecondSpatialFrequency_cpd + thisSpatialFrequencyDifference, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

            elseif staircase.allChangeDirections == 2
                
                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), stimulus.orientation_deg, [], [], [], [], kPsychDontDoRotation, [180, staircase.SecondSpatialFrequency_cpd - thisSpatialFrequencyDifference, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);
        else

            for thisGabor = 1:stimulus.nItems

                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), staircase.SecondOrientation_deg(thisStaircaseTrial)+90, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0])
            
            end

        end

        secondIntervalTime = Screen('Flip', ptbWindow, blankTime + timing.intervalDuration);

        % Flip to blank

        finishTime = Screen('Flip', ptbWindow, secondIntervalTime + timing.displayDuration);

        % Get keypress response 'first' or 'second'

        % Get response and store

        waitResponse = 1;

        while waitResponse

            [keySecs, keyCode] = KbWait(-1,2);
            pressedKey = find(keyCode);

            if (pressedKey==79)||(pressedKey==80)
                waitResponse = 0;
            end

        end

        allResponses(thisStaircase,thisStaircaseTrial) = pressedKey;

        % Determine if response is correct

        if pressedKey + whichOne == 81; % Correct repsonses will always add up to this number

            ifCorrect = 1;

        else

            ifCorrect = 0;

        end

        allCorrect(thisStaircase,thisStaircaseTrial) = ifCorrect;

        % Update staircase

        allData(thisStaircase) = QuestUpdate(allData(thisStaircase), thisLogOrientationDifference, ifCorrect);

    end

    % Save staircase parameters
    
    allFirstOrientations(thisStaircase,:) = staircase.FirstOrientation_deg;
    allSecondOrientations(thisStaircase,:) = staircase.SecondOrientation_deg;
    
    % Save data file
    
end

