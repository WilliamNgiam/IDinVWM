% Individual differences in perceptual ability and visual working memory

% This study aims to examine whether an individual's ability to
% discriminate between different orientation and spatial frequency of Gabor
% patches predicts the individual's visual working memory (VWM) ability on 
% a change-detection task.

% Run_IDinVWM.m should be run before this
% This code is for the orientation discrimination task in Experiment 1 of 
% WN's PhD

% WN started writing this June 2015

% -------------------------------------------------------------------------

% Participant parameters received from Run_IDinVWM.m

% Create mask textures prior to experiment
% 
% MaskForGaborStimuli;

% Set up experiment parameters

experiment.nCalibrationTrials = 20;
experiment.nPracticeTrials = 10;
experiment.thisStaircase = 1;

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
audio.sampleRate = 48000;           % Default audio sample rate
defaultTone = audio.toneAmplitude*sin(linspace(0,2*pi*audio.toneFreq*audio.toneLength,audio.sampleRate*audio.toneLength));

% Set up stimulus parameters

stimulus.size_dva = 4;            % Item size in degrees of visual angle
stimulus.eccentricity_dva = 4;      % Eccentricity in degrees of visual angle
stimulus.nArrayItems = 2;           % Number of items shown in array
stimulus.nItems = 2;                % Number of items in discrimination task
stimulus.maxDiff = 90;              % Maximum difference (degrees)

stimulus.fixationOn = 1;            % Fixation on (1) or off (0)
stimulus.fixationSize_dva = .25;    % Fixation size in degrees of visual angle

stimulus.frequency_cpd = 5;         % Spatial frequency of Gabor patches (in cycles per degree)
stimulus.sc = 1/2.5;                % Spatial constant of Gabor patches (in degrees)
stimulus.contrast = 1;              % Contrast of Gabor patches
stimulus.aspectratio = 1;           % Aspect ratio of Gabor patches

% Set up temporal parameters

timing.fixationDuration = 0.5;                      % Duration fixation point is displayed
timing.displayDuration = 1;                         % Duration Gabors are displayed for each interval
timing.maskDuration = .5;                           % Duration of mask following Gabor display
timing.intervalDuration = 1.5;                        % Duration between intervals
timing.blankDuration = 2;                           % Duration of inter-trial interval (blank)

% Calculate temporal parameters
%nMaskFrames = round(timing.maskArrayDuration/timing.maskFrameDuration);

% Set up staircase parameters

staircase.pThreshold = .82;     % Optimal for 2-AFC
staircase.tGuess = 0;           % log10(PAL_Weibull(PFparams,staircase.pThreshold, 'Inverse')); % I think determine this from average of pilot
staircase.tGuessSD = 2;         % Not sure if this is appropriate
staircase.beta = 2;             % From simulation
staircase.delta = .01;          % lapseValues(sd); % Lapse rate
staircase.gamma = .5;           % Guess rate, 2-AFC 
staircase.grain = .01;          % Step size
staircase.range = 4;            % Range is tguess +(-range/2:grain:range/2)

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

% % Enable alpha blending for typical drawing of masked textures
% Screen('BlendFunction', ptbWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Calculate equipment parameters

equipment.mpd = (equipment.viewDist)*tan(deg2rad(2*stimulus.eccentricity_dva))/stimulus.eccentricity_dva; % Calculate mm per degree of visual angle to the ecccentricity of the stimuli
equipment.ppd = equipment.ppm*equipment.mpd;

% Calculate spatial parameters

stimulus.size_pix = round(stimulus.size_dva*equipment.ppd);                     % Item size in pixels
stimulus.eccentricity_pix = round(stimulus.eccentricity_dva*equipment.ppd);     % Eccentricity of stimulus in pixels
stimulus.fixationSize_pix = stimulus.fixationSize_dva*equipment.ppd;            % Fixation cross size in pixels

stimulus.frequency_pix = stimulus.frequency_cpd/equipment.ppd;                  % Spatial frequency of Gabor patches in cycles per pixel
stimulus.sc_pix = round(stimulus.sc*equipment.ppd);                             % Spatial constant of Gabor patches in pixels

% Instruction Text

taskExplanationText = ['For this block, you will be shown a pair of stimuli at two different intervals.\n' ...
    'In one pair, the orientation of the stimuli will be different to each other.\n' ...
    'While the other pair will have the same orientation.\n' ...
    'You need to respond with which pair you think was different.\n' ...
    'Press any key for the next set of instructions'];
    
responseInstructionText = ['If you think the first pair was different, respond with a left arrow keypress.\n' ...
    'If you think the second pair was different, respond with a right arrow keypress.\n' ...
    'Press a response key to continue.'];

startCalibrationText = ['You will begin with some calibration trials.\n' ...
    'Press a response key to start the calibration trials.'];

startBlockText = ['Press any key to start the block'];

finishText = ['You have completed the experiment.\n' ...
     'Please inform the experimenter.'];
% ----------------------------------------------------------------------- %

%                           Begin Experiment

% ----------------------------------------------------------------------- %

%                       Orientation Discrimination                        %

% Practice trials?

allDifferences = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase); % For saving the orientation difference used on each trial

% Create staircases

for thisStaircase = 1:experiment.nStaircases
    
    allData(thisStaircase) = QuestCreate(staircase.tGuess,staircase.tGuessSD,staircase.pThreshold,staircase.beta,staircase.delta,staircase.gamma,staircase.grain,staircase.range);
    
end

% Estimate staircase parameters by testing across range of differences

% Set up "calibration" parameters

differencesToTest_deg = linspace(90,0.1,experiment.nCalibrationTrials);
differencesToTest_rad = deg2rad(differencesToTest_deg);

whichEstimateInterval = mod(randperm(experiment.nCalibrationTrials),2)+1;  % Determines which pair will have different Gabors
whichChangeDirection = mod(randperm(experiment.nCalibrationTrials),2)+1;   % Determines which direction the orientation will change

calibration.FirstOrientation_deg = round(rand(1,experiment.nCalibrationTrials)*180);
calibration.SecondOrientation_deg = round(rand(1,experiment.nCalibrationTrials)*180);

calibration.FirstOrientation_rad = deg2rad(calibration.FirstOrientation_deg);
calibration.SecondOrientation_rad = deg2rad(calibration.SecondOrientation_deg);

% Set up location rects of Gabors for every trial

calibration.itemBaseTheta_rad = repmat(linspace(0,2*pi-(2*pi/stimulus.nArrayItems), stimulus.nArrayItems)',1,experiment.nCalibrationTrials);
calibration.arrayJitter_rad = ((2*pi)/stimulus.nArrayItems)*repmat(rand(1,experiment.nCalibrationTrials), stimulus.nArrayItems, 1);
calibration.allItemTheta_rad = calibration.itemBaseTheta_rad+calibration.arrayJitter_rad;
calibration.allItemTheta_deg = rad2deg(calibration.allItemTheta_rad);

% Instruction Text

HideCursor;

DrawFormattedText(ptbWindow,taskExplanationText, 'center', 'center', colour.textVal);
Screen('Flip',ptbWindow);
waitResponse = 1;

while waitResponse

    [time, keyCode] = KbWait(-1,2);
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

% Prevent Vernier discrimination task

for thisCalibrationTrial = 1:experiment.nCalibrationTrials

    minimumDifference_deg = 45;                                 % Minimum difference between orientation and position in degrees
    minimumDifference_rad = deg2rad(minimumDifference_deg);     % Minimum difference between orientation and position in radians
    
    if calibration.FirstOrientation_rad(thisCalibrationTrial) < pi/2
    
        if deg2rad(90)-minimumDifference_rad < calibration.FirstOrientation_rad(thisCalibrationTrial) + calibration.allItemTheta_rad(1,thisCalibrationTrial) ...
                && calibration.FirstOrientation_rad(thisCalibrationTrial) + calibration.allItemTheta_rad(1,thisCalibrationTrial) < deg2rad(90)+minimumDifference_rad

            calibration.FirstOrientation_rad(thisCalibrationTrial) = calibration.FirstOrientation_rad(thisCalibrationTrial) + (minimumDifference_rad*2);
            calibration.FirstOrientation_deg(thisCalibrationTrial) = rad2deg(calibration.FirstOrientation_rad(thisCalibrationTrial));

        end
        
    elseif calibration.FirstOrientation_rad(thisCalibrationTrial) >= pi/2
        
        if deg2rad(270)-minimumDifference_rad < calibration.FirstOrientation_rad(thisCalibrationTrial) + calibration.allItemTheta_rad(1,thisCalibrationTrial) ...
                && calibration.FirstOrientation_rad(thisCalibrationTrial) + calibration.allItemTheta_rad(1,thisCalibrationTrial) < deg2rad(270)+minimumDifference_rad
        
            calibration.FirstOrientation_rad(thisCalibrationTrial) = calibration.FirstOrientation_rad(thisCalibrationTrial) + (minimumDifference_rad*2);
            calibration.FirstOrientation_deg(thisCalibrationTrial) = rad2deg(calibration.FirstOrientation_rad(thisCalibrationTrial));
        
        end 
    
    end
    
    if calibration.SecondOrientation_rad(thisCalibrationTrial) < pi/2
    
        if deg2rad(90)-minimumDifference_rad < calibration.SecondOrientation_rad(thisCalibrationTrial) + calibration.allItemTheta_rad(2,thisCalibrationTrial) ...
                && calibration.SecondOrientation_rad(thisCalibrationTrial) + calibration.allItemTheta_rad(2,thisCalibrationTrial) < deg2rad(90)+minimumDifference_rad

            calibration.SecondOrientation_rad(thisCalibrationTrial) = calibration.SecondOrientation_rad(thisCalibrationTrial) + (minimumDifference_rad*2);
            calibration.SecondOrientation_deg(thisCalibrationTrial) = rad2deg(calibration.SecondOrientation_rad(thisCalibrationTrial));

        end
        
    elseif calibration.SecondOrientation_rad(thisCalibrationTrial) >= pi/2
        
        if deg2rad(270)-minimumDifference_rad < calibration.SecondOrientation_rad(thisCalibrationTrial) + calibration.allItemTheta_rad(2,thisCalibrationTrial) ...
                && calibration.SecondOrientation_rad(thisCalibrationTrial) + calibration.allItemTheta_rad(2,thisCalibrationTrial) < deg2rad(270)+minimumDifference_rad
        
            calibration.SecondOrientation_rad(thisCalibrationTrial) = calibration.SecondOrientation_rad(thisCalibrationTrial) + (minimumDifference_rad*2);
            calibration.SecondOrientation_deg(thisCalibrationTrial) = rad2deg(calibration.SecondOrientation_rad(thisCalibrationTrial));
        
        end 
    
    end

end

% Create Gabor textures

[gaborid, gaborrect] = CreateProceduralGabor(ptbWindow,stimulus.size_pix,stimulus.size_pix, 0, [0.5 0.5 0.5 0], 1, 0.5);
       
calibration.allCalibrationResponses = NaN(1,experiment.nCalibrationTrials);
calibration.allCalibrationCorrect = NaN(1,experiment.nCalibrationTrials);
calibration.allCalibrationChange = whichEstimateInterval;

 
for thisCalibrationTrial = 1:experiment.nCalibrationTrials
    
    % On each trial, create a new noise array to change mask texture
    
    % Wipe previous mask textures to clear memory
    
    if thisCalibrationTrial ~= 1
        
        Screen('Close',maskTexture);
        
    end
    
    % Create a noise image the same size as your Gabor texture, randomising each
    % pixel's luminance value (uniform over -1,1).

    a = 0;
    b = 1;
    noiseArray = a + (b-a).*rand(size(imageArray));

    % Take the Fourier transform of the noise image and set aside the phase
    % spectrum.

    fftnoiseImage = fft2(noiseArray,SFT_Size,SFT_Size);
    frnoiseImage = abs(fftnoiseImage);
    phnoiseImage = angle(fftnoiseImage);
    
    % Take the inverse transform using the amplitude spectrum and phase
    % spectrum.
    maskednoise = ifft2(maskwhatever .* (phnoiseImage*exp(1i)), SFT_Size, SFT_Size);

    % Crop
    middleRow = floor(SFT_Size/2);
    nRows = round(stimulus.size_pix/2);
    startRow = middleRow - nRows + 1;
    endRow = middleRow + nRows;

    maskednoise = real(maskednoise(startRow:endRow,startRow:endRow));

    % Normalise
    maskednoise = (maskednoise-min(maskednoise(:))) / (max(maskednoise(:)) - min(maskednoise(:)));
    maskednoise = (maskednoise*2)-1; % Balance
    mask = blobArray.*maskednoise;
    mask = (mask-min(mask(:)))/(max(mask(:))-min(mask(:)));
    rgbImage = repmat(uint8(255.*mask),[1 1 3]);
    maskTexture = Screen('MakeTexture',ptbWindow,rgbImage);
    
    % Create rects
    
    fixRect = [0 0 stimulus.fixationSize_pix stimulus.fixationSize_pix];
    fixRect = CenterRectOnPoint(fixRect, screenCentreX, screenCentreY);

    itemRect = [0 0 stimulus.size_pix stimulus.size_pix];

    theseItemTheta = calibration.allItemTheta_rad(:,thisCalibrationTrial);
    [theseItemX, theseItemY] = pol2cart(theseItemTheta,stimulus.eccentricity_pix*ones(stimulus.nArrayItems,1));

    itemRects = NaN(4,stimulus.nArrayItems);
    maskRects = NaN(4,stimulus.nArrayItems);
    
    for thisItem = 1:stimulus.nArrayItems
    
        itemRects(:,thisItem) = CenterRectOnPoint(itemRect,theseItemX(thisItem)+screenCentreX,theseItemY(thisItem)+screenCentreY)';
        
    end
    
    maskRects = itemRects;
    
    % Determine which interval will have the different pair

    whichOne = whichEstimateInterval(thisCalibrationTrial);

    % Determine what the size of the difference in orientation will be
    
    thisOrientationDifference = differencesToTest_deg(thisCalibrationTrial);
    thisLogOrientationDifference = log10(thisOrientationDifference);
    
    % Display fixation cross

    if stimulus.fixationOn

        Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

    end
    
    % Display fixation with tone
    
    startTime = Screen('Flip', ptbWindow, responseTime);
    
    % Play tone to signify start of trial

     PsychPortAudio('Start', ptbAudioPort, 1, startTime, 0);
     
     WaitSecs(.5);
            
    % Draw first set of Gabors
    
    if stimulus.fixationOn

        Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

    end
    
    if whichOne == 1 % First interval is different

        for thisGabor = 1:stimulus.nArrayItems-1

            Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), calibration.FirstOrientation_deg(thisCalibrationTrial), [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);
                
        end
        
        if whichChangeDirection(thisCalibrationTrial) == 1 % Determines which direction the orientation change occurs
            
            Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), calibration.FirstOrientation_deg(thisCalibrationTrial)+thisOrientationDifference, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);
         
        elseif whichChangeDirection(thisCalibrationTrial) == 2
            
            Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), calibration.FirstOrientation_deg(thisCalibrationTrial)-thisOrientationDifference, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

        end

    else
        
        for thisGabor = 1:stimulus.nItems
            
            Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), calibration.FirstOrientation_deg(thisCalibrationTrial), [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);
                
        end
        
    end
    
    % Display first intervals of Gabor
    
    firstIntervalTime = Screen('Flip', ptbWindow, startTime + timing.blankDuration);
    
    % Display masks
    
    if stimulus.fixationOn

        Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

    end
    
    for thisGabor = 1:stimulus.nItems
    
        Screen('DrawTexture',ptbWindow,maskTexture, [],maskRects(:,thisGabor));
        
    end
    
    firstMaskTime = Screen('Flip',ptbWindow, firstIntervalTime + timing.displayDuration);
    
    % Display blank (with fixation) for interval

    if stimulus.fixationOn

        Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

    end
    
    blankTime = Screen('Flip', ptbWindow, firstMaskTime + timing.maskDuration);
    
    % Draw second set of Gabors
    
    if stimulus.fixationOn

        Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

    end
    
    if whichOne == 2 % Second interval is different

        for thisGabor = 1:stimulus.nItems-1

            Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), calibration.SecondOrientation_deg(thisCalibrationTrial), [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0])
        
        end

        if whichChangeDirection(thisCalibrationTrial) == 1
            
            Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), calibration.SecondOrientation_deg(thisCalibrationTrial)+thisOrientationDifference, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);
        
        elseif whichChangeDirection(thisCalibrationTrial) == 2
            
            Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), calibration.SecondOrientation_deg(thisCalibrationTrial)-thisOrientationDifference, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);
        
        end
        
    else

        for thisGabor = 1:stimulus.nItems

            Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), calibration.SecondOrientation_deg(thisCalibrationTrial), [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0])
        end

    end

    secondIntervalTime = Screen('Flip', ptbWindow, blankTime + timing.intervalDuration);
    
    % Draw mask
    
    if stimulus.fixationOn

        Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

    end
    
    for thisGabor = 1:stimulus.nItems
    
        Screen('DrawTexture',ptbWindow,maskTexture,[], maskRects(:,thisGabor));
        
    end
    
    secondMaskTime = Screen('Flip',ptbWindow,secondIntervalTime + timing.displayDuration);
    
    % Flip to blank
    
    finishTime = Screen('Flip', ptbWindow, secondMaskTime + timing.maskDuration);
    
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
    
    calibration.allCalibrationResponses(thisCalibrationTrial) = pressedKey;
    
    % Determine if response is correct
    
    if pressedKey + whichOne == 81; % Correct repsonses will always add up to this number
    
        ifCorrect = 1;

    else
        
        ifCorrect = 0;
        
    end
    
    calibration.allCalibrationCorrect(thisCalibrationTrial) = ifCorrect;
    
    % Update all staircases
    
    for thisStaircase = 1:experiment.nStaircases

        allData(thisStaircase) = QuestUpdate(allData(thisStaircase), thisLogOrientationDifference, ifCorrect);

    end

end

staircase.allFirstOrientations = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);
staircase.allSecondOrientations = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);
staircase.allDifferences = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);
staircase.allResponses = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);
staircase.allCorrect = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);

for thisStaircase = 1:experiment.nStaircases
    
    startStaircaseText = ['You will now begin block ' num2str(experiment.thisStaircase) ' of ' num2str(experiment.nStaircases) ' blocks.\n' ...
    'Press a response key to start the block.'];

    % Set up staircase parameters
    
    staircase.allChanges = mod(randperm(experiment.nTrialsPerStaircase),2)+1;           % Determines which pair will have different Gabors
    staircase.allChangeDirection = mod(randperm(experiment.nTrialsPerStaircase),2)+1;   % Determines which change direction

    staircase.FirstOrientation_deg = round(rand(1,experiment.nTrialsPerStaircase)*180);
    staircase.SecondOrientation_deg = round(rand(1,experiment.nTrialsPerStaircase)*180);

    staircase.FirstOrientation_rad = deg2rad(staircase.FirstOrientation_deg);
    staircase.SecondOrientation_rad = deg2rad(staircase.SecondOrientation_deg);

    % Set up location rects of Gabors for every trial

    staircase.itemBaseTheta_rad = repmat(linspace(0,2*pi-(2*pi/stimulus.nArrayItems), stimulus.nArrayItems)',1,experiment.nTrialsPerStaircase);
    staircase.arrayJitter_rad = ((2*pi)/stimulus.nArrayItems)*repmat(rand(1,experiment.nTrialsPerStaircase), stimulus.nArrayItems, 1);
    staircase.allItemTheta_rad = staircase.itemBaseTheta_rad+staircase.arrayJitter_rad;
    staircase.allItemTheta_deg = rad2deg(staircase.allItemTheta_rad);
    
    % Prevent Vernier discrimination task

    minimumDifference_deg = 45;                                 % Minimum difference between orientation and position in degrees
    minimumDifference_rad = deg2rad(minimumDifference_deg);     % Minimum difference between orientation and position in radians

    for thisStaircaseTrial = 1:experiment.nTrialsPerStaircase
        
        if staircase.FirstOrientation_rad(thisStaircaseTrial) < pi/2
    
            if deg2rad(90)-minimumDifference_rad < staircase.FirstOrientation_rad(thisStaircaseTrial) + staircase.allItemTheta_rad(1,thisStaircaseTrial) ...
                    && staircase.FirstOrientation_rad(thisStaircaseTrial) + staircase.allItemTheta_rad(1,thisStaircaseTrial) < deg2rad(90)+minimumDifference_rad

                staircase.FirstOrientation_rad(thisStaircaseTrial) = staircase.FirstOrientation_rad(thisCalibrationTrial) + (minimumDifference_rad*2);
                staircase.FirstOrientation_deg(thisStaircaseTrial) = rad2deg(staircase.FirstOrientation_rad(thisStaircaseTrial));

            end

        elseif staircase.FirstOrientation_rad(thisStaircaseTrial) >= pi/2

            if deg2rad(270)-minimumDifference_rad < staircase.FirstOrientation_rad(thisStaircaseTrial) + staircase.allItemTheta_rad(1,thisStaircaseTrial) ...
                    && staircase.FirstOrientation_rad(thisStaircaseTrial) + staircase.allItemTheta_rad(1,thisStaircaseTrial) < deg2rad(270)+minimumDifference_rad

                staircase.FirstOrientation_rad(thisStaircaseTrial) = staircase.FirstOrientation_rad(thisStaircaseTrial) + (minimumDifference_rad*2);
                staircase.FirstOrientation_deg(thisStaircaseTrial) = rad2deg(staircase.FirstOrientation_rad(thisStaircaseTrial));

            end 

        end

        if staircase.SecondOrientation_rad(thisStaircaseTrial) < pi/2

            if deg2rad(90)-minimumDifference_rad < staircase.SecondOrientation_rad(thisStaircaseTrial) + staircase.allItemTheta_rad(2,thisStaircaseTrial) ...
                    && staircase.SecondOrientation_rad(thisStaircaseTrial) + staircase.allItemTheta_rad(2,thisStaircaseTrial) < deg2rad(90)+minimumDifference_rad

                staircase.SecondOrientation_rad(thisStaircaseTrial) = staircase.SecondOrientation_rad(thisStaircaseTrial) + (minimumDifference_rad*2);
                staircase.SecondOrientation_deg(thisStaircaseTrial) = rad2deg(staircase.SecondOrientation_rad(thisStaircaseTrial));

            end

        elseif staircase.SecondOrientation_rad(thisStaircaseTrial) >= pi/2

            if deg2rad(270)-minimumDifference_rad < staircase.SecondOrientation_rad(thisStaircaseTrial) + staircase.allItemTheta_rad(2,thisStaircaseTrial) ...
                    && staircase.SecondOrientation_rad(thisStaircaseTrial) + staircase.allItemTheta_rad(2,thisStaircaseTrial) < deg2rad(270)+minimumDifference_rad

                staircase.SecondOrientation_rad(thisStaircaseTrial) = staircase.SecondOrientation_rad(thisStaircaseTrial) + (minimumDifference_rad*2);
                staircase.SecondOrientation_deg(thisCalibrationTrial) = rad2deg(staircase.SecondOrientation_rad(thisStaircaseTrial));

            end 

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

    for thisStaircaseTrial = 1:experiment.nTrialsPerStaircase

        % On each trial, create a new noise array to change mask texture
    
        % Wipe previous mask textures to clear memory

        if thisStaircaseTrial ~= 1

            Screen('Close',maskTexture);

        end

        % Create a noise image the same size as your Gabor texture, randomising each
        % pixel's luminance value (uniform over -1,1).

        a = 0;
        b = 1;
        noiseArray = a + (b-a).*rand(size(imageArray));

        % Take the Fourier transform of the noise image and set aside the phase
        % spectrum.

        fftnoiseImage = fft2(noiseArray,SFT_Size,SFT_Size);
        frnoiseImage = abs(fftnoiseImage);
        phnoiseImage = angle(fftnoiseImage);

        % Take the inverse transform using the amplitude spectrum and phase
        % spectrum.
        maskednoise = ifft2(maskwhatever .* (phnoiseImage*exp(1i)), SFT_Size, SFT_Size);

        % Crop
        middleRow = floor(SFT_Size/2);
        nRows = round(stimulus.size_pix/2);
        startRow = middleRow - nRows + 1;
        endRow = middleRow + nRows;

        maskednoise = real(maskednoise(startRow:endRow,startRow:endRow));

        % Normalise
        maskednoise = (maskednoise-min(maskednoise(:))) / (max(maskednoise(:)) - min(maskednoise(:)));
        maskednoise = (maskednoise*2)-1; % Balance
        mask = blobArray.*maskednoise;
        mask = (mask-min(mask(:)))/(max(mask(:))-min(mask(:)));
        rgbImage = repmat(uint8(255.*mask),[1 1 3]);
        maskTexture = Screen('MakeTexture',ptbWindow,rgbImage);
        
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

        whichOne = staircase.allChanges(thisStaircaseTrial);

        % Determine what the size of the difference in orientation will be

        thisLogOrientationDifference = QuestMean(allData(thisStaircase));
        thisOrientationDifference = 10.^thisLogOrientationDifference;
        staircase.allDifferences(thisStaircase,thisStaircaseTrial) = thisOrientationDifference;
        staircase.allLogDifferences(thisStaircase,thisStaircaseTrial) = thisLogOrientationDifference;

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

                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), staircase.FirstOrientation_deg(thisStaircaseTrial), [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

            end

            if staircase.allChangeDirection(thisStaircaseTrial) == 1
                
                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), staircase.FirstOrientation_deg(thisStaircaseTrial)+thisOrientationDifference, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

            elseif staircase.allChangeDirection(thisStaircaseTrial) == 2
                
                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), staircase.FirstOrientation_deg(thisStaircaseTrial)-thisOrientationDifference, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

            end
                
        else

            for thisGabor = 1:stimulus.nItems

                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), staircase.FirstOrientation_deg(thisStaircaseTrial), [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

            end

        end

        % Display first intervals of Gabor

        firstIntervalTime = Screen('Flip', ptbWindow, startTime + timing.blankDuration);
        
        % Draw mask
        
        if stimulus.fixationOn

            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

        end
    
        for thisGabor = 1:stimulus.nItems
    
            Screen('DrawTexture',ptbWindow,maskTexture, [],itemRects(:,thisGabor));
        
        end
    
        firstMaskTime = Screen('Flip', ptbWindow, firstIntervalTime + timing.displayDuration);
        
        % Display blank (with fixation) for interval

        if stimulus.fixationOn

            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

        end

        blankTime = Screen('Flip', ptbWindow, firstMaskTime + timing.maskDuration);

        % Draw second set of Gabors

        if stimulus.fixationOn

            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

        end

        if whichOne == 2 % Second interval is different

           for thisGabor = 1:stimulus.nItems-1

                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), staircase.SecondOrientation_deg(thisStaircaseTrial), [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0])
            
           end

           if staircase.allChangeDirection(thisStaircaseTrial) == 1
               
               Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), staircase.SecondOrientation_deg(thisStaircaseTrial)+thisOrientationDifference, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

           elseif staircase.allChangeDirection(thisStaircaseTrial) == 2
               
               Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,stimulus.nItems), staircase.SecondOrientation_deg(thisStaircaseTrial)-thisOrientationDifference, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

           end
           
        else

            for thisGabor = 1:stimulus.nItems

                Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), staircase.SecondOrientation_deg(thisStaircaseTrial), [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0])
            end

        end

        secondIntervalTime = Screen('Flip', ptbWindow, blankTime + timing.intervalDuration);
        
        % Draw mask
        
         if stimulus.fixationOn

            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

        end
    
        for thisGabor = 1:stimulus.nItems
    
            Screen('DrawTexture',ptbWindow,maskTexture, [],itemRects(:,thisGabor));
        
        end
        
        secondMaskTime = Screen('Flip', ptbWindow, secondIntervalTime + timing.displayDuration);
        
        % Flip to blank

        finishTime = Screen('Flip', ptbWindow, secondMaskTime + timing.maskDuration);

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

        staircase.allResponses(thisStaircase,thisStaircaseTrial) = pressedKey;

        % Determine if response is correct

        if pressedKey + whichOne == 81; % Correct repsonses will always add up to this number

            ifCorrect = 1;

        else

            ifCorrect = 0;

        end

        staircase.allCorrect(thisStaircase,thisStaircaseTrial) = ifCorrect;

        % Update staircase

        allData(thisStaircase) = QuestUpdate(allData(thisStaircase), thisLogOrientationDifference, ifCorrect);

    end

    % Save staircase parameters
    
    staircase.allFirstOrientations(thisStaircase,:) = staircase.FirstOrientation_deg;
    staircase.allSecondOrientations(thisStaircase,:) = staircase.SecondOrientation_deg;
    experiment.thisStaircase = experiment.thisStaircase+1;
    
end

% Show end text

DrawFormattedText(ptbWindow,finishText, 'center', 'center', colour.textVal);
Screen('Flip',ptbWindow);

% Save data file
cd(saveDirectory);
fileNumber = 1;
newFileName = [participant.ID '_IDinVWM_' num2str(fileNumber) '.mat'];

while exist([saveDirectory newFileName],'file')
fileNumber = fileNumber+1;
newFileName = [participant.ID '_IDinVWM_' num2str(fileNumber) '.mat'];
end

save(newFileName, 'colour', 'experiment', 'stimulus', 'participant', 'staircase', 'stimulus', 'calibration', 'startTime', 'finishTime');

% Save user data file
cd(userDirectory);
save(participant.userFile, 'participant', 'allData');

% Plot each staircase
for thisStaircase = 1:experiment.nStaircases
    
    staircaseFigureName = ['Staircase:_' num2str(thisStaircase)];
    figure('Color', 'white', 'Name', staircaseFigureName);
    subplot(1,2,1);
    plot(1:experiment.nTrialsPerStaircase, staircase.allDifferences(thisStaircase,:),'k-');
    subplot(1,2,2);
    plot(1:experiment.nTrialsPerStaircase, staircase.allLogDifferences(thisStaircase,:),'k-');
    
end
