% Individual differences in perceptual ability and visual working memory

% This study aims to examine whether an individual's ability to
% discriminate between different orientation and spatial frequency of Gabor
% patches predicts the individual's visual working memory (VWM) ability on 
% a change-detection task.

% Run_IDinVWM.m should be run before this
% This code is for the orientation change-detection task in Experiment 1 of 
% WN's PhD

% WN started writing this July 2015

% -------------------------------------------------------------------------

% Get participant parameters from Run_IDinVWM.

% AssertOpenGL;

% Set up experiment parameters

experiment.nCalibrationTrials = 20;
experiment.nPracticeTrials = 10;

% Set up equipment parameters

equipment.viewDist = 500;                           % Viewing distance in mm
equipment.ppm = 2.7;                                % Pixels per mm (HP P1120, 1024 x 768, 120 Hz) % Remeasure
equipment.gammaVals = 1.0./[2.6434 2.2312 2.171];   % Gamma values for CRT in GT519

colour.fixVal = 1;
colour.textVal = 0;

% Set up colour parameters

colour.blackVal = 0;
colour.greyVal = 0.5;
colour.whiteVal = 1;

% Set up stimulus parameters

stimulus.size_dva = 4;            % Item size in degrees of visual angle
stimulus.eccentricity_dva = 4;      % Eccentricity in degrees of visual angle
stimulus.nArrayItems = 2;           % Number of items shown in array
stimulus.nItems = 2;                % Number of items in discrimination task
stimulus.maxDiff = 90;              % Maximum difference (degrees)

stimulus.fixationOn = 1;            % Fixation on (1) or off (0)
stimulus.fixationSize_dva = .25;    % Fixation size in degrees of visual angle

stimulus.frequency_cpd = 5;     % Spatial frequency of Gabor patches (in cycles per degree)
stimulus.sc = 1/2.5;                % Spatial constant of Gabor patches (in degrees)
stimulus.contrast = 1;              % Contrast of Gabor patches
stimulus.aspectratio = 1;           % Aspect ratio of Gabor patches

% Set up audio parameters

audio.toneLength = .1;
audio.toneFreq = 880;
audio.toneAmplitude = 1.0;
audio.sampleRate = 48000;           % Default audio sample rate
defaultTone = audio.toneAmplitude*sin(linspace(0,2*pi*audio.toneFreq*audio.toneLength,audio.sampleRate*audio.toneLength));

% Set up temporal parameters

timing.blankDuration = 1;         % Duration of inter-trial interval (blank) (after response?)
timing.fixationDuration = 0.5;    % Duration of fixation point
timing.memoryDuration = 1;        % Duration of memory array in change-detection task
timing.maskDuration = 0.5;        % Duration of mask in change-detection task
timing.SOADuration = 0.5;         % Duration of stimulus onset asynchrony between memory arrays

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

taskExplanationText = ['For this block, you will be shown an array of stimuli at two different intervals.\n' ...
    'The first array will contain ' num2str(stimulus.nArrayItems) ' gratings, followed by a mask.\n' ...
    'A second array will follow that is identical except for one grating.\n' ...
    'Press any key for the next set of instructions'];
    
responseInstructionText = ['You need to respond with which grating you think has changed.\n' ...
    'Click on the grating that you think has changed using the mouse.\n' ...
    'Press any key to continue.'];

startCalibrationText = ['You will begin with some calibration trials.\n' ...
    'Press any key to start the calibration trials.'];

startBlockText = ['Press any key to start the block'];

% ----------------------------------------------------------------------- %

%                           Begin Experiment

% ----------------------------------------------------------------------- %

%                        Orientation Change-Detection                     %

% Create staircases

for thisStaircase = 1:experiment.nStaircases
    
    allData(thisStaircase) = QuestCreate(staircase.tGuess,staircase.tGuessSD,staircase.pThreshold,staircase.beta,staircase.delta,staircase.gamma,staircase.grain,staircase.range);
    
end

% Estimate staircase parameters

% Set up calibration parameters

differencesToTest_deg = linspace(90,0.1,experiment.nCalibrationTrials);
differencesToTest_rad = deg2rad(differencesToTest_deg);

% Might need a change on every trial for the staircase to work.
% changeOrNoChange = mod(randperm(experiment.nCalibrationTrials),2)+1;      % Determines whether the trial will be same or different
allChanges = mod(randperm(experiment.nCalibrationTrials),stimulus.nArrayItems)+1;              % Randomises which Gabor changes on every changed trial
allChangeDirection = mod(randperm(experiment.nCalibrationTrials),2)+1;      % Determines the direction of the change on each trial (clockwise or anticlockwise)

calibration.allOrientations_deg = NaN(stimulus.nArrayItems, experiment.nCalibrationTrials);
calibration.allOrientations_rad = NaN(stimulus.nArrayItems, experiment.nCalibrationTrials);

% Randomise the orientation of all items in the array (for each trial)
    
calibration.allOrientations_deg = rand(stimulus.nArrayItems,experiment.nCalibrationTrials)*180;

% Prevent orientations of all items in the array being too similar to each
% other

for thisTrial = 1:experiment.nCalibrationTrials
    
    minimumDifference = 30;
    
    % Create a matrix of the differences
    
    differenceMatrix = zeros(stimulus.nArrayItems);
    checking = 1;
    
    while checking
    
        for thisReferenceItem = 1:stimulus.nArrayItems

            for thisComparisonItem = 1:stimulus.nArrayItems

                differenceMatrix(thisReferenceItem,thisComparisonItem) = calibration.allOrientations_deg(thisComparisonItem,thisTrial) - calibration.allOrientations_deg(thisReferenceItem,thisTrial);

            end

        end

        % Check the values of the matrix

        tooSimilar = find(differenceMatrix > 0 & differenceMatrix < minimumDifference);
        isEmpty = isempty(tooSimilar);
        
        if isEmpty == 0
        
            howManyAreSimilar = numel(tooSimilar);

            for thisDifference = 1:howManyAreSimilar

                theReferenceOne = ceil(tooSimilar(thisDifference)/stimulus.nArrayItems);
                theComparisonOne = mod(tooSimilar(thisDifference),stimulus.nArrayItems);
                
                if theComparisonOne == 0
                    
                    theComparisonOne = stimulus.nArrayItems;
                    
                end

                if calibration.allOrientations_deg(theReferenceOne,thisTrial) >= calibration.allOrientations_deg(theComparisonOne,thisTrial)

                    calibration.allOrientations_deg(theReferenceOne,thisTrial) = calibration.allOrientations_deg(theReferenceOne,thisTrial) + minimumDifference/2;
                    calibration.allOrientations_deg(theComparisonOne,thisTrial) = calibration.allOrientations_deg(theComparisonOne,thisTrial) - minimumDifference/2;
                   
                    if calibration.allOrientations_deg(theReferenceOne,thisTrial) > 180
                        
                        calibration.allOrientations_deg(theReferenceOne,thisTrial) = calibration.allOrientations_deg(theReferenceOne,thisTrial) - 180;
                        
                    end
                    
                    if calibration.allOrientations_deg(theComparisonOne,thisTrial) < 0
                        
                        calibration.allOrientations_deg(theComparisonOne,thisTrial) = calibration.allOrientations_deg(theComparisonOne,thisTrial) + 180;
                        
                    end
                    
                elseif calibration.allOrientations_deg(theReferenceOne,thisTrial) < calibration.allOrientations_deg(theComparisonOne,thisTrial)

                    calibration.allOrientations_deg(theReferenceOne,thisTrial) = calibration.allOrientations_deg(theReferenceOne,thisTrial) - minimumDifference/2;
                    calibration.allOrientations_deg(theComparisonOne,thisTrial) = calibration.allOrientations_deg(theComparisonOne,thisTrial) + minimumDifference/2;
                    
                    if calibration.allOrientations_deg(theReferenceOne,thisTrial) < 0
                        
                        calibration.allOrientations_deg(theReferenceOne,thisTrial) = calibration.allOrientations_deg(theReferenceOne,thisTrial) + 180;
                        
                    end
                    
                    if calibration.allOrientations_deg(theComparisonOne,thisTrial) > 180
                        
                        calibration.allOrientations_deg(theComparisonOne,thisTrial) = calibration.allOrientations_deg(theComparisonOne,thisTrial) - 180;
                        
                    end
                    
                end

            end
    
        elseif isEmpty == 1
            
            checking = 0;
            
        end
        
    end    
    
end

calibration.allOrientations_rad = deg2rad(calibration.allOrientations_deg);

% Set up location parameters

calibration.itemBaseTheta_rad = repmat(linspace(0,2*pi-(2*pi/stimulus.nArrayItems), stimulus.nArrayItems)',1,experiment.nCalibrationTrials);
calibration.arrayJitter_rad = ((2*pi)/stimulus.nArrayItems)*repmat(rand(1,experiment.nCalibrationTrials), stimulus.nArrayItems, 1);
calibration.allItemTheta_rad = calibration.itemBaseTheta_rad+calibration.arrayJitter_rad;
calibration.allItemTheta_deg = rad2deg(calibration.allItemTheta_rad);

calibration.allCorrect = NaN(1,experiment.nCalibrationTrials);

% % Prevent Vernier discrimination
% 
% minimumDifference_deg = 45;
% minimumDifference_rad = deg2rad(minimumDifference_deg);
% 
% for thisTrial = 1:experiment.nCalibrationTrials
%     
%     for thisItem = 1:stimulus.nArrayItems
%         
%         if calibration.allOrientations_rad(thisItem,thisTrial) < pi/2
%             
%             if deg2rad(90)-minimumDifference_rad < calibration.allOrientations_rad(thisItem,thisTrial) + calibration.allItemTheta_rad(thisItem,thisTrial) ...
%                 && calibration.allOrientations_rad(thisItem,thisTrial) + calibration.allItemTheta_rad(thisItem,thisTrial) < deg2rad(90)+minimumDifference_rad
% 
%                 calibration.allOrientations_rad(thisItem,thisTrial) = calibration.allOrientations_rad(thisItem,thisTrial) + (minimumDifference_rad*2);
%                 calibration.allOrientations_rad(thisItem,thisTrial) = rad2deg(calibration.allOrientations_rad(thisItem,thisTrial));
% 
%             end
%             
%         elseif calibration.allOrientations_rad(thisItem,thisTrial) > pi/2
%             
%             if deg2rad(270)-minimumDifference_rad < calibration.allOrientations_rad(thisItem,thisTrial) + calibration.allItemTheta_rad(thisItem,thisTrial) ...
%                 && calibration.allOrientations_rad(thisItem,thisTrial) + calibration.allItemTheta_rad(thisItem,thisTrial) < deg2rad(270)+minimumDifference_rad
%         
%                 calibration.allOrientations_rad(thisItem,thisTrial) = calibration.allOrientations_rad(thisItem,thisTrial) + (minimumDifference_rad*2);
%                 calibration.allOrientations_rad(thisItem,thisTrial) = rad2deg(calibration.allOrientations_rad(thisItem,thisTrial));
%         
%             end 
%         
%         end
%         
%     end
%     
% end

% Wait for response

% Instruction Text

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
    waitResponse = 0;

end

% Wait for keypress response to start calibration

DrawFormattedText(ptbWindow,startCalibrationText, 'center', 'center', colour.textVal);
Screen('Flip',ptbWindow);
waitResponse = 1;

while waitResponse

    [startCalibrationTime, keyCode] = KbWait(-1,2);
    waitResponse = 0;
    
end
    
responseTime = startCalibrationTime;

% Begin change-detection task

% Begin calibration

% Create Gabor textures

[gaborid, gaborrect] = CreateProceduralGabor(ptbWindow,stimulus.size_pix,stimulus.size_pix, 0, [.5 .5 .5 0], 1, .5);
   
thisChange = 1;

for thisCalibrationTrial = 1:experiment.nCalibrationTrials
    
    % On each trial, create a new noise array to change mask texture
    
    % Wipe previous mask textures to save memory
    
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
    
    HideCursor;
    
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
            
    % Draw memory array
    
    if stimulus.fixationOn

        Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

    end
    
    for thisGabor = 1:stimulus.nArrayItems
        
        Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), calibration.allOrientations_deg(thisGabor,thisCalibrationTrial), [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);
            
    end
    
    memoryTime = Screen('Flip', ptbWindow, startTime + timing.fixationDuration);
    
    % Draw mask
    
    if stimulus.fixationOn

        Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

    end

    for thisGabor = 1:stimulus.nArrayItems

        Screen('DrawTexture',ptbWindow,maskTexture, [],itemRects(:,thisGabor));

    end
    
    maskTime = Screen('Flip',ptbWindow, memoryTime + timing.memoryDuration);
    
    % Blank until test array
    
    blankTime = Screen('Flip',ptbWindow, maskTime + timing.maskDuration);
    
    % Draw test array (will always be different)
    
    whichChanged = allChanges(thisCalibrationTrial);                    % Tells you which Gabor is changed
    whichChangeDirection = allChangeDirection(thisCalibrationTrial);   % Tells you which direction the Gabor is changed  
    
    for thisGabor = 1:stimulus.nArrayItems

        if thisGabor ~= whichChanged

            Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), calibration.allOrientations_deg(thisGabor,thisCalibrationTrial), [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

        elseif thisGabor == whichChanged
            
            if whichChangeDirection == 1; % Add orientation difference on (anticlockwise change)

            Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), calibration.allOrientations_deg(thisGabor,thisCalibrationTrial)+thisOrientationDifference, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);
            
            elseif whichChangeDirection == 2; % Subtract orientation difference (clockwise change)
                
            Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), calibration.allOrientations_deg(thisGabor,thisCalibrationTrial)-thisOrientationDifference, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);
            
            end
            
        end
        
    end
        
    testTime = Screen('Flip', ptbWindow, blankTime + timing.blankDuration);
    
    % Flip to blank
    
%     waitResponseTime = Screen('Flip', ptbWindow, testTime + timing.memoryDuration);
%     
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

    % Flip to blank.
    
    endTrialTime = Screen('Flip',ptbWindow);
    
    ClickedTarget = find(CheckResponse);
    calibration.allResponses(thisCalibrationTrial) = ClickedTarget;

    if ClickedTarget == allChanges(thisCalibrationTrial)
        
        ifCorrect = 1;
        
    elseif ClickedTarget ~= allChanges(thisCalibrationTrial)
        
        ifCorrect = 0;
        
    end

    calibration.allCorrect(thisCalibrationTrial) = ifCorrect;
    
    % Update all staircases
    
    allData(thisStaircase) = QuestUpdate(allData(thisStaircase), thisLogOrientationDifference, ifCorrect);

end

% Begin staircase

% Set up staircase parameters

staircase.allOrientations = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);
staircase.allDifferences = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);
staircase.allResponses = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);
staircase.allCorrect = NaN(experiment.nStaircases,experiment.nTrialsPerStaircase);

for thisStaircase = 1:experiment.nStaircases

    startStaircaseText = ['You will now begin block ' num2str(thisStaircase) ' of ' num2str(experiment.nStaircases) ' blocks.\n' ...
    'Press any key to start the block.'];
    
    % Set up staircase parameters
    
    staircase.allChanges = mod(randperm(experiment.nTrialsPerStaircase),stimulus.nArrayItems)+1;        % Determines which of the Gabors in each array changes
    staircase.allChangeDirection = mod(randperm(experiment.nTrialsPerStaircase),2)+1;                   % Determines which change direction

    % Determine the orientations of each Gabor in the array
        
    staircase.allOrientations_deg = round(rand(stimulus.nArrayItems,experiment.nTrialsPerStaircase)* 180);
    
    % Prevent orientations of all items in the array being too similar to each
    % other

    for thisTrial = 1:experiment.nTrialsPerStaircase

        minimumDifference = 30;

        % Create a matrix of the differences

        differenceMatrix = zeros(stimulus.nArrayItems);
        checking = 1;

        while checking

            for thisReferenceItem = 1:stimulus.nArrayItems

                for thisComparisonItem = 1:stimulus.nArrayItems

                    differenceMatrix(thisReferenceItem,thisComparisonItem) = staircase.allOrientations_deg(thisComparisonItem,thisTrial) - staircase.allOrientations_deg(thisReferenceItem,thisTrial);

                end

            end

            % Check the values of the matrix

            tooSimilar = find(differenceMatrix > 0 & differenceMatrix < minimumDifference);
            isEmpty = isempty(tooSimilar);

            if isEmpty == 0

                howManyAreSimilar = numel(tooSimilar);

                for thisDifference = 1:howManyAreSimilar

                    theReferenceOne = ceil(tooSimilar(thisDifference)/stimulus.nArrayItems);
                    theComparisonOne = mod(tooSimilar(thisDifference),stimulus.nArrayItems);

                    if theComparisonOne == 0

                        theComparisonOne = stimulus.nArrayItems;

                    end

                    if staircase.allOrientations_deg(theReferenceOne,thisTrial) >= staircase.allOrientations_deg(theComparisonOne,thisTrial)

                        staircase.allOrientations_deg(theReferenceOne,thisTrial) = staircase.allOrientations_deg(theReferenceOne,thisTrial) + minimumDifference/2;
                        staircase.allOrientations_deg(theComparisonOne,thisTrial) = staircase.allOrientations_deg(theComparisonOne,thisTrial) - minimumDifference/2;

                        if staircase.allOrientations_deg(theReferenceOne,thisTrial) > 180

                            staircase.allOrientations_deg(theReferenceOne,thisTrial) = staircase.allOrientations_deg(theReferenceOne,thisTrial) - 180;

                        end

                        if staircase.allOrientations_deg(theComparisonOne,thisTrial) < 0

                            staircase.allOrientations_deg(theComparisonOne,thisTrial) = staircase.allOrientations_deg(theComparisonOne,thisTrial) + 180;

                        end

                    elseif staircase.allOrientations_deg(theReferenceOne,thisTrial) < staircase.allOrientations_deg(theComparisonOne,thisTrial)

                        staircase.allOrientations_deg(theReferenceOne,thisTrial) = staircase.allOrientations_deg(theReferenceOne,thisTrial) - minimumDifference/2;
                        staircase.allOrientations_deg(theComparisonOne,thisTrial) = staircase.allOrientations_deg(theComparisonOne,thisTrial) + minimumDifference/2;

                        if staircase.allOrientations_deg(theReferenceOne,thisTrial) < 0

                            staircase.allOrientations_deg(theReferenceOne,thisTrial) = staircase.allOrientations_deg(theReferenceOne,thisTrial) + 180;

                        end

                        if staircase.allOrientations_deg(theComparisonOne,thisTrial) > 180

                            staircase.allOrientations_deg(theComparisonOne,thisTrial) = staircase.allOrientations_deg(theComparisonOne,thisTrial) - 180;

                        end

                    end

                end

            elseif isEmpty == 1

                checking = 0;

            end

        end    

    end
    
    staircase.allOrientations_rad = deg2rad(staircase.allOrientations_deg);

    % Set up location rects of Gabors for every trial

    staircase.itemBaseTheta_rad = repmat(linspace(0,2*pi-(2*pi/stimulus.nArrayItems), stimulus.nArrayItems)',1,experiment.nTrialsPerStaircase);
    staircase.arrayJitter_rad = ((2*pi)/stimulus.nArrayItems)*repmat(rand(1,experiment.nTrialsPerStaircase), stimulus.nArrayItems, 1);
    staircase.allItemTheta_rad = staircase.itemBaseTheta_rad+staircase.arrayJitter_rad;
    staircase.allItemTheta_deg = rad2deg(staircase.allItemTheta_rad);
    
%     % Prevent Vernier discrimination task
% 
%     minimumDifference_deg = 45;                                 % Minimum difference between orientation and position in degrees
%     minimumDifference_rad = deg2rad(minimumDifference_deg);     % Minimum difference between orientation and position in radians
% 
%     for thisTrial = 1:experiment.nTrialsPerStaircase
%     
%         for thisItem = 1:stimulus.nArrayItems
% 
%             if staircase.allOrientations_rad(thisItem,thisTrial) < pi/2
% 
%                 if staircase.allItemTheta_rad(thisItem,thisTrial) > 0 && staircase.allItemTheta_rad(thisItem,thisTrial) < pi/2
% 
%                     staircase.allOrientations_rad(thisItem,thisTrial) = staircase.allOrientations_rad(thisItem,thisTrial) + pi/2;
%                     staircase.allOrientations_deg(thisItem,thisTrial) = rad2deg(staircase.allOrientations_rad(thisItem,thisTrial));
% 
%                 elseif staircase.allItemTheta_rad(thisItem,thisTrial) > pi && staircase.allItemTheta_rad(thisItem,thisTrial) < 3*pi/2
% 
%                     staircase.allOrientations_rad(thisItem,thisTrial) = staircase.allOrientations_rad(thisItem,thisTrial) + pi/2;
%                     staircase.allOrientations_deg(thisItem,thisTrial) = rad2deg(staircase.allOrientations_rad(thisItem,thisTrial));
% 
%                 end
% 
%             elseif staircase.allOrientations_rad(thisItem,thisTrial) > pi/2
% 
%                 if staircase.allItemTheta_rad(thisItem,thisTrial) > pi/2 && staircase.allItemTheta_rad(thisItem,thisTrial) < pi
% 
%                     staircase.allOrientations_rad(thisItem,thisTrial) = staircase.allOrientations_rad(thisItem,thisTrial) + pi/2;
%                     staircase.allOrientations_deg(thisItem,thisTrial) = rad2deg(staircase.allOrientations_rad(thisItem,thisTrial));
% 
%                 elseif staircase.allItemTheta_rad(thisItem,thisTrial) > 3*pi/2 && staircase.allItemTheta_rad(thisItem,thisTrial) < 2*pi
% 
%                     staircase.allOrientations_rad(thisItem,thisTrial) = staircase.allOrientations_rad(thisItem,thisTrial) + pi/2;
%                     staircase.allOrientations_deg(thisItem,thisTrial) = rad2deg(staircase.allOrientations_rad(thisItem,thisTrial));
% 
%                 end
% 
%             end
% 
%         end
% 
%     end

    
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

        % Wipe previous mask textures
        
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

        for thisGabor = 1:stimulus.nItems

            Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), staircase.allOrientations_deg(thisGabor,thisStaircaseTrial), [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

        end

        % Display first intervals of Gabor

        firstIntervalTime = Screen('Flip', ptbWindow, startTime + timing.blankDuration);

        % Display mask
        
        if stimulus.fixationOn

            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

        end
    
        for thisGabor = 1:stimulus.nItems
    
            Screen('DrawTexture',ptbWindow,maskTexture, [],itemRects(:,thisGabor));
        
        end
        
        maskTime = Screen('Flip', ptbWindow, firstIntervalTime + timing.memoryDuration);
        
        % Display blank (with fixation) for interval

        if stimulus.fixationOn

            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

        end

        blankTime = Screen('Flip', ptbWindow, maskTime + timing.maskDuration);

        % Draw test Array of Gabors

        if stimulus.fixationOn

            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);

        end

        for thisGabor = 1:stimulus.nArrayItems
            
            if thisGabor ~= staircase.allChanges(thisStaircaseTrial)
                
               Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), staircase.allOrientations_deg(thisGabor,thisStaircaseTrial), [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

            elseif thisGabor == staircase.allChanges(thisStaircaseTrial)

               if staircase.allChangeDirection(thisStaircaseTrial) == 1

                   Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), staircase.allOrientations_deg(thisGabor,thisStaircaseTrial)+thisOrientationDifference, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

               elseif staircase.allChangeDirection(thisStaircaseTrial) == 2

                   Screen('DrawTexture', ptbWindow, gaborid, [], itemRects(:,thisGabor), staircase.allOrientations_deg(thisGabor,thisStaircaseTrial)-thisOrientationDifference, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);

               end
           
            end
            
        end

        testTime = Screen('Flip', ptbWindow, blankTime + timing.blankDuration);

%         % Flip to blank
% 
%         finishTime = Screen('Flip', ptbWindow, secondIntervalTime + timing.displayDuration);

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

        HideCursor;

        endTrialTime = Screen('Flip',ptbWindow);
        
        staircase.allResponses(thisStaircase,thisStaircaseTrial) = ClickedTarget;

        % Determine if response is correct

        if ClickedTarget == staircase.allChanges(thisStaircaseTrial)

            ifCorrect = 1;

        elseif ClickedTarget ~= staircase.allChanges(thisStaircaseTrial)

            ifCorrect = 0;

        end

        staircase.allCorrect(thisStaircase,thisStaircaseTrial) = ifCorrect;

        % Update staircase

        allData(thisStaircase) = QuestUpdate(allData(thisStaircase), thisLogOrientationDifference, ifCorrect);

    end

    experiment.thisStaircase = experiment.thisStaircase+1;
    
end

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

    

