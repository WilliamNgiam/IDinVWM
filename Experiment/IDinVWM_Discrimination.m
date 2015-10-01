% Individual differences in perceptual ability and visual working memory

% This study aims to examine whether an individual's ability to
% discriminate between different orientation and spatial frequency of Gabor
% patches predicts the individual's visual working memory (VWM) ability on 
% a change-detection task.

% This code is for the orientation discrimination task in Experiment 1 of 
% WN's PhD

% WN started writing this June 2015

% -------------------------------------------------------------------------

% Set up equipment parameters

equipment.viewDist = 500;                           % Viewing distance in mm
equipment.ppm = 2.7;                                % Pixels per mm (HP P1120, 1024 x 768, 120 Hz) % Remeasure
equipment.gammaVals = 1.0./[2.6434 2.2312 2.171];    % Gamma values for CRT in GT519

% Set up colour parameters

colour.blackVal = 0;
colour.greyVal = 0.5;
colour.whiteVal = 1;

% Set up stimulus parameters

stimulus.size_dva = 2.5;            % Item size in degrees of visual angle
stimulus.eccentricity_dva = 4;      % Eccentricity in degrees of visual angle
stimulus.nArrayItems = 2;           % Number of items shown in array

stimulus.nItems = 2;                % Number of items in discrimination task

% Set up temporal parameters

stimulus.fixationDuration = 0.5;    % Duration fixation point is displayed
stimulus.displayDuration = 0.5;     % Duration Gabors are displayed for each interval
stimulus.intervalDuration = 1;      % Duration between intervals
stimulus.blankDuration = 1;         % Duration of inter-trial interval (blank)

% Set up staircase parameters

staircase.pThreshold = .82; % Optimal for 2-AFC
staircase.tGuess = log10(PAL_Weibull(PFparams,staircase.pThreshold, 'Inverse')); % Gives the threshold from the Weibull PF fit
staircase.tGuessSD = 0.5; % Not sure if this is appropriate
staircase.beta = 3.5; % Standard
staircase.delta = lapseValues(sd); % Lapse rate
staircase.gamma = .5; % Guess rate, 2-AFC 

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
screenCentreX = round(ScreenWidth/2);
screenCentreY = round(ScreenHeigh/2);
flipInterval = Screen('GetFlipInterval', ptbWindow);

% Calculate equipment parameters

equipment.mpd = (equipment.viewDist)*tan(deg2rad(2*stimulus.eccentricity))/stimulus.eccentricity; % Calculate mm per degree of visual angle to the ecccentricity of the stimuli
equipment.ppd = equipment.ppm*equipment.mpd;

% Calculate spatial parameters

stimulus.size_pix = stimulus.size*equipment.ppd;                     % Item size in pixels
stimulus.eccentricity_pix = stimulus.eccentricity*equipment.ppd;     % Eccentricity of stimulus in pixels

% ----------------------------------------------------------------------- %

%                           Begin Experiment

% ----------------------------------------------------------------------- %

%                       Orientation Discrimination                        %

% Create staircases

for thisStaircase = 1:experiment.nStaircases
    
    allData(thisparticipantNum,thisID) = QuestCreate(staircase.tGuess,staircase.tGuessSD,staircase.beta,staircase.delta,staircase.gamma);
    
% Instruction Text

odExplanationText = ['For this block, you will be shown a pair of stimuli at two different intervals./n/n' ...
    'In one pair, the orientation of the stimuli will be different to each other./n/n' ...
    'While the other pair will have the same orientation./n/n' ...
    'You need to respond with which pair you think was different.../n/n'];
    
responseInstructionText = ['If you think the first pair was different, respond with the "F" keypress./n/n' ...
    'If you think the second pair was different, respond with the "J" keypress./n/n'];

practiceResponseText = ['Press any key to start the practice trials'];

experimentResponseText = ['Press any key to start the block'];
        

% Wait for response

for thisStaircase = 1:experiment.nStaircases

    % Set up trial parameters

    whichInterval = mod(randperm(experiment.nTrialsPerStaircase),2);    % Determines which pair will have different Gabors

    % Determine orientation of Gabors for each trial

    gabor.FirstOrientation_rad = (round(rand(1,experiment.nTrialsPerStaircase)*180))*2*pi/180;  % Quick fix to get orientations to display
    gabor.SecondOrientation_rad = (round(rand(1,experiment.nTrialsPerStaircase)*180))*2*pi/180;
    
    gabor.OrientationDifference = 15*2*pi/180;   % Set for now, should be determined by staircase

    % Set up location rects of Gabors for every trial

    itemBaseTheta_rad = repmat(linspace(0,2*pi-(2*pi/stimulus.nArrayItems), stimulus.nArrayItems)',1,experiment.nTrialsPerStaircase);
    arrayJitter_rad = ((2*pi)/stimulus.nArrayItems)*repmat(rand(1,experiment.nTrialsPerStaircase), stimulus.nArrayItems, 1);
    allItemTheta_rad = itemBaseTheta_rad+arrayJitter_rad;

    % Restrict location of Gabors to prevent Vernier discrimination
    
    for thisTrial = 1:experiment.nTrialsPerStaircase
        
        minimumDifference_deg = 45;                                 % Minimum difference between orientation and position in degrees
        minimumDifference_rad = deg2rad(minimumDifference_Deg);     % Minimum difference between orientation and position in radians

        if abs(GaborFirstOrientation_rad(thisTrial) - allItemTheta_rad(thisTrial)) < minimumDifference_rad
        
            GaborFirstOrientation_rad(thisTrial) = GaborOrientation_rad(thisTrial) + minimumDifference_rad
        
        else 
        
            % Do nothing
        
        end
        
        if abs(GaborSecondOrientation_rad(thisTrial) - allItemTheta_rad(thisTrial)) < minimumDifference_rad
            
            GaborSecondOrientation_rad(thisTrial) = GaborOrientation_rad(thisTrial) + minimumDifference_rad
            
        else
            
            % Do nothing
            
        end

    % Create Gabor textures

    [gaborid, gaborrect] = CreateProceduralGabor(ptbWindow,GaborWidth,GaborHeight);

        for thisTrial = 1:nTrialsPerStaircase

            % Determine which interval will have the different pair
            
            whichIsDifferent = whichInterval(thisTrial);
            
            % Display fixation cross
            
            

            % Draw first set of Gabors

            if whichIsDifferent = 1
            
                for thisGabor = 1:stimulus.nItems

                    Screen('DrawTexture', ptbWindow, gaborid, [], itemRects, gaborFirstOrientation_rad(thisGabor), [], [], modulateColor, [], kPsychDontDoRotation, [phase+180, thisGaborFrequency, sc, contrast, aspectratio, 0, 0, 0]);

                end

            
            % Display blank for interval

            % Draw second set of Gabors

            % Get keypress response 'first' or 'second'






