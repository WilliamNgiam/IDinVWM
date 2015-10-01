% Individual differences in perceptual ability and visual working memory

% This study aims to examine whether an individual's ability to
% discriminate between different orientation and spatial frequency of Gabor
% patches predicts the individual's visual working memory (VWM) ability on 
% a change-detection task.

% This code is for the change-detection tasks in Experiment 1 of WN's PhD

% WN started writing this June 2015

% -------------------------------------------------------------------------

clear all;
Screen('CloseAll');

addpath(genpath('/Users/wngi5916/Documents/MATLAB/'));

AssertOpenGL;

% Set up experiment parameters

% Set up equipment parameters

equipment.viewDist = 500;                           % Viewing distance in mm
equipment.ppm = 2.7;                                % Pixels per mm (HP P1120, 1024 x 768, 120 Hz) % Remeasure
equipment.gammaVals = 1.0./2.6434 2.2312 2.171];    % Gamma values for CRT in GT519

% Set up colour parameters

colour.blackVal = 0;
colour.greyVal = 0.5;
colour.whiteVal = 1;

% Set up stimulus parameters

stimulus.size_dva = 2.5;            % Item size in degrees of visual angle
stimulus.eccentricity_dva = 4;      % Eccentricity in degrees of visual angle

stimulus.nItems = 4;                % Number of items in array for change-detection task

% Set up temporal parameters

stimulus.fixationDuration = 0.5;    % Duration of fixation point
stimulus.blankDuration = 1;         % Duration of inter-trial interval (blank)

stimulus.memoryDuration = 0.5;      % Duration of memory array in change-detection task
stimulus.maskDuration = 0.5;        % Duration of mask in change-detection task

% Set up staircase parameters

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

itemBaseTheta_rad = repmat(linspace(0,2*pi-(2*pi/stimulus.nArrayItems), stimulus.nArrayItems)',1,nTrialsPerSession);
arrayJitter_rad = ((2*pi)/stimulus.nArrayItems)*repmat(rand(1,nTrialsPerSession), stimulus.nArrayItems, 1);
allItemTheta = itembaseTheta_rad+arrayJitter_rad;

% Get participant ID

% ----------------------------------------------------------------------- %

%                           Begin Experiment

% ----------------------------------------------------------------------- %

% Create Gabor textures

%                Condition 1 - Orientation Change-Detection               %

% Instruction Text

% Wait for response






