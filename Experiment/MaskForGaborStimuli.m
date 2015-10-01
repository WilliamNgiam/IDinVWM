% Set up equipment parameters
close all;

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

stimulus.size_dva = 4;              % Item size in degrees of visual angle
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

stimulus.truncate =         4;          % Standard deviations
% stimulus.truncate_p =       ceil(stimulus.truncate*stimulus.SD_p);
% Calculate equipment parameters

equipment.mpd = (equipment.viewDist)*tan(deg2rad(2*stimulus.eccentricity_dva))/stimulus.eccentricity_dva; % Calculate mm per degree of visual angle to the ecccentricity of the stimuli
equipment.ppd = equipment.ppm*equipment.mpd;

% Calculate spatial parameters

stimulus.size_pix = round(stimulus.size_dva*equipment.ppd);                     % Item size in pixels
stimulus.eccentricity_pix = round(stimulus.eccentricity_dva*equipment.ppd);     % Eccentricity of stimulus in pixels
stimulus.fixationSize_pix = stimulus.fixationSize_dva*equipment.ppd;            % Fixation cross size in pixels

stimulus.frequency_pix = stimulus.frequency_cpd/equipment.ppd;                  % Spatial frequency of Gabor patches in cycles per pixel
stimulus.sc_pix = round(stimulus.sc*equipment.ppd);                             % Spatial constant of Gabor patches in pixels

% Set up Psychtoolbox

screenID = max(Screen('Screens'));
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'EnablePseudoGrayOutput');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
Screen('Preference','SkipSyncTests',2);
[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, colour.greyVal);
% PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
[screenWidth, screenHeight] = RectSize(winRect);
screenCentreX = round(screenWidth/2);
screenCentreY = round(screenHeight/2);
flipInterval = Screen('GetFlipInterval', ptbWindow);
% stimulus.size_dva = 10; %2.5            % Item size in degrees of visual angle
% stimulus.size_pix = round(stimulus.size_dva*equipment.ppd);                     % Item size in pixels

% Create mask textures

% Create a horizontal or vertical Gabor
[gaborid, gaborrect] = CreateProceduralGabor(ptbWindow,stimulus.size_pix,stimulus.size_pix, 0, [0.5 0.5 0.5 0], 1, .5);
Screen('DrawTexture', ptbWindow, gaborid, [], gaborrect, 0, [], [], [], [], kPsychDontDoRotation, [180, stimulus.frequency_pix, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);
Screen('Flip',ptbWindow);
imageArray = Screen('GetImage', ptbWindow, gaborrect); % This is uint8
imageArray = double(imageArray)/255;
imageArray = imageArray(:,:,1); % Use one layer as grayscale

% Normalise
imageArray = (imageArray-min(imageArray(:)))/(max(imageArray(:))-min(imageArray(:)));
imageArray = (imageArray*2)-1; % Balance

% Fourier transformation
SFT_Size = 2^12;
fftGaborImage = fftshift(fft2(imageArray, SFT_Size, SFT_Size));

% Amplitude spectrum
fr = abs(fftGaborImage);
% Phase spectrum
ph = angle(fftGaborImage);
sca;

spatialFT =                 abs(fftshift(fft2(imageArray,SFT_Size,SFT_Size))).^2;

spatialFigure =             figure('Color', 'white');
xlabel('x (degrees)');
ylabel('y (degrees)');
hold on;
imageArray =              (imageArray-min(min(imageArray)))./(max(max(imageArray))-min(min(imageArray)));
imagesc([-8 8],[-8 8],imageArray);
colormap(gca,gray(2^12));
axis square;
axis image;
set(gca,'tickdir','out','xtick',-8:4:8,'ytick',-8:4:8,'xminortick','on','yminortick','on');
box on;
sFigCB = colorbar;

nyqF = equipment.ppd/2;

sftFigure =                 figure('Color', 'white');
xlabel('x (cpd)');
ylabel('y (cpd)');
hold on;
spatialFT =              (spatialFT-min(min(spatialFT)))./(max(max(spatialFT))-min(min(spatialFT)));
imagesc([-nyqF nyqF],[-nyqF nyqF],spatialFT);
colormap(gca,jet(2^12));
axis square;
axis image;
axis([-10 10 -10 10]);
set(gca,'tickdir','out','xtick',-10:1:10,'ytick',-10:1:10,'xminortick','on','yminortick','on');
box on;
sFTCB = colorbar;

% Isolate the DC (0 cpd) line across the relevant dimension
figure('Color', 'white');
xlabel('radians');
frequencyAxis = linspace(-nyqF,nyqF,SFT_Size);
plot(frequencyAxis,fr(SFT_Size/2,:));

% Fit a Gaussian to that 2-D function, note the standard deviation
startBit = ceil(SFT_Size/2);
endBit = SFT_Size;


f = fit(frequencyAxis(startBit:endBit).',fr(2048,(startBit:endBit)).','gauss1');
plot(f,frequencyAxis(startBit:endBit).',fr(2048,(startBit:endBit)).');

% Get standard deviation

sd = f.c1;

% Create an image of a Gaussian annulus the same way you'd make a gaussian
% blob, but using radial coordinates.

imSize = stimulus.size_pix;                           % image size: n X n
sigma = sd;

% make linear ramp
X = linspace(-nyqF,nyqF,SFT_Size);                           % X is a vector from 1 to imageSize
%X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5
[Xm, Ym] = meshgrid(X,X);
[th, r] = cart2pol(Xm,Ym);
%s = sigma / imSize;  
figure;
%gauss = exp( -(((R.^2)+(R.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
gauss = exp(-((r-stimulus.frequency_cpd)./(sd/sqrt(2))).^2);


imagesc([-nyqF nyqF],[-nyqF nyqF], gauss);                        % display
colormap(gca,jet(2^12));
axis square;
axis image;
axis([-10 10 -10 10]);
set(gca,'tickdir','out','xtick',-10:1:10,'ytick',-10:1:10,'xminortick','on','yminortick','on');
box on;
sFTCB = colorbar;
figure;

% Confirm that annulus is matched in SF spectrum to the Gabor
imagesc([-nyqF nyqF],[-nyqF nyqF], gauss-spatialFT);                        % display
colormap(gca,jet(2^12));
axis square;
axis image;
axis([-10 10 -10 10]);
set(gca,'tickdir','out','xtick',-10:1:10,'ytick',-10:1:10,'xminortick','on','yminortick','on');
box on;
sFTCB = colorbar;
figure;


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

% Multiply the amplitude spectrum by the Gaussian annulus.

maskwhatever = fftshift(gauss);

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

% When you draw this, draw it masked by a procedural gabor with SF = 0, to
% give you a Gaussian envelope.
imagesc([-8 8],[-8 8],maskednoise);
colormap(gca,gray(2^12));
axis square;
axis image;
set(gca,'tickdir','out','xtick',-8:4:8,'ytick',-8:4:8,'xminortick','on','yminortick','on');
box on;
sFigCB = colorbar;

% Now check the amplitude spectrum of the resulting image
figure;
fftOfnoiseImage = fft2(maskednoise,SFT_Size,SFT_Size);
frOfnoiseImage = abs(fftshift(fftOfnoiseImage));
imagesc([-nyqF nyqF],[-nyqF nyqF], frOfnoiseImage);                        % display
colormap(gca,jet(2^12));
axis square;
axis image;
axis([-10 10 -10 10]);
set(gca,'tickdir','out','xtick',-10:1:10,'ytick',-10:1:10,'xminortick','on','yminortick','on');
box on;
sFTCB = colorbar;
% 
% Set up Psychtoolbox

screenID = max(Screen('Screens'));
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'EnablePseudoGrayOutput');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
Screen('Preference','SkipSyncTests',2);
[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, colour.greyVal);
%PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
[screenWidth, screenHeight] = RectSize(winRect);
screenCentreX = round(screenWidth/2);
screenCentreY = round(screenHeight/2);
flipInterval = Screen('GetFlipInterval', ptbWindow);
% stimulus.size_dva = 10; %2.5            % Item size in degrees of visual angle
% stimulus.size_pix = round(stimulus.size_dva*equipment.ppd);                     % Item size in pixels

% Create mask textures

% Create a gaussian blob
[gaborid, gaborrect] = CreateProceduralGabor(ptbWindow,stimulus.size_pix,stimulus.size_pix, 0, [0.5 0.5 0.5 0], 1, .5);
Screen('DrawTexture', ptbWindow, gaborid, [], gaborrect, 0, [], [], [], [], kPsychDontDoRotation, [90, 0, stimulus.sc_pix, stimulus.contrast, stimulus.aspectratio, 0, 0, 0]);
Screen('Flip',ptbWindow);
blobArray = Screen('GetImage', ptbWindow, gaborrect); % This is uint8
blobArray = double(blobArray)/255;
blobArray = blobArray(:,:,1); % Use one layer as grayscale

% Normalise
blobArray = (blobArray-min(blobArray(:)))/(max(blobArray(:))-min(blobArray(:)));

% Normalise
maskednoise = (maskednoise-min(maskednoise(:))) / (max(maskednoise(:)) - min(maskednoise(:)));
maskednoise = (2*maskednoise)-1;
mask = blobArray.*maskednoise;
mask = (mask-min(mask(:)))/(max(mask(:))-min(mask(:)));

% mask = (mask+1)/2; % Normalise to a range of 0.5 to 1
rgbImage = repmat(uint8(255.*mask),[1 1 3]);
maskTexture = Screen('MakeTexture',ptbWindow,rgbImage);
Screen('DrawTexture',ptbWindow,maskTexture);
Screen('Flip',ptbWindow);