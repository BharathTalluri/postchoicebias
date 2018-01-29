%% Commitment RDK code
% example script to practice coarse motion discimination (which will be
% used for thresholding before they do any fine task)
% ----------------------------------------------------------------

close all; clear all; clc;
dbstop if error;
addpath('D:\USERS\AnneUrai\Commitment');
addpath('D:\USERS\AnneUrai\Commitment\stats');

%% SETUP

setup.training  = true;
% ask for subject number, default: 0
setup.participant       = input('Participant number? ');
if isempty(setup.participant)
    setup.participant   = 20; %test
end

% INITIALIZE PSYCHTOOLBOX
window.dist            = 50;
window.width           = 40;
window.skipChecks      = 1;
[window, audio]        = SetupPTB(window); %load all the psychtoolbox things

% present loading text on the screen
Screen('TextSize', window.h, 15);
Screen('DrawText', window.h, 'Loading...', window.center(1)*0.60, window.center(2)*0.75 , [255 255 255] );
Screen('Flip', window.h);

% INITIALIZE EXPERIMENTAL PARAMETERS
% make sure that you load the correct config whether you're in the
% threshold definition or main experiment...

setup.cancel = false; %this turns true if escape is pressed

% example - easy
setup.cohlevel          = .4;
setup.ntrials           = 37; % number of trials per coherence level for this whole session

% timing
setup.fixtime           = (.5 + .5*rand(1, setup.ntrials)); %generate the fixation time 0.6 - 0.8s (36-48 frames)
setup.viewingtime       = .75; % maximum viewing duration in seconds
setup.nframes           = ceil(setup.viewingtime * window.frameRate); %number of frames the stimulus is windowed

%% (output) dots: set stimuli (dots) parameters

[dots, fix]                 = setupDots(window, setup);
switch mod(setup.participant,2), % counterbalance for every two subjects
    case 0
       dots.direction = 0:5:180;
    case 1
       dots.direction = 180:5:360;
end

dots.direction = [0 180 360 0 180 360];
% randomize instead
dots.direction = Shuffle(dots.direction);

% TURN SCREEN BLACK & window INTRODUCTORY TEXT
Screen('FillRect', window.h, window.black);
Screen('DrawText',window.h, 'Are you ready?', window.center(1)*0.60, window.center(2)*0.75 , [255 255 255] );
Screen('Flip', window.h);

WaitSecs(.1);
KbWait(); % after windowing the introductory text, wait for a keypress

%% LOOP OVER TRIALS
block = 1;
for trial = 1 : setup.ntrials, % only loop over trials
    
    clear keyIsDown secs keyCode;
    WaitSecs(.5);
    
    % variable fixation interval before stim onset
    window      = drawFixation(window, fix, dots); % fixation
    Screen('Flip', window.h); % dontclear
    
    % preload all the dot coordinates before starting the trial
    stim(:, :, :)      = dots_limitedlifetime(setup, window, dots, block, trial, '');
    
    % PRESENT STIMULUS
    TimingCnt = GetSecs + window.frameDur - window.slack;
    
    % LOOP OVER FRAMES
    for frameNum = 1:setup.nframes,
        
        % draws the dots on the screen for this frame
        Screen('DrawDots',window.h,squeeze(stim(frameNum, :, :)), ...
            dots.size, dots.color, window.center, 2);
        
        %Draw the fixation and Flip
        window      = drawFixation(window, fix, dots); % fixation
        
        [flip.stim.VBL(block, trial, frameNum), ...
            flip.stim.StimOns(block, trial, frameNum), ...
            flip.stim.FlipTS(block, trial, frameNum), ...
            flip.stim.Missed(block, trial, frameNum), ...
            flip.stim.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
        TimingCnt = flip.stim.VBL(block, trial, frameNum) + window.frameDur - window.slack;
    end % end loop over frames
        
    % draws the dots on the screen for this frame
    Screen('DrawDots',window.h,squeeze(stim(frameNum, :, :)), ...
        dots.size, dots.color, window.center, 2);
    window      = drawFixation(window, fix, dots); % fixation
    
    % calculate the coordinates of the reference mark line
    xstart  = (dots.radius-(deg2pix(window, 3))) .*cos(dots.direction(block, trial)*pi/180);
    xend    = (dots.radius+(deg2pix(window, 1))) .*cos(dots.direction(block, trial)*pi/180);
    
    ystart  = (dots.radius-(deg2pix(window, 3))) .*sin(dots.direction(block, trial)*pi/180);
    yend    = (dots.radius+(deg2pix(window, 1))) .*sin(dots.direction(block, trial)*pi/180);
    
    % also add a dial in the reference direction for this subject
    Screen('DrawLines', window.h, [xstart xend; ystart yend], 2, dots.color, window.center);
    Screen('Flip', window.h); % dontclear
    
    [secs, keyCode, deltaSecs] = KbWait();
    try
        switch KbName(keyCode),
            case 'ESCAPE',
                setup.cancel = true;
            case 'esc'
                setup.cancel = true;
            otherwise
                disp('did not recognize key pressed');
        end
    catch
        disp('did not recognize key pressed');
    end
    
    % break out of all trials if ESC was pressed
    if setup.cancel,
        break
        warning('experiment was manually terminated');
    end
    
end % END LOOP OVER TRIALS

% exit gracefully
disp('done!');
Screen('CloseAll'); ShowCursor;
PsychPortAudio('Stop', audio.h);
sca;