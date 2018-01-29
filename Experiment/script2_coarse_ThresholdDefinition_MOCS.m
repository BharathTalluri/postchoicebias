%% Code to determine a subject's threshold before starting the real training
% Anne Urai, March 2014
% -----------------------------------------------------------------

clear all; close all; clc;
addpath('D:\USERS\AnneUrai\Commitment');
addpath('D:\USERS\AnneUrai\Commitment\stats');

% general setup
% ask for subject number, default: 0
setup.participant       = input('Participant number? ');
if isempty(setup.participant)
    setup.participant   = 20; %test
end

window.dist            = 50;
window.width           = 40;
window.skipChecks      = 1;
[window, audio]        = SetupPTB(window); %load all the psychtoolbox things

Screen('TextSize', window.h, 15);
Screen('DrawText',window.h, 'Loading...', window.center(1)*0.60, window.center(2)*0.75 , [255 255 255] );
Screen('Flip', window.h);

%% Initialize setup
% the function configuration contains all the important experimental
% design parameters

[setup, dots, fix, results, flip] = configExperiment_threshdef(window, setup);

% preload all the dot coordinates
stim = nan(setup.nblocks, setup.ntrials, setup.nframes, 2, dots.nDots);
% only one
for block = 1:setup.nblocks,
    for trial = 1:setup.ntrials,
        % preload all the dot coordinates before starting the trial
        stim(block, trial, :, :, :)      = dots_limitedlifetime(setup, window, dots, block, trial, '');
        Screen('DrawText',window.h, sprintf('Loading %.1f percent', round(sub2ind([setup.ntrials setup.nblocks ], trial, block)/setup.totalntrials*100)), window.center(1)*0.60, window.center(2)*0.75 , [255 255 255] );
        Screen('Flip', window.h);
    end
end

for block = 1:setup.nblocks,
    
    Screen('FillRect', window.h, window.black);
    Screen('DrawText',window.h, 'Are you ready?', window.center(1)*0.60, window.center(2)*0.75 , [255 255 255] );
    switch setup.counterbalancing
        case 1 %odd participants
            Screen('DrawText',window.h, 'Press Z for upwards', window.center(1)*0.60, window.center(2)*0.6 , [255 255 255] );
            Screen('DrawText',window.h, 'and M for downwards motion.', window.center(1)*0.60, window.center(2)*0.7 , [255 255 255] );
        case 0 %even participants
            Screen('DrawText',window.h, 'Press M for upwards', window.center(1)*0.60, window.center(2)*0.6 , [255 255 255] );
            Screen('DrawText',window.h, 'and Z for downwards motion.', window.center(1)*0.60, window.center(2)*0.7 , [255 255 255] );
    end
    Screen('Flip', window.h);
    
    WaitSecs(.1);
    KbWait(); %after windowing the introductory text, wait for a keypress
    
    %% start the loop over trials
    for trial = 1:setup.ntrials,
        
        % ISI of 1 second;
        Screen('Flip', window.h);
        WaitSecs(1);
        
        % variable fixation before stim onset
        for frameNum = 1:ceil(setup.fixtime(block, trial)*window.frameRate),
            %Draw the fixation and Flip
            window      = drawFixation(window, fix, dots); % fixation
            Screen('Flip', window.h); % don't clear so the same stimuli are drawn
        end
        
        %% Loop through the stimulus
        TimingCnt = GetSecs + window.frameDur - window.slack;
        
        for frameNum = 1:setup.nframes,
            %Draw all dots at once
            window  = drawAllDots(window, dots, block, trial, stim, frameNum);
            
            %Draw the fixation and Flip
            window      = drawFixation(window, fix, dots); % fixation
            Screen('DrawingFinished', window.h);
            
            [flip.stim.VBL(block, trial, frameNum), ...
                flip.stim.StimOns(block, trial, frameNum), ...
                flip.stim.FlipTS(block, trial, frameNum), ...
                flip.stim.Missed(block, trial, frameNum), ...
                flip.stim.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
            TimingCnt = flip.stim.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            if frameNum == 1,
                fprintf('\n block%d_trial%d_stim%d_dir%d \n', block, trial, dots.direction(block, trial));
            end
        end
        
        %% collect response
        while GetSecs-flip.stim.VBL(block,trial,frameNum) < setup.resptime && isnan(results.response(block, trial)),
            % when no response has been given, and the maximum
            % response time hasnt been reached
            [keyIsDown, secs, keyCode]  = KbCheck();
            
            if keyIsDown,
                try
                    switch KbName(keyCode)
                        case 'z', % left target 1, right target 2
                            switch setup.counterbalancing
                                case 1 %odd participants
                                    results.response(block, trial) = 270;
                                case 0
                                    results.response(block, trial) = 90;
                            end
                        case 'm',
                            switch setup.counterbalancing,
                                case 1
                                    results.response(block, trial) = 90;
                                case 0
                                    results.response(block, trial) = 270;
                            end
                        case 'ESCAPE', % if escape is pressed, exit the experiment
                            setup.cancel = true;
                            results.response(block, trial) = NaN;
                        case 'esc', % if escape is pressed, exit the experiment
                            setup.cancel = true;
                            results.response(block, trial) = NaN;
                        otherwise % if any other key was pressed, fill in a NaN
                            results.response(block, trial) = NaN;
                    end
                catch % in case it doesnt understand the input. for example when two keys are pressed at the same time
                    results.response(block, trial) = NaN;
                end
            else
                results.response(block, trial) = NaN; %if no key was pressed, NaN
            end
        end %button pressed
        
        % calculate reaction time (even if no button was pressed)
        results.key{block, trial}   = KbName(keyCode);
        results.RT(block, trial)    = secs - flip.stim.VBL(block,trial,frameNum); %frameNum is the last frame that got presented
        
        %% code for correct responses
        % this bit of the code needs to be edited so that the
        % correct flag is applied when response matches stimulus
        % (eg. up down)
        if results.response(block, trial) == dots.direction(block,trial), %whether motion is stronger than 50% or not
            results.correct(block,trial) = true;
        elseif isnan(results.response(block, trial)),
            results.correct(block,trial) = NaN;
            results.RT(block,trial) = NaN; %set RT to NaN to easily filter out trials without a response
        else results.correct(block,trial) = false;
        end
        
        fprintf('\n block%d_trial%d_resp%d_correct%d_RT%.3f \n', block, trial, results.response(block, trial), results.correct(block, trial), results.RT(block, trial));
        
        % break out of all trials if ESC was pressed
        if setup.cancel,
            break
            warning('experiment was manually terminated');
        end
        
    end %end trial loop
    
    % break out of all trials if ESC was pressed
    if setup.cancel,
        break
        warning('experiment was manually terminated');
    end
    
    %% break text
    Screen('DrawText',window.h, sprintf('Finished block %d.', block) , ...
        window.center(1)*0.60, window.center(2)*0.40 , [255 255 255] );
    
    if block < setup.nblocks,
        Screen('DrawText',window.h, 'Take a break!', window.center(1)*0.6, window.center(2), [255 255 255] );
    else % finish
        Screen('DrawText',window.h, 'Done! Please call the experimenter.', window.center(1)*0.60, window.center(2) , [255 255 255] );
    end
    Screen('Flip', window.h);
    WaitSecs(1);
    KbWait();
    
    % present some useful into to the experimenter
    fprintf('\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n');
    fprintf('\n Finished block %d out of %d! \n', block, setup.nblocks)
    % compute the overall binary accuracy
    fprintf('\n Binary choice accuracy = %.2f percent \n \n', 100*nanmean(results.correct(block, :)));
    fprintf('\n Binary choice RT = %.2f s \n \n', nanmean(results.RT(block, :)));
    fprintf('\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n');
    
    % break out of all 1s if ESC was pressed
    if setup.cancel == true,
        warning('experiment was manually terminated');
    end
end

% wrap up and save
setup.datetime = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
% create subject specific file and save - add unique datestring to avoid any overwriting of files
setup.filename = sprintf('Data/P%d_thresholdMOCS_%s.mat', setup.participant, setup.datetime);
save(setup.filename, '-mat', 'setup', 'window', 'stim', 'dots', 'fix', 'results', 'flip');
% !! important, make sure all the variables of interest are saved to
% possibly recreate the exact experiment later.

% exit gracefully
%KbWait();
disp('done!');
Screen('CloseAll'); ShowCursor;
PsychPortAudio('Stop', audio.h);
sca;

% take all results and plot the full psychometric extravaganza
% make bootstrapping work!s
[datapoints1, datapoints2, datapoints3, fit1, fit2, fit3] = FitThreshold(setup.filename, 'fit', 'bootstrap');

save(sprintf('Data/P%d_threshold.mat', setup.participant), 'fit1');