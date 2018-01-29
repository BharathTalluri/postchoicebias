%% Commitment RDK code
% example script to practice coarse motion discimination (which will be
% used for thresholding before they do any fine task)
% ----------------------------------------------------------------

close all; clear all; clc;
dbstop if error;
addpath('D:\USERS\AnneUrai\Commitment');
addpath('D:\USERS\AnneUrai\Commitment\stats');

%% SETUP

setup.Eye       = false;
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
[setup, dots, fix, results, sound] = configExperiment_example(setup, window, audio);

% INITIALIZE DOTS COORDINATES
stim = nan(setup.nblocks, setup.ntrials, setup.nframes, 2, dots.nDots);
for block = 1 : setup.startatblock : setup.nblocks,
    for trial = 1 : setup.ntrials,
        % preload all the dot coordinates before starting the trial
        stim(block, trial, :, :, :)      = dots_limitedlifetime(setup, window, dots, block, trial, '');
        Screen('DrawText',window.h, sprintf('Loading %.1f percent', round(sub2ind([setup.ntrials setup.nblocks ], trial, block))/setup.totalntrials*100), window.center(1)*0.60, window.center(2)*0.75 , [255 255 255] );
        Screen('Flip', window.h);
    end % end loop over trials
end % end loop over blocks

%% LOOP OVER BLOCKS
for block = setup.startatblock : setup.nblocks,
    
    % initialize EyeLink for this block
    if setup.Eye == true,
        edfFile = ELconfig(window, setup, block);
    end
    
    % TURN SCREEN BLACK & window INTRODUCTORY TEXT
    Screen('FillRect', window.h, window.black);
    Screen('DrawText',window.h, 'Are you ready?', window.center(1)*0.60, window.center(2)*0.75 , [255 255 255] );
    switch setup.counterbalancing
        case 1 %odd participants
            Screen('DrawText',window.h, 'Press Z for upwards and M for downwards motion', window.center(1)*0.60, window.center(2)*0.6 , [255 255 255] );
        case 0 %even participants
            Screen('DrawText',window.h, 'Press M for upwards and Z for downwards motion', window.center(1)*0.60, window.center(2)*0.6 , [255 255 255] );
    end
    Screen('Flip', window.h);
    
    WaitSecs(.1);
    KbWait(); % after windowing the introductory text, wait for a keypress

    %% LOOP OVER TRIALS
    for trial = 1 : setup.ntrials
        
        clear keyIsDown secs keyCode;

        WaitSecs(.5);
        KbWait(); % wait for a keypress: start the trial as soon as a key is pressed
        WaitSecs(.1);

        % variable fixation interval before stim onset
        for frameNum = 1:ceil(setup.fixtime(block, trial)*window.frameRate),
            
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.fix.VBL(block, trial, frameNum), ...
                flip.fix.StimOns(block, trial, frameNum), ...
                flip.fix.FlipTS(block, trial, frameNum), ...
                flip.fix.Missed(block, trial, frameNum), ...
                flip.fix.beampos(block, trial, frameNum)] = Screen('Flip', window.h, [], 1); % dontclear
        end
        
        if setup.Eye,
            Eyelink ('Message', sprintf('block%d_trial%d_fix', block, trial));
        else
            disp(sprintf('block%d_trial%d_fix', block, trial));
        end
        
        % play reference stimulus onset tone
        % PsychPortAudio('SetLoop',audio.h, sound.tonepos(3,1), sound.tonepos(3,2));
        % PsychPortAudio('Start', audio.h); %like flip
        
        % PRESENT STIMULUS
        TimingCnt = flip.fix.VBL(block, trial, frameNum) + window.frameDur - window.slack;
        
        % LOOP OVER FRAMES
        for frameNum = 1:setup.nframes,
            
            window = drawAllDots(window, dots, block, trial, stim, frameNum);
            
            %Draw the fixation and Flip
            window      = drawFixation(window, fix, dots); % fixation
            
            [flip.stim.VBL(block, trial, frameNum), ...
                flip.stim.StimOns(block, trial, frameNum), ...
                flip.stim.FlipTS(block, trial, frameNum), ...
                flip.stim.Missed(block, trial, frameNum), ...
                flip.stim.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
            TimingCnt = flip.stim.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            
            % send stim onset to EL
            if frameNum == 1,
                if setup.Eye == true ,
                    Eyelink ('Message', sprintf('block%d_trial%d_stim_dir%d', block, trial, dots.direction(block, trial)));
                else %otherwise, output in command line
                    disp(sprintf('block%d_trial%d_stim_dir%d', block, trial, dots.direction(block, trial)));
                end
            end
            
            % check for response given
            [keyIsDown, secs, keyCode] = KbCheck();
            
            if keyIsDown,
                break % exit the loop over frames if a response is made
            end
        end % end loop over frames
        
        switch KbName(keyCode)
            case 'z',
                switch setup.counterbalancing
                    case 1 %odd participants
                        results.response(1, trial) = 270;
                    case 0
                        results.response(1, trial) = 90;
                end
            case 'm',
                switch setup.counterbalancing,
                    case 1
                        results.response(1, trial) = 90;
                    case 0
                        results.response(1, trial) = 270;
                end
            case 'ESCAPE', % if escape is pressed, exit the experiment
                disp('ESCAPE key pressed')
                setup.cancel = true;
                results.response(block, trial) = NaN;
            case 'esc', % if escape is pressed, exit the experiment
                disp('esc key pressed')
                setup.cancel = true;
                results.response(block, trial) = NaN;
            otherwise % if any other key was pressed, fill in a NaN
                disp('other key pressed')
                results.response(block, trial) = NaN;
        end % end switch
        
        % calculate reaction time (even if no button was pressed)
        results.key{block, trial}   = KbName(keyCode);
        results.RT(block, trial)    = secs - flip.stim.VBL(block,trial, 1);
        
        %% CODE CORRECT RESPONSE
        
        if results.response(block, trial) == dots.direction(block,trial),
            results.correct(block,trial) = true;
        elseif isnan(results.response(block, trial)),
            results.correct(block,trial) = NaN;
            results.RT(block,trial) = NaN; %set RT to NaN to easily filter out trials without a response
        else
            results.correct(block,trial) = false; % wrong trial
        end
        
        % send response  to EL
        if setup.Eye == true,
            Eyelink ('Message', sprintf('block%d_trial%d_stimdir%d_respkey%d_correct%d_RT%.4f', ...
                block, trial, dots.direction(block, trial), results.response(block, trial), results.correct(block, trial), results.RT(block, trial)));
        else
            disp(sprintf('block%d_trial%d_resp_key%d_correct%d_RT%.4f', block, trial, results.response(block, trial), results.correct(block, trial), results.RT(block, trial)));
        end
        
       % present auditory feedback
        %correct
        if results.correct(block,trial) == true, %correct
            PsychPortAudio('SetLoop',audio.h, sound.tonepos(1,1), sound.tonepos(1,2));
            if setup.training, Screen('DrawText', window.h, 'correct', window.center(1), window.center(2)*.9 , [255 255 255] ); end
            
            %incorrect
        elseif results.correct(block,trial) == false,
            PsychPortAudio('SetLoop',audio.h, sound.tonepos(2,1), sound.tonepos(2,2));
            if setup.training, Screen('DrawText', window.h, 'wrong', window.center(1), window.center(2)*.9 , [255 255 255] ); end
            
            %no response given
        elseif isnan(results.correct(block,trial)),
            PsychPortAudio('SetLoop',audio.h, sound.tonepos(2,1), sound.tonepos(2,2));
            if setup.training, Screen('DrawText', window.h, 'no resp', window.center(1), window.center(2)*.9 , [255 255 255] ); end
            
            % if there was no deviation from the reference, select random feedback
        elseif sign(setup.direction.stimone(block,trial)-setup.refdirec(block, trial)) == 0 && ~isnan(results.correct(block,trial)),
            
            throwcoin = randi(2);
            switch throwcoin,
                case 1
                    PsychPortAudio('SetLoop',audio.h, sound.tonepos(1,1), sound.tonepos(1,2));
                    if setup.training, Screen('DrawText', window.h, 'correct', window.center(1), window.center(2)*.9 , [255 255 255] ); end
                case 2
                    PsychPortAudio('SetLoop',audio.h, sound.tonepos(2,1), sound.tonepos(2,2));
                    if setup.training, Screen('DrawText', window.h, 'wrong', window.center(1), window.center(2)*.9 , [255 255 255] ); end
            end
        else
            %unrecognized response
            setup.cancel = true;
            warning('could not determine which feedback to give');
        end
        results.feedbackonset(block, trial) = GetSecs;
        PsychPortAudio('Start', audio.h);
        Screen('Flip', window.h);

        
        %send trigger to EL
        if setup.Eye,
            Eyelink ('Message', sprintf('block%d_trial%d_feedback_correct%d', block, trial, results.correct(block, trial)));
        else
            disp(sprintf('block%d_trial%d_feedback_correct%d', block, trial, results.correct(block, trial)));
        end
       
        
        % break out of all trials if ESC was pressed
        if setup.cancel,
            break
            warning('experiment was manually terminated');
        end
        
    end % END LOOP OVER TRIALS
    
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
    
end % END LOOP OVER BLOCKS

nanmean(results.correct)
% exit gracefully
disp('done!');
Screen('CloseAll'); ShowCursor;
PsychPortAudio('Stop', audio.h);
sca;