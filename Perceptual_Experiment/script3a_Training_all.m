%% Commitment RDK code
% this models the methodology of Bronfman, Brezis & Usher with random dot
% stimuli. To be used by Ana Vojvodic.

% Anne Urai, Jan 2015

% this experiment does the training in stages:
% 1. one interval, choice-only (with feedback)
% 2. one interval, choice and no-choice (with feedback, correct feedback after valid nochoice resp)
% 3. two intervals, with estimation (choice0 and choice1)
% 4. two intervals, choice-1 choice1 choice0
% ----------------------------------------------------------------

close all; clear all; clc;
dbstop if error;
%addpath('D:\USERS\AnneUrai\Commitment');
%addpath('D:\USERS\AnneUrai\Commitment\stats');

% general setup
setup.Eye           = false; % true if using Eyelink
setup.cancel        = false; % becomes true if escape is pressed, will abort experiment (but save data)
setup.training      = true;

% ask for subject number, default: 0
setup.participant       = input('Participant number? ');
if isempty(setup.participant)
    setup.participant   = 20; %test
end

% ask which training we want to do
whichtraining           = input('Which training? 1-4 ');

% Setup the PsychToolbox
window.dist             = 40; % viewing distance in cm ( this is fixed in the chinrest, D1.33)
window.width            = 30; % physical width of the screen in cm, 40 cm in D1.33
window.skipChecks       = 1; % set to 1 to skip VBL tests and avoid warnings
[window, audio]         = SetupPTB(window); %load all the psychtoolbox things

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% INITIALIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Screen('TextSize', window.h, 15);
Screen('DrawText',window.h, 'Loading...', window.center(1)*0.60, window.center(2)*0.75 , [255 255 255] );
Screen('Flip', window.h);

% the function configuration contains all the important experimental design parameters
[setup, dots, fix, results, sound, stim] = configuration_training(window, audio, setup);

% training - specifics
switch whichtraining
    case 1
        % for this binary training, give feedback on every trial (always choice only);
        setup.choice = -1*ones(setup.nblocks, setup.ntrials);
    case 2
        setup.choice(find(setup.choice==1)) = -1;
    case 3,
        % only choice1 and choice0, no choice-1 trials anymore
        setup.choice = rand(setup.nblocks, setup.ntrials) > .5;
       % setup.choice(find(setup.choice==1))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP OVER BLOCKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while setup.cancel == false,
    
    for block = 1:setup.nblocks,

        % turn the screen black after eyelink setup
        Screen('FillRect', window.h, window.black);
        Screen('DrawText',window.h, 'Are you ready?', window.center(1)*0.60, window.center(2)*0.75 , [255 255 255] );
        Screen('DrawText',window.h, 'Use the mouse buttons to indicate your choice', window.center(1)*0.60, window.center(2)*0.6 , [255 255 255] );
        Screen('Flip', window.h);
        
        WaitSecs(.1);
        MouseWait; %after displaying the introductory text, wait for a mouse click
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LOOP OVER TRIALS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for trial = 1:setup.ntrials,
            
            Screen('FillRect', window.h, window.black);
            Screen('Flip', window.h); % flip once, so stationary dots
            WaitSecs(.1);
            MouseWait; %wait for a keypress
            
            % variable fixation interval
            for frameNum = 1:ceil(setup.fixtime(block, trial)*window.frameRate),
                window      = drawFixation(window, fix, dots); % fixation
                
                [flip.fix.VBL(block, trial, frameNum), ...
                    flip.fix.StimOns(block, trial, frameNum), ...
                    flip.fix.FlipTS(block, trial, frameNum), ...
                    flip.fix.Missed(block, trial, frameNum), ...
                    flip.fix.beampos(block, trial, frameNum)] = Screen('Flip', window.h);
            end
            
            if setup.Eye == true,
                Eyelink ('Message', sprintf('block%d_trial%d_fix', block, trial));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % FIRST INTERVAL
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            TimingCnt = flip.fix.VBL(block, trial, frameNum) + window.frameDur - window.slack;
            
            for frameNum = 1:setup.nframes,
                window          = drawAllDots_refmark(setup, window, dots, block, trial, stim.one, frameNum);
                window          = drawFixation(window, fix, dots); % fixation
                Screen('DrawingFinished', window.h);
                
                [flip.int1.VBL(block, trial, frameNum), ...
                    flip.int1.StimOns(block, trial, frameNum), ...
                    flip.int1.FlipTS(block, trial, frameNum), ...
                    flip.int1.Missed(block, trial, frameNum), ...
                    flip.int1.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
                
                if setup.Eye == true && frameNum == 1,
                    Eyelink ('Message', sprintf('block%d_trial%d_int1_ref%d_dir%d', block, trial, setup.refdirec(block, trial), setup.direction.stimone(block, trial)));
                elseif setup.Eye == false && frameNum == 1,
                    disp(sprintf('block%d_trial%d_int1_ref%d_dir%d', block, trial, setup.refdirec(block, trial), setup.direction.stimone(block, trial)));
                end
                TimingCnt = flip.int1.VBL(block, trial, frameNum) + window.frameDur - window.slack;
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BINARY CHOICE OR NO CHOICE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            window      = drawAllDots_refmark(setup, window, dots, block, trial, stim.one, frameNum);
            window      = drawFixation(window, fix, dots); % fixation
            
            % present the sound and visual cue for prompt
            switch abs(setup.choice(block, trial)),
                case 1 % make a binary choice
                    PsychPortAudio('SetLoop',audio.h, sound.tonepos(4,1), sound.tonepos(4,2));
                    if whichtraining < 3, Screen('DrawText', window.h, 'choice', window.center(1)*0.9, window.center(2)*.9 , [255 255 255] ); end
                case 0 % no choice
                    PsychPortAudio('SetLoop',audio.h, sound.tonepos(5,1), sound.tonepos(5,2));
                    if whichtraining < 3, Screen('DrawText', window.h, 'wheel', window.center(1)*0.9, window.center(2)*.9 , [255 255 255] ); end
            end
            
            responset = Screen('Flip', window.h);
            PsychPortAudio('Start', audio.h); %play the tone
            
            % now collect the right type of response
            while GetSecs-responset < setup.choicetime,
                
                % always wait for a specific time, so all intervals have the same length
                [x,y,clicks] = GetMouse;
                
                if any(clicks) && isnan(results.binary.response(block, trial)); % if a click was made and it's the first one on this trial
                    secs = GetSecs; % reaction time
                    results.binary.clicks{block, trial} = clicks; %save the actual click
                    
                    % counterbalance
                    switch mod(setup.participant,2),
                        case 0 % 0-180deg, lower half
                            if clicks(1) == 1,
                                results.binary.response(block, trial) = 1; % clockwise
                            elseif clicks(3) == 1,
                                results.binary.response(block, trial) = -1; % counterclockwise
                            elseif clicks(2) == 1,
                                results.binary.response(block, trial) = 0; % wheel
                            else % if any other key was pressed, fill in a NaN
                                results.binary.response(block, trial) = NaN;
                            end
                        case 1 % 180-360deg, upper half
                            if clicks(1) == 1,
                                results.binary.response(block, trial) = -1; % counterclockwise
                            elseif clicks(3) == 1,
                                results.binary.response(block, trial) = 1; % clockwise
                            elseif clicks(2) == 1,
                                results.binary.response(block, trial) = 0; % wheel
                            else % if any other key was pressed, fill in a NaN
                                results.binary.response(block, trial) = NaN;
                            end
                    end
                    % calculate RT (only if a response was given)
                    results.binary.RT(block, trial)     = secs - flip.int1.VBL(block,trial, frameNum);
                end
                
                % exit the whole experiment if escape is pressed
                [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
                if keyIsDown,
                    try
                        switch KbName(keyCode),
                            case 'ESCAPE',
                                setup.cancel = true;
                            case 'esc'
                                setup.cancel = true;
                        end
                    catch
                        warning('could not understand keyboard input');
                    end
                end
            end % interval
            
            % CALCULATE CORRECTNESS OF THE RESPONSE
            if abs(setup.choice(block, trial)) == 1, % choice and choice-only trials
                if results.binary.response(block, trial) == sign(setup.direction.stimone(block,trial)-setup.refdirec(block, trial)),
                    % correct answer
                    results.binary.correct(block,trial) = 1;
                elseif results.binary.response(block, trial) == -sign(setup.direction.stimone(block,trial)-setup.refdirec(block, trial)),
                    % incorrect answer
                    results.binary.correct(block,trial) = 0;
                elseif sign(setup.direction.stimone(block,trial)-setup.refdirec(block, trial)) == 0 && abs(results.binary.response(block, trial))==1,
                    % first interval had zero directional evidence
                    results.binary.correct(block,trial) = 0;
                elseif results.binary.response(block, trial) == 0, % wheel pressed
                    results.binary.correct(block,trial) = -1;
                else % unrecognized or no resp
                    results.binary.correct(block,trial) = -1;
                end
                
            elseif abs(setup.choice(block, trial)) == 0, % no-choice trials, make sure the wheel was pressed
                if results.binary.response(block, trial) == 0,
                    results.binary.correct(block, trial) = 1;
                else % if anything other than the wheel was pressed
                    results.binary.correct(block,trial) = -1;
                end
            end
            
            % send response to EL
            if setup.Eye == true,
                Eyelink ('Message', sprintf('block%d_trial%d_choice%d_ref%d_dir%d_resp%d_correct%d_RT%.3f', block, trial, setup.choice(block, trial), setup.refdirec(block, trial), setup.direction.stimone(block, trial), results.binary.response(block, trial), results.binary.correct(block, trial), results.binary.RT(block, trial)));
                fprintf('\n block%d_trial%d_choice%d_ref%d_dir%d_resp%d_correct%d_RT%.3f \n', block, trial, setup.choice(block, trial), setup.refdirec(block, trial), setup.direction.stimone(block, trial), results.binary.response(block, trial), results.binary.correct(block, trial), results.binary.RT(block, trial));
            end
            
            if results.binary.correct(block,trial) == -1;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % invalid feedback
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                PsychPortAudio('SetLoop',audio.h, sound.tonepos(3,1), sound.tonepos(3,2));
                PsychPortAudio('Start', audio.h);
                
                Screen('DrawText', window.h, 'INVALID', window.center(1), window.center(2)*.9 , [255 255 255] );
                Screen('Flip', window.h);
                disp('INVALID RESPONSE');
                
            elseif setup.choice(block, trial) == -1, % choice only
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % normal feedback
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % present auditory feedback
                if results.binary.correct(block,trial) == 1, % correct
                    PsychPortAudio('SetLoop',audio.h, sound.tonepos(1,1), sound.tonepos(1,2));
                    if whichtraining < 4, Screen('DrawText', window.h, 'correct', window.center(1), window.center(2)*.9 , [255 255 255] ); end
                    
                elseif results.binary.correct(block,trial) == 0,    % incorrect
                    PsychPortAudio('SetLoop',audio.h, sound.tonepos(2,1), sound.tonepos(2,2));
                    if whichtraining < 4, Screen('DrawText', window.h, 'error', window.center(1), window.center(2)*.9 , [255 255 255] ); end
                    
                elseif sign(setup.direction.stimone(block,trial)-setup.refdirec(block, trial)) == 0 && abs(results.binary.response(block, trial))== 1,
                    % if there was no deviation from the reference, select random feedback
                    throwcoin = randi(2);
                    switch throwcoin,
                        case 1
                            PsychPortAudio('SetLoop',audio.h, sound.tonepos(1,1), sound.tonepos(1,2));
                            if whichtraining < 4, Screen('DrawText', window.h, 'correct', window.center(1), window.center(2)*.9 , [255 255 255] ); end
                        case 2
                            PsychPortAudio('SetLoop',audio.h, sound.tonepos(2,1), sound.tonepos(2,2));
                            if whichtraining < 4, Screen('DrawText', window.h, 'error', window.center(1), window.center(2)*.9 , [255 255 255] ); end
                    end
                else
                    % unrecognized response
                    setup.cancel = true;
                    warning('could not determine which feedback to give');
                end
                results.feedbackonset(block, trial) = GetSecs;
                PsychPortAudio('Start', audio.h);
                Screen('Flip', window.h);
                
                %send trigger to EL
                if setup.Eye,
                    Eyelink ('Message', sprintf('block%d_trial%d_feedback_correct%d', block, trial, results.binary.correct(block, trial)));
                    fprintf('\n block%d_trial%d_feedback_correct%d \n', block, trial, results.binary.correct(block, trial));
                end
                
                WaitSecs(.5); % leave visual feedback on the screen for a bit, and wait a bit after feedback
                
                
            elseif setup.choice(block, trial)>-1, % if this trial has a second interval
                
                if whichtraining > 2,
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % SECOND MOTION STIMULUS
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    TimingCnt = GetSecs +  window.frameDur - window.slack;
                    
                    for frameNum = 1:setup.nframes,
                        
                        window      = drawAllDots_refmark(setup, window, dots, block, trial, stim.two, frameNum);
                        window      = drawFixation(window, fix, dots); % fixation
                        Screen('DrawingFinished', window.h);
                        
                        [flip.int2.VBL(block, trial, frameNum), ...
                            flip.int2.StimOns(block, trial, frameNum), ...
                            flip.int2.FlipTS(block, trial, frameNum), ...
                            flip.int2.Missed(block, trial, frameNum), ...
                            flip.int2.beampos(block, trial, frameNum)] = Screen('Flip', window.h, TimingCnt);
                        
                        % send stim onset to EL
                        if setup.Eye == true && frameNum == 1,
                            Eyelink ('Message', sprintf('block%d_trial%d_int2_ref%d_dir%d', block, trial, setup.refdirec(block, trial), setup.direction.stimtwo(block, trial)));
                            fprintf('\n block%d_trial%d_int2_ref%d_dir%d \n', block, trial, setup.refdirec(block, trial), setup.direction.stimtwo(block, trial));
                        end
                        
                        TimingCnt = flip.int2.VBL(block, trial, frameNum) + window.frameDur - window.slack;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % CONTINUOUS ESTIMATION
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    results.estimation.response(block, trial) = setup.refdirec(block, trial);
                    choicemade = false;
                    
                    % put the mouse in the right place
                    SetMouse(round((dots.radius+(deg2pix(window, 1))).*cos(results.estimation.response(block, trial)*pi/180))+ window.center(1), ...
                        round((dots.radius+(deg2pix(window, 1))).*sin(results.estimation.response(block, trial)*pi/180))+window.center(2));
                    ShowCursor('CrossHair');
                    
                    while GetSecs - flip.int2.VBL(block,trial,frameNum) < setup.estimationmaxtime && choicemade == false,
                        
                        %  add a white dial in the reference direction for this subject
                        window = drawAllDots_refmark(setup, window, dots, block, trial, stim.two, frameNum);
                        
                        % get the mouse position on this frame
                        [x, y, clicks]  = GetMouse();
                        
                        % convert into the angle and save
                        results.estimation.response(block, trial) = atan2(y-window.center(2), x-window.center(1)) * (180/pi);
                        
                        % display the current line in redx =
                        window      = drawFixation(window, fix, dots); % fixation
                        Screen('DrawLines', window.h, ...
                            [(dots.radius-(deg2pix(window, 3))).*cos(results.estimation.response(block, trial)*pi/180), ...
                            (dots.radius+(deg2pix(window, 1))).*cos(results.estimation.response(block, trial)*pi/180); (dots.radius-(deg2pix(window, 3))).*sin(results.estimation.response(block, trial)*pi/180), ...
                            (dots.radius+(deg2pix(window, 1))).*sin(results.estimation.response(block, trial)*pi/180)], ...
                            2, fix.color, window.center);
                        Screen('Flip', window.h);
                        
                        % get response
                        if any(clicks),
                            choicemade = true;
                            results.estimation.RT(block, trial)    = GetSecs - flip.int2.VBL(block,trial, frameNum);
                            % the last results.estimation.response is the last one that was saved (the same as the one that returned 'clicked'
                        else
                            % continue the while loop
                        end
                        
                        % exit the whole experiment if escape is pressed
                        [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
                        if keyIsDown,
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
                        end
                        
                    end % estimation phase end
                    
                    if setup.Eye,
                        Eyelink ('Message', sprintf('block%d_trial%d_estimation_ref%d_meandir%.2f_resp%d_RT%.3f', block, trial, setup.refdirec(block, trial), setup.direction.meanestimate(block, trial), results.estimation.response(block, trial),  results.estimation.RT(block, trial)));
                        fprintf('\n block%d_trial%d_estimation_ref%d_meandir%.2f_resp%f_RT%.3f \n', block, trial, setup.refdirec(block, trial), setup.direction.meanestimate(block, trial), results.estimation.response(block, trial),  results.estimation.RT(block, trial));
                    end
                    HideCursor;
                end
            end % choices
            
            % break out of all trials if ESC was pressed
            if setup.cancel == true,
                break
                warning('experiment was manually terminated');
            end
            
        end %end trial loop for this block
        
        % break text
        Screen('DrawText',window.h, sprintf('Finished block %d!', block) , ...
            window.center(1)*0.60, window.center(2)*0.40 , [255 255 255] );
        if block < setup.nblocks,
            Screen('DrawText',window.h, 'Take a break!', window.center(1)*0.6, window.center(2)*.6, [255 255 255] );
        else % finish
            Screen('DrawText',window.h, 'Done! Please call the experimenter.', window.center(1)*0.60, window.center(2)*.6 , [255 255 255] );
        end
        Screen('Flip', window.h);

        % break out of all blocks if ESC was pressed
        if setup.cancel == true,
            break
            warning('experiment was manually terminated');
        end
        
        WaitSecs(.1);
        % wait for keypress to start with the next block
        MouseWait;
        WaitSecs(.1);
        
    end % end this block
end

% !!!!!!!!!!!!!!!!!!!!!!!!! SAVE THE FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% wrap up and save
setup.datetime = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
% create subject specific file and save - add unique datestring to avoid any overwriting of files
setup.filename = sprintf('Data/P%d_training_%s.mat', setup.participant, setup.datetime);

save(setup.filename, '-mat', 'setup', 'window', 'stim', 'dots', 'fix', 'results', 'audio', 'sound', 'flip');
disp('SAVED FILE TO DISK');
% !! important, make sure all the variables of interest are saved to
% possibly recreate the exact experiment later.

disp('done!'); Screen('CloseAll'); ShowCursor;
PsychPortAudio('Stop', audio.h);
sca;

