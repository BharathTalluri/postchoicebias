function [] = script6_Main_CommitmentRDK_mouse()
%% Commitment RDK code
% this models the methodology of Bronfman, Brezis & Usher with random dot
% stimuli. To be used by Ana Vojvodic.

% Anne Urai, March 2014
% ----------------------------------------------------------------

clear all; close all; clc;
addpath('D:\USERS\AnneUrai\Commitment');
addpath('D:\USERS\AnneUrai\Commitment\stats');
%cd('D:\USERS\AnneUrai\Commitment\');
pwd

% general setup
setup.startatblock  = 1; % if anything went wrong during the session, can easily restart at a later block
setup.Eye           = false; % true if using Eyelink
setup.cancel        = false; % becomes true if escape is pressed, will abort experiment (but save data)

% ask for subject number, default: 0
setup.participant       = input('Participant number? ');
if isempty(setup.participant)
    setup.participant   = 20; %test
end

% ask for session number, default: 0
setup.session        = input('Session? ');
if isempty(setup.session)
    setup.session    = 0; %test
end

if setup.startatblock > 1, %
    % load the file that was generated at block 1
    cd Data/
    allfiles = dir(fullfile('*.mat'));
    for file = 1:length(allfiles),
        % find the file that matches the participant and session nr
        if strncmp(sprintf('P%d_s%d_', setup.participant, setup.session), allfiles(file).name, 6);
            thisfile = allfiles(file).name;
        end
    end
    load(thisfile);
    disp('PREVIOUSLY CREATED FILE LOADED IN');
end

% Setup the PsychToolbox
window.dist             = 50; % viewing distance in cm ( this is fixed in the chinrest, D1.09)
window.width            = 40; % physical width of the CRT screen in cm, 40 cm in D1.09
window.skipChecks       = 1; % set to 1 to skip VBL tests and avoid warnings
[window, audio]         = SetupPTB(window); %load all the psychtoolbox things
HideCursor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Screen('TextSize', window.h, 15);
Screen('DrawText',window.h, 'Loading...', window.center(1)*0.60, window.center(2)*0.75 , [255 255 255] );
Screen('Flip', window.h);

% if restarting the experiment after some failure (eg. unable to move
% the EyeLink file), the script will use the previously generated config
% and continue at the next block

if setup.startatblock == 1, %
    % create the design and load stimuli
    % the function configuration contains all the important experimental design parameters
    [setup, dots, fix, results, sound, stim] = configuration_main(window, audio, setup);
else
    % fill the sound buffer again
    PsychPortAudio('FillBuffer', audio.h, sound.tonebuf);
end
cd ..

setup.filename = sprintf('~/Data/MotionEnergy/Test_90coh.mat');
save(setup.filename, '-mat', 'dots','setup', 'stim', 'window' );

sca;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP OVER BLOCKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for block = setup.startatblock:setup.nblocks,
    
    if setup.Eye == true,
        edfFile = ELconfig(window, setup, block);
    end
    
    % turn the screen black after eyelink setup
    Screen('FillRect', window.h, window.black);
    Screen('DrawText', window.h, 'Are you ready?', window.center(1)*0.60, window.center(2)*0.75 , [255 255 255] );
    Screen('DrawText', window.h, 'Click the mouse to start.', window.center(1)*0.60, window.center(2)*0.6 , [255 255 255] );
    
    Screen('Flip', window.h);
    
    WaitSecs(.1);
    MouseWait; % after displaying the introductory text, wait for a mouse click
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOOP OVER TRIALS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for trial = 1:setup.ntrials,
        
        % allow for blink breaks every setup.blinkinterval trials
        % the subjects can continue by pressing any key
        if setup.Eye,
            Eyelink ('Message', 'blinkbreak_start');
        end
        
        Screen('FillRect', window.h, window.black);
        Screen('Flip', window.h); % flip once, so stationary dots
        %WaitSecs(1); % 1s ISI instead of manual!
        MouseWait; %wait for a keypress
        
        if setup.Eye,
            Eyelink ('Message', 'blinkbreak_end');
        end
        
        % variable fixation interval
        for frameNum = 1:ceil(setup.fixtime(block, trial)*window.frameRate),
            
            window      = drawFixation(window, fix, dots); % fixation
            [flip.fix.VBL(block, trial, frameNum), ...
                flip.fix.StimOns(block, trial, frameNum), ...
                flip.fix.FlipTS(block, trial, frameNum), ...
                flip.fix.Missed(block, trial, frameNum), ...
                flip.fix.beampos(block, trial, frameNum)] = Screen('Flip', window.h);
            
            if setup.Eye == true && frameNum == 1,
                Eyelink ('Message', sprintf('block%d_trial%d_fix', block, trial));
                fprintf('\n block%d_trial%d_fix \n', block, trial)
            end
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
                fprintf('\n block%d_trial%d_int1_ref%d_dir%d \n', block, trial, setup.refdirec(block, trial), setup.direction.stimone(block, trial));
            end
            TimingCnt = flip.int1.VBL(block, trial, frameNum) + window.frameDur - window.slack;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BINARY CHOICE OR NO CHOICE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        window      = drawAllDots_refmark(setup, window, dots, block, trial, stim.one, frameNum);
        window      = drawFixation(window, fix, dots); % fixation
        
        % present the choice or no-choice prompt
        switch abs(setup.choice(block, trial)),
            case 1 % make a binary choice
                PsychPortAudio('SetLoop',audio.h, sound.tonepos(4,1), sound.tonepos(4,2));
            case 0 %no choice
                PsychPortAudio('SetLoop',audio.h, sound.tonepos(5,1), sound.tonepos(5,2));
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
                        otherwise
                            disp('did not recognize key pressed');
                    end
                catch
                    disp('did not recognize key pressed');
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
                results.binary.correct(block, trial) = NaN;
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
            disp('INVALID RESPONSE');
            
        elseif setup.choice(block, trial) == -1, % choice only
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % normal feedback
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % present auditory feedback
            if results.binary.correct(block,trial) == 1, % correct
                PsychPortAudio('SetLoop',audio.h, sound.tonepos(1,1), sound.tonepos(1,2));
                
            elseif results.binary.correct(block,trial) == 0,    %incorrect
                PsychPortAudio('SetLoop',audio.h, sound.tonepos(2,1), sound.tonepos(2,2));
                
            elseif sign(setup.direction.stimone(block,trial)-setup.refdirec(block, trial)) == 0 && abs(results.binary.response(block, trial))== 1,
                % if there was no deviation from the reference, select random feedback
                throwcoin = randi(2);
                switch throwcoin,
                    case 1
                        PsychPortAudio('SetLoop',audio.h, sound.tonepos(1,1), sound.tonepos(1,2));
                    case 2
                        PsychPortAudio('SetLoop',audio.h, sound.tonepos(2,1), sound.tonepos(2,2));
                end
            else
                % unrecognized response
                setup.cancel = true;
                warning('could not determine which feedback to give');
            end
            results.feedbackonset(block, trial) = GetSecs;
            PsychPortAudio('Start', audio.h);
            
            %send trigger to EL
            if setup.Eye,
                Eyelink ('Message', sprintf('block%d_trial%d_feedback_correct%d', block, trial, results.binary.correct(block, trial)));
                fprintf('\n block%d_trial%d_feedback_correct%d \n', block, trial, results.binary.correct(block, trial));
            end
            
            Screen('Flip', window.h);
            WaitSecs(.5); % leave visual feedback on the screen for a bit, and wait a bit after feedback
            
        elseif setup.choice(block, trial)>-1, % if this trial has a second interval
            
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
            % ShowCursor('CrossHair');
            
            while GetSecs - flip.int2.VBL(block,trial,frameNum) < setup.estimationmaxtime && choicemade == false,
                
                %  add a white dial in the reference direction for this subject
                window = drawAllDots_refmark(setup, window, dots, block, trial, stim.two, frameNum);
                
                % get the mouse position on this frame
                [x, y, clicks]  = GetMouse();
                
                % convert into the angle and save
                results.estimation.response(block, trial) = atan2(y-window.center(2), x-window.center(1)) * (180/pi);
                
                % display the current line in red
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
            %HideCursor;
            
        end % choices
        
        % break out of all trials if ESC was pressed
        if setup.cancel == true,
            break
            warning('experiment was manually terminated');
        end
        
    end %end trial loop for this block
    
    % break text
    Screen('DrawText',window.h, sprintf('Finished block %d.', block) , ...
        window.center(1)*0.60, window.center(2)*0.40 , [255 255 255] );
    % calculate the vector with response estimates
    switch mod(setup.participant,2), % counterbalance for every two subjects
        case 0
            esterror = abs(results.estimation.response(block, :) - setup.refdirec(block, :) - setup.direction.meanestimate(block, :));
        case 1
            esterror = abs(results.estimation.response(block, :)+ 360 - setup.refdirec(block, :) - setup.direction.meanestimate(block, :));
    end
    esterror(isnan(esterror)) = [];
    
    % draw
    Screen('DrawText',window.h, sprintf('Your estimation was on average '), ...
        window.center(1)*0.60, window.center(2)*0.60 , [255 255 255] );
    Screen('DrawText',window.h, sprintf('%.3f degrees away ', median(esterror)), ...
        window.center(1)*0.60, window.center(2)*0.70 , [255 255 255] );
    Screen('DrawText',window.h, sprintf('from the true direction of motion.'), ...
        window.center(1)*0.60, window.center(2)*0.80 , [255 255 255] );
    Screen('Flip', window.h);
    
    % draw some text for the experimenter
    fprintf('\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n');
    fprintf('\n Finished block %d out of %d! \n', block, setup.nblocks)
    fprintf('\n Estimation was median %.3f degrees away from the true direction of motion \n', median(esterror));
    % compute the overall binary accuracy
    trials = find((setup.direction.stimone(block, :)-setup.refdirec(block, :))~=0 & abs(setup.choice(block, :))==1 & ~isnan(results.binary.response(block, :)) & results.binary.correct(block, :) > -1);
    fprintf('\n Binary choice accuracy = %.2f percent \n \n', 100*mean(results.binary.correct(block, trials)));
    fprintf('\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n');
    
    %% save the EL file for this block
    setup.datetime      = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    setup.eyefilename = sprintf('D:/USERS/AnneUrai/Commitment/Data/P%d_s%d_b%d_%s.edf', setup.participant, setup.session, block, setup.datetime);
    Eyelink('CloseFile');
    Eyelink('WaitForModeReady', 500);
    try
        status              = Eyelink('ReceiveFile',edfFile, setup.eyefilename); %this collects the file from the eyelink
        disp(status);
        disp(['File ' setup.eyefilename ' saved to disk']);
    catch
        warning(['File ' setup.eyefilename ' not saved to disk']);
    end
    
    % break out of all blocks if ESC was pressed
    if setup.cancel == true,
        break
        warning('experiment was manually terminated');
    end
    
    % continue to the next block
    if block < setup.nblocks,
        Screen('DrawText',window.h, 'Take a break! Click to continue.', window.center(1)*0.4, window.center(2)*.9, [255 255 255] );
    else % finish
        Screen('DrawText',window.h, 'Done! Please call the experimenter.', window.center(1)*0.60, window.center(2)*.9 , [255 255 255] );
    end
    Screen('Flip', window.h);
    
    WaitSecs(.1);
    % wait for keypress to start with the next block
    MouseWait;
    WaitSecs(.1);
    
end % end this block

% !!!!!!!!!!!!!!!!!!!!!!!!! SAVE THE FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

setup.datetime = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
% create subject specific file and save - add unique datestring to avoid any overwriting of files
setup.filename = sprintf('D:/USERS/AnneUrai/Commitment/Data/P%d_s%d_%s.mat', setup.participant, setup.session, setup.datetime);
save(setup.filename, '-mat', 'audio','dots', 'fix', 'flip','results','setup',  'sound', 'stim', 'window' );
disp(['SAVED FILE TO DISK ' setup.filename]);

% close the eyelink
if setup.Eye == true,
    Eyelink('StopRecording');
end

disp('done!'); Screen('CloseAll'); ShowCursor;
PsychPortAudio('Stop', audio.h);
sca;

if setup.Eye,
    msgbox('MEASURE DUMMY PUPIL NOW!');
end

%% show the correlation between their estimate and the actual mean direction
% for all the blocks in this session
estimates = results.estimation.response-setup.refdirec;
switch mod(setup.participant,2), % recompute trials just across the boundary of the unit circle
    case 0
        estimates(estimates<-200) = estimates(estimates<-200)+360;
    case 1
        estimates(estimates<200) = estimates(estimates<200)+360;
end
figure; scatter(estimates(:), setup.direction.meanestimate(:));
ylim([-25 25]);
ylabel('actual direction'); xlabel('estimation');

%% also show a measure of flip performance
figure;
for b = 1:setup.nblocks,
    for t = 1:setup.ntrials,
        try
            subplot(221);
            plot(diff(squeeze(flip.int1.StimOns(b,t,:))));
            if any(diff(squeeze(flip.int1.StimOns(b,t,:))) < 0.01 | diff(squeeze(flip.int1.StimOns(b,t,:))) > .02),
                fprintf('weird fliptime block %d, trial %d \n', b, t);
            end
            hold on;
        end
        try
            subplot(222);
            plot(diff(squeeze(flip.int2.StimOns(b,t,:))));
            hold on;
        end
        drawnow; %WaitSecs(0.05);
    end
end
subplot(221);
axis tight; title('int1'); xlabel('frames'); ylabel('frameDur');
axis tight; title('int2'); xlabel('frames'); ylabel('frameDur');

end
