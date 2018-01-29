function [setup, dots, fix, results, sound] = configExperiment_example(setup, window, audio)
% function for the configuration of the experimental design

setup.cancel = false; %this turns true if escape is pressed

% example - easy
setup.cohlevel         = .25;
setup.directions        = [90 270]; % up or down

setup.totalntrials      = 20; % number of trials per coherence level for this whole session
setup.startatblock      = 1;
setup.nblocks           = 1; % number of blocks
setup.ntrials           = setup.totalntrials/setup.nblocks; % total number of trials per block
setup.counterbalancing  = mod(setup.participant,2);

% timing
setup.fixtime           = (.5 + .5*rand(setup.nblocks, setup.ntrials)); %generate the fixation time 0.6 - 0.8s (36-48 frames)
setup.viewingtime       = 10; % maximum viewing duration in seconds
setup.nframes           = ceil(setup.viewingtime * window.frameRate); %number of frames the stimulus is windowed

%% (output) dots: set stimuli (dots) parameters

[dots, fix]                 = setupDots(window, setup);
dots.direction              = Shuffle(repmat(setup.directions, 1, setup.totalntrials));

%% (output) results: preallocate results structures

results.response            = NaN(setup.nblocks,setup.ntrials);
results.correct             = NaN(setup.nblocks,setup.ntrials);
results.RT                  = NaN(setup.nblocks,setup.ntrials);
results.stimonset           = NaN(setup.nblocks,setup.ntrials,setup.nframes);

% auditory feedback, using some ENS functions (probably from Valentin)
sound = setupSounds(setup, audio);

end % end of the configExperiment function