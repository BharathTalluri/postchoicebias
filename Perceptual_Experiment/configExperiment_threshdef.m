function [setup, dots, fix, results, flip] = configExperiment_threshdef(window, setup)

% general experimental design
setup.cohlevel          = [0 2.5, 5, 10, 20, 40] / 100;
setup.trialrep          = 100; %nr of trials per coherence level (per block) - should be 100
setup.totalntrials      = round(setup.trialrep * length(setup.cohlevel));
setup.nblocks           = 6;
setup.ntrials           = setup.totalntrials / setup.nblocks;
assert(~mod(setup.ntrials,1), 'Number of trials per block is not integer');
setup.counterbalancing  = mod(setup.participant,2);
% use only one block if at all possible
setup.cancel = false;

% timing
setup.fixtime           = (.5 + .5*rand(setup.nblocks, setup.ntrials)); %generate the fixation time 0.6 - 0.8s (36-48 frames)
setup.viewingtime       = .75; % viewing duration in seconds (fixed in this script, or maximum viewing duration in RT paradigm
setup.nframes           = ceil(setup.viewingtime*window.frameRate); %number of frames the stimulus is displayed
setup.resptime          = 3; % maximum time for response

%% design

[dots, fix]             = setupDots(window, setup);

% do some randomization
dots.direction          = Shuffle(repmat([90 270]', [setup.ntrials/2 setup.nblocks]))';
dots.coherence          = Shuffle(repmat(setup.cohlevel', [ceil(setup.ntrials/length(setup.cohlevel)) setup.nblocks ]))';
dots.coherence          = dots.coherence(1:setup.nblocks, 1:setup.ntrials);

%% preallocate results and stimuli structure
% preallocation is a good habit to make sure that Matlab knows how big your
% output structures will be. You might run into memory problems if your
% structures grow on each loop - Matlab will have to find a new chunk of
% memory each time which costs significant time.

results.response            = NaN(setup.nblocks, setup.ntrials);
results.correct             = NaN(setup.nblocks, setup.ntrials);
results.RT                  = NaN(setup.nblocks, setup.ntrials);

% preallocate a full flip structure to store the output of every dynamic
% flip
flip.fix.VBL            = nan(setup.nblocks, setup.ntrials, ceil(max(max(setup.fixtime))/window.frameDur));
flip.fix.StimOns        = flip.fix.VBL;
flip.fix.FlipTS         = flip.fix.VBL;
flip.fix.Missed         = flip.fix.VBL;
flip.fix.beampos        = flip.fix.VBL;

flip.stim.VBL           = nan(setup.nblocks, setup.ntrials, setup.nframes);
flip.stim.StimOns       = flip.stim.VBL;
flip.stim.FlipTS        = flip.stim.VBL;
flip.stim.Missed        = flip.stim.VBL;
flip.stim.beampos       = flip.stim.VBL;

flip.resptime.VBL        = nan(setup.nblocks, setup.ntrials, ceil(setup.resptime/window.frameDur));
flip.resptime.StimOns    = flip.resptime.VBL;
flip.resptime.FlipTS     = flip.resptime.VBL;
flip.resptime.Missed     = flip.resptime.VBL;
flip.resptime.beampos    = flip.resptime.VBL;

end