function [setup, dots, fix, results, sound, stim] = configuration_main(window, audio, setup)
% this function contains all the important information about the
% experimental design, and pregenerates the stimuli
%
% Anne Urai, University of Amsterdam, 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% permutations of the possible directions
perm        = unique(nchoosek(repmat([-20 -10 0 10 20], 1,5), 2), 'rows');

% remove the permutations that are too extreme
perm(5,:) = [];
perm(20, :) = [];

% general experimental design
try
    load(sprintf('Data/P%d_threshold.mat', setup.participant), 'fit1');
    setup.cohlevel         = fit1.x(dsearchn(fit1.y', 0.80)); % 80% correct point of the 'raw' Weibull fit
    fprintf('using individual threshold, %.3f for this sj \n', setup.cohlevel);
catch
    warning('No threshold file found, using canonical 90% threshold');
    setup.cohlevel = .9;
end

setup.trialrep          = 5; % total nr of trials per permutation - within choice only
setup.percentagechoiceonlytrials = 3; % one third choice-1
setup.totalntrials      = setup.trialrep*setup.percentagechoiceonlytrials*length(perm);
setup.nblocks           = 5;
setup.ntrials           = setup.totalntrials / setup.nblocks; %nr of trials per block
%make sure this leaves an integer nr of trials per block!
assert(~mod(setup.ntrials,1), 'Number of trials per block is not integer');

% trial repetition (within choice trials, 1/3rd of all trials)
fulldesign              = vertcat(repmat(perm, [setup.trialrep 1]));
meanperm                = mean(fulldesign')'; % do hist(meanperm, 20) to plot the distribution
switch mod(setup.participant,2), % counterbalance for every two subjects
    case 0
        refdirec    = randi([1 180], [length(fulldesign) 1]); % generate random reference directions for each trial
    case 1
        refdirec    = randi([180 360], [length(fulldesign) 1]);
end
fulldesign        = [fulldesign meanperm refdirec]; %add the mean in the third column

% add this subject's reference direction to it to get the actual motion
% direction
fulldesign(:, 1)      = fulldesign(:,1) + fulldesign(:, 4);
fulldesign(:, 2)      = fulldesign(:,2) + fulldesign(:, 4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIMING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setup.fixtime           = (.5 + .5*rand(setup.nblocks, setup.ntrials)); %generate the fixation time 0.6 - 0.8s (36-48 frames)
setup.stimtime          = .75; % viewing duration in seconds (fixed in this script, or maximum viewing duration in RT paradigm
setup.nframes           = round(window.frameRate*setup.stimtime);
setup.choicetime        = 2;  % identical timing for response with or without choice
setup.estimationmaxtime = 10; % maximum time for estimation task
setup.blinkbreakntrials = 1;  % blink break each trial
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dots, fix] = setupDots(window, setup);

% put the directions into the dots structure
dots.direction.stimone = fulldesign(:,1)';
dots.direction.stimtwo = fulldesign(:,2)';

% preload  the dot coordinates with nans
% dimensions: trial (1-115), frames, x/y, number of dots.
tmpstim.one = nan(length(fulldesign), setup.nframes, 2, dots.nDots);
tmpstim.two = nan(length(fulldesign), setup.nframes, 2, dots.nDots);

for t = 1:setup.nblocks*setup.ntrials/setup.percentagechoiceonlytrials,
    
    % preload all the dot coordinates (for choice)
    tmpstim.one(t, :, :, :)      = dots_limitedlifetime(setup, window, dots, 1, t, 'one');
    tmpstim.two(t, :, :, :)      = dots_limitedlifetime(setup, window, dots, 1, t, 'two');
    
    % display the progress of loading
    if ~ isempty(Screen('Windows')), %if there is a window open
        Screen('DrawText',window.h, sprintf('Loading %d percent', ...
            round(t/(setup.nblocks*setup.ntrials/setup.percentagechoiceonlytrials)*100)), ...
            window.center(1)*0.60, window.center(2)*0.75 , [255 255 255] );
        Screen('Flip', window.h);
    end
end

% duplicate the exact dots coordinates so that they appear (i) once in choice
% (ii) once in nochoice and (iii) once in choice only - x3
tmpstim2.one = repmat(tmpstim.one, [3 1 1 1]);
tmpstim2.two = repmat(tmpstim.two, [3 1 1 1]);

% also duplicate the design matrix
% assign all these design options either choice (1), no choice (0) or
% choice only (-1)
fulldesign2 = vertcat([fulldesign ones(length(fulldesign), 1)], ...
    [fulldesign zeros(length(fulldesign), 1)], ...
    [fulldesign(1:end, :) -ones(length(fulldesign), 1)]);

% now randomly order all those trials within today's session
setup.randomize     = randperm(length(fulldesign2)); % keep the seed to be able to reconstruct the randomization!
%setup.randomize     = 1:length(fulldesign2); % only for testing!
setup.fulldesign    = fulldesign2(setup.randomize, :);

% crucial: shuffle the actual stimuli in the same way to be the same as
% the design matrix
stim.one = tmpstim2.one(setup.randomize, :, :, :);
stim.two = tmpstim2.two(setup.randomize, :, :, :);

%!!!! for each stimulus combination (first and second direction), the
%choice and no-choice versions now have identical dots but are all randomly
%shuffled throughout today's session.

% put everything back nicely into separate design components
setup.direction.stimone          = reshape(setup.fulldesign(:,1), [setup.ntrials, setup.nblocks])';
setup.direction.stimtwo          = reshape(setup.fulldesign(:,2), [setup.ntrials, setup.nblocks])';
setup.direction.meanestimate     = reshape(setup.fulldesign(:,3), [setup.ntrials, setup.nblocks])';
setup.refdirec                   = reshape(setup.fulldesign(:,4), [setup.ntrials, setup.nblocks])';
setup.choice                     = reshape(setup.fulldesign(:,5), [setup.ntrials, setup.nblocks])';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOUNDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sound = setupSounds(setup, audio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREALLOCATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preallocate results and stimuli structure
% preallocation is a good habit to make sure that Matlab knows how big your
% output structures will be. You might run into memory problems if your
% structures grow on each loop - Matlab will have to find a new chunk of
% memory each time which costs significant time.

results.binary.response            = NaN(setup.nblocks,setup.ntrials);
results.binary.correct             = NaN(setup.nblocks,setup.ntrials);
results.binary.RT                  = NaN(setup.nblocks,setup.ntrials);

% create results structures for both the binary and continuous judgments
% they will make on each trial
results.estimation.response        = NaN(setup.nblocks,setup.ntrials);
results.estimation.correct         = NaN(setup.nblocks,setup.ntrials);
results.estimation.RT              = NaN(setup.nblocks,setup.ntrials);

% preallocate a full flip structure to store the output of every dynamic
% flip
flip.fix.VBL            = nan(setup.nblocks, setup.ntrials, ceil(max(max(setup.fixtime))/window.frameDur));
flip.fix.StimOns        = flip.fix.VBL;
flip.fix.FlipTS         = flip.fix.VBL;
flip.fix.Missed         = flip.fix.VBL;
flip.fix.beampos        = flip.fix.VBL;

flip.int1.VBL           = nan(setup.nblocks, setup.ntrials, setup.nframes);
flip.int1.StimOns       = flip.int1.VBL;
flip.int1.FlipTS        = flip.int1.VBL;
flip.int1.Missed        = flip.int1.VBL;
flip.int1.beampos       = flip.int1.VBL;

flip.interval.VBL        = nan(setup.nblocks, setup.ntrials, ceil(max(max(setup.choicetime))/window.frameDur));
flip.interval.StimOns    = flip.interval.VBL;
flip.interval.FlipTS     = flip.interval.VBL;
flip.interval.Missed     = flip.interval.VBL;
flip.interval.beampos    = flip.interval.VBL;

flip.int2.VBL           = nan(setup.nblocks, setup.ntrials, setup.nframes);
flip.int2.int2Ons       = flip.int2.VBL;
flip.int2.FlipTS        = flip.int2.VBL;
flip.int2.Missed        = flip.int2.VBL;
flip.int2.beampos       = flip.int2.VBL;

end
