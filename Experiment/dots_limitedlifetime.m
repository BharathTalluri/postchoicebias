function [ stimuli] = dots_limitedlifetime(setup, window, dots, block, trial, whichstim)

% Generates a dot patch according to a limited lifetime
% algorithm, see:
% Pilly, P. K., & Seitz, A. R. (2009).
% What a difference a parameter makes: a psychophysical comparison
% of random dot motion algorithms. Vision Research, 49(13), 1599?612.
%--------------------------------------------------------------------------

NVAR        = dots.nvar; % interleaved sequence
if length(setup.cohlevel)> 1,
    COHERENCE   = dots.coherence(block, trial);
else
    COHERENCE   = dots.coherence;
end

LIFETIME    = dots.lifetime;
NDOTS       = dots.nDots; % 774, 3097
SPEED       = dots.speed; % (in degrees per second)

switch whichstim
    case 'one'
        DIRECTION   = dots.direction.stimone(block, trial);
    case 'two'
        DIRECTION   = dots.direction.stimtwo(block, trial);
    otherwise
        DIRECTION   = dots.direction(block, trial);
end

% calculate the displacement in x and y directions (per frame)
STEPX = SPEED/window.frameRate*cos(DIRECTION*pi/180);
STEPY = SPEED/window.frameRate*sin(DIRECTION*pi/180);

NFR         = ceil(setup.nframes/NVAR);
RADIUS      = dots.radius;
INNER       = dots.innerspace; %inner circle doesn't contain dots

for var = 1:NVAR, %create separate variants
    
    % generate random starting points within a circular aperture
    rad     = RADIUS*sqrt(rand(NDOTS,1));
    theta   = 2*pi*rand(NDOTS,1);
    pos     = [rad.*cos(theta) rad.*sin(theta)]; %convert to cartesian coordinates
    
    % assign each dot a scalar that indicates how long it has been
    % 'alive' (that is, a signal dot)
    % Each dot will have a integer value 'life' which is how many frames the
    % dot has been going.  The starting 'life' of each dot will be a random
    % number between 0 and dots.lifetime-1 so that they don't all 'die' on the
    % same frame:
    
    life    = ceil(rand(1,NDOTS)*LIFETIME);
    
    for frameNum = 1:NFR,
        
        % define the noise dots
        noisedots           = round(NDOTS*(1-COHERENCE)); %nr of noisedots based on the coherence
        t1                  = rand(NDOTS,1);
        [t1,t2]             = sort(t1);
        noiseindex          = t2(1:noisedots); %random subset of dots
        % define signal dots
        signalindex         = t2(noisedots+1:end); %the dots that are signal
        
        % move signal dots with a certain speed in the right direction
        pos(signalindex,:)  = [STEPX+pos(signalindex,1) STEPY+pos(signalindex,2)];
        
        % replot the noisedots somewhere in the aperture
        rad                 = RADIUS*sqrt(rand(noisedots,1));
        theta               = 2*pi*rand(noisedots,1);
        pos(noiseindex,:)   = [rad.*cos(theta) rad.*sin(theta)]; %convert to cartesian coordinates
        
        %increment the 'life' of each dot
        life                = life+1;
        
        %find the 'dead' dots
        deadindex           = mod(life,LIFETIME)==0;
        deaddots            = length(find(deadindex==1));
        
        %replace the positions of the dead dots to a random location
        rad                 = RADIUS*sqrt(rand(deaddots,1));
        theta               = 2*pi*rand(deaddots,1);
        pos(deadindex,:)    = [rad.*cos(theta) rad.*sin(theta)]; %convert to cartesian coordinates
        
        %find the dots that have left the aperture
        outindex            = find(sqrt((pos(:,1).^2)+(pos(:,2).^2)) >= RADIUS); %index dots that are outside the aperture
        % wrap them around in the direction where they came from
        [theta, rad]        = cart2pol(pos(outindex, 1), pos(outindex,2));
        theta               = theta + pi; %move to the other side of the circle
        rad                 = RADIUS*ones(length(rad),1);
        pos(outindex, :)    = [rad.*cos(theta) rad.*sin(theta)]; %move back to the top, only change y coordinate
        
        %find the dots that are too close to the fixation
        innerindex          = find(sqrt((pos(:,1).^2)+(pos(:,2).^2)) <= INNER);
        [theta, rad]        = cart2pol(pos(innerindex, 1), pos(innerindex, 2));
        rad                 = INNER + (RADIUS - INNER)*sqrt(rand(length(innerindex),1)); %random radius
        theta               = theta + pi;
        pos(innerindex, :)  = [rad.*cos(theta) rad.*sin(theta)]; %move back to the top, only change y coordinate
        
        % save the positions per variant
        temp{var,frameNum}  = pos';
    end
end

%% this is not the most elegant code, but is used to interleave the multiple variants of the stimulus created
stimuli = nan(NFR, 2, NDOTS);

if NVAR > 1,
    % interleave the cell arrays with position information
    Interleaved = reshape(temp, 1, NFR*NVAR);
    
    for f = 1:length(Interleaved),
        % allocate all the dot positions from this trial to the stimulus structure
        stimuli(f, :, :) = Interleaved{f};
    end
    
    % possibly cut the last frames off to have the correct dims
    stimuli = stimuli(1:setup.nframes, :, :);
    
else %one variant, no interleaving
    for f = 1:NFR,
        stimuli(f, :, :) = temp{f};
    end
end


end