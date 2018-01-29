function [ stimuli] = dots_brownian(setup, window, dots, block, trial, whichstim)

% Generates a dot patch according to a brownian motion
% algorithm, see:
% Pilly, P. K., & Seitz, A. R. (2009).
% What a difference a parameter makes: a psychophysical comparison
% of random dot motion algorithms. Vision Research, 49(13), 1599?612.
%--------------------------------------------------------------------------

COHERENCE   = dots.coherence;
NDOTS       = dots.nDots; % 774, 3097
SPEED       = dots.speed/3; % (in degrees per second)

switch whichstim
    case 'one'
        DIRECTION   = dots.direction.stimone(block, trial);
    case 'two'
        DIRECTION   = dots.direction.stimtwo(block, trial);
end

% calculate the displacement in x and y directions (per frame)
STEPX = SPEED/window.frameRate*cos(DIRECTION*pi/180);
STEPY = SPEED/window.frameRate*sin(DIRECTION*pi/180);

NFR         = ceil(setup.nframes);
RADIUS      = dots.radius;
INNER       = dots.innerspace; %inner circle doesn't contain dots

% generate random starting points within a circular aperture
rad     = RADIUS*sqrt(rand(NDOTS,1));
theta   = 2*pi*rand(NDOTS,1);
pos     = [rad.*cos(theta) rad.*sin(theta)]; %convert to cartesian coordinates

for frameNum = 1:NFR,
    
    % define the noise dots
    noisedots           = round(NDOTS*(1-COHERENCE)); %nr of noisedots based on the coherence
    t1                  = rand(NDOTS,1);
    [t1,t2]             = sort(t1);
    noiseindex          = t2(1:noisedots); %random subset of dots
    % define signal dots
    signalindex         = t2(noisedots+1:end); %the dots that are signal
    
    % define the motion direction of the noise dots
    noisedir            = randi(360,[1, noisedots]); %random direction between 0 and 180 degrees in degrees
    
    % move the noise dots in a noise direction
    noisemotion         = [(SPEED.*cos(noisedir*pi/180)/window.frameRate); (SPEED.*sin(noisedir*pi/180)/window.frameRate)]';
    pos(noiseindex,:)   =  pos(noiseindex,:) + noisemotion; %update
    
    % move the signal dots in the signal direction
    signalmotion        = [(SPEED*cos(DIRECTION*pi/180)/window.frameRate); (SPEED*sin(DIRECTION*pi/180)/window.frameRate)]';
    pos(signalindex,:)  =  [pos(signalindex,1) + signalmotion(1) pos(signalindex,2) + signalmotion(2)]; %convert to cartesian coordinates
    
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
    % move them to the other side of the inner annulus?
    rad                 = INNER; % border of the inner space
    theta               = theta + pi; % other side
    pos(innerindex, :)  = [rad.*cos(theta) rad.*sin(theta)]; %move back to the top, only change y coordinate
    
    % save the positions per variant
    stimuli(frameNum, :, :)  = pos';
end

end