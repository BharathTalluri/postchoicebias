function [ window ] = drawAllDots_refmark(setup, window, dots, block, trial, stimuli, frameNum)

% find the index of this stimulus from the block and trial
trialindex = sub2ind(size(setup.refdirec'), trial, block);
%disp(trialindex);

% draws the dots on the screen for this frame
Screen('DrawDots',window.h,squeeze(stimuli(trialindex, frameNum, :, :)), ...
    dots.size, dots.color, window.center, 2);

% calculate the coordinates of the reference mark line
xstart  = (dots.radius-(deg2pix(window, 3))) .*cos(setup.refdirec(block, trial)*pi/180);
xend    = (dots.radius+(deg2pix(window, 1))) .*cos(setup.refdirec(block, trial)*pi/180);

ystart  = (dots.radius-(deg2pix(window, 3))) .*sin(setup.refdirec(block, trial)*pi/180);
yend    = (dots.radius+(deg2pix(window, 1))) .*sin(setup.refdirec(block, trial)*pi/180);

% also add a dial in the reference direction for this subject
Screen('DrawLines', window.h, [xstart xend; ystart yend], 2, dots.color, window.center);

end