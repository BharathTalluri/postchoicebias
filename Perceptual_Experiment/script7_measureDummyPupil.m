
%% Commitment RDK code
% this models the methodology of Bronfman, Brezis & Usher with random dot
% stimuli. To be used by Ana Vojvodic.

% Anne Urai, March 2014
% ----------------------------------------------------------------

clear all; close all; clc;

% general setup
setup.Eye           = true; % true if using Eyelink

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

% Setup the PsychToolbox
window.dist             = 50; % viewing distance in cm ( this is fixed in the chinrest, D1.09)
window.width            = 40; % physical width of the CRT screen in cm, 40 cm in D1.09
window.skipChecks       = 1; % set to 1 to skip VBL tests and avoid warnings
[window, audio]         = SetupPTB(window); %load all the psychtoolbox things

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
block = 0;
Screen('Flip', window.h);
% setup eyelink, here make sure that the dummy eye is seen
edfFile = ELconfig(window, setup, block);

% workaround
setup.cohlevel = 0;
[dots, fix] = setupDots(window, setup);

% show fixation
Screen('FillRect', window.h, window.black);
window      = drawFixation(window, fix, dots); % fixation
Screen('Flip', window.h);

% record for 1.5 minute to see how much camera jitter there is?
WaitSecs(30);

%% save the EL file for this block

setup.datetime = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
%  fprintf('Receiving data file ''%s''\n', edfFile );
setup.eyefilename = sprintf('Data/P%d_s%d_b%d_%s_dumnyPupil.edf', setup.participant, setup.session, block, setup.datetime);

status = Eyelink('ReceiveFile', edfFile, setup.eyefilename); %this collects the file from the eyelink
disp(status);
disp(['File ' setup.eyefilename ' saved to disk']);

% close the eyelink
Eyelink('StopRecording');
Eyelink('CloseFile');

disp('done!'); Screen('CloseAll'); ShowCursor;
PsychPortAudio('Stop', audio.h);
sca;