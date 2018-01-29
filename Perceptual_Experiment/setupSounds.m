function sound = setupSounds(setup, audio)

% auditory feedback, using some ENS functions (probably from Valentin)
a = [1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 ]; % counterbalance

switch a(setup.participant)
    case 1
        sound.feedback.correct      = CreateTone(880, 0.150, audio.freq); % 150 ms, 880 Hz
        sound.feedback.incorrect    = CreateTone(200 ,0.150, audio.freq); % 150 ms, 200 Hz
        sound.feedback.invalid      = CreateTone(200 ,0.500, audio.freq); % longer error tone
        sound.prompt.choice         = wavread(char('Sosumi.wav'))';
        sound.prompt.nochoice       = wavread(char('Hero.wav'))';
    case 2
        sound.feedback.correct      = CreateTone(200 ,0.150, audio.freq); % 150 ms, 200 Hz
        sound.feedback.incorrect    = CreateTone(880, 0.150, audio.freq); % 150  ms, 880 Hz
        sound.feedback.invalid      = CreateTone(880, 0.500, audio.freq); % longer error tone
        sound.prompt.nochoice       = wavread(char('Sosumi.wav'))';
        sound.prompt.choice         = wavread(char('Hero.wav'))';
end

% concatenate and fill the tone buffer
[sound.tonebuf, sound.tonepos] = CreateAudioBuffer(sound.feedback.correct, ...
    sound.feedback.incorrect, sound.feedback.invalid, sound.prompt.choice, sound.prompt.nochoice);

PsychPortAudio('FillBuffer', audio.h, sound.tonebuf);

end