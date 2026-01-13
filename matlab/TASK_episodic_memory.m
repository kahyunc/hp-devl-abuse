% Episodic memory task for adolescents
% :: TASK - presented using MATLAB R2020b (MathWorks, USA) with Psychophysics Toolbox Version 3 
% :: code written by Kahyun Choi


%%
sca;
close all;


%% input sbj information

% sbj information
sbjDate = input('Date(YYMMDD): ','s');
sbjNum = input('Sbj Number: ','s');
sbjAge = input('Age: ','s');
sbjSex = input('M/F: ','s');

% check if restart
txt_name = './results/*_WCDWWW_diary.txt';
d = dir(txt_name); [~, idx] = max([d.datenum]); recent_txt = d(idx).name;
underbar_idx = strfind(recent_txt, '_');
recent_num = recent_txt(strfind(recent_txt, 'sbj')+3:underbar_idx(2)-1);

isRestart = false;
if str2double(recent_num) == str2double(sbjNum)
    restartQ = input('\nsame subject done before, continue? (1 = yes): ');
    if restartQ == 1
        isRestart = true;
    end
end
if isRestart
    fprintf('\n\n***** Task continues! sbj A%s *****\n', sbjNum);
    fprintf('\nPLEASE MAKE SURE ALL THE LAST RESULTS ARE SAVED\n');
    fprintf(' - if not, STOP running and run the last section of this script now!\n\n');
    
    mat_file_name = ['./results/', sbjDate, '_sbj', sbjNum,'_WCDWWW.mat'];
    if ~exist(mat_file_name, 'file')
        error('no mat file exists');
    end
end

ready = input('\nready? (enter) : ', 's');
if isRestart
    n_restart = n_restart + 1;
    fprintf('\n\n***** sbj A%s - task #%d *****\n\n\n', sbjNum, n_restart);
else
    n_restart = 1;
    clearvars -except sbjDate sbjNum sbjAge sbjSex isRestart n_restart
end


%% preference feature setting
% feature setting

% load subject preference data
file_name='./data_2_learning/WCDI_learning_data_subjID=';

d = dir('./data_2_learning/*.mat');
[~,idx] = max([d.datenum]); % find the last modified file
filename = d(idx).name; % name of file
load(['./data_2_learning/' filename], 'feature_animals_rank');

preferAnimal=feature_animals_rank(1:5);

% randomize animal order and type
preferAnimalNum = preferAnimal;
learnedAnimalNum = [preferAnimalNum(1) preferAnimalNum(3) preferAnimalNum(5)];


%% screen setting 

% ------- screen setup -------
PsychDefaultSetup(2);

Screen('Preference', 'SkipSyncTests', 1);

screenNumber = max(Screen('Screens'));
% screenNumber = 1;

% PsychDebugWindowConfiguration      % debug mode

% color
white = WhiteIndex(screenNumber);   
black = BlackIndex(screenNumber);
grey = white/2;

% open the screen   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black, [], 32, 2);  
% [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black, [0 0 1920 1200], 32, 2);

% window size
[screenXpixels, screenYpixels] = Screen('WindowSize', window); 

% display size
[width_mm, height_mm] = Screen('DisplaySize', screenNumber);

% frame duration
ifi = Screen('GetFlipInterval', window);    

% center coordinate
[xCenter, yCenter] = RectCenter(windowRect);    

% blend function for the screen
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');  

topPriorityLevel = MaxPriority(window);

% number of frames to wait before re-drawing
waitframes = 1;



%% experiment setting

% ------- fixation setup -------
% screen size ratio (standard: 1920*1080)
xRatio = xCenter/960;
yRatio = yCenter/540;

% fixation gif size & position
fixGifCoords = [xCenter-(50*xRatio) yCenter-(50*yRatio) xCenter+(50*xRatio) yCenter+(50*yRatio)];


% ------- interference setup -------
% 8 times in one block, short / long

i_trials = 8;
i_1currBlock = 0;
i_2currBlock = 0;

% 1: left, 2: right
answerfile = './images/times/time_answer.mat';
load(answerfile);   % interference answer load (i_short & i_long)


% ------- keyboard info -------
KbName('UnifyKeyNames');

escapeKey = KbName('ESCAPE');
spaceKey = KbName('space');
signalKey = KbName('5%');

leftKey = KbName('1!');
rightKey = KbName('4$');
selectKey = KbName('3#');
% num2Key = KbName('2@');

RestrictKeysForKbCheck([escapeKey spaceKey signalKey leftKey rightKey selectKey]);


% ------- image & text -------
imageFolder = [cd '/images/'];
textFolder = [cd '/texts/'];

encInstText = imread([textFolder 'encInstText.jpg']);
simInstText = imread([textFolder 'simInstText.jpg']);
retInstText = imread([textFolder 'retInstText.jpg']);
encCInstText = imread([textFolder 'encCInstText.jpg']);
simCInstText = imread([textFolder 'simCInstText.jpg']);
runStartText = imread([textFolder 'runStartText.jpg']);
expEndText = imread([textFolder 'expEndText.jpg']);
interInstText_S = imread([textFolder 'interInstText_S.jpg']);
interInstText_L = imread([textFolder 'interInstText_L.jpg']);
interText_S = imread([textFolder 'interText_S.jpg']);
interText_L = imread([textFolder 'interText_L.jpg']);

[animal1, ~, alpha] = imread([imageFolder 'animal1_fox.png']);
animal1(:, :, 4) = alpha;
[animal2, ~, alpha] = imread([imageFolder 'animal2_rabbit.png']);
animal2(:, :, 4) = alpha;
[animal3, ~, alpha] = imread([imageFolder 'animal3_pig.png']);
animal3(:, :, 4) = alpha;
[animal4, ~, alpha] = imread([imageFolder 'animal4_tiger.png']);
animal4(:, :, 4) = alpha;
[animal5, ~, alpha] = imread([imageFolder 'animal5_dog.png']);
animal5(:, :, 4) = alpha;
[animal6, ~, alpha] = imread([imageFolder 'animal6_panda.png']); % control
animal6(:, :, 4) = alpha;
animalImg = {animal1 animal2 animal3 animal4 animal5};
[selectShade, ~, alpha] = imread([imageFolder 'selected_shade.png']);
selectShade(:, :, 4) = alpha;
    
contextA = imread([imageFolder 'bgnd1_green_far.png']);
contextB = imread([imageFolder 'bgnd2_khaki_far.png']);
contextC = imread([imageFolder 'bgnd3_brown.png']);
contextA_what = imread([imageFolder 'bgnd1_what.png']);
contextB_what = imread([imageFolder 'bgnd2_what.png']);
contextC_what = imread([imageFolder 'bgnd3_what.png']);
sel_L = imread([imageFolder 'select_L.png']);
sel_R = imread([imageFolder 'select_R.png']);

fixImg = imread([imageFolder 'fix.gif'], 'frames', 'all');

for i = 1:64
    tmpImg = imread([imageFolder 'times/' num2str(i) '.jpg']);
    timeImgs{i} = tmpImg;
end


% ------- XY coordinates -------
% ratio for 1920*1200
y1Ratio = yCenter/600;

% encoding animal location coordinates (animal size 200*200 for 1920*1200)
encLocAXY = [xCenter-(672*xRatio) yCenter-(445*y1Ratio) xCenter-(492*xRatio) yCenter-(265*y1Ratio);     xCenter-(422*xRatio) yCenter-(217*y1Ratio) xCenter-(232*xRatio) yCenter-(27*y1Ratio); ...
             xCenter-(531*xRatio) yCenter+(153*y1Ratio) xCenter-(331*xRatio) yCenter+(353*y1Ratio);     xCenter+(90*xRatio) yCenter-(72*y1Ratio) xCenter+(300*xRatio) yCenter+(138*y1Ratio); ...
             xCenter+(456*xRatio) yCenter+(289*y1Ratio) xCenter+(656*xRatio) yCenter+(489*y1Ratio);     xCenter+(463*xRatio) yCenter-(321*y1Ratio) xCenter+(653*xRatio) yCenter-(131*y1Ratio)];   

encLocBXY = [xCenter-(762*xRatio) yCenter-(269*y1Ratio) xCenter-(572*xRatio) yCenter-(79*y1Ratio);      xCenter-(493*xRatio) yCenter+(149*y1Ratio) xCenter-(283*xRatio) yCenter+(359*y1Ratio); ...
             xCenter-(130*xRatio) yCenter+(25*y1Ratio) xCenter+(70*xRatio) yCenter+(225*y1Ratio);       xCenter+(36*xRatio) yCenter-(333*y1Ratio) xCenter+(226*xRatio) yCenter-(143*y1Ratio); ...
             xCenter+(368*xRatio) yCenter-(418*y1Ratio) xCenter+(548*xRatio) yCenter-(238*y1Ratio);     xCenter+(602*xRatio) yCenter-(100*y1Ratio) xCenter+(802*xRatio) yCenter+(100*y1Ratio)];

encLocCXY = [xCenter-(819*xRatio) yCenter-(126*y1Ratio) xCenter-(619*xRatio) yCenter+(74*y1Ratio);      xCenter-(533*xRatio) yCenter-(126*y1Ratio) xCenter-(333*xRatio) yCenter+(74*y1Ratio); ...
             xCenter+(-247*xRatio) yCenter-(126*y1Ratio) xCenter-(47*xRatio) yCenter+(74*y1Ratio);      xCenter+(39*xRatio) yCenter-(126*y1Ratio) xCenter+(239*xRatio) yCenter+(74*y1Ratio)];

% simulation & retrieval animal coordinates (animal size 200*200 for 1920*1080)
retAnimalXY = [xCenter-(710*xRatio) yCenter+(310*yRatio) xCenter-(510*xRatio) yCenter+(510*yRatio); xCenter-(414*xRatio) yCenter+(310*yRatio) xCenter-(214*xRatio) yCenter+(510*yRatio); ...
               xCenter-(100*xRatio) yCenter+(310*yRatio) xCenter+(100*xRatio) yCenter+(510*yRatio); xCenter+(196*xRatio) yCenter+(310*yRatio) xCenter+(396*xRatio) yCenter+(510*yRatio); ...
               xCenter+(510*xRatio) yCenter+(310*yRatio) xCenter+(710*xRatio) yCenter+(510*yRatio)];

% shade location coordinates (shade size 188*47 for 1920*1080)
retLocAXY = [xCenter-(545*xRatio) yCenter-(323*yRatio) xCenter-(365*xRatio) yCenter-(276*yRatio); xCenter-(348*xRatio) yCenter-(153*yRatio) xCenter-(168*xRatio) yCenter-(106*yRatio); ...
             xCenter-(429*xRatio) yCenter+(125*yRatio) xCenter-(249*xRatio) yCenter+(172*yRatio); xCenter+(61*xRatio) yCenter-(28*yRatio) xCenter+(241*xRatio) yCenter+(19*yRatio); ...
             xCenter+(343*xRatio) yCenter+(223*yRatio) xCenter+(523*xRatio) yCenter+(270*yRatio); xCenter+(343*xRatio) yCenter-(226*yRatio) xCenter+(523*xRatio) yCenter-(179*yRatio)];
retLocBXY = [xCenter-(611*xRatio) yCenter-(187*yRatio) xCenter-(431*xRatio) yCenter-(140*yRatio); xCenter-(392*xRatio) yCenter+(133*yRatio) xCenter-(212*xRatio) yCenter+(180*yRatio); ...
             xCenter-(114*xRatio) yCenter+(32*yRatio)  xCenter+(66*xRatio) yCenter+(79*yRatio);   xCenter+(10*xRatio) yCenter-(232*yRatio) xCenter+(190*xRatio) yCenter-(185*yRatio); ...
             xCenter+(267*xRatio) yCenter-(302*yRatio) xCenter+(447*xRatio) yCenter-(255*yRatio); xCenter+(455*xRatio) yCenter-(59*yRatio) xCenter+(635*xRatio) yCenter-(12*yRatio)];
    
% background
encBGCoords = [50*xRatio 50*y1Ratio (2*xCenter)-(50*xRatio) (2*yCenter)-(50*y1Ratio)];
simBGCoords = [50*xRatio 50*y1Ratio (2*xCenter)-(50*xRatio) (2*yCenter)-(50*y1Ratio)];
retBGCoords = [xCenter-(710*xRatio) yCenter-(500*yRatio) xCenter+(710*xRatio) yCenter+(300*yRatio)];

% interference
i_imgCoords = [xCenter-(680*xRatio) yCenter-(510*y1Ratio) xCenter+(680*xRatio) yCenter+(510*y1Ratio)];
i_selImgCoords = [xCenter-(680*xRatio) yCenter+(30*y1Ratio) xCenter+(680*xRatio) yCenter+(540*y1Ratio)];
i_textImgCoords = [xCenter-(960*xRatio) yCenter-(600*y1Ratio) xCenter+(960*xRatio) yCenter-(376*y1Ratio)];


% ------- design matrix -------
if ~isRestart
    rng('shuffle');
    % randomize context/www order
    conOrder = [1 1 2 2];            % 1: A, 2: B
    conOrder = conOrder(randperm(4));

    contextOrder = [];                      % 1: A, 2: B
    for i = 1:size(conOrder,2)
        if conOrder(i) == 1         % A-B
            contextOrder = [contextOrder 1 2];
        elseif conOrder(i) == 2     % B-A
            contextOrder = [contextOrder 2 1];
        end
    end

    % randomize interference order
    i_order = randperm(64);

    i_typeOrder = [1 1 2 2]; % 1: short-long, 2: long-short
    i_typeOrder = i_typeOrder(randperm(4));
    i_type = [];
    for i = 1:size(i_typeOrder,2) % 1: short, 2: long
        if i_typeOrder(i) == 1
            i_type = [i_type 1 2];
        elseif i_typeOrder(i) == 2
            i_type = [i_type 2 1];
        end
    end

    % what case: exclude 1, 2, 3, 4, 5, 1, 3, 5 -> order randomize
    excludeAnimal = [1 2 3 4 5 1 3 5];
    excludeAnimal = excludeAnimal(randperm(8));
    
% ------- save subject information and trial information -------

    Data = []; 
    Data.sbj.sbjDate = sbjDate;
    Data.sbj.sbjNum = sbjNum;
    Data.sbj.sbjAge = sbjAge;
    Data.sbj.sbjSex = sbjSex;
    Data.sbj.preferAnimal = preferAnimalNum;
    Data.sbj.learnedAnimal = learnedAnimalNum;
    Data.trial.trialinfo.context = contextOrder;
    Data.trial.trialinfo.animal = [];
    Data.trial.trialinfo.location = [];
    Data.trial.trialinfo.interference.type = i_type;
    Data.trial.trialinfo.interference.stim = i_order;
    Data.trial.encoding = [];
    Data.trial.retrieval = [];
    Data.trial.interference = [];
    Data.trial.control = [];
    Data.time = [];
    Data.shot = [];
    
else
    mat_file_name = ['./results/', sbjDate, '_sbj', sbjNum,'_WCDWWW.mat'];
    load(mat_file_name);
    
    contextOrder = Data.trial.trialinfo.context;
    i_type = Data.trial.trialinfo.interference.type;
    i_order = Data.trial.trialinfo.interference.stim;
end


if isRestart
    diary_www = ['./results/', sbjDate, '_sbj', sbjNum, '_', num2str(n_restart), '_WCDWWW_diary.txt'];
else
    diary_www = ['./results/', sbjDate, '_sbj', sbjNum,'_WCDWWW_diary.txt'];
end
diary(diary_www);


%% task

if isRestart
    % restart from the last run
    currRun = currRun - 1;
    currTrial = currRun * 2;
    if currTrial == 0
        isControl = 0;
    else
        isControl = 1;
    end
else
    % start default
    currRun = 0;
    currTrial = 0;
    isControl = 0;
end

aniAcc = []; locAcc = []; interAcc = []; interSum99 = [];

isRunning = true;
expStart = tic;

while isRunning

% check for control trials (2 trials -> 1 control trial)
if mod(currTrial,2)     % odd
    isControl = 0;
else                    % even
    if currTrial > 0
        if isControl
            isControl = 0;
        else
            isControl = 1;
        end
    end
end

% trial & run ++
if ~isControl
    currTrial = currTrial + 1;

    if mod(currTrial,2)     % odd
        currRun = currRun + 1;
    end
end
   

%% encoding
% total 8 trials, 2 trials for each run

% ------- encoding setting -------

if ~isControl
    % animal select
    A = [1 2 3 4 5]; B = excludeAnimal(currTrial);
    selAnimalNum = setdiff(A, B);
    selAnimalNum = selAnimalNum(randperm(4));

    % location random select
    C = [1 2 3 4 5 6]; 
    C = C(randperm(6));
    selLocationNum = C(1:4);
    
    % currTrial context www order
    currContext = contextOrder(currTrial);
    
    % encoding setting save    
    Data.trial.trialinfo.animal = [Data.trial.trialinfo.animal; selAnimalNum];
    Data.trial.trialinfo.location = [Data.trial.trialinfo.location; selLocationNum];
end


% ------- instruction text -------

% Run starts: check sbj's condition (press space bar to start)

if mod(currTrial,2) && ~isControl
    runStartTextT = Screen('MakeTexture', window, runStartText);
    Screen('DrawTexture', window, runStartTextT, [], [0 0 2*xCenter 2*yCenter], 0);    
    Screen('Flip', window);

    waitPress = true; waitShot = true; 
    while waitPress
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(spaceKey)
            waitPress = false;
            fprintf('\nrun starts... waiting for the signal input...\n\n');
            while waitShot
                [keyIsDown, secs, keyCode] = KbCheck;
                if keyCode(signalKey)   % start run after shot signal
                    runStart = tic;
                    shotTime = GetSecs;
                    Data.time(end+1).run = currRun;
                    Data.time(end).trial = currTrial;
                    Data.time(end).control = isControl;
                    Data.time(end).label = 'run start';
                    Data.time(end).time = shotTime;
                    Data.time(end).Time = toc(expStart);
                    Data.time(end).runTime = toc(runStart);
                    
                    Data.shot(end+1).run = currRun;
                    Data.shot(end).label = 'start';
                    Data.shot(end).time = shotTime;
                    
                    fprintf('signal in!\n\n');
                    waitShot = false;
                elseif keyCode(escapeKey)
                    sca;
                    waitPress = false; isRunning = false;
                    fprintf('\n*** Experiment terminated ***\n');
                    fprintf('Current run: %d, Current trial: %d, Control: %d\n', currRun, currTrial, isControl);
                    fprintf('Pause timing: Before run starts\n\n');
                    return
                end
            end
        elseif keyCode(escapeKey)
            isPaused = true;
            while keyIsDown
                [keyIsDown, ~, ~] = KbCheck;
            end
            fprintf('\n* Press space bar to continue ** Press ESC to stop *\n\n');
            while isPaused
                [keyIsDown, secs, keyCode] = KbCheck;
                if keyCode(spaceKey)
                    isPaused = false; waitPress = false;
                    fprintf('continue...\n\n');
                elseif keyCode(escapeKey)
                    sca;
                    waitPress = false; isRunning = false;
                    fprintf('\n*** Experiment terminated ***\n');
                    fprintf('Current run: %d, Current trial: %d, Control: %d\n', currRun, currTrial, isControl);
                    fprintf('Pause timing: Before run starts\n\n');
                    return
                end
            end
        elseif keyCode(signalKey)
            fprintf('S\n');
            while keyIsDown
                [keyIsDown, ~, ~] = KbCheck;
            end
            isFinal = false;
            while ~isFinal
                shotTime = GetSecs;
                Data.shot(end+1).run = currRun;
                Data.shot(end).label = 'end';
                Data.shot(end).time = shotTime;
                isFinal = true;
            end
        end
    end
end

if ~isControl
    fprintf('\n- Trial %d Start -\n\n', currTrial);
    fprintf('** Encoding **\n');
    fprintf('Current run: %d, Current trial: %d, Control: %d\n', currRun, currTrial, isControl);
    fprintf('Animal: %s, Location: %s, Context: %d\n\n', mat2str(selAnimalNum), mat2str(selLocationNum), currContext);
elseif isControl
    fprintf('** Encoding **\n');
    fprintf('Current run: %d, Current trial: %d, Control: %d\n\n', currRun, currTrial, isControl);
end

% instruction text

if ~isControl
    encInstTextT = Screen('MakeTexture', window, encInstText);
    Screen('DrawTexture', window, encInstTextT, [], [0 0 2*xCenter 2*yCenter], 0);
    Screen('Flip', window);
    
elseif isControl
    encCInstTextT = Screen('MakeTexture', window, encCInstText);
    Screen('DrawTexture', window, encCInstTextT, [], [0 0 2*xCenter 2*yCenter], 0);
    Screen('Flip', window);
end

Data.time(end+1).run = currRun;
Data.time(end).trial = currTrial;
Data.time(end).control = isControl;
Data.time(end).label = 'enc instruction';
Data.time(end).time = GetSecs;
Data.time(end).Time = toc(expStart);
Data.time(end).runTime = toc(runStart);

% 3s
waitTime = 0; tic;
while waitTime < 3
    waitTime = toc;
    
    [keyIsDown, secs, keyCode] = KbCheck;
    if keyCode(escapeKey)
        isPaused = true;
        while keyIsDown
            [keyIsDown, ~, ~] = KbCheck;
        end
        
        Data.time(end+1).run = currRun;
        Data.time(end).trial = currTrial;
        Data.time(end).control = isControl;
        Data.time(end).label = 'PAUSE';
        Data.time(end).time = GetSecs;
        Data.time(end).Time = toc(expStart);
        Data.time(end).runTime = toc(runStart);
        fprintf('\n*** Paused ***\nPRESS SPACE BAR TO CONTINUE, PRESS ESC TO STOP\n\n');
        
        while isPaused
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(spaceKey)
                isPaused = false; tic; % restart 3s
                fprintf('continue...\n\n');
            elseif keyCode(escapeKey)
                sca;
                isRunning = false;
                fprintf('\n*** Experiment terminated ***\n');
                fprintf('Current run: %d, Current trial: %d, Control: %d\n', currRun, currTrial, isControl);
                fprintf('Pause timing: Before encoding\n\n');
                return
            end
        end
    end
end


% ------- fixation -------

% 7s
fixImgT = Screen('MakeTexture', window, fixImg(:,:,:,1));
Screen('DrawTexture', window, fixImgT, [], fixGifCoords, 0);
Screen('Flip', window);

Data.time(end+1).run = currRun;
Data.time(end).trial = currTrial;
Data.time(end).control = isControl;
Data.time(end).label = 'enc fixation';
Data.time(end).time = GetSecs;
Data.time(end).Time = toc(expStart);
Data.time(end).runTime = toc(runStart);

waitTime = 0; tStart = tic; fixFrame = 1;
while waitTime < 7
    waitTime = toc(tStart);
    waitImgUpdate = true; tic;
    while waitImgUpdate
        waitImgTime = toc;
        waitTime = toc(tStart);
        if waitTime > 7
            break
        end
        if waitImgTime > 0.05
            if fixFrame > size(fixImg,4)-1
                fixFrame = 1;
            else
                fixFrame = fixFrame + 1;
            end
            fixImgT = Screen('MakeTexture', window, fixImg(:,:,:,fixFrame));
            Screen('DrawTexture', window, fixImgT, [], fixGifCoords, 0);
            Screen('Flip', window);
            waitImgUpdate = false;
        end
    end
end


% ------- encoding animation ------- 

Data.time(end+1).run = currRun;
Data.time(end).trial = currTrial;
Data.time(end).control = isControl;
Data.time(end).label = 'enc start';
Data.time(end).time = GetSecs;
Data.time(end).Time = toc(expStart);
Data.time(end).runTime = toc(runStart);


if ~isControl
    % background 1.5s
    if currContext == 1         % A
        % animal location of context A
        locationXY = encLocAXY; 

        % background only (1.5s)
        encBGImgT = Screen('MakeTexture', window, contextA);
        Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0);

    elseif currContext == 2     % B
        % animal location of context B
        locationXY = encLocBXY;

        % background only (1.5s)
        encBGImgT = Screen('MakeTexture', window, contextB);
        Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0);
    end
    
    Screen('Flip', window);

    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'enc start bg';
    Data.time(end).time = GetSecs;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

    WaitSecs(1.5);
    
    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'enc animation';
    Data.time(end).time = GetSecs;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

    % animal on (0.05s) -> on (1.9s) -> on (0.05s) -> off (0.5s)
    for i = 1:size(selAnimalNum,2)
        encCurrAnimal = animalImg{selAnimalNum(i)};
        encCurrAnimalT = Screen('MakeTexture', window, encCurrAnimal);

        % animal location & location for animation
        encCurrLoc = locationXY(selLocationNum(i),:);
        encCurrLoc_pre = encCurrLoc;
        encCurrLoc_pre(2) = encCurrLoc_pre(2)+(10*yRatio);
        encCurrLoc_pre(4) = encCurrLoc_pre(4)+(10*yRatio);
        
        % 0.05s
        Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0); % background image
        Screen('DrawTexture', window, encCurrAnimalT, [], encCurrLoc_pre, 0); % curr animal image
        Screen('Flip', window);
        
        Data.trial.encoding(end+1).run = currRun;
        Data.trial.encoding(end).trial = currTrial;
        Data.trial.encoding(end).control = isControl;
        Data.trial.encoding(end).animal = selAnimalNum(i);
        Data.trial.encoding(end).location = selLocationNum(i);
        Data.trial.encoding(end).time = GetSecs;
        Data.trial.encoding(end).label = 'on';
        WaitSecs(0.05);
        
        % 1.9s
        Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0); % background image
        Screen('DrawTexture', window, encCurrAnimalT, [], encCurrLoc, 0); % curr animal image
        Screen('Flip', window);
        WaitSecs(1.9);

        % 0.05s
        Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0); % background image
        Screen('DrawTexture', window, encCurrAnimalT, [], encCurrLoc_pre, 0); % curr animal image
        Screen('Flip', window);
        WaitSecs(0.05);
        
        % 0.5s
        Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0); % background image
        Screen('Flip', window);
        
        Data.trial.encoding(end+1).run = currRun;
        Data.trial.encoding(end).trial = currTrial;
        Data.trial.encoding(end).control = isControl;
        Data.trial.encoding(end).animal = selAnimalNum(i);
        Data.trial.encoding(end).location = selLocationNum(i);
        Data.trial.encoding(end).time = GetSecs;
        Data.trial.encoding(end).label = 'off';
        WaitSecs(0.5);
    end

    % background 0.5s
    Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0); % background image
    Screen('Flip', window);
    WaitSecs(0.5);

elseif isControl
    
    % animal location of control context
    locationXY = encLocCXY;
    
    % background only (1.5s)
    encBGImgT = Screen('MakeTexture', window, contextC);
    Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0);
    Screen('Flip', window);
    
    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'enc start bg';
    Data.time(end).time = GetSecs;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

    WaitSecs(1.5);

    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'enc animation';
    Data.time(end).time = GetSecs;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

    % animal on (0.05s) -> on (1.9s) -> on (0.05s) -> off (0.5s)
    for i = 1:size(locationXY,1)
        ctrlAnimal = animal6;
        ctrlAnimalT = Screen('MakeTexture', window, ctrlAnimal);

        encCurrLoc = locationXY(i,:);
        encCurrLoc_pre = encCurrLoc;
        encCurrLoc_pre(2) = encCurrLoc_pre(2)+(10*yRatio);
        encCurrLoc_pre(4) = encCurrLoc_pre(4)+(10*yRatio);
        
        % 0.05s
        Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0); % background image
        Screen('DrawTexture', window, ctrlAnimalT, [], encCurrLoc_pre, 0); % curr animal image
        Screen('Flip', window);
        
        Data.trial.control(end+1).run = currRun;
        Data.trial.control(end).time = GetSecs;
        Data.trial.control(end).label = 'enc animal on';
        WaitSecs(0.05);
        
        % 1.9s
        Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0); % background image
        Screen('DrawTexture', window, ctrlAnimalT, [], encCurrLoc, 0); % curr animal image
        Screen('Flip', window);
        WaitSecs(1.9);

        % 0.05s
        Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0); % background image
        Screen('DrawTexture', window, ctrlAnimalT, [], encCurrLoc_pre, 0); % curr animal image
        Screen('Flip', window);
        WaitSecs(0.05);
        
        % 0.5s
        Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0); % background image
        Screen('Flip', window);
        
        Data.trial.control(end+1).run = currRun;
        Data.trial.control(end).time = GetSecs;
        Data.trial.control(end).label = 'enc animal off';
        WaitSecs(0.5);
    end

    % background 0.5s
    Screen('DrawTexture', window, encBGImgT, [], encBGCoords, 0); % background image
    Screen('Flip', window);
    WaitSecs(0.5);

end

Data.time(end+1).run = currRun;
Data.time(end).trial = currTrial;
Data.time(end).control = isControl;
Data.time(end).label = 'enc end';
Data.time(end).time = GetSecs;
Data.time(end).Time = toc(expStart);
Data.time(end).runTime = toc(runStart);


%% interference

% short & long

if ~isControl
    
    % current trial info
    if i_type(currTrial) == 1
        currType = 1;   % short
    elseif i_type(currTrial) == 2
        currType = 2;   % long
    end
    
    if currType == 1
        currTypeT = 'short';
        i_answer = i_short;
        i_1currBlock = i_1currBlock + 1;
        i_currBlock = i_1currBlock;

    elseif currType == 2
        currTypeT = 'long';
        i_answer = i_long;
        i_2currBlock = i_2currBlock + 1;
        i_currBlock = i_2currBlock;
    end

        
    % ------- instruction text -------

    if currType == 1        % short
        interInstTextT = Screen('MakeTexture', window, interInstText_S);
        interTextT = Screen('MakeTexture', window, interText_S);
    elseif currType == 2    % long
        interInstTextT = Screen('MakeTexture', window, interInstText_L);
        interTextT = Screen('MakeTexture', window, interText_L);
    end
    Screen('DrawTexture', window, interInstTextT, [], [0 0 2*xCenter 2*yCenter], 0);
    Screen('Flip', window);
       
    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'int instruction';
    Data.time(end).time = GetSecs;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

    waitTime = 0; tic;
    while waitTime < 3      % inst text 3s
        waitTime = toc;

        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            isPaused = true;
            while keyIsDown
                [keyIsDown, ~, ~] = KbCheck;
            end

            Data.time(end+1).run = currRun;
            Data.time(end).trial = currTrial;
            Data.time(end).control = isControl;
            Data.time(end).label = 'PAUSE';
            Data.time(end).time = GetSecs;
            Data.time(end).Time = toc(expStart);
            Data.time(end).runTime = toc(runStart);
            fprintf('\n*** Paused ***\nPRESS SPACE BAR TO CONTINUE, PRESS ESC TO STOP\n\n');

            while isPaused
                [keyIsDown, secs, keyCode] = KbCheck;
                if keyCode(spaceKey)
                    isPaused = false; tic; % restart 4s
                    fprintf('continue...\n\n');
                elseif keyCode(escapeKey)
                    sca;
                    isRunning = false;
                    fprintf('\n*** Experiment terminated ***\n');
                    fprintf('Current run: %d, Current trial: %d, Control: %d\n', currRun, currTrial, isControl);
                    fprintf('Pause timing: Before interference\n\n');
                    return
                end
            end
        end           
    end

    fprintf('** interference **\n');
    fprintf('Current run: %d, Current trial: %d\n', currRun, currTrial);
    fprintf('Condition: %s #%d\n\n', currTypeT, i_currBlock);
    
    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'int start';
    Data.time(end).time = GetSecs;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

% ------- fixation -------
    
    Screen('DrawTexture', window, fixImgT, [], fixGifCoords, 0);
    Screen('Flip', window);

    waitTime = 0; tStart = tic; fixFrame = 1;
    while waitTime < 1
        waitTime = toc(tStart);
        waitImgUpdate = true; tic;
        while waitImgUpdate
            waitImgTime = toc;
            waitTime = toc(tStart);
            if waitTime > 1
                break
            end
            if waitImgTime > 0.05
                if fixFrame > size(fixImg,4)-1
                    fixFrame = 1;
                else
                    fixFrame = fixFrame + 1;
                end
                fixImgT = Screen('MakeTexture', window, fixImg(:,:,:,fixFrame));
                Screen('DrawTexture', window, fixImgT, [], fixGifCoords, 0);
                Screen('Flip', window);
                waitImgUpdate = false;
            end
        end
    end

% ------- interference task -------

    % stimulus 4s, answer 0.5s -> total 8 trials in one block

    i_response = [];
    
    for i_currTrial = 1:i_trials
        
        btnPress = 99; currRT = 99; currCorrect = 99;
        
        % stimulus setting
        i_currStim = i_order(8*(currTrial - 1) + i_currTrial);   
        i_currImg = rgb2gray(timeImgs{i_currStim});
        i_currImgT = Screen('MakeTexture', window, i_currImg);
        
        i_sel_L_T = Screen('MakeTexture', window, sel_L);
        i_sel_R_T = Screen('MakeTexture', window, sel_R);
        
        % - stimulus 4s -
        Screen('DrawTexture', window, i_currImgT, [], i_imgCoords, 0);
        Screen('DrawTexture', window, interTextT, [], i_textImgCoords, 0);
        Screen('Flip', window);
        
        Data.trial.interference(end+1).trial = currTrial;
        Data.trial.interference(end).order = i_currTrial;
        Data.trial.interference(end).type = currType;
        Data.trial.interference(end).stimulus = i_currStim;
        Data.trial.interference(end).onTime = GetSecs;
        
        % key input
        tStart = GetSecs; inputEnd = 0; keyPressed = 0;
        
        waitTime = 0; tic;
        while waitTime < 4
            waitTime = toc;
            [keyIsDown, secs, keyCode] = KbCheck;
            if ~keyPressed
                if keyCode(leftKey)
                    currRT = secs - tStart;
                    Screen('DrawTexture', window, i_currImgT, [], i_imgCoords, 0);
                    Screen('DrawTexture', window, interTextT, [], i_textImgCoords, 0);
                    Screen('DrawTexture', window, i_sel_L_T, [], i_selImgCoords, 0);
                    Screen('Flip', window);
                    WaitSecs(0.5);

                    btnPress = 1;
                    keyPressed = 1;
                    break;
                    
                elseif keyCode(rightKey)
                    currRT = secs - tStart;
                    Screen('DrawTexture', window, i_currImgT, [], i_imgCoords, 0);
                    Screen('DrawTexture', window, interTextT, [], i_textImgCoords, 0);
                    Screen('DrawTexture', window, i_sel_R_T, [], i_selImgCoords, 0);
                    Screen('Flip', window);
                    WaitSecs(0.5);
                    
                    btnPress = 2;
                    keyPressed = 1;
                    break;
                end
            end
        end
        
        if btnPress == 99
            currCorrect = 99;
        else
            if i_answer(i_currStim) == 1 && btnPress == 1
                currCorrect = 1;
            elseif i_answer(i_currStim) == 2 && btnPress == 2
                currCorrect = 1;
            else
                currCorrect = 0;
            end
        end
        
        Data.trial.interference(end).offTime = GetSecs;
        Data.trial.interference(end).response = btnPress;
        Data.trial.interference(end).RT = currRT;
        Data.trial.interference(end).correct = currCorrect;
        
        fprintf('Current trial: %d, Word number: %d\n', i_currTrial, i_currStim);
        fprintf('Response: %d, RT: %.4f, Correct: %d\n\n', btnPress, currRT, currCorrect);
    end
    
elseif isControl
    
    % fixation 3s
    fixImgT = Screen('MakeTexture', window, fixImg(:,:,:,1));
    Screen('DrawTexture', window, fixImgT, [], fixGifCoords, 0);
    Screen('Flip', window);

    waitTime = 0; tStart = tic; fixFrame = 1;
    while waitTime < 3
        waitTime = toc(tStart);
        waitImgUpdate = true; tic;
        while waitImgUpdate
            waitImgTime = toc;
            waitTime = toc(tStart);
            if waitTime > 3
                break
            end
            if waitImgTime > 0.05
                if fixFrame > size(fixImg,4)-1
                    fixFrame = 1;
                else
                    fixFrame = fixFrame + 1;
                end
                fixImgT = Screen('MakeTexture', window, fixImg(:,:,:,fixFrame));
                Screen('DrawTexture', window, fixImgT, [], fixGifCoords, 0);
                Screen('Flip', window);
                waitImgUpdate = false;
            end
        end
    end
end


%% retrieval

% ------- instruction text (simulation) -------

if ~isControl
    simInstTextT = Screen('MakeTexture', window, simInstText);
    Screen('DrawTexture', window, simInstTextT, [], [0 0 2*xCenter 2*yCenter], 0);
    Screen('Flip', window);
    
elseif isControl
    simCInstTextT = Screen('MakeTexture', window, simCInstText);
    Screen('DrawTexture', window, simCInstTextT, [], [0 0 2*xCenter 2*yCenter], 0);
    Screen('Flip', window);
end

Data.time(end+1).run = currRun;
Data.time(end).trial = currTrial;
Data.time(end).control = isControl;
Data.time(end).label = 'sim instruction';
Data.time(end).time = GetSecs;
Data.time(end).Time = toc(expStart);
Data.time(end).runTime = toc(runStart);

waitTime = 0; tic;
while waitTime < 3
    waitTime = toc;
    
    [keyIsDown, secs, keyCode] = KbCheck;
    if keyCode(escapeKey)
        isPaused = true;
        while keyIsDown
            [keyIsDown, ~, ~] = KbCheck;
        end

        Data.time(end+1).run = currRun;
        Data.time(end).trial = currTrial;
        Data.time(end).control = isControl;
        Data.time(end).label = 'PAUSE';
        Data.time(end).time = GetSecs;
        Data.time(end).Time = toc(expStart);
        Data.time(end).runTime = toc(runStart);
        fprintf('\n*** Paused ***\nPRESS SPACE BAR TO CONTINUE, PRESS ESC TO STOP\n\n');

        while isPaused
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(spaceKey)
                isPaused = false; tic; % restart 3s
                fprintf('continue...\n\n');
            elseif keyCode(escapeKey)
                sca;
                isRunning = false;
                fprintf('\n*** Experiment terminated ***\n');
                fprintf('Current run: %d, Current trial: %d, Control: %d\n', currRun, currTrial, isControl);
                fprintf('Pause timing: Before retrieval\n\n');
                return
            end
        end
    end
end
fprintf('** retrieval **\n');
fprintf('Current run: %d, Current trial: %d, Control: %d\n\n', currRun, currTrial, isControl);


% ------- fixation -------

% 7s
fixImgT = Screen('MakeTexture', window, fixImg(:,:,:,1));
Screen('DrawTexture', window, fixImgT, [], fixGifCoords, 0);
Screen('Flip', window);

Data.time(end+1).run = currRun;
Data.time(end).trial = currTrial;
Data.time(end).control = isControl;
Data.time(end).label = 'sim fixation';
Data.time(end).time = GetSecs;
Data.time(end).Time = toc(expStart);
Data.time(end).runTime = toc(runStart);

waitTime = 0; tStart = tic; fixFrame = 1;
while waitTime < 7
    waitTime = toc(tStart);
    waitImgUpdate = true; tic;
    while waitImgUpdate
        waitImgTime = toc;
        waitTime = toc(tStart);
        if waitTime > 7
            break
        end
        if waitImgTime > 0.05
            if fixFrame > size(fixImg,4)-1
                fixFrame = 1;
            else
                fixFrame = fixFrame + 1;
            end
            fixImgT = Screen('MakeTexture', window, fixImg(:,:,:,fixFrame));
            Screen('DrawTexture', window, fixImgT, [], fixGifCoords, 0);
            Screen('Flip', window);
            waitImgUpdate = false;
        end
    end
end


% ------- simulation -------

% 10s

if ~isControl   
    fprintf('** retrieval **\n');
    fprintf('Context: %d\n\n', currContext);
    
    A = [1 2 3 4 5];
    whatAnimalOrder = A(randperm(5));
            
    if currContext == 1         % context A
        simImgT = Screen('MakeTexture', window, contextA_what);
        Screen('DrawTexture', window, simImgT, [], simBGCoords, 0);

    elseif currContext == 2     % context B
        simImgT = Screen('MakeTexture', window, contextB_what);
        Screen('DrawTexture', window, simImgT, [], simBGCoords, 0);
    end
    Screen('Flip', window);

    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'sim start';
    Data.time(end).time = GetSecs;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

    WaitSecs(10);
    
    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'sim end';
    Data.time(end).time = GetSecs;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

elseif isControl
    
    fprintf('** retrieval **\n\n');
    
    simImgT = Screen('MakeTexture', window, contextC_what);
    Screen('DrawTexture', window, simImgT, [], simBGCoords, 0);    
    Screen('Flip', window);
    
    tmpTime = GetSecs;
    
    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'sim start';
    Data.time(end).time = tmpTime;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

    Data.trial.control(end+1).run = currRun;
    Data.trial.control(end).time = tmpTime;
    Data.trial.control(end).label = 'sim start';
    
    WaitSecs(10);
    
    tmpTime = GetSecs;
    
    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'sim end';
    Data.time(end).time = tmpTime;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

    Data.trial.control(end+1).run = currRun;
    Data.trial.control(end).time = tmpTime;
    Data.trial.control(end).label = 'sim end';
end


% ------- selection (except for control) -------

if ~isControl
    
    % - instruction text -
    % 3s
    retInstTextT = Screen('MakeTexture', window, retInstText);
    Screen('DrawTexture', window, retInstTextT, [], [0 0 2*xCenter 2*yCenter], 0);
    Screen('Flip', window);
    
    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'ret instruction';
    Data.time(end).time = GetSecs;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

    waitTime = 0; tic;
    while waitTime < 3
        waitTime = toc;

        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            isPaused = true;
            while keyIsDown
                [keyIsDown, ~, ~] = KbCheck;
            end

            Data.time(end+1).run = currRun;
            Data.time(end).trial = currTrial;
            Data.time(end).control = isControl;
            Data.time(end).label = 'PAUSE';
            Data.time(end).time = GetSecs;
            Data.time(end).Time = toc(expStart);
            Data.time(end).runTime = toc(runStart);
            fprintf('\n*** Paused ***\nPRESS SPACE BAR TO CONTINUE, PRESS ESC TO STOP\n\n');

            while isPaused
                [keyIsDown, secs, keyCode] = KbCheck;
                if keyCode(spaceKey)
                    isPaused = false; tic; % restart 3s
                    fprintf('continue...\n\n');
                elseif keyCode(escapeKey)
                    sca;
                    isRunning = false;
                    fprintf('\n*** Experiment terminated ***\n');
                    fprintf('Current run: %d, Current trial: %d, Control: %d\n', currRun, currTrial, isControl);
                    fprintf('Pause timing: Before reenactment\n\n');
                    return
                end
            end
        end
    end
    
    % - selection setting - 
    
    % shade image
    selShadeImg = selectShade;
    selShadeImgT = Screen('MakeTexture', window, selShadeImg);

    % draw background
    if currContext == 1         % A
        retLocationXY = retLocAXY;
        retBGImgT = Screen('MakeTexture', window, contextA);
        Screen('DrawTexture', window, retBGImgT, [], retBGCoords, 0);
        
    elseif currContext == 2     % B
        retLocationXY = retLocBXY;
        retBGImgT = Screen('MakeTexture', window, contextB);
        Screen('DrawTexture', window, retBGImgT, [], retBGCoords, 0);
    end
    
    % draw animal image (whatAnimalOrder)
    for i = 1:length(whatAnimalOrder)
        whatAnimalImg = animalImg{whatAnimalOrder(i)};
        whatAnimalImgT = Screen('MakeTexture', window, whatAnimalImg);
        Screen('DrawTexture', window, whatAnimalImgT, [], retAnimalXY(i,:), 0);
    end
    
    Screen('Flip', window);
    
    fprintf('** reenactment **/n');

    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'ret start';
    Data.time(end).time = GetSecs;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

    selectOrder = 1;
    while selectOrder < 5
        
        % - animal select -
        
        currSelAnimal = randi(5,1);
        
        Data.time(end+1).run = currRun;
        Data.time(end).trial = currTrial;
        Data.time(end).control = isControl;
        Data.time(end).label = 'select start';
        Data.time(end).time = GetSecs;
        Data.time(end).Time = toc(expStart);
        Data.time(end).runTime = toc(runStart);

        selectKeyPress = false;
        while ~selectKeyPress
            
            % current selected animal
            if currSelAnimal < 1
                currSelAnimal = currSelAnimal + 5;
            elseif currSelAnimal > 5
                currSelAnimal = currSelAnimal - 5;
            end
            
            % background & animals
            Screen('DrawTexture', window, retBGImgT, [], retBGCoords, 0);
            for i = 1:length(whatAnimalOrder)
                whatAnimalImg = animalImg{whatAnimalOrder(i)};
                whatAnimalImgT = Screen('MakeTexture', window, whatAnimalImg);
                
                if i == currSelAnimal
                    selAnimalXY = retAnimalXY(i, :);
                    selAnimalXY(1) = selAnimalXY(1)-20;
                    selAnimalXY(2) = selAnimalXY(2)-20;
                    selAnimalXY(3) = selAnimalXY(3)+20;
                    selAnimalXY(4) = selAnimalXY(4)+20;                 
                else
                    selAnimalXY = retAnimalXY(i, :);
                end
                Screen('DrawTexture', window, whatAnimalImgT, [], selAnimalXY, 0);
            end

            Screen('Flip', window);
            
            % key press
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(selectKey)
                while keyIsDown
                    [keyIsDown, ~, ~] = KbCheck;
                end
                
                selectedAnimal = whatAnimalOrder(currSelAnimal);
                
                Data.trial.retrieval(end+1).run = currRun;
                Data.trial.retrieval(end).trial = currTrial;
                Data.trial.retrieval(end).order = selectOrder;
                Data.trial.retrieval(end).animal = selectedAnimal;
                Data.trial.retrieval(end).animalT = GetSecs;
                
                selectKeyPress = true;

            elseif keyCode(leftKey)
                while keyIsDown
                    [keyIsDown, ~, ~] = KbCheck;
                end
                
                currSelAnimal = currSelAnimal - 1;
                
            elseif keyCode(rightKey)
                while keyIsDown
                    [keyIsDown, ~, ~] = KbCheck;
                end
                
                currSelAnimal = currSelAnimal + 1;
            end
        end
        
        if selectedAnimal == selAnimalNum(selectOrder)
            selAnimalAns = 1;
        else
            selAnimalAns = 0;
        end
      
        Data.trial.retrieval(end).animalEnc = selAnimalNum(selectOrder);
        Data.trial.retrieval(end).animalAcc = selAnimalAns;

        fprintf('Order: %d\nAnimal: %d, Response: %d, Correct: %d\n', selectOrder, selAnimalNum(selectOrder), selectedAnimal, selAnimalAns);
        
        % - location select -

        currSelLoc = randi(6,1);
        
        selectKeyPress = false;
        while ~selectKeyPress
            
            % current selected shade
            if currSelLoc < 1
                currSelLoc = currSelLoc + 6;
            elseif currSelLoc > 6
                currSelLoc = currSelLoc - 6;
            end
            
            % background & animals
            Screen('DrawTexture', window, retBGImgT, [], retBGCoords, 0);
            
            for i = 1:length(whatAnimalOrder)
                whatAnimalImg = animalImg{whatAnimalOrder(i)};
                whatAnimalImgT = Screen('MakeTexture', window, whatAnimalImg);
                Screen('DrawTexture', window, whatAnimalImgT, [], retAnimalXY(i,:), 0);
            end
            
            Screen('DrawTexture', window, selShadeImgT, [], retLocationXY(currSelLoc, :), 0); % red shade
            Screen('Flip', window);
            
            % key press
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(selectKey)
                while keyIsDown
                    [keyIsDown, ~, ~] = KbCheck;
                end
                
                selectedLoc = currSelLoc;
                
                Data.trial.retrieval(end).location = selectedLoc;
                Data.trial.retrieval(end).locationT = GetSecs;
                
                selectKeyPress = true;

            elseif keyCode(leftKey)
                while keyIsDown
                    [keyIsDown, ~, ~] = KbCheck;
                end
                
                currSelLoc = currSelLoc - 1;
                
            elseif keyCode(rightKey)
                while keyIsDown
                    [keyIsDown, ~, ~] = KbCheck;
                end
                
                currSelLoc = currSelLoc + 1;
            end            
        end
        
        if selectedLoc == selLocationNum(selectOrder)
            selLocationAns = 1;
        else
            selLocationAns = 0;
        end
        
        Data.trial.retrieval(end).locationEnc = selLocationNum(selectOrder);
        Data.trial.retrieval(end).locationAcc = selLocationAns;

        fprintf('Location: %d, Response: %d, Correct: %d\n\n', selLocationNum(selectOrder), selectedLoc, selLocationAns);

        Data.time(end+1).run = currRun;
        Data.time(end).trial = currTrial;
        Data.time(end).control = isControl;
        Data.time(end).label = 'select end';
        Data.time(end).time = GetSecs;
        Data.time(end).Time = toc(expStart);
        Data.time(end).runTime = toc(runStart);

        % - show selected motion -
        
        % wait (0.5s)
        Screen('DrawTexture', window, retBGImgT, [], retBGCoords, 0);

        for i = 1:length(whatAnimalOrder)
            whatAnimalImg = animalImg{whatAnimalOrder(i)};
            whatAnimalImgT = Screen('MakeTexture', window, whatAnimalImg);
            Screen('DrawTexture', window, whatAnimalImgT, [], retAnimalXY(i,:), 0);
        end
        Screen('Flip', window);
        WaitSecs(0.5);
        
        % on (1s)
        Screen('DrawTexture', window, retBGImgT, [], retBGCoords, 0);

        for i = 1:length(whatAnimalOrder)
            whatAnimalImg = animalImg{whatAnimalOrder(i)};
            whatAnimalImgT = Screen('MakeTexture', window, whatAnimalImg);
            Screen('DrawTexture', window, whatAnimalImgT, [], retAnimalXY(i,:), 0);
        end
        
        selAnimalImg = animalImg{selectedAnimal};
        selAnimalImgT = Screen('MakeTexture', window, selAnimalImg);
        
        selLocationXY = retLocationXY(selectedLoc, :);
        selLocationXY(1) = selLocationXY(1)+(20*xRatio);
        selLocationXY(2) = selLocationXY(2)-(120*yRatio);
        selLocationXY(3) = selLocationXY(1)+(148*xRatio);
        selLocationXY(4) = selLocationXY(2)+(148*yRatio);
        
        Screen('DrawTexture', window, selAnimalImgT, [], selLocationXY, 0);
        Screen('Flip', window);

        Data.trial.retrieval(end).motionOn = GetSecs;
        WaitSecs(1);
        
        % off (0.5s)
        Screen('DrawTexture', window, retBGImgT, [], retBGCoords, 0);

        for i = 1:length(whatAnimalOrder)
            whatAnimalImg = animalImg{whatAnimalOrder(i)};
            whatAnimalImgT = Screen('MakeTexture', window, whatAnimalImg);
            Screen('DrawTexture', window, whatAnimalImgT, [], retAnimalXY(i,:), 0);
        end
        Screen('Flip', window);
        
        Data.time(end+1).run = currRun;
        Data.time(end).trial = currTrial;
        Data.time(end).control = isControl;
        Data.time(end).label = 'motion end';
        Data.time(end).time = GetSecs;
        Data.time(end).Time = toc(expStart);
        Data.time(end).runTime = toc(runStart);
        Data.trial.retrieval(end).motionOff = GetSecs;

        WaitSecs(0.5);
        
        % select next order
        selectOrder = selectOrder + 1;
    end
    
    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'ret end';
    Data.time(end).time = GetSecs;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);
end   

% trial end

Data.time(end+1).run = currRun;
Data.time(end).trial = currTrial;
Data.time(end).control = isControl;
Data.time(end).label = 'trial end';
Data.time(end).time = GetSecs;
Data.time(end).Time = toc(expStart);
Data.time(end).runTime = toc(runStart);

if ~isControl
    % trial accuracy calculation
    anians = []; 
    locans = []; 
    interans = [];
    for trialNum = 1:size(Data.trial.retrieval,2)
        if Data.trial.retrieval(trialNum).trial == currTrial
            anians = [anians Data.trial.retrieval(trialNum).animalAcc];
            locans = [locans Data.trial.retrieval(trialNum).locationAcc];
        end
    end
    for trialNum = 1:size(Data.trial.interference,2)
        if Data.trial.interference(trialNum).trial == currTrial
            interans = [interans Data.trial.interference(trialNum).correct];
        end
    end
    aniAcc = [aniAcc mean(anians)]; locAcc = [locAcc mean(locans)]; 

    interans = []; 
    inter99 = 0; 
    for i = size(Data.trial.interference,2)-7:size(Data.trial.interference,2)
        if Data.trial.interference(i).correct == 99
        	inter99 = inter99 + 1;
        else
        	interans = [interans Data.trial.interference(i).correct];
        end
    end
    interAcc = [interAcc mean(interans)]; 
    interSum99 = [interSum99 inter99];

    fprintf('- Trial %d End -\n\n', currTrial);
    fprintf('Trial retrieval accuracy: animal - %.3f, location - %.3f\n', mean(anians), mean(locans));
    fprintf('Trial interference accuracy: interference - %.3f, no response - %d\n\n\n', mean(interans), inter99);
elseif isControl
    fprintf('- Control Trial %d End -\n\n\n', currRun);
    fprintf('*** RUN %d END ***\n\n\n', currRun);
end


%%
% experiment end (4th run, control finished)

if isControl && currRun == 4
    expEndTextT = Screen('MakeTexture', window, expEndText);
    Screen('DrawTexture', window, expEndTextT, [], [0 0 2*xCenter 2*yCenter], 0);
    Screen('Flip', window);
    
    Data.time(end+1).run = currRun;
    Data.time(end).trial = currTrial;
    Data.time(end).control = isControl;
    Data.time(end).label = 'experiment end';
    Data.time(end).time = GetSecs;
    Data.time(end).Time = toc(expStart);
    Data.time(end).runTime = toc(runStart);

    allaniAcc = mean(aniAcc); alllocAcc = mean(locAcc); allinterAcc = nanmean(interAcc);
    
    fprintf('*** EXPERIMENT END ***\n\n');
    fprintf('Sbj retrieval accuracy: animal - %.3f, location - %.3f\n', allaniAcc, alllocAcc);
    fprintf('Sbj interference accuracy: interference - %.3f, no response - %d\n\n', allinterAcc, sum(interSum99));
    
    WaitSecs(5);
    
    isRunning = false;
    sca;
end

end

%% save

if isRestart
    fileName = ['./results/', sbjDate, '_sbj', sbjNum, '_', num2str(n_restart), '_WCDWWW.mat'];
else
    fileName = ['./results/', sbjDate, '_sbj', sbjNum,'_WCDWWW.mat'];
end
save(fileName, 'Data');

diary(diary_www);
diary off


