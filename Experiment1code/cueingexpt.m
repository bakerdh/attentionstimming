function cueingexpt

% Posner cueing EEG experiment
% DHB 23/1/17

close all;
W = what; E.exptpath = strcat(W.path,'/'); clear W;
useVP = 1;

pID = inputdlg('Enter participant ID number');
E.subj = pID{1,1};

ST.SF = 2;
ST.npixelsperdegree = 36;       % at 57cm
ST.gratingsize = 4*ST.npixelsperdegree;
ST.ncycles = ST.SF*ST.gratingsize/ST.npixelsperdegree;

ST.duration = 0.2;
ST.cuedur = 0.2;
ST.postcuedurmin = 0.4;        % variable ISI from 0.4 - 0.6s
ST.postcuevarrange = 0.2;
ST.IBI = 10;                    % enforced break of 10s between blocks to stop and start EEG recording
ST.ITI = 1;
ST.ITIstd = 0.2;        % also variable post response interval between trials
ST.basecontrast = 0.5;
ST.inccontrast = 0.6;  % small contrast increment to one half of grating

E.blockstorun = 4;
E.ntrialsperblock = 200;

WaitSecs(0.01);                % important to load in the MEX file before the expt starts
GetSecs;
InitializePsychSound;
tr = PsychPortAudio('Open',[],[],[],[],1);
PsychPortAudio('FillBuffer', tr, MakeBeep(440*sqrt(2),0.05,44100).*0.5);

fname = strcat('Results/',E.subj,'results.mat');
if exist(fname,'file')
    load(fname);
else
    E.cuesidelist(1:400) = 1;
    E.cuesidelist(401:800) = 2;
    E.cuecongruentlist(1:800) = 1;
    E.cuecongruentlist(301:400) = 2;
    E.cuecongruentlist(701:800) = 2;
    E.trialorderlist = randperm(800);
    E.trialcounter = 0;
    R.allresponses = zeros(1,800);
    R.responsetime = zeros(1,800);
    R.trialtime = zeros(1,800);
    E.currentblock = 0;
    save(fname,'E','R');
end


try                 % start the 'try/catch' loop
    
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    PsychGPUControl('SetDitheringEnabled', 0);
    
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    if useVP        % using a ViewPixx or ProPixx
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        Datapixx('Open');
        Datapixx('EnableVideoScanningBacklight');       % Only required if a VIEWPixx.
        Datapixx('RegWr');
        
        PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
        PsychImaging('AddTask', 'General', 'EnableDataPixxM16OutputWithOverlay');
        [w, winRect] = PsychImaging('OpenWindow', screenNumber, 0, [], [], [], 0, 0, kPsychNeedFastBackingStore);
        Screen('LoadNormalizedGammaTable', w, linspace(0,1,256)'*ones(1,3), 0); % THIS IS THE IMPORTANT THING TO DO, NOTE THE LAST ARGUMENT IS 0.
        
        HideCursor;
        ST.greylevel = doimagegamma(0.5);
    else
        rect = [0 0 1024 1024];
        [w, winRect] = Screen('OpenWindow',screenNumber,128,rect);
        ST.greylevel = round(255*doimagegamma(0.5));
    end
    
    [ST.width, ST.height] = Screen('WindowSize', w);
    ifi=Screen('GetFlipInterval', w);
    ifims = ifi*1000;
    
    if useVP
        if ifims>10
            disp('Refresh rate is not 120Hz');
        end
    else
        ifims = 1000/60;           % fix as though 120Hz for testing
    end
    
    Screen('TextSize', w, 24);
    
    Screen('FillRect', w, ST.greylevel);
    drawfixation(w, ST.width, ST.height);
    lastflip = Screen('Flip', w);
    
    window = make_soft_window(ST.gratingsize,ST.gratingsize,0.9);
    
    gratingL = mkgrating(ST.gratingsize, ST.ncycles, 90, 0, ST.basecontrast) .* window;
    gratingH = mkgrating(ST.gratingsize, ST.ncycles, 90, 0, ST.inccontrast) .* window;
    grating(1:(ST.gratingsize/2),:) = gratingL(1:(ST.gratingsize/2),:);
    grating((((ST.gratingsize/2)+1):ST.gratingsize),:) = gratingH((((ST.gratingsize/2)+1):ST.gratingsize),:);
    comp = (1+grating)/2; comp = doimagegamma(comp);
    gtextures(1) = Screen('MakeTexture', w, comp, [], [], 2);
    grating(1:(ST.gratingsize/2),:) = gratingH(1:(ST.gratingsize/2),:);
    grating((((ST.gratingsize/2)+1):ST.gratingsize),:) = gratingL((((ST.gratingsize/2)+1):ST.gratingsize),:);
    comp = (1+grating)/2; comp = doimagegamma(comp);
    gtextures(2) = Screen('MakeTexture', w, comp, [], [], 2);
    
    facestim = imread('facel.bmp');
    facestim = double(facestim)/255;
    facestim = doimagegamma(facestim);
    facetextures(1) = Screen('MakeTexture', w, facestim*255);
    facestim = imread('facer.bmp');
    facestim = double(facestim)/255;
    facestim = doimagegamma(facestim);
    facetextures(2) = Screen('MakeTexture', w, facestim*255);
    facesize = size(facestim);
    r1 = [1 1 facesize(1)/2 facesize(2)/2];
    faceRect = CenterRectOnPoint(r1, ST.width*0.5, ST.height*0.5);
    
    ST.nframes = round(1000/ifims);       % no of frames in 1 sec
    targetframes = round(ST.duration*ST.nframes);
    
    r1 = [1 1 ST.gratingsize ST.gratingsize];
    destRect(:,1) = CenterRectOnPoint(r1, ST.width*0.35, ST.height*0.5);
    destRect(:,2) = CenterRectOnPoint(r1, ST.width*0.65, ST.height*0.5);
    
    
    breakcode = 0;
    
    while E.currentblock < E.blockstorun
        
        E.currentblock = E.currentblock + 1;
        block = E.currentblock;
        cuetype = 2;        % always face cues

        R.reportedifi(E.currentblock) = ifi;
        
        PsychPortAudio('Start', tr);
        
        Screen('FillRect', w, ST.greylevel);
        drawfixation(w, ST.width, ST.height);
        Screen('DrawText', w, 'Start recording, then click mouse to begin next block',0,0,0);
        lastflip = Screen('Flip', w);
        
        mouseloop;
        
        Screen('FillRect', w, ST.greylevel);
        drawfixation(w, ST.width, ST.height);
        lastflip = Screen('Flip', w);
        WaitSecs(2);
        
        if useVP
            for n = 1:3                         % 3 trigger pulses to indicate start of block
                Datapixx('SetDoutValues', transformindex(1));
                Datapixx('RegWrRd');
                WaitSecs(0.05);
                Datapixx('SetDoutValues', 0);
                Datapixx('RegWrRd');
                WaitSecs(0.05);
            end
        end
        WaitSecs(1);
        
        Screen('FillRect', w, ST.greylevel);
        drawfixation(w, ST.width, ST.height);
        lastflip = Screen('Flip', w);
        
        trial = 0;
        delayuntilnextstim = 0;
        responseoffset = lastflip;
        startofblock = lastflip;
        
        while trial < E.ntrialsperblock
            trial = trial + 1;
            E.trialcounter = E.trialcounter + 1;
            cuecongruent = E.cuecongruentlist(E.trialorderlist(E.trialcounter));
            cueside = E.cuesidelist(E.trialorderlist(E.trialcounter));
            if cuecongruent==1
                targetside = cueside;       % congruent
            else
                targetside = 3-cueside;     % incongruent
            end
            
            
            Screen('FillRect', w, ST.greylevel);
            drawfixation(w, ST.width, ST.height);
            lastflip = Screen('Flip', w, responseoffset+delayuntilnextstim);
            
            linelength = 25;
            
            Screen('FillRect', w, ST.greylevel);
            % show the cue first
            if cuetype==1       % arrow cue
                Screen('DrawLine', w, [0], (ST.width/2)-(linelength+1), ST.height/2, (ST.width/2)+linelength, ST.height/2, 5);
                if cueside==1
                    Screen('DrawLine', w, [0], (ST.width/2)-(linelength+1), ST.height/2, (ST.width/2)-(linelength*0.25), ST.height/2+(linelength*0.25), 5);
                    Screen('DrawLine', w, [0], (ST.width/2)-(linelength+1), ST.height/2, (ST.width/2)-(linelength*0.25), ST.height/2-(linelength*0.25+1), 5);
                else
                    Screen('DrawLine', w, [0], (ST.width/2)+linelength, ST.height/2, (ST.width/2)+(linelength*0.25), ST.height/2+(linelength*0.25), 5);
                    Screen('DrawLine', w, [0], (ST.width/2)+linelength, ST.height/2, (ST.width/2)+(linelength*0.25), ST.height/2-(linelength*0.25+1), 5);
                end
            else                % face cue
                Screen('DrawTexture', w, facetextures(cueside), [], faceRect);
            end
            lastflip = Screen('Flip', w);
            
            if useVP
                Datapixx('SetDoutValues', transformindex(100*cuetype + cueside));     % cue type and direction encoded
                Datapixx('RegWrRd');
            end
            
            Screen('FillRect', w, ST.greylevel);
            drawfixation(w, ST.width, ST.height);
            lastflip = Screen('Flip', w, lastflip + ST.cuedur);
            
            if useVP
                Datapixx('SetDoutValues', 0);       % turn off trigger encoding cue
                Datapixx('RegWrRd');
            end
            
            Screen('FillRect', w, ST.greylevel);
            drawfixation(w, ST.width, ST.height);
            lastflip = Screen('Flip', w, lastflip + ST.postcuedurmin + rand*ST.postcuevarrange);
            
            % now show the stimulus
            incside = ceil(rand*2);      % target increment is on left or right of the grating
            n = 0;
            while n < targetframes
                n = n+1;
                
                Screen('FillRect', w, ST.greylevel);
                Screen('DrawTexture', w, gtextures(incside), [], destRect(:,targetside));
                drawfixation(w, ST.width, ST.height);
                lastflip = Screen('Flip', w);
                
                if useVP
                    if n<3
                        Datapixx('SetDoutValues', transformindex(targetside*10+incside));     % first trigger also contains condition code
                    else
                        Datapixx('SetDoutValues', 0);
                    end
                    Datapixx('RegWrRd');
                    
                end
                if n==1
                    trialstart = lastflip-startofblock;
                end
            end
            Screen('FillRect', w, ST.greylevel);
            drawfixation(w, ST.width, ST.height);
            trialoffset = Screen('Flip', w);
            
            delayuntilnextstim = ST.ITI + ST.ITIstd*randn;
            
            exitcode = 0;
            while exitcode==0
                
                [x,y,buttons] = GetMouse;
                [keyIsDown, responseoffset, keyCode] = KbCheck;
                
                if buttons(1) || keyCode(KbName('UpArrow'))
                    resp = 1;
                    exitcode = 1;
                elseif sum(buttons(2:end)) || keyCode(KbName('DownArrow'))
                    resp = 2;
                    exitcode = 1;
                elseif keyCode(KbName('Escape'))
                    resp = 0;
                    breakcode = 1;
                    exitcode = 1;
                    trial = 1000;
                    E.currentblock = 100;
                end
            end
            if useVP        % send a trigger to record the response in the EEG file
                Datapixx('SetDoutValues', transformindex(30+resp));
                Datapixx('RegWrRd');
                WaitSecs(0.05);
                Datapixx('SetDoutValues', 0);
                Datapixx('RegWrRd');
                WaitSecs(0.05);
            end
            
            R.allresponses(E.trialorderlist(E.trialcounter)) = resp;
            R.responsetime(E.trialorderlist(E.trialcounter)) = responseoffset-trialoffset;
            R.trialtime(E.trialorderlist(E.trialcounter)) = trialstart;
            
            Screen('FillRect', w, ST.greylevel);
            drawfixation(w, ST.width, ST.height);
            Screen('Flip', w);
            
        end
        
        if ~breakcode
            save(fname, 'R', 'E');       % save every block of trials
            if E.currentblock < E.blockstorun
                
                Screen('FillRect', w, ST.greylevel);
                drawfixation(w, ST.width, ST.height);
                Screen('DrawText', w, 'Stop recording and take a break',0,0,0);
                lastflip = Screen('Flip', w);
                
                WaitSecs(ST.IBI);
            end
        end
        
    end
catch
    lasterr
end

Screen('Close',gtextures);

Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);

Screen('Flip', w);
ShowCursor;

if useVP        % using a ViewPixx or ProPixx
    Datapixx('DisableVideoScanningBacklight');
    WaitSecs(0.1);
    Datapixx('RegWr');
end

WaitSecs(0.1);
Screen('CloseAll');
WaitSecs(0.1);
if useVP
    Datapixx('Close');
end
WaitSecs(0.1);
PsychPortAudio('Close', tr);

end
%--------------------------------------------------------------------------
function drawfixation(w, width, height)

Screen('DrawLine', w, [0], (width/2-4)-1, (height/2)-1, (width/2+4)-1, (height/2)-1, 2);
Screen('DrawLine', w, [0], (width/2)-1, (height/2+4)-1, (width/2)-1, (height/2-4)-1, 2);

end
%--------------------------------------------------------------------------
function output = doimagegamma(i)

% gamma corrects the stimuli before sending to Bits++
% adapted from Mark's code, DHB 29.01.08

% parameters from last gamma correct, 12/8/14 on VPixx in M16 mode
k = 0.5344;
Lmax = 101.8721;
j0 = 0.1292;
gamma = 1.9252;
%%%%%

i0 = 0;
imax = 1;                                               % Bits++ always scaled between 0 and 1
imean = (i0+imax)/2;
jmax = 1;

Lmin = k + (Lmax-k)*(max(-j0,0)/(jmax-j0) ).^gamma;     % Eqn 2, with j set to 0, to get Lmin
Lmin = max(Lmin,0);                                     % ensure Lmin not <0
Lmean = (Lmin+Lmax)/2;
L = Lmean + (Lmax-Lmean)*(i-imean)/(imax-imean);        % desired luminance values Eqn 4
j = ((L - k)/(Lmax-k)).^(1/gamma)*(jmax - j0) + j0;     % These are the gamma-corrected lut values, j: Eqn 3
output = max(j,j0);                                     % Eqn 3 conditional
output = double(output);

end
%--------------------------------------------------------------------------------------------------
function output = transformindex(input)

% fixes the binary inputs for the EEG amplifier because the pins are in a different order from the ViewPixx
% desired numbers must be <256
% DHB 18/8/14

truebits = 2.^(2:2:24);
dn = dec2bin(input,length(truebits));
output = 0;
for m = 1:length(truebits)
    output = output + truebits(m)*str2num(dn(end-m+1));
end

end
%--------------------------------------------------------------------------
function mouseloop

exitcode = 0;

while exitcode==0
    [x,y,buttons] = GetMouse;
    
    if sum(buttons)>0
        exitcode = 1;
    end
end

end
%--------------------------------------------------------------------------------------------------
function imag1 = mkgrating(Regionsize, f, o, p, c)

%TSM; 26.6.03
% modified by DHB to make single component gratings only, scaled from -1 to 1
% f is spatial frequency, scaled as cycles per image
% o is orientation (degrees), p is phase (degrees relative to centre), c is contrast

p = p*pi/180;
o = o*2*pi/360;		% convert from degrees to radians
f = f/Regionsize;
x0 = ((Regionsize+1)/2);
y0 = x0;

u = f .* cos(o) * 2 * pi;
v = f .* sin(o) * 2 * pi;

imag1 = zeros(Regionsize, Regionsize);
[xx, yy] = meshgrid(1:Regionsize, 1:Regionsize);

imag1(:,:) = (c .* sin(u .*(xx-x0) + v.*(yy-y0) + p));

end
%--------------------------------------------------------------------------
function mask = make_soft_window(W,H,D)

% Mark's code for making a raised cosine window

% SYNTAX: mask = make_soft_window(W,H,[D])
% ends an array 'mask' that is 1 inside the circular window, shading to zero outside
% W, H are the width and height of the whole (rectangular or square) array, in pixels
% Diameter of the soft window at half-height defaults to 0.90 units
%    where 1 unit = image width or height (whichever is the smaller)
% Smoothing is by convolution with a cosine half-cycle of width 0.1 units
% Optional parameter D specifies this diameter (in relative units, range 0 -> 1)
% MAG, 27.2.04

%soft window parameters
if nargin<3, D = 0.9; end % sets default diameter to 0.9
radius = min(W*D/2,H*D/2);% radius in pixels
blur = 2*(min(W/2,H/2) - radius);  % blur half-cycle
L = blur;
X1 = [-L/2:L/2];

% 1-D blur function (applied twice, in x and y)
WinKernel = cos(X1*pi/L); % half-cycle cosine
%image coordinates - X and Y arrays
X = [1:W] - W/2;
Y = [1:H] - H/2;
xx = repmat(X,H,1);
yy = repmat(Y',1,W);

% make circular soft window
mask = single((xx.*xx + yy.*yy) < radius^2); % logical 0 outside the circle,1 inside it
mask = conv2(WinKernel,WinKernel,mask,'same'); 	% smooth the mask
mask = mask/max(max(mask));						% scale the mask 0-1
% figure(2);plot(X,mask(H/2,:),'r-',Y,mask(:,W/2));
mask = double(mask);
end
%--------------------------------------------------------------------------------------------------