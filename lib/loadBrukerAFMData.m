function [z,d,t,dt,varargout] = loadBrukerAFMData(filepath,medium)
%loadBrukerAFMData Extract AFM Force Curve Data from Bruker Files
%   This function will take in the filepath to a bruker file and then load
%   the data inside it. The standard outputs are the z-extension (z),
%   deflection (d), time array (t), and timestep (dt). In addition, the 
%   function can provide the cantilever stiffness (k), the retract-phase
%   deflection, retract-phase z-extension, deflection scaling factor, and
%   z-extension scaling factor. These last two factors are used to convert
%   raw data in Volts/Least-Significant-Bits (or V/LSB) to distance (nm).

% Open an instance of the given file
fid = fopen(filepath,'r');

% Find the offsets in the binary
sHSDC_datalength = di_header_find(filepath,'\Data length:');
sHSDC_dataoffset = di_header_find(filepath,'\Data offset:');

% Search through the file and get some useful positions
for i = 1:length(sHSDC_datalength)
    fseek(fid,sHSDC_datalength(i),-1);
    line = fgets(fid);
    datalength(i) = extract_num(line);

    if i > 1
        fseek(fid,sHSDC_dataoffset(i-1),-1);
        line = fgets(fid);
        dataoffset(i-1) = extract_num(line);
    end
end

% Rewind the file and start reading the data
frewind(fid);
fseek(fid,datalength(1),-1);
rawData = (fread(fid, datalength(1)/2, 'int16'));

% Store our data in a temporary variable to avoid overwriting anything
% accidentally
usefulData = rawData;

% Create start and end indices for our dataset based on the standard size
% of the binary files
a_ind = [1 datalength(2)/4+1 datalength(2)/2+1 datalength(2)/2+datalength(3)/4+1];
b_ind = [datalength(2)/4 datalength(2)/2 (datalength(2)/2)+(datalength(3)/4) (datalength(2)/2+datalength(3)/2)];

if b_ind(end) > length(usefulData)
    b_ind(end) = length(usefulData);
end

% Grab the four phases of data we care about: Z-raw and Z-sensor
approachcurvesraw = flip(usefulData(a_ind(1):b_ind(1)));
approachcurvesraw(approachcurvesraw==0) = [];
retractcurvesraw = flip(usefulData(a_ind(2):b_ind(2)));
retractcurvesraw(retractcurvesraw==0) = [];
zapproachraw = flip(usefulData(a_ind(3):b_ind(3)));
zapproachraw(zapproachraw==0) = [];
zretractraw = flip(usefulData(a_ind(4):b_ind(4)));
zretractraw(zretractraw==0) = [];

% Find the spring constant
k_temp = di_header_find(filepath,'\Spring Constant:');
fseek(fid,k_temp(1),-1);
line = fgets(fid);
k_cantilever = extract_num(line);

% Grab the deflection InVOLS
InVOLS_temp = di_header_find(filepath,'\@Sens. DeflSens:');
fseek(fid,InVOLS_temp(1),-1);
line = fgets(fid);
defl_InVOLS = extract_num(line); % Value in nm/V

% Grab the zSensor InVOLS
InVOLS2_temp = di_header_find(filepath,'\@Sens. ZsensSens:');
fseek(fid,InVOLS2_temp(1),-1);
line = fgets(fid);
zSens_InVOLS = extract_num(line); % Value in nm/V

% Grab the deflection scaling so we can convert deflection
deflScale_temp = di_header_find(filepath,'\@4:Z scale: V [Sens. DeflSens]');
fseek(fid,deflScale_temp,-1);
line = fgets(fid);
SensDeflSens_InVOLS = extract_num(line); % Value in V/LSB

% Grab the Z-sensor scaling so we can convert z-sensor
zScale_temp = di_header_find(filepath,'\@4:Z scale: V [Sens. ZsensSens]');
fseek(fid,zScale_temp,-1);
line = fgets(fid);
SensZSens_InVOLS = extract_num(line); % Value in V/LSB

% Calculate the scaling factors (these are in nm/LSB)
deflScaling = defl_InVOLS*SensDeflSens_InVOLS;
zPosnScaling = zSens_InVOLS*SensZSens_InVOLS;

% Scale our data
approachDefl = (approachcurvesraw(:)-(approachcurvesraw(1))).*deflScaling;
retractDefl = (retractcurvesraw(:)-(approachcurvesraw(1))).*deflScaling;
zApproach = (zapproachraw(:)-(zapproachraw(1))).*zPosnScaling.*(1e6).*(1e-9);
zRetract = (zretractraw(:)-(zapproachraw(1))).*zPosnScaling.*(1e6).*(1e-9);

% Make the time array
frewind(fid);
dt_temp = di_header_find(filepath,'\@Sample period: V (0.1000000 us/LSB)');
fseek(fid,dt_temp(2),-1);
line = fgets(fid);
line = erase(line,'\@Sample period: V (0.1000000 us/LSB)');
dt = extract_num(line).*1e-6; % Value in us

dataError_pos = ( length(zApproach) - find(abs(flip(gradient(zApproach)))>10*median(abs(flip(gradient(zApproach)))),1) );    % Find last position above 0
dataError_pos2 = ( length(approachDefl) - find(abs(flip(gradient(approachDefl)))>10*median(abs(gradient(approachDefl))),1) );    % Find where deflection blows up
dataError_pos = max([dataError_pos dataError_pos2]);   % Choose worst cutoff

if isempty(dataError_pos) || dataError_pos == 0
    z = zApproach*(1e-6);
    d = approachDefl*(1e-9);
else
    z = zApproach(dataError_pos:end)*(1e-6);
    d = approachDefl(dataError_pos:end)*(1e-9);
end
if length(z) ~= length(d)
    idx = min([length(z),length(d)]);
    z = z((end-idx+1):end);
    d = d((end-idx+1):end);
end
t = ((1:size(d,1)).*dt)';

fclose(fid);

% Now perform the basic AFM offset work to correct everything
% First, make some smooth data to find the place to use for correction
ysmooth = smoothdata(d);
[~,idx] = min(ysmooth);

% Perform the correction for measurement in liquid, if applicable
if strcmpi(medium,'liquid')
    % Do a quick fit
    ydata = d(1:idx);
    xdata = (1:length(ydata))';
    P = polyfit(xdata,ydata,1);
    yshift = P(1)*xdata+P(2);
    yshift_all = [yshift;yshift(end).*ones(length(d)-length(yshift),1)];
    
    % Visualize the correction
%     figure
%     hold on
%     plot(z,d,'bo')
%     plot(z(xdata),P(1)*xdata+P(2),'g-','linewidth',3)
%     plot(z,d-yshift_all,'r-')
%     hold off
    
    % Correct the approach data
    d = d - yshift_all;
end

% Calculate Deflection Offset
d_0_mean = mean(d(1:idx));
d = d - d_0_mean;
z = z - z(1) + d;

% Handle optional output arguments
for k = 1:nargout
    switch k
        case 1
            varargout{k} = k_cantilever;
        case 2
            varargout{k} = retractDefl;
        case 3
            varargout{k} = zRetract;
        case 4
            varargout{k} = deflScaling;
        case 5
            varargout{k} = zPosnScaling;            
        otherwise
            error('Too many output arguments requested!');
    end
end

end

