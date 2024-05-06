%
%
%  Randomizes Position of Spots for Imaris 8.4.1
%
%  Copyright Francesco Marangoni 2017
%
%
%  Installation:
%
%  - Copy this file into the XTensions folder in the Imaris installation directory
%  - You will find this function in the Image Processing menu, submenu
%  Spots Functions
%
%    <CustomTools>
%      <Menu>
%      <Submenu name="Spots Functions">
%        <Item name="Randomize Spot Distribution" icon="Matlab">
%          <Command>MatlabXT::XTRandomizeSpotDistribution(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%    </CustomTools>
%
%  Description:
%
%   The user has to create a Spots component. 
%   This XTension reads the position of the center of every Spot
%   and translates it to a new, RANDOM position within the Volume. 
%   This is useful to calculate whether the original distribution 
%   of Spots was random or not.
%
%   Usage Example:
%       - Load "yourfile.ims".
%       - Create Spots (default values) and select it
%       - Select tools tab from spots component
%       - Start "Randomize Spot Distribution"
%
function XTRandomizeSpotDistribution(aImarisApplicationID)

if isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
    vImarisApplication = aImarisApplicationID;
else
    % connect to Imaris interface
    javaaddpath ImarisLib.jar
    vImarisLib = ImarisLib;
    if ischar(aImarisApplicationID)
        aImarisApplicationID = round(str2double(aImarisApplicationID));
    end
    vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
end

% Connects to selected Spots
vSpots = vImarisApplication.GetFactory.ToSpots(vImarisApplication.GetSurpassSelection);

% Gives error message if no Spots are selected
if (vImarisApplication.GetFactory.IsSpots(vImarisApplication.GetSurpassSelection) == 0)
    msgbox('Please select some Spots!');
    return;
end

% Asks user the Channel depicting where cells are
vAns = inputdlg({'Please enter the NUMBER of the Channel showing where cells are:'}, ...
    'Randomize cells within allowed space',1,{'8'});
if isempty(vAns), return, end
vChID =abs(str2double(vAns{1}));

vSpotsTest = vImarisApplication.GetFactory.CreateSpots; % initializes spots to retrieve stats

vXmin = vImarisApplication.GetDataSet.GetExtendMinX;
vYmin = vImarisApplication.GetDataSet.GetExtendMinY;
vZmin = vImarisApplication.GetDataSet.GetExtendMinZ;
vXmax = vImarisApplication.GetDataSet.GetExtendMaxX;
vYmax = vImarisApplication.GetDataSet.GetExtendMaxY;
vZmax = vImarisApplication.GetDataSet.GetExtendMaxZ;

% Retrieve Spot objects info
 vSpotsRadius = vSpots.GetRadiiXYZ;
 vSpotsPosT = vSpots.GetIndicesT;

% Calculates Spot objects info one by one

for vIndex = 0:length(vSpotsRadius)-1
n =  0;
vSpotIntSumValues = 0;
while vSpotIntSumValues < 6.6  % MODIFY HERE THRESHOLD TO ENTER RANDOMIZATION 
     XCenter = random('unif', vXmin, vXmax); % Randomizes coordinates of new center of mass
     YCenter = random('unif', vYmin, vYmax);
     ZCenter = random('unif', vZmin, vZmax);

     vSpotsTest.Set([[XCenter, YCenter, ZCenter]], [0], [10]);

     vAllSpotStatistics = vSpotsTest.GetStatistics;
     vNamesSpotStat = cell(vAllSpotStatistics.mNames);  %creates cell array with Stat names
     vFactorsSpotStat = cell(vAllSpotStatistics.mFactors);  %creates cell array with rows being Category, Channel, Collection, Time
     vValuesSpotStat = vAllSpotStatistics.mValues;  %creates cell array with Values
     vIntSumSpotIndex = strmatch('Intensity Mean', vNamesSpotStat);  % gets indices of values belonging to Collection "Intensity Sum"
     vIntSumSpotIndex(:,2) = str2double(vFactorsSpotStat(2,(min(vIntSumSpotIndex):max(vIntSumSpotIndex))))';  % combines the previous info and Channel info
     vSpotIntSumIndexCh = vIntSumSpotIndex(find(vIntSumSpotIndex(:,2)== vChID));  % get indices of values belonging to "vChID" channel
     vSpotIntSumValues = vValuesSpotStat(vSpotIntSumIndexCh,:);  % retrieves fluorescence values
     n = n+1;
end
     vRandPosXYZ(vIndex+1,1) = XCenter; % Puts coordinates passing test to array of new coordinates
     vRandPosXYZ(vIndex+1,2) = YCenter;
     vRandPosXYZ(vIndex+1,3) = ZCenter;
     vIndex
     n
     
end

vSpotsNew = vImarisApplication.GetFactory.CreateSpots; % initializes new spots
vSpotsNew.Set(vRandPosXYZ, vSpotsPosT, vSpotsRadius(:,1));
vSpotsNew.SetRadiiXYZ(vSpotsRadius);
vSpotsNew.SetName(sprintf([char(vSpots.GetName) ' RANDOMIZED']));
vImarisApplication.GetSurpassScene.AddChild(vSpotsNew, -1); %Adds Spots to scene

