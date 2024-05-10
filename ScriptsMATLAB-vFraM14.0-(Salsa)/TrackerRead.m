% EDITED BY FRANCESCO MARANGONI 

function [varargout] = TrackerRead_March(infile, cyc_time, max_vel, min_track_length, min_disp, min_inst_vel, varargin)
% defines the file as a function file. Output variables in [], input variables in().  
% infile is the input file, cyc_time, max_vel, min_track_length, min_disp, and other variables (varargin) come from the CellTracker GUI.

wh = waitbar(0, 'Reading and Analyzing file...');

fid = fopen(infile);
% opens the variable "infile" together with the file identifier (fuction fopen, -1 = not allowed, 1 = allowed, 2 = error)
if (fid <= 0),
    error('Cannot open data file provided.  Check name and/or path.');
end  % "End" is used to terminate a cycle or a conditional statement (if)

all_names = {};
all_data = {};
data = [];
num_tracks = 1;
numlines = 1;
current_track = NaN;
%{} is an emty cell of an array of arrays, [] is an empty array

%READS THE INPUT FILE AS ARRAY OF ARRAYS
while ~feof(fid),   % while (~) not (feof) end of file named "fid" 
    tline = fgetl(fid); %reads a line without its terminator (fgetl) and assignes it to the variable tline
    if (length(tline) < 4), % if length of characters in tline is less than 4 
        %ignore this bad bad line
    elseif (strcmp(upper(tline(1:4)), 'NAME')),  % or if elements of tline are NAME (uppercase)
        % We are in the first line.
    else  %in all other cases
        lineofdata = {};
        datacount = 1;
        while (length(tline) > 0),  % while the read line (tline) has elements (is >0)
            [token, tline] = strtok(tline, char(9));  % assign to the variable "token" the value of the line until (char(9), and save the remainder in tline
            lineofdata{datacount} = token; %assign to the array of array lineofdata, position 1 (then 2, 3, and so on until the last field - they are 10 in the input file) the value contained in "token".
            datacount = datacount+1; %increments datacount (field number)
        end  % so this cycle memorizes a line of data subdividing it in fields based on tab. 
        
        if (strcmp(lineofdata{2}, 'N/A')),  % if the position 2 of lineofdata is N/A
            %ignore
        elseif (str2double(lineofdata{2}) ~= current_track),  % if the position 2 of lineofdata (converted in double-precision number) is different (~=) from the "current_track" value
            % We are at the beginning of a cell track, because we have detected a new track number
            if ~isempty(data),  % if the array "data" is not empty
                all_data{num_tracks} = data;  % assign to the array of array all_data (at position "num_tracks", initially 1) the value of "data"
                num_tracks = num_tracks + 1;  %increments num-tracks
            end
            % this cycle subdivides the input in tracks, adding the result to the array of arrays line by line

            current_track = str2double(lineofdata{2});  % establishes the number of the current track
            numlines = 1;  % resets numlines to 1
            data = [];  %resets data
            for i = 1:10,  %This line defines the number of fields to handle.  MODIFY IT IN CASE YOU WANNA ADD MORE FIELDS.
                data(numlines, i) = str2num(lineofdata{i+1});  % creates array "data" from the combined array "lineofdata" without position 1 (i.e.erases information about the name of each object) and starts a count (numlines) of objects composing each track. 
            end
            numlines = numlines + 1; %increments the number of the line
            
            lineofdata = {};  %resets lineofdata, which is now unuseful since its data (but the object name) have been copied to "data". 

        else
            %Data past first line.
            for i = 1:10,  %This line defines the number of fields to handle.  MODIFY IT IN CASE YOU WANNA ADD MORE FIELDS.
                data(numlines, i) = str2num(lineofdata{i+1});
            end
            numlines = numlines + 1;
        end
	% the same thing happens to subsequent lines.  Until now, tracks are
	% defined only using their number
    end
end
% This block reads the input file, subdividing it into fields (except for the object name, which is erased) and tracks

if ~isempty(data),
    all_data{num_tracks} = data;
    num_tracks = num_tracks + 1;
end
% saves last track

fclose(fid);    %closes input file
waitbar(0.25);  % waitbar proceeds to one quarter

% open output files
analysis_fid = fopen(sprintf('%s_analysis.txt', strtok(infile, '.')), 'w');
if (analysis_fid <= 0),
    error('Cannot create analysis output file.');
end

summary_fid = fopen(sprintf('%s_summary.txt', strtok(infile, '.')), 'w');
if (summary_fid <= 0),
    error('Cannot create summary output file.');
end

msd_fid = fopen(sprintf('%s_MSD.txt', strtok(infile, '.')), 'w');
if (msd_fid <= 0),
    error('Cannot create MSD output file.');
end

%try
    % print header with tabs (\t) and carriage returns (\n)
    fprintf(analysis_fid, 'Track #\tStep #\tCentroid X\tCentroid Y\tCentroid Z\tDelta 3D Dist (microns)\tDelta x\tDelta y\tDelta z\tInst. 3D Vel (microns/min)\tInst. Vel. X\tInst. Vel. Y\tInst. Vel. Z\t3D Disp (micron)\tDisp. X\tDisp. Y\tDisp. Z\tAngle Change\tInst. Vel. 2D\tMotile Angle Change\tTimepoint\tDelta Time\tSignaling Index\tTarget 1 Distance\tTarget 2 Distance\tQualitative interaction Tumor\tQualitative interaction 2\tn-point displacement\t Segment ID\t Segment Duration (min)\t Segment AUC\t Peak Max magnitude\t Peak Average Magnitude\t\n');
    fprintf(summary_fid,'Track #\tTotal Dist.\tTotal X\tTotal Y\tTotal Z\tMean Vel. 3D (microns/min)\tAvg. Vel X\tAvg. Vel Y\tAvg. Vel Z\tTotal Disp\tTotal Disp. X\tTotal Disp. Y\tTotal Disp. Z\tAvg. Angle Change\tConfinement Ratio\tAvg. Vel 2D\tTrack Duration (min)\tArrest Coefficient\tAvg. Signaling Index\tTrack position\tPercent time signaling\tExceeded Max Allowed Velocity\n');
    fprintf(msd_fid,'Track #\tStep Size\tD(t)\tMSD\n');
    
    totsign = 0;
    totsignpar = 0;
    totsignstr = 0;
    tpnb = 0;
    
    %smoothwindow=input('How much do you want velocity smoothed (ODD values only)? ');
    
    minutes=input('Minutes for n-minute displacement ');
    npoint = minutes*60/cyc_time; 
    %minutes2=input('Minutes before and after. (Minimum 1, for neighbor calculation) ');
    %npoint2 = minutes2*60/cyc_time; 
    npoint2 = 1;
    % target = input('Insert X Y Z coordinates of target [between brackets] ');
    
    for i=1:length(all_data),  %i becomes the track number
        cur_data = all_data{i};
        
        if (size(cur_data, 1) < (min_track_length - 1)),  
            continue;
        end 
	% if track duration (cur_data, field 1) is smaller or equal to the minimum allowed, it skips the analysis
        
        % distances
        del_pos = diff(cur_data(:, 3:5));  % assigns to the array "del_pos" the difference (diff) between two consecutive positions (:) of the array "cur_data" at positions 3,4,5 (x,y,z)
        del_dist = sqrt(del_pos(:,1).^2 + del_pos(:,2).^2 + del_pos(:,3).^2); % assigns distances to array "del_dist" .
    	del_dist_2D = sqrt(del_pos(:,1).^2 + del_pos(:,2).^2); 
	
	
        del_time = diff(cur_data(:, 2));   % computation of time distance between points ADDED 1/16/2011
        del_time_minusfirst = del_time(2:end)*cyc_time/60;  %deletes first element of del_time, useful to calculate AUC of velocity and NFAT
        
        
        % velocities
        inst_vel = del_pos./(cyc_time/60);  % ORIGINAL LINE array "inst_vel" stores instantaneous velocities on x,y,z SCREWS UP ID 2 POINTS ARE NOT CONSECUTIVE.
        inst_vel_3D = (del_dist./del_time)/(cyc_time/60);  % MODIFIED LINE 1/16/2011
        inst_vel_2D = (del_dist_2D./del_time)/(cyc_time/60);   % MODIFIED LINE  2/9/2011
        vel_norm_bkgsub = (inst_vel_3D-4)/10.4;  %normalizes velocities using upper bound 14.4 and lower bound 4.  Formula:  xnorm = (x-LB)/(UB-LB)
        vel_norm_bkgsub(vel_norm_bkgsub<0) = 0;  % sets negative values to zero
        velareainc = zeros(length(vel_norm_bkgsub)-1, 1);
        for b = 2:length(inst_vel_3D),  %useful to compute AUC of normalized, background subtracted 3D instantaneous velocities 3/6/2012
            velareainc(b-1) = ((vel_norm_bkgsub(b-1)+vel_norm_bkgsub(b))/2)*del_time_minusfirst(b-1);
        end
        
        % displacements
        del_disp = cur_data(2:end, 3:5) - repmat(cur_data(1, 3:5), [size(cur_data, 1)-1 1]);
        del_disp_3D = sqrt(del_disp(:,1).^2+del_disp(:,2).^2 + del_disp(:,3).^2);
        
        
        sign_index = cur_data(:,6); % modified 3/2/11, signaling info is in the 7th column of input (6 here because name column is cropped)
        sign_index_minusfirst = sign_index(2:end);  %useful to compute AUC of velocity and NFAT 3/6/2012
        NFAT_norm_bkgsub = (sign_index_minusfirst-65)/35;  %normalizes NFAT indexes using upper bound 100 and lower bound 65.  Formula:  xnorm = (x-LB)/(UB-LB)
        NFAT_norm_bkgsub(NFAT_norm_bkgsub<0) = 0;  % sets negative values to zero
        NFATareainc = zeros(length(NFAT_norm_bkgsub)-1, 1);
        for a = 2:length(NFAT_norm_bkgsub),  %useful to compute AUC of velocity and NFAT 3/6/2012
            NFATareainc(a-1) = ((NFAT_norm_bkgsub(a-1)+NFAT_norm_bkgsub(a))/2)*del_time_minusfirst(a-1);
        end
        
        AUCsign = sum(NFATareainc);  %added FraM 3/6/2012 to compute AUC of NFAT
        AUCvel = sum(velareainc);   %added FraM 3/6/2012 to compute AUC of Vel
        
        
        
        
        interaction1 = cur_data(:,7);  % modified 3/2/11, interaction info is in the 8th column of input (7 here because name column is cropped)
        interaction2 = cur_data(:,8);  % modified 4/11/11, interaction info is in the 9th column of input (8 here because name column is cropped)
        tum_interaction = cur_data(:,9);  % modified 5/27/11, tum_interaction info is in the 10th column of input (9 here because name column is cropped)
        DC_interaction = cur_data(:,10);  %modified 2/9/14, DC_interaction info is in the 11th column of input (10 here because name column is cropped)
        
        
        %-------- SEGMENT SPLITTING MODULE  2/9/14
        %If you need smoothing of data, implemet it here such as %RA_inst_vel_3D=smooth(inst_vel_3D, smoothwindow);  
        
        % insert constants here, if needed
        %score = [];
        %for j = 1:size(cur_data,1),
        %    if %(CRITERIA FOR GOOD SCORE), 
        %        score(j) = 1;
        %    else
        %        score(j) = 0;
        %    end
        %end % this block gives a good score to timepoints fulfilling criteria
        
        %  Detection of high quality segments
        %start = [];
        %finish = [];
        %del = [];
        %minlength = 5;  % minimum length of segments (timepoints)
        %gap = 1;        % max tolerated gap (timepoints)
        %start = find((diff([0,score,0])) ==1 );  % finds indexes of segment starting points
        %finish = find((diff([0,score,0])) ==-1 ) -1;  %  finds indexes of segment ending points
        
        %for z = 1:(length(finish)-1),  % this for cycle looks for indices of finish and start points closer than gap, and stores them in an array to be used to delete them
        %    if finish(z)+gap+1 >= start(z+1), 
        %        del = [del, z];
        %    end
        %end
        %start(del+1) = [];
        %finish(del) = [];
        %del = [];  % del reinitialization
        
        %for y = 1:length(start),  % this cycle calculates the length of each segment and rejects the ones below minimum lenght
        %    if (finish(y)-start(y)) <= minlength,
        %        del = [del, y];
        %    end
        %end
        %start(del) = [];
        %finish(del) = [];
        
        %displ2D = NaN(size(inst_vel_2D));  % prepares vectors for display change name as needed
        %displanglechange = NaN(size(anglevessel));  % prepares vectors for display change name as needed
        %segmentID = NaN(size(cur_data));  % this is probably to be kept as-is
        
        %for z = 1:length(start), % puts values into coherent vectors
        %    displ2D(start(z):finish(z)-1) = inst_vel_2D(start(z):finish(z)-1);  %change this vector as needed
        %    displanglechange(start(z):finish(z)-1) = anglevessel(start(z):finish(z)-1);%change this vector as needed
        %    segmentID(start(z):finish(z)) = ((z)*0.01)+cur_data(1,1); %leave this vector as-is
        %end
        %-----------------------------------------------
        
        %-------- SEGMENT SPLITTING BASED ON AVG SI 9/18/18
        % transform sign_index to fill spaces up with NaN
        full = cur_data(end, 2);
        full_sign_index = NaN(full,1);
        j=cur_data(1, 2);
        full_sign_index(j,1) = sign_index(1,1);
        for i = 1:size(del_time),
            j=j+del_time(i);
            full_sign_index(j,1) = sign_index((i+1),1);
        end
           
        MA_sign_index=movmean(full_sign_index, 1, 'omitnan');  % second term in smooth is the smoothing window and must be odd. Right now, no smoothing applied. substitute 2*npoint2+1 is smoothing is wanted.
        nanindex = find(isnan(full_sign_index));
        for i = 1:length(nanindex),
            MA_sign_index(nanindex(i)) = NaN;
        end
                
        score = zeros(size(MA_sign_index));  % this block gives a good score to timepoints fulfilling criterion1 (existence)
        for j = 1:size(MA_sign_index),
            if ~isnan(MA_sign_index(j))  %(GOOD SCORE = exists), 
                score(j) = 1;
            else
                score(j) = 0;
            end
        end 
        
        max_nonexistent = 3; % max tolerated nonexistent timepoint (minutes) - rehabilitates timepoints < max nonexistent
        max_nonex_npoint = max_nonexistent*60/cyc_time;  % converts minutes to timepoints
        count = 1 ;
        while count < size(score, 1),
            if score(count) == 0,
                count2 = count;
                while score(count2) == 0,
                    count2=count2+1;
                end    
                if (count2-count) <= max_nonex_npoint,
                     score(count:count2) = 1;
                end
                count = count2;
            else
                count = count+1;
            end
        end
                                   
        
        score1 = zeros(size(MA_sign_index));  % this block gives a good score to timepoints fulfilling criterion2 (signaling)
        baseline = prctile(MA_sign_index, 30) + (prctile(MA_sign_index, 30)-min(MA_sign_index));
        for j = 1:size(MA_sign_index),
                if ((MA_sign_index(j) >baseline) | (isnan(MA_sign_index(j)))),  %(GOOD SCORE = signaling or undetermined), 
                score1(j) = 1;
            else
                score1(j) = 0;
            end
        end 
                                       
        % gap = 1;        % max tolerated gap (minutes) 
        % gap_npoint = gap*60/cyc_time;  % converts minutes to timepoints
        gap_npoint = 1;
        count = 1 ;
        count2 = count;
        while count < size(score1, 1),
            if score1(count, 1) == 0,
                count2 = count;
                while ((score1(count2, 1) == 0) & (count2 < size(score1, 1))),
                    count2=count2+1;
                end    
                if (count2-count) <= gap_npoint,
                     score1(count:count2, 1) = 1;
                end
                count = count2;
            else
                count = count+1;
            end
        end
                
        score(:,1) = score(:,1).*score1(:,1); % puts the two conditions together
        
        %  Detection of high quality segments
        start = [];
        finish = [];
        del = [];
        % minlength = 1*60/cyc_time;  % minimum length of segments (minutes)
        minlength = 4; % minlength defined as timepoints, here it corresponds to 10 sec
        start = find((diff([0;score;0])) ==1 );  % finds indexes of segment starting points
        finish = find((diff([0;score;0])) ==-1 ) -1;  %  finds indexes of segment ending points
        for y = 1:length(start),  % this cycle calculates the length and AUC of each segment and rejects the ones not passing QC
            if (finish(y)-start(y)) < minlength,
                del = [del, y];
            else
                AUC = (trapz(MA_sign_index(start(y):finish(y)))-(baseline*(finish(y)-start(y))))*(cyc_time/60);
                Dur = (finish(y)-start(y))*(cyc_time/60);
                if Dur <= 1 & AUC <=1000
                    del = [del, y];
                elseif  Dur > 1 & (AUC/Dur)<=800
                    del = [del, y];
                end
            end
        end
        start(del) = [];
        finish(del) = [];
        
        full_segmentID = NaN(size(MA_sign_index)); 
        full_seg_duration = NaN(size(MA_sign_index));
        full_seg_AUC = NaN(size(MA_sign_index));
        full_seg_Peak = NaN(size(MA_sign_index));
        full_seg_AvgPeak = NaN(size(MA_sign_index));
        
        for z = 1:length(start), % puts values into coherent vectors
            full_segmentID(start(z):finish(z)) = ((z)*0.01)+cur_data(1,1);  %gives a name to segments
            full_seg_duration(start(z)) = (finish(z)-start(z))*(cyc_time/60);  %gives segment duration pairing it with the starting timepoint
            full_seg_AUC(start(z)) = (trapz(MA_sign_index(start(z):finish(z)))-(baseline*(finish(z)-start(z))))*(cyc_time/60); % computation of baseline subtracted AUC
            full_seg_Peak(start(z)) = max(MA_sign_index(start(z):finish(z)))-baseline; % computation of baseline subtracted max peak
            full_seg_AvgPeak(start(z)) = mean(MA_sign_index(start(z):finish(z)))-baseline; % computation of baseline subtracted average magnitude
        end
                
        % This block erases nonexistent timepoints
        segmentID = [];
        seg_duration = [];
        seg_AUC = [];
        seg_Peak = [];
        seg_AvgPeak = [];
        for w = 1:size(cur_data(:,2)),
            segmentID = [segmentID; full_segmentID(cur_data(w,2))];
            seg_duration = [seg_duration; full_seg_duration(cur_data(w,2))];
            seg_AUC = [seg_AUC; full_seg_AUC(cur_data(w,2))];
            seg_Peak = [seg_Peak; full_seg_Peak(cur_data(w,2))];
            seg_AvgPeak = [seg_AvgPeak; full_seg_AvgPeak(cur_data(w,2))];
        end
                              
                
        %-----------------------------------------------
        
        %percentage of signaling events in custom interval npoint2  9/18/18
        %THIS IS THE NEIGHBORHOOD ANALYSIS
        sign_index_pad = [NaN((npoint2),1); full_sign_index(:,1); NaN((npoint2),1)];  %pads sign_index array with 0 in order to calculate moving window on extremes
        sign_neighbors_temp = zeros(size(full_sign_index));
        for countsign = (npoint2:(size(sign_index_pad)-npoint2)),  
            numsig=0;
            percsig=0;
            countwindow=0;
            if sign_index_pad(countsign,1) > 0.45
                for neighbor = (countsign-npoint2):(countsign+npoint2)
                    if ~isnan(sign_index_pad(neighbor,1)),
                        countwindow = countwindow +1;
                    end
                    if sign_index_pad(neighbor,1) > 0.45,
                        numsig=numsig+1;
                    end
                end
                percsig = numsig/countwindow*100;
                sign_neighbors_temp((countsign-npoint2),1) = percsig;
            end
        end
        
        % This block erases nonexistent timepoints
        sign_neighbors = [];
        for w = 1:size(cur_data(:,2)),
            sign_neighbors = [sign_neighbors; sign_neighbors_temp(cur_data(w,2))];
        end
        
        
        %--------------------------------------------------
        
        % n-point displacement considering nonexistent timepoints 9/27/2018 
        % fill missing timepoints with NaN
        full_positions = NaN(full,3);
        j=cur_data(1,2);
        full_positions(j,1:3) = cur_data(1,3:5);
        for i = 1:size(del_time),
            j=j+del_time(i);
            full_positions(j,1:3) = cur_data((i+1),3:5);
        end
        % calculates n-minutes displacement
        del_disp_n = full_positions((1+npoint):(end), 1:3) - full_positions(1:(end-npoint), 1:3);% Calculates n-point displacement in each dimension
        del_disp_3D_n_full = sqrt(del_disp_n(:,1).^2+del_disp_n(:,2).^2 + del_disp_n(:,3).^2);  % Calculates n-point displacement in 3D
        %del_disp_3D_n(del_disp_3D_n>15)=16;  % poses upper limit of displacement to 15
        del_disp_3D_n_full(((end+1):(end+npoint)),:) = NaN;  % Fills the n-point 3D displacement arrays
        % removes nonexisting timepoints
        del_disp_3D_n = [];
        for w = 1:size(cur_data(:,2)),
            del_disp_3D_n = [del_disp_3D_n; del_disp_3D_n_full(cur_data(w,2))];
        end
        %------------------------------------------------

        % Average signaling in parenchyma and stroma, and calculation of NFAT efficiency 2/1/12
        sign_1 = 0;
        par_sign_1 = 0;
        str_sign_1 = 0;
        maxtimepoint = 0;
        avgsign = 0;
        avgsignpar = 0;
        avgsignstr = 0;
        fastNFAT = 0;
        NFATeff = '';
        
        maxtimepoint = max (cur_data(:,2));
        
        for conta = 1:size(cur_data,1),  
            
            if cur_data(conta,6) > 0.6,
                sign_1 = sign_1 +1;
                if cur_data(conta,9) == 0,
                    str_sign_1 = str_sign_1 + 1;
                else
                    par_sign_1 = par_sign_1 + 1;
                end
                %if inst_vel_3D(conta-1) > 4,  % fast steps are defined by velocity above 4 microns/min
                %    fastNFAT = fastNFAT + 1;
                %end
            end
        end
        
        %if sign_1>4,  %excludes calculation in tracks with too few signaling events
        %    NFATeff = fastNFAT/sign_1;
        %end
        
        
        avgsign = (sign_1/size(cur_data,1));
        %avgsignpar = par_sign_1;
        %avgsignstr = str_sign_1;
     
        % Calculation of the number of APCs scanned during track  2/1/2012
        %APCnum = size(find(unique(cur_data(:,10))));  %counts nonzero unique elements
        
        
       
        %Classification of tracks based on location 5/31/11
        if (cur_data(1,9) == 0) & (cur_data(end,9) ==0),
            pos_classifier = 'str-str';
        elseif (cur_data(1,9) ~= 0) & (cur_data(end,9) ==0),
            pos_classifier = 'par-str';
        elseif (cur_data(1,9) == 0) & (cur_data(end,9) ~=0),
            pos_classifier = 'str-par';
        else
            pos_classifier = 'par-par';
        end
        
        % directionality to target (Francesco 11/6/2022)
      % direction = [];
      % target_array = [];
      % dir_vel = [];
      % for k = 1:(size(cur_data, 1)-1),
      %     target_array(k, 1:3) = target - cur_data(k, 3:5);
         
      % end
      
      % for j = 1:(size(del_pos,1)),
      %      if (norm(del_pos(j,:)) == 0) || (norm(target_array(j,:)) == 0),
      %          direction(j) = 0;
      %      else
      %          direction(j) = dot(del_pos(j,:), target_array(j,:))/(norm(del_pos(j,:))*norm(target_array(j,:)));
      %      end          
      % end
       
      % dir_vel = inst_vel_3D .* direction'; % calculates directed speed
       
       %------------
      
      
        
        % angle change
        del_angle = [];
        for j = 1:(size(del_pos,1)-1),
            if (norm(del_pos(j,:)) == 0) || (norm(del_pos(j+1,:)) == 0),
                del_angle(j) = 0;
            else
                del_angle(j) = acos(dot(del_pos(j,:), del_pos(j+1,:))/(norm(del_pos(j,:))*norm(del_pos(j+1,:))))*180/pi;
            end
            % A positive cross product along z is counter-clockwise
            %  We shall call this a negative angle
            temp = cross(del_pos(j,:), del_pos(j+1, :));
            del_angle(j) = -sign(temp(3))*del_angle(j);           
        end
   
        % make the column of angle changes without non-motile cell angles
        % (defined by comparing inst. disp. to user's entered min)
        motile = [];
        for c = 1:length(del_disp_3D)-1,
            if (del_disp_3D(c+1) >= min_disp),
                motile(c) = del_angle(c);
            else
                motile(c) = NaN;
            end
        end

	% ARREST COEFFICIENT COMPUTATION (CONSIDERS DISTANT JOINING)
	% FRANCESCO MARANGONI 1/17/2011
	nonmotile_arr_coeff = [];
	for d = 1:length(inst_vel_3D),
	   if (inst_vel_3D(d) <= min_inst_vel),
		 nonmotile_arr_coeff(d) = del_time(d);
           else
             nonmotile_arr_coeff(d) = 0;
           end
        end
        
        
        
        exceeded = find(abs(inst_vel_3D) >= max_vel);  % finds indexes of instantaneous velocities above the limit
        exceeded = [0; exceeded; size(inst_vel_3D,1)+1];  % THIS MAY CAUSE PROBLEMS: cuts tracks based on overcoming the max allowed velocity.  So we now set the max allowed velocity in the GUI to 10000um/min
        
        for j = 1:(length(exceeded)-1),
            start_pos = exceeded(j) + 1;
            end_pos = exceeded(j+1) - 1;
            if ((end_pos - start_pos) < (min_track_length - 1)),
                continue;
            end
            
            %make a column of angle changes without non-motile cell angles as
            %defined by comparing the average dist for a given cell to the
            %user's entered value
            
            avg_dist = mean(del_dist(start_pos:end_pos));
            motile2 = [];
            
            if (avg_dist >= min_disp),
                motile2 = del_angle(start_pos:end_pos-1);
            else
                motile2 = NaN*ones(size(del_angle(start_pos:end_pos-1)));
            end
            
            %print first line of analysis output
            fprintf(analysis_fid, '%d\t%d\t%5.2f\t%5.2f\t%5.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%5.2f\t%s\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n', cur_data(1,1), 1, cur_data(start_pos, 3:5), ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',cur_data(1,2), ' ', sign_index(1), interaction1(1), interaction2(1), tum_interaction(1), DC_interaction(1), del_disp_3D_n(1), segmentID(1), seg_duration(1), seg_AUC(1), seg_Peak(1), seg_AvgPeak(1));
            %print second line of analysis output (no angle)
            fprintf(analysis_fid, '%d\t%d\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%s\t%5.2f\t%s\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t\n', ...
                cur_data(1,1), 2, cur_data(start_pos+1,3:5), del_dist(start_pos), del_pos(start_pos,:), inst_vel_3D(start_pos), inst_vel(start_pos,:), del_disp_3D(start_pos), del_disp(start_pos,:), ' ', inst_vel_2D(start_pos), ' ', cur_data(start_pos+1,2), del_time(start_pos), sign_index(start_pos+1), interaction1(start_pos+1), interaction2(start_pos+1), tum_interaction(start_pos+1), DC_interaction(start_pos+1), del_disp_3D_n(start_pos+1), segmentID(start_pos+1), seg_duration(start_pos+1), seg_AUC(start_pos+1), seg_Peak(start_pos+1), seg_AvgPeak(start_pos+1)); %modified 1/27/2011 to make MATLAB consider empty cells as such
            %print the rest of the analysis file
            for t = (start_pos + 1):end_pos,
                fprintf(analysis_fid, '%d\t%d\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t\n', ...
                    cur_data(1,1), t+1, cur_data(t+1,3:5), del_dist(t), del_pos(t,:), inst_vel_3D(t), inst_vel(t,:), del_disp_3D(t), del_disp(t,:), del_angle(t-1), inst_vel_2D(t), motile(t-1), cur_data(t+1,2), del_time(t), sign_index(t+1), interaction1(t+1), interaction2(t+1), tum_interaction(t+1), DC_interaction(t+1), del_disp_3D_n(t+1), segmentID(t+1), seg_duration(t+1), seg_AUC(t+1), seg_Peak(t+1), seg_AvgPeak(t+1));
            end    
            fprintf(analysis_fid,'---\t\n');
            
            xyz_del_disp = cur_data(end_pos+1,3:5) - cur_data(start_pos,3:5);
            ttl_del_disp = sqrt(sum(xyz_del_disp.^2));
            xyz_del_dist = sum(abs(del_pos(start_pos:end_pos,:)),1);
            ttl_del_dist = sum(del_dist(start_pos:end_pos));
	    ttl_del_dist_2D = sum(del_dist_2D(start_pos:end_pos));
            CI = ttl_del_disp/ttl_del_dist;
            xyz_avg_vel = mean(inst_vel(start_pos:end_pos),1); 
            ttl_avg_vel = mean (inst_vel_3D(start_pos:end_pos));
            avg_del_angle = my_nanmean(del_angle(start_pos:end_pos-1));
            %stdev_abs_angle=nanstd(abs(del_angle(start_pos:end_pos-1))); %added 6/27/2011
            % avg_2d_vel = mean(inst_vel_2D(start_pos:end_pos));
	        duration = sum(del_time(start_pos:end_pos)*cyc_time/60);  % LINE ADDED 1/16/2011
            mean_2D_vel = ttl_del_dist_2D/duration; %Line added 2/9/2011
	        mean_vel_3D = ttl_del_dist/duration;  % LINE ADDED 1/17/2011
	        arr_coeff = sum(nonmotile_arr_coeff(start_pos:end_pos)) / sum(del_time(start_pos:end_pos)); % LINE ADDED 1/17/2011
            mean_sign = mean(sign_index(start_pos:end_pos));  %added 5/31/11
            %motility_coeff = (ttl_del_disp^2)/(6*duration);    %added 5/31/11
            PercTimeSignaling = sum(seg_duration,'omitnan') / duration * 100;
            duration2 = sum(del_time(start_pos+1:end_pos)*cyc_time/60);
            EfficiencyAUC = AUCvel*AUCsign;  % added FraM 3/6/2012 
            EfficiencyAUC_dur = EfficiencyAUC/(duration2^2);
            
            % calculate and print out msd
            Thresh=10;
            CT = 8;
            alltime=[];
            alldist=[];
            coeffs=[];
            for step = 1:(end_pos-start_pos+1),
                dist = [];
                for cur_pos = start_pos:(end_pos+1-step),                
                    dist(cur_pos-start_pos+1) = sqrt(sum((cur_data(cur_pos+step, 3:5) - cur_data(cur_pos, 3:5)).^2));
                end
                fprintf(msd_fid, '%d\t%d\t%5.2f\t%5.2f\n', cur_data(1,1), sum(del_time(1:step))*cyc_time, mean(dist), mean(dist)^2);  
                alltime(step)= sqrt(sum(del_time(1:step))*cyc_time/60);
                alldist(step)= mean(dist);
                if (mean(dist)>Thresh & CT==8),
                    CT = sqrt(sum(del_time(1:step))*cyc_time/60);
                end
            end
            
            [valmax, imax] = max(alldist);  %this finds the maximum value of mean displacement
            pers = 0;
            for i=1:length(alltime),  %this finds the index of the first timepoint bigger than the persistence time
                if (alltime(i)>1 & pers==0),
                    pers = i;
                end
            end
            alltime = alltime(pers:imax);  %keeps elements from 3 (immediately after persistence time) to max
            alldist = alldist(pers:imax);
            coeffs=polyfit(alltime,alldist,1);
            fprintf(msd_fid, '-----\t\n');
           
            %print out summary line
            fprintf(summary_fid,'%d\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%s\t%5.2f', ...
                cur_data(1,1),ttl_del_dist, xyz_del_dist(1:3), mean_vel_3D, xyz_avg_vel(1:3), ttl_del_disp, xyz_del_disp(1:3), avg_del_angle, CI, mean_2D_vel, duration, arr_coeff, mean_sign, pos_classifier, PercTimeSignaling);
            if length(exceeded)>2,
                fprintf(summary_fid,'\t*');
            end
            fprintf(summary_fid,'\n');
            
            
        end
        waitbar(0.25 + 0.75*(i/length(all_data)));
        
    totsign = totsign + sign_1;
    totsignpar = totsignpar + par_sign_1;
    totsignstr = totsignstr + str_sign_1;
    
    if maxtimepoint>tpnb,
        tpnb = maxtimepoint;
    end
        
    end

    
    nc = totsign/tpnb;
    nc_PARENCHYMA = totsignpar/tpnb;
    nc_STROMA = totsignstr/tpnb;
    T = (tpnb-1) *cyc_time/60;


    varargout{1} = all_data;
    
    fclose(analysis_fid);
    fclose(summary_fid);
    fclose(msd_fid);
    close(wh);
% catch
%     fprintf('An error occurred during tracker file processing.\nLast Error: %s\n', lasterr);
%     fclose(analysis_fid);
%     fclose(summary_fid);
%     fclose(msd_fid);
%     close(wh);
% end

