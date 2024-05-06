% MODIFIED BY FRANCESCO MARANGONI 


function [varargout] = GraphAnalysis(min_inst_vel, cyc_time, varargin)

% [varargout] = GraphAnalysis(infile)
% Inputs:
%  infile - Can either be a single string of a filename to plot, or a cell
%       array of filenames.

warning off all;
num_files = length(varargin);
filelist = varargin;

% Create the figures for plotting
distrib_fig = figure; hold on;
axis square;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
distrib_fig_two = figure; hold on;
axis square;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
cell_mot_chg_fig = figure; hold on;
axis square;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
scatter_fig = figure; hold on;
axis square;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
cell_fft_fig = figure; hold on;
axis square;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
% scatter_fig_two = figure; hold on;
% set(gcf, 'color', 'white');
% set(gca, 'color', 'white');
% set(gcf, 'InvertHardCopy', 'off');
contour_plot = figure; hold on;
axis square;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
change_hist_fig = figure; hold on;
axis square;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
sample_fig = figure; hold on;
axis square;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
fft_fig = figure; hold on;
axis square;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
anglechg_fig = figure; 
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
% plots added by Francesco 1/25/2011
sign_time_plot = figure; hold on;
%axis ([0 61 0 100]);
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');  %  Ensures that hardcopies will have the same color as figures on screen
sign_vel_plot = figure; hold on;
%axis ([0 50 0 109]);
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
% plots added by Francesco 2/5/2011
sign_scheme = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
sign_vel_scheme = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
% plots added by Francesco 2/13/2011
sign_vel_scheme_2D = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
% plot added 3/2/11 Francesco Marangoni
sign_inter_scheme_2D = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
% plot added 4/11/11 Francesco Marangoni
sign_inter_scheme_CTRL = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
% plot XY tracks + clor coding FraM 25/1/2012
sample_fig_ccode = figure; hold on;
axis square;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
% Velocity as function of time, color coding tumor position
vel_time_cctumpos = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
% XY color code speed
XYspeed = figure; hold on;
axis square;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
% Velocity as function of time, color coding distance to target 1
vel_time_cctarget1 = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');

count = 0;
colors = {'r', 'g', 'b', 'm', 'k', 'y', 'c'};
cont_colors = {'r', 'g', 'b', 'c', 'm', 'k', 'y'};
hue0 = [0, 0, 1];  % palette added by Francesco Marangoni 6/13/2011 - blue 
hue1 = [0, 0.5, 1];
hue2 = [0, 1, 1];  % teal
hue3 = [0, 1, 0.5];  
hue4 = [0, 1, 0];   % green
hue5 = [0.5, 1, 0]; 
hue6 = [1, 1, 0];   % - yellow
hue7 = [1, 0.5, 0];
hue8 = [1, 0, 0];   % - red

Xmin = input('Insert Xmin');
Xmax = input('Insert Xmax');
Ymin = input('Insert Ymin');
Ymax = input('Insert Ymax');
Zmin = input('Insert Zmin');
Zmax = input('Insert Zmax');


wh = waitbar(0, 'Graphing file(s)...');

for cur_file = 1:num_files,   % cycles every open file one at a time
    clear infile;

    infile = filelist{cur_file};
    
    fid = fopen(infile);
    
    % init variables
    all_data = {};
    num_tracks = 1;
    numlines = 1;
    data = [];
    
    % ignore first line
    tline = fgetl(fid);
    
    % read data out of file
    while ~feof(fid),  
        % get the whole line
        tline = fgetl(fid); % this works fine
        if (strcmp(upper(tline(1:3)), '---')),
            % so, we are at the end of a cell's track
            all_data{num_tracks} = data; % all data is a line array in which each position is the matrix of data for a given track (timepoints, parameters)
            numlines = 1;
            data = [];
            lineofdata = {};
            num_tracks = num_tracks + 1;
            
        else
            % read in lines of data for a given track
            datacount = 1;
            lineofdata = {};
            % read in the chunks of the line
            
            while (length(tline) > 0),
                [token, tline] = strtok(tline, char(9));  %token is chunk of tline till tab (char9). Puts the remain of line into tline.
                lineofdata{datacount} = token;  % attibutes tokens (line chunks) to column arrays. HERE IS THE PROBLEM! IT DOES NOT CONSIDER EMPTY CELLS!!!  SOLVED THROUGH PUTTING A SPACE (Trackerread .m file)
                datacount = datacount+1;
            end
            % transfer the chunks into numbers - does it fine
            for i = 1:35,  % modified 5/2/17 to handle more input parameters.  MODIFY HERE TO ADD MORE PARAMETERS
                try
                    data(numlines, i) = str2num(lineofdata{i});  %translates each column in numbers. If try statement is not verified, does catch (attributes NaN)  PROBLEM IS UPSTREAM, WHEN NUMBERS ARE SEPARATED IN CHUNKS
                catch
                    data(numlines, i) = nan;  
                end
            end
            numlines = numlines + 1;
        end
    end
    
    % Plot of Distance to target vs time. Good for neutrophils 10/17/2022 ----------------------
    
    figure % draws figure of Distance as function of Time, color coding Directionality 
    hold on;
    axis([0 180 0 400]); % set axis limits
    daspect ([1 3 1]); % sets x axis 3 times as y
    title('Distance to target vs time')
    xlabel('Time (min)')
    ylabel('Distance (micron)')
   
  
    % Retrieves relevant data from all_data and plots them
    for i = 1:length(all_data),
        temp = all_data{i};
        Time3 = [];
        Dist3 = [];
        Dir3= [];
        Time3 = [Time3; temp(find(~isnan(temp(:,22))),22)]; % creates column array appending values from all tracks that are not empty in the corresponding column - col 22 is timepoint.
        Time3 = (Time3 -1) * (cyc_time/60); % Translates to actual time
        Dist3 = [Dist3; temp(find(~isnan(temp(:,25))),25)];  % col 25 is distance to target
        Dir3 = [Dir3; temp(find(~isnan(temp(:,34))),34)];    % col 34 is directionality
        Time3(1) = []; % removes first position from Time 3 and Dist 3
        Dist3(1) = [];
        
            for step = 1:(length(Time3)-1),
                
                 if Dir3(step+1)<-0.9,
                     hue = [0 0 1];  %blue for min
                 elseif ((Dir3(step+1)<-0.7) & (Dir3(step+1)>=-0.9)),
                     hue = [0 0.25 1];
                 elseif ((Dir3(step+1)<-0.5) & (Dir3(step+1)>=-0.7)),
                     hue = [0 0.5 1];
                % elseif ((Dir3(step+1)<-0.3) & (Dir3(step+1)>=-0.5)),
                %     hue=[0 0.75 1];
                 elseif ((Dir3(step+1)<0.5) & (Dir3(step+1)>=-0.5)),
                     hue = [0.5 0.5 0.5]; % gray
                % elseif ((Dir3(step+1)<0.5) & (Dir3(step+1)>=0.3)),
                %    hue=[1 0.75 0];  
                 elseif ((Dir3(step+1)<0.7) & (Dir3(step+1)>=0.5)),
                    hue= [1 0.5 0];  
                 elseif ((Dir3(step+1)<0.9) & (Dir3(step+1)>=0.7)),
                    hue=[1 0.25 0];
                 elseif ((Dir3(step+1)>=0.9)),
                    hue = [1 0 0];  % red
                 end
                 
                 if Dir3(step+1)>0.5,
                     thick = 1;
                 else
                     thick = 0.1;
                 end
                
                 plot(Time3(step:step+1), Dist3(step:step+1),'Color', hue, 'Linewidth', thick);
            end
    end
    
    hold off;
 
    % 3D plot of velocity 02/11/2021 ---------------------
    
    Xcoord2 = [];
    Ycoord2 = [];
    Zcoord2 = [];
    Speed2 = [];
    % Retrieves relevant data from all_data
    for i = 1:length(all_data),
        temp = all_data{i};
        Xcoord2 = [Xcoord2; temp(find(~isnan(temp(:,10))),3)];  % creates column array appending values from all tracks that are not empty in col 29
        Ycoord2 = [Ycoord2; temp(find(~isnan(temp(:,10))),4)];
        Zcoord2 = [Zcoord2; temp(find(~isnan(temp(:,10))),5)];
        Speed2 = [Speed2; temp(find(~isnan(temp(:,10))),10)];
    end
    % Creates a 3D scatter plot using the scatter3 function
    figure
    scatter3(Xcoord2, Ycoord2, Zcoord2, 15, Speed2, 'filled')
    axis equal; % set identical lenght for data units
    axis([Xmin Xmax Ymin Ymax Zmin Zmax]); % set axis limits
    view(2)  % another nice view is view(-34, 14). if orthogonal view is preferred, put view(2)
    % Add title and axis labels
    title('3D Velocity')
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid off
    % Add a colorbar with tick labels
    colormap(flipud(jet));  %sets colormap LUT to "hot", red first then yellow.  "Jet"  Gives the whole spectrum from cold to hot. flipud flips the colormap.
    caxis([0 10]);
    colorbar;
    %colorbar('Location', 'EastOutside', 'YTickLabel',...
    %{'4 micron', '8 micron', '12 micron', '16 micron', ...
    %'20 micron', '24 micron', '28 micron'});
    % END OF 3D plot of 3D Velocity ----------------
    
    % 3D plot of n-point displacement 5/2/2017 ---------------------
    
    Xcoord = [];
    Ycoord = [];
    Zcoord = [];
    Magnet = [];
    % Retrieves relevant data from all_data
    for i = 1:length(all_data),
        temp = all_data{i};
        Xcoord = [Xcoord; temp(find(~isnan(temp(:,29))),3)];  % creates column array appending values from all tracks that are not empty in col 29
        Ycoord = [Ycoord; temp(find(~isnan(temp(:,29))),4)];
        Zcoord = [Zcoord; temp(find(~isnan(temp(:,29))),5)];
        Magnet = [Magnet; temp(find(~isnan(temp(:,29))),29)];
    end
    % Creates a 3D scatter plot using the scatter3 function
    figure
    scatter3(Xcoord, Ycoord, Zcoord, 25, Magnet, 'filled')
    axis equal; % set identical lenght for data units
    axis([Xmin Xmax Ymin Ymax Zmin Zmax]); % set axis limits
    view(2)  % another nice view is view(-34, 14). if orthogonal view is preferred, put view(2)
    % Add title and axis labels
    title('n-minutes displacement')
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid off
    % Add a colorbar with tick labels
    colormap(flipud(jet));  %sets colormap LUT to "hot", red first then yellow.  "Jet"  Gives the whole spectrum from cold to hot. flipud flips the colormap.
    caxis([0 50]);
    colorbar;
    %colorbar('Location', 'EastOutside', 'YTickLabel',...
    %{'4 micron', '8 micron', '12 micron', '16 micron', ...
    %'20 micron', '24 micron', '28 micron'});
    % END OF 3D plot of n-point displacement ----------------
    
    
    % Scheme of signaling/interaction vs. time by tracks in 2D Francesco Marangoni 3/2/2011
      
    grey = [0.9 0.9 0.9];
    smooth_sig = input('SI Smoothing ');
    smooth_vel = input('3D Velocity Smoothing ');
    max_distance = input('Insert the radius for target scanning (microns)');
    thr = input('Insert the threshold for interaction (microns)');
    gap = 0.3; %this is the distance between 2 consecutive traces FraM 4/11/11
    for i = 1:length(all_data),
        temp = all_data{i};
        time = [];
        time = [time; temp(find(~isnan(temp(:,23))),23)];
        
        signaling = [];
        interaction = [];
        interaction2 = [];
        el_time = zeros(length(time)+1); 
        ranking = zeros(length(time)+1);
        threshold = zeros(length(time)+1);
        signhue = zeros(1,3);
        
        signaling = [signaling; temp(find(~isnan(temp(:,24))),24)];
        interaction = [interaction; temp(find(~isnan(temp(:,25))),25)];
        interaction2 = [interaction2; temp(find(~isnan(temp(:,26))),26)];
        
        ranking(:) = (1+gap)*i; %modified FraM 4/11/11
        threshold(:) = ((1+gap)*i)-(thr/max_distance);
        interaction(:) = (((1+gap)*i)-(interaction(:)./max_distance));
        interaction2(:) = (((1+gap)*i)-(interaction2(:)./max_distance));
        
       
        %mainteraction = NaN(length(interaction),1);
        %mainteraction2 = NaN(length(interaction2),1);
        %masignaling = NaN(length(signaling),1);
        %for k = (1+smooth):(length(mainteraction)-smooth),
        %    mainteraction(k) = mean(interaction((k-smooth):(k+smooth)));
        %    mainteraction2(k) = mean(interaction2((k-smooth):(k+smooth)));
        %    masignaling(k) = nanmean(signaling((k-smooth):(k+smooth)));
        %end
       
        mainteraction = movmean(interaction, 1);
        mainteraction2 = movmean(interaction2, 1);
        masignaling = movmean(signaling, (2*smooth_sig+1));
        
        el_time(1) = 0;
        for j = 2:length(el_time),
            el_time(j) = sum(time(1:(j-1)))*cyc_time/60;
        end
        
        %plot first figure signaling/interaction with Ag
        figure(sign_inter_scheme_2D);
        rectangle('Position', [0, ((1+gap)*i)-1, el_time(length(el_time)), 1], 'EdgeColor', grey, 'FaceColor', grey);
        
        for step = 1:(length(el_time)-1),
            if (interaction(step) < ((1+gap)*i)-1) | (interaction(step+1) < ((1+gap)*i)-1),
                continue
            else
            plot(el_time(step:step+1), mainteraction(step:step+1),'Color', 'k', 'Linewidth', 1);
            end
        end
        
        for step = 1:(length(el_time)),
            if interaction(step) < ((1+gap)*i)-1,
                continue
            else
            if masignaling(step)<= 0.45,
                 signhue = hue0;
            elseif ((masignaling(step)>0.45) & (masignaling(step)<=0.493)),
                 signhue = hue1;
            elseif ((masignaling(step)>0.493) & (masignaling(step)<=0.536)),
                 signhue = hue2;
            elseif ((masignaling(step)>0.536) & (masignaling(step)<=0.579)),
                 signhue = hue3; 
            elseif ((masignaling(step)>0.579) & (masignaling(step)<=0.622)),
                 signhue = hue4;
            elseif ((masignaling(step)>0.622) & (masignaling(step)<=0.665)),
                 signhue = hue5;
            elseif ((masignaling(step)>0.665) & (masignaling(step)<=0.708)),
                 signhue = hue6;
            elseif ((masignaling(step)>0.708) & (masignaling(step)<=0.75)),
                 signhue = hue7;
            else
                 signhue = hue8;
            end
            end
            plot(el_time(step), mainteraction(step), 'Marker', 'o', 'MarkerFaceColor', signhue, 'MarkerEdgeColor', signhue, 'MarkerSize', 4);
        end
        
        
        
        line([0 el_time(length(el_time))], [((1+gap)*i)-(thr/max_distance) ((1+gap)*i)-(thr/max_distance)], [0 0], 'Linestyle', '--', 'Color', 'k');
        %rectangle('Position', [0.2, ((1+gap)*(i-1)+0.2),
        %el_time(length(el_time)), gap-0.2], 'EdgeColor', 'w', 'FaceColor', 'w');
        
        
        %plot first figure signaling/interaction with CTRL
        figure(sign_inter_scheme_CTRL);
        rectangle('Position', [0, ((1+gap)*i)-1, el_time(length(el_time)), 1], 'EdgeColor', grey, 'FaceColor', grey);
        
        for step = 1:(length(el_time)-1),
             if (interaction2(step) < ((1+gap)*i)-1) | (interaction2(step+1) < ((1+gap)*i)-1),
                continue
             else
                plot(el_time(step:step+1), mainteraction2(step:step+1),'Color', 'k', 'Linewidth', 1);
             end
        end
        
        for step = 1:(length(el_time)),
            if interaction2(step) < ((1+gap)*i)-1,
                continue
            elseif (interaction(step)>=threshold(step)),  
                signhue = [0.5 0.5 0.5];
            else    
           if masignaling(step)<= 0.45,
                 signhue = hue0;
            elseif ((masignaling(step)>0.45) & (masignaling(step)<=0.493)),
                 signhue = hue1;
            elseif ((masignaling(step)>0.493) & (masignaling(step)<=0.536)),
                 signhue = hue2;
            elseif ((masignaling(step)>0.536) & (masignaling(step)<=0.579)),
                 signhue = hue3; 
            elseif ((masignaling(step)>0.579) & (masignaling(step)<=0.622)),
                 signhue = hue4;
            elseif ((masignaling(step)>0.622) & (masignaling(step)<=0.665)),
                 signhue = hue5;
            elseif ((masignaling(step)>0.665) & (masignaling(step)<=0.708)),
                 signhue = hue6;
            elseif ((masignaling(step)>0.708) & (masignaling(step)<=0.75)),
                 signhue = hue7;
            else
                 signhue = hue8;
            end
              
            end
            
            plot(el_time(step), mainteraction2(step), 'Marker', 'o', 'MarkerFaceColor', signhue, 'MarkerEdgeColor', signhue, 'MarkerSize', 4);
           
        end
        
        
        
        line([0 el_time(length(el_time))], [((1+gap)*i)-(thr/max_distance) ((1+gap)*i)-(thr/max_distance)], [0 0], 'Linestyle', '--', 'Color', 'k');
        %rectangle('Position', [0.2, ((1+gap)*(i-1)+0.2),
        %el_time(length(el_time)), gap-0.2], 'EdgeColor', 'w', 'FaceColor', 'w');
        
    end
    
    figure(sign_inter_scheme_2D);
    set(gca, 'YTickLabel', [-(max_distance); 0]);
    set(gca, 'YTick', [gap:1:(gap+1)]);
    set(gca, 'Tickdir', 'out');
    axis ([0, inf, 0, ((1+gap)*length(all_data)+1)]);
    xlabel('Time (min)');
    ylabel('Distance to target (\mum)');
    title ('Signaling/Distance vs. Time in 2D');
    
    
    figure(sign_inter_scheme_CTRL);
    set(gca, 'YTickLabel', [-(max_distance); 0]);
    set(gca, 'YTick', [gap:1:(gap+1)]);
    set(gca, 'Tickdir', 'out');
    axis ([0, inf, 0, ((1+gap)*length(all_data)+1)]);
    xlabel('Time (min)');
    ylabel('Distance to target (\mum)');
    title ('Signaling/Distance vs. Time - APC2');
   
    % Scheme of signaling in all tracks by color coding  Francesco Marangoni 5/27/2011
    figure(sign_scheme);  % opens signaling scheme figure
    
    
    for i = 1:length(all_data),
        time = [];
        signaling = [];
        tumor = [];
        signhue = zeros(1,3);
        temp = all_data{i};
        time = [time; temp(find(~isnan(temp(:,22))),22)];
        signaling = [signaling; temp(find(~isnan(temp(:,24))),24)];
        tumor = [tumor; temp(find(~isnan(temp(:,27))),27)];
        ranking = zeros(length(time)+1);
        ranking(:) = i-1; 
        
        masignaling = movmean (signaling, (2*smooth_sig+1));
        
        for step = 1:((length(masignaling)-1)),
            plot(time(step:step+1), ranking(step:step+1),'Color', 'k', 'Linewidth', 1);
        end
        
        for step = 1:((length(masignaling)-1)),
          if masignaling(step)<= 0.45,
                 signhue = hue0;
            elseif ((masignaling(step)>0.45) & (masignaling(step)<=0.493)),
                 signhue = hue1;
            elseif ((masignaling(step)>0.493) & (masignaling(step)<=0.536)),
                 signhue = hue2;
            elseif ((masignaling(step)>0.536) & (masignaling(step)<=0.579)),
                 signhue = hue3; 
            elseif ((masignaling(step)>0.579) & (masignaling(step)<=0.622)),
                 signhue = hue4;
            elseif ((masignaling(step)>0.622) & (masignaling(step)<=0.665)),
                 signhue = hue5;
            elseif ((masignaling(step)>0.665) & (masignaling(step)<=0.708)),
                 signhue = hue6;
            elseif ((masignaling(step)>0.708) & (masignaling(step)<=0.75)),
                 signhue = hue7;
            else
                 signhue = hue8;
            end
            
            if tumor(step)==0,
                symbol = '^';
            else
                symbol = 'o';
            end
                
            plot(time(step), ranking(step), 'Marker', symbol, 'MarkerFaceColor', signhue, 'MarkerEdgeColor', signhue, 'MarkerSize', 5);
            
        end
        
    end
    set(gca, 'YTick', [0:length(all_data)]);
    set(gca, 'Tickdir', 'out');
    xlabel('Timepoint');
    ylabel('Track number');
    title ('Schematic view signaling vs. time');
    
    
    
    %---------------------
    % Scheme of signaling/velocity vs. time by tracks in 2D Francesco Marangoni 2/13/2011 CHANGE HERE
    figure(sign_vel_scheme_2D);  % opens signaling scheme 2D figure
   
    grey = [0.9 0.9 0.9];
    max_velocity = input('Insert the highest instantaneous velocity (microns/min)');
    for i = 1:length(all_data),
        time = [];
        signaling = [];
        mod_velocity = [];
        tum1 = [];
        DC1 = [];
        signhue = zeros(1,3);
        temp = all_data{i};
        time = [time; temp(find(~isnan(temp(:,23))),23)];
        signaling = [signaling; temp(find(~isnan(temp(:,24))),24)];
        signaling(1) = [];  %removes first element of vector
        mod_velocity = [mod_velocity; temp(find(~isnan(temp(:,10))),10)];
        mod_velocity(:) = ((1+gap)*i+(mod_velocity(:)./max_velocity)-1);
        tum1 = [tum1; temp(find(~isnan(temp(:,27))),27)];
        DC1 = [DC1; temp(find(~isnan(temp(:,28))),28)];
        el_time = zeros(length(time)); 
        ranking = zeros(length(time));
        ranking(:) = 2*i; 
        
        mamodvelocity = movmean (mod_velocity, (2*smooth_vel+1));
        masignaling = movmean (signaling, (2*smooth_sig+1));
               
        for j = 1:length(time),
            el_time(j) = sum(time(1:j))*cyc_time/60;
        end
        
        rectangle('Position', [0, ((1+gap)*i)-1, el_time(length(el_time)), 1], 'EdgeColor', grey, 'FaceColor', grey);
        
        for step = 1:(length(el_time)-1),
            plot(el_time(step:step+1), mamodvelocity(step:step+1),'Color', 'k', 'Linewidth', 1);
        end
        
        for step = 1:(length(el_time)),
            if masignaling(step) == -1,
                signhue = [0.4 0.4 0.4];  %greys out errors
            elseif ((masignaling(step)>0) & (masignaling(step)<=0.45)),
                 signhue = hue0;
            elseif ((masignaling(step)>0.45) & (masignaling(step)<=0.493)),
                 signhue = hue1;
            elseif ((masignaling(step)>0.493) & (masignaling(step)<=0.536)),
                 signhue = hue2;
            elseif ((masignaling(step)>0.536) & (masignaling(step)<=0.579)),
                 signhue = hue3; 
            elseif ((masignaling(step)>0.579) & (masignaling(step)<=0.622)),
                 signhue = hue4;
            elseif ((masignaling(step)>0.622) & (masignaling(step)<=0.665)),
                 signhue = hue5;
            elseif ((masignaling(step)>0.665) & (masignaling(step)<=0.708)),
                 signhue = hue6;
            elseif ((masignaling(step)>0.708) & (masignaling(step)<=0.75)),
                 signhue = hue7;
            else
                 signhue = hue8;
            end
            
                     
            
            if tum1(step)==0,
                symbol = '^';
            else
                symbol = 'o';
            end
            
            if DC1(step)==0,
                fill = [1 1 1]; %if no DC, fill is white
            else
                fill = signhue;
            end    
            
            plot(el_time(step), mamodvelocity(step), 'Marker', symbol, 'MarkerFaceColor', fill, 'MarkerEdgeColor', signhue, 'MarkerSize', 5);
            
        end
        
        rectangle('Position', [0.2, ((1+gap)*(i-1)+0.1), el_time(length(el_time)), gap-0.2], 'EdgeColor', 'w', 'FaceColor', 'w');
    end
    set(gca, 'YTickLabel', [0; max_velocity]);
    set(gca, 'YTick', [gap:1:(1+gap)]);
    set(gca, 'Tickdir', 'out');
    axis ([0, inf, 0, ((1+gap)*length(all_data)+1)]);
    xlabel('Time (min)');
    ylabel('Instantaneous velocity (\mum/sec)');
    title ('Schematic view signaling/Velocity vs. Time in 2D');
    
    %---------------------
    
    
    % Graph of CTL velocity vs. time, color coding position to tumor 4/4/2012
    
    figure(vel_time_cctumpos);
   
    grey = [0.9 0.9 0.9];
    for i = 1:length(all_data),
        time = [];
        mod_velocity = [];
        tum1 = [];
        hue = zeros(1,3);
        temp = all_data{i};
        time = [time; temp(find(~isnan(temp(:,23))),23)];
        mod_velocity = [mod_velocity; temp(find(~isnan(temp(:,10))),10)];
        mod_velocity(:) = ((1+gap)*i+(mod_velocity(:)./max_velocity)-1);
        tum1 = [tum1; temp(find(~isnan(temp(:,27))),27)];
        tum1(1) = [];  % gets rid of first element of tum array
        el_time = zeros(length(time)); 
        ranking = zeros(length(time));
        ranking(:) = 2*i; 
        
        mamodvelocity = movmean (mod_velocity, (2*smooth_vel+1));
              
        
        for j = 1:length(time),
            el_time(j) = sum(time(1:j))*cyc_time/60;
        end
        
        rectangle('Position', [0, ((1+gap)*i)-1, el_time(length(el_time)), 1], 'EdgeColor', grey, 'FaceColor', grey);
        
        for step = 1:(length(el_time)-1),
            plot(el_time(step:step+1), mamodvelocity(step:step+1),'Color', 'k', 'Linewidth', 1);
        end
        
        for step = 1:(length(el_time)),
            
            if tum1(step)==0,
                hue = hue0;  %stroma is blue
            else
                hue = hue8;  %parenchyma is red
            end
            
            plot(el_time(step), mamodvelocity(step), 'Marker', 'o', 'MarkerFaceColor', hue, 'MarkerEdgeColor', hue, 'MarkerSize', 4);
            
        end
        
        rectangle('Position', [0.2, ((1+gap)*(i-1)+0.1), el_time(length(el_time)), gap-0.2], 'EdgeColor', 'w', 'FaceColor', 'w');
    end
    set(gca, 'YTickLabel', [0; max_velocity]);
    set(gca, 'YTick', [gap:1:(1+gap)]);
    set(gca, 'Tickdir', 'out');
    axis ([0, inf, 0, ((1+gap)*length(all_data)+1)]);
    xlabel('Time (min)');
    ylabel('Instantaneous velocity (\mum/sec)');
    title ('Schematic view Velocity vs. Time color coding tumor position in 2D');

    %----------------------
    %----------------------  Velocity vs. Time, color coding distance to target FraM 8/2/2012
    figure(vel_time_cctarget1);
   
    grey = [0.9 0.9 0.9];
    for i = 1:length(all_data),
        time = [];
        mod_velocity = [];
        targ1 = [];
        hue = zeros(1,3);
        temp = all_data{i};
        time = [time; temp(find(~isnan(temp(:,23))),23)];
        mod_velocity = [mod_velocity; temp(find(~isnan(temp(:,10))),10)];
        mod_velocity(:) = ((1+gap)*i+(mod_velocity(:)./max_velocity)-1);
        targ1 = [targ1; temp(find(~isnan(temp(:,25))),25)];
        targ1(1) = [];  % gets rid of first element of targ1 array
        el_time = zeros(length(time)); 
        ranking = zeros(length(time));
        ranking(:) = 2*i; 
        
        mamodvelocity = movmean (mod_velocity, (2*smooth_vel+1));
       
        for j = 1:length(time),
            el_time(j) = sum(time(1:j))*cyc_time/60;
        end
        
        rectangle('Position', [0, ((1+gap)*i)-1, el_time(length(el_time)), 1], 'EdgeColor', grey, 'FaceColor', grey);
        
        for step = 1:(length(el_time)-1),
            plot(el_time(step:step+1), mamodvelocity(step:step+1),'Color', 'k', 'Linewidth', 1);
        end
        
        for step = 1:(length(el_time)),
            
            if targ1(step)>26,
                hue = hue0;  %dist>10um is blue
            elseif ((targ1(step)>24) & (targ1(step)<=26)),
                hue = hue1;
            elseif ((targ1(step)>22) & (targ1(step)<=24)),
                hue = hue2;
            elseif ((targ1(step)>20) & (targ1(step)<=22)),
                hue=hue3;
            elseif ((targ1(step)>18) & (targ1(step)<=20)),
                hue=hue4;
            elseif ((targ1(step)>16) & (targ1(step)<=18)),
                hue=hue5;  
            elseif ((targ1(step)>14) & (targ1(step)<=16)),
                hue=hue6;  
            elseif ((targ1(step)>12) & (targ1(step)<=14)),
                hue=hue7;
            else
                hue = hue8;  %dist<12um is red
            end
            
            plot(el_time(step), mamodvelocity(step), 'Marker', 'o', 'MarkerFaceColor', hue, 'MarkerEdgeColor', hue, 'MarkerSize', 4);
            
        end
        
        rectangle('Position', [0.2, ((1+gap)*(i-1)+0.1), el_time(length(el_time)), gap-0.2], 'EdgeColor', 'w', 'FaceColor', 'w');
    end
    set(gca, 'YTickLabel', [0; max_velocity]);
    set(gca, 'YTick', [gap:1:(1+gap)]);
    set(gca, 'Tickdir', 'out');
    axis ([0, inf, 0, ((1+gap)*length(all_data)+1)]);
    xlabel('Time (min)');
    ylabel('Instantaneous velocity (\mum/sec)');
    title ('Schematic view Velocity vs. Time color coding Distance to target1 in 2D');

    %----------------------
        
    % Scheme of signaling vs. velocity vs. time 3D Francesco Marangoni
    % 2/6/2011
    figure(sign_vel_scheme);  % opens signaling velocity scheme figure
    
    
    for i = 1:length(all_data),
        time1 = [];
        signaling1 = [];
        velocity1 = [];
        signhue1 = zeros(1,3);
        temp = all_data{i};
        velocity1 = [velocity1; temp(find(~isnan(temp(:,10))),10)];
        time1 = [time1; temp(find(~isnan(temp(:,23))),23)];
        signaling1 = [signaling1; temp(find(~isnan(temp(:,24))),24)];
        el_time1 = zeros(length(time1));
        ranking1 = zeros(length(time1));
        ranking1(:) = i; 
        
        
        for j = 1:length(time1),
            el_time1(j) = sum(time1(1:j))*cyc_time/60;
        end
        
        for step = 1:(length(el_time1)-1),
            if signaling1(step)<= 0.45,
                 signhue1 = hue0;
            elseif ((signaling1(step)>0.45) & (signaling1(step)<=0.493)),
                 signhue1 = hue1;
            elseif ((signaling1(step)>0.493) & (signaling1(step)<=0.536)),
                 signhue1 = hue2;
            elseif ((signaling1(step)>0.536) & (signaling1(step)<=0.579)),
                 signhue1 = hue3; 
            elseif ((signaling1(step)>0.579) & (signaling1(step)<=0.622)),
                 signhue1 = hue4;
            elseif ((signaling1(step)>0.622) & (signaling1(step)<=0.665)),
                 signhue1 = hue5;
            elseif ((signaling1(step)>0.665) & (signaling1(step)<=0.708)),
                 signhue1 = hue6;
            elseif ((signaling1(step)>0.708) & (signaling1(step)<=0.75)),
                 signhue1 = hue7;
            else
                 signhue1 = hue8;
            end
            %plot(el_time1(step:(step+1)), velocity1(step:(step+1)), 'Color', signhue1, 'Linewidth', 1);
            plot3(ranking1(step:(step+1)), el_time1(step:(step+1)), velocity1(step:(step+1)), 'Color', signhue1, 'Linewidth', 1);
            view(-189, 20);
        end
    end
    %axis ([0, 60, 0, 60, 0, (length(all_data)+1)]);
    axis square;
    ylabel('Time (min)');
    zlabel('Instantaneous velocity (\mum/s)');
    xlabel('Track Number');
    title ('Schematic view signaling/velocity vs. time');







    % 3D velocity histogram
    x = 0:2:40;   %defines bins
    y = [];
    for i = 1:length(all_data),  % i is a track indicator
        temp = all_data{i};
        y = [y; temp(find(~isnan(temp(:,10))),10)]; % puts temp(:,10) in temp(:,10) if it is a number. Concatenates an array y with that info.
    end
%     [phat, pci] = gamfit(y);
%     [m,v] = gamstat(phat(1), phat(2));
%     [lambdahat,lambdaci] = poissfit(y);
%     [weibp, weibci] = weibfit(y);
%     fprintf('Distribution %d: %f (%f - %f); %f (%f - %f)\n', cur_file, phat(1), pci(1,:), phat(2), pci(2,:));
%     fprintf('\tMean - %f; Variance - %f\n', m, v);
%     fprintf('\tLambda - %f (%f - %f)\n', lambdahat, lambdaci);
%     fprintf('\tWeibull - %f (%f - %f); %f (%f - %f)\n', weibp(1), weibci(1,:), weibp(2), weibci(2,:));
    for_ttest{cur_file} = y;
    n = hist(y,x);  %designs histogram
    distrib_all(cur_file, :) = n/length(y);
        
    % Scatterplot instantaneous velocity x versus instantaneous velocity y
    figure(scatter_fig);
    x = [];
    y = [];
    for i = 1:length(all_data),
        temp = all_data{i};
        y = [y; temp(find(~isnan(temp(:,12))),12)];
        x = [x; temp(find(~isnan(temp(:,11))),11)];
    end
    plot(x,y, '.', 'Color', colors{mod(cur_file-1,7)+1});
    xlabel('Inst. velocity in x (\mum/s)');
    ylabel('Inst. velocity in y (\mum/s)');
    %title(sprintf('Instantaneous velocity in x vs. Instantaneous velocity in y - %s', strrep(strtok(infile,'.'), '_', '\_')));
    title('Instantaneous velocity in x vs. Instantaneous velocity in y');
    axis([-40 40 -40 40]);
    grid on;

    % Scatterplot Inst Vel 3D vs. Signaling Index  Francesco Marangoni
    figure(sign_vel_plot);
    x1 = [];
    y1 = [];
    for i = 1:length(all_data),  % puts together velocity with signaling 
        temp = all_data{i};
        y1 = [y1; temp(find(~isnan(temp(:,24))),24)];
        y1(1) = [];  %removes first element of vector
        x1 = [x1; temp(find(~isnan(temp(:,10))),10)];
    end
    
    max1 = movmean (x1, (2*smooth_vel+1));
    may1 = movmean (y1, (2*smooth_sig+1));
    
    
    plot (max1,may1,'bo','MarkerFaceColor','b','MarkerSize',4);
    xlabel('Inst. 3D velocity (\mum/s)');
    ylabel('Signaling Index');
    title ('Instantaneous 3D velocity vs. Signaling');
    grid on;
    
%     figure(scatter_fig_two);
%    [xmu, xsig, xmuci, xsigci] = normfit(x);
%   [ymu, ysig, ymuci, ysigci] = normfit(y);
    xplot = -40:2:40;
    yplot = -40:2:40;
    clear zz;
    clear z;
%     for g = 1:length(xplot),
%         for h = 1:length(yplot),
%             z(g,h) = normpdf(xplot(g), xmu, xsig) * normpdf(yplot(h), ymu, ysig);
%         end
%     end
    for g = 1:length(xplot)-1,
        for h = 1:length(yplot)-1,
            tempx = logical(x >= xplot(g)) & logical(x <= xplot(g+1));
            tempy = logical(y >= yplot(h)) & logical(y <= yplot(h+1));
            zz(g,h) = sum(ismember(find(tempx), find(tempy)));
        end
    end
%     [tx, ty] = meshgrid(xplot, yplot);
%     contour3(tx, ty, z);
%     xlabel('Inst. velocity in x (\mum/s)');
%     ylabel('Inst. velocity in y (\mum/s)');
%     %title(sprintf('Instantaneous velocity in x vs. Instantaneous velocity in y - %s', strrep(strtok(infile,'.'), '_', '\_')));
%     title('Instantaneous velocity in x vs. Instantaneous velocity in y -- Gaussian Fit');
%     axis([-40 40 -40 40]);
%     grid on;
    
    %contour plot
    figure(contour_plot);
    [tx, ty] = meshgrid(-39:2:39, -39:2:39);
    contour3(tx, ty, zz, 20, cont_colors{mod(cur_file-1, 7)+1});
    xlabel('Inst. velocity in x (\mum/s)');
    ylabel('Inst. velocity in y (\mum/s)');
    %title(sprintf('Instantaneous velocity in x vs. Instantaneous velocity in y - %s', strrep(strtok(infile,'.'), '_', '\_')));
    title('Instantaneous velocity in x vs. Instantaneous velocity in y');
    axis([-40 40 -40 40]);
    grid on;
    
    
    %2D vel hist data generation
    hx = 0:2:40;
    hist_data_2d{cur_file} = sqrt(x.^2 + y.^2);
    n = hist(sqrt(x.^2 + y.^2), hx);
    distrib_2d(cur_file, :) = n/length(y);    
    
    %Absolute Angle Change histogram (in progress)
    %figure;
    angle = [];
    for i = 1:length(all_data),
        temp = all_data{i};
        for j = 1:size(temp,1),
            if (((temp(j,7)) == 0) || (temp(j,8)) == 0),
                angle = [angle; 0];
            else
                temp_t = atan(temp(j,8)/temp(j,7));
                temp_t = (temp_t/pi)*180;
                if (~isnan(temp_t)),
                    angle = [angle; temp_t];
                end
            end
        end
        
    end
    x =-180:10:180;
    n = hist(angle, x);
    %bar(n)
    %title('Absolute Angle Histogram');
    %output absolute angle array
    %filename = strcat(infile,'absangle.txt');
    %dlmwrite(filename, n, '\t');
    
    
    % Angle Change histogram
    figure(change_hist_fig);
    x = -180:10:180;
    y = [];
    for i = 1:length(all_data),
        temp = all_data{i};
        for j = 1:size(temp,1),
            if ~isnan(temp(j,20)),
                y = [y; temp(j,20)];
            end
        end
    end
    n = hist(y, x);
    all_change(cur_file, :) = n/length(y);
    
 


    
    %polar plot of angle change
    figure(anglechg_fig);
    dir = (x * pi/180);
    polar(dir, n, colors{mod(cur_file-1, 7)+1})
    hold on;
    title('Angle Change');
    
    % Representative Tracks (need to plot w/all starting from 0)
    figure(sample_fig);
    hold on;
    x = [];
    y = [];
    axis ([-100 100 -100 100]);
    set (gca, 'FontSize', 18, 'LineWidth', 2);
    for i = 1:length(all_data),
        temp = all_data{i};
        y = [0; temp(:,7)];  
        x = [0; temp(:,8)];
        track = zeros(size(temp,1), 2);
        for j = 2:size(temp,1),
            track(j,1) = track(j-1,1)+temp(j,7);  %NOTE x and y were swapped in original script
            track(j,2) = track(j-1,2)+temp(j,8);
        end    
        if (num_files == 1),
            plot(track(:,1), track(:,2), 'Color', colors{mod(count, 7)+1}, 'LineWidth', 1);
        else
            plot(track(:,1), track(:,2), 'Color', colors{mod(cur_file-1, 7)+1}, 'LineWidth', 1);
        end
        count = count + 1;

        xlabel('\Deltax (\mum)');
        ylabel('\Deltay (\mum)');
        title('Cell Tracks');


    end
    [legend_h,object_h,plot_h,text_strings] = legend(strrep(filelist,'_','\_'), 'Location', 'SouthOutside');
    for i = 1:length(filelist),
        set(object_h(i), 'FontSize', 12);
    end
    for i = 1:length(plot_h),
        set(plot_h(i), 'Color', colors{mod(i-1, 7)+1});
    end
   
    % XY SCHEME + COLOR CODING (FraM 25/1/2012)
    figure(sample_fig_ccode);
    hold on;
    axis ([Xmin Xmax Ymin Ymax]);
    set (gca, 'FontSize', 10, 'LineWidth', 2, 'Tickdir', 'out');
    for i = 1:length(all_data),
        x = [];
        y = [];
        signaling = [];
        signhue = zeros(1,3);
        temp = all_data{i};
        signaling = [signaling; temp(find(~isnan(temp(:,24))),24)];
        y = [y; temp(find(~isnan(temp(:,4))),4)];
        x = [x; temp(find(~isnan(temp(:,3))),3)];
        %y = Ymax-y;
        masignaling1 = movmean(signaling, (2*smooth_sig+1));
    
    for step = 1:(length(masignaling1)-1),
            if masignaling1(step) == -1,
                signhue = [0.4 0.4 0.4];  %greys out errors
            elseif ((masignaling1(step)>0) & (masignaling1(step)<=0.45)),
                 signhue = hue0;
            elseif ((masignaling1(step)>0.45) & (masignaling1(step)<=0.493)),
                 signhue = hue1;
            elseif ((masignaling1(step)>0.493) & (masignaling1(step)<=0.536)),
                 signhue = hue2;
            elseif ((masignaling1(step)>0.536) & (masignaling1(step)<=0.579)),
                 signhue = hue3; 
            elseif ((masignaling1(step)>0.579) & (masignaling1(step)<=0.622)),
                 signhue = hue4;
            elseif ((masignaling1(step)>0.622) & (masignaling1(step)<=0.665)),
                 signhue = hue5;
            elseif ((masignaling1(step)>0.665) & (masignaling1(step)<=0.708)),
                 signhue = hue6;
            elseif ((masignaling1(step)>0.708) & (masignaling1(step)<=0.75)),
                 signhue = hue7;
            else
                 signhue = hue8;
            end
            
            plot(x(step:step+1), y(step:step+1), 'Color', signhue, 'Linewidth', 2);
        end
   
        title('Cell Tracks with signaling');
        xlabel('x (\mum)');
        ylabel('y (\mum)');

    end
    
    %----------------
    %XY scheme color coding speed FraM 8/1/2012
    figure(XYspeed);
    hold on;
    axis ([Xmin Xmax Ymin Ymax]);
    set (gca, 'FontSize', 10, 'LineWidth', 2, 'Tickdir', 'out');
    
       for i = 1:length(all_data),
        x = [];
        y = [];
        speed = [];
        signhue = zeros(1,3);
        temp = all_data{i};
        speed = [speed; temp(find(~isnan(temp(:,10))),10)];
        y = [y; temp(find(~isnan(temp(:,4))),4)];
        x = [x; temp(find(~isnan(temp(:,3))),3)];
        y(1) = [];
        x(1) = [];
        %y = Ymax-y;
    
        for step = 1:(length(speed)-1),
            if speed(step) == -1,
                signhue = [0.4 0.4 0.4];  %greys out errors
            elseif ((speed(step)>9) & (speed(step)<=100)),
                 signhue = hue0;
            elseif ((speed(step)>8) & (speed(step)<=9)),
                 signhue = hue1;
            elseif ((speed(step)>7) & (speed(step)<=8)),
                 signhue = hue2;
            elseif ((speed(step)>6) & (speed(step)<=7)),
                 signhue = hue3; 
            elseif ((speed(step)>5) & (speed(step)<=6)),
                 signhue = hue4;
            elseif ((speed(step)>4) & (speed(step)<=5)),
                 signhue = hue5;
            elseif ((speed(step)>3) & (speed(step)<=4)),
                 signhue = hue6;
            elseif ((speed(step)>2) & (speed(step)<=3)),
                 signhue = hue7;
            else
                 signhue = hue8;
            end
            
            plot(x(step:step+1), y(step:step+1), 'Color', signhue, 'Linewidth', 2);
        end
    
        title('Cell Tracks with speed');
        xlabel('x (\mum)');
        ylabel('y (\mum)');

    end
    
    
    %----------------
    
    
    
    
    %Calculate fft for each cell
    figure(cell_fft_fig);
    max = 0;
    for i = 1:length(all_data),
        if (max<size(all_data{i},1)),
            max = size(all_data{i},1);
        end
    end
    max = ceil(max/2)*2;
    count = 0;
    clear fft_cell_array fft_cell_toplot_array fft_array fft_toplot_array;
    for i = 1:length(all_data),
        if (size(all_data{i},1) >= 10),
            count = count + 1;
            temp = all_data{i};
            vel_temp = sum(temp(2:end,11:12).^2,2)';
            fft_array(count, :) = zeros(1, max);
            fft_array(count, 1:length(vel_temp)) = detrend(vel_temp);
        end
    end
    for i = 1:size(fft_array,1),
        Y = fft(fft_array(i,:), max);
        Pyy = abs(Y)/max;
        fft_toplot_array(i,:) = Pyy;
    end    
    % THIS ASSUMES YOUR CYCLE TIME IS 15s (so, a minute requires *4)
    % ASSUMPTION NOW WRONG, TO BE FIXED IF THIS ANALYSIS BECOMES VALUABLE
    % (Francesco Marangoni 1/27/2011
    x = (0:max/2)/max*4;
    plot(x,fft_toplot_array(:,1:(max/2)+1), 'Color', colors{mod(cur_file-1, 7)+1}, 'LineWidth', 1.25);   
    xlabel('Frequency (cycles/min)');
    ylabel('Power');
    title('Power Spectrum of Each Cell''s Velocity Change over Time');
    
    %output single cell fft
    k = findstr(infile,'/');
    if isempty(k),
        infile = infile(1:end-12);
    else
        infile = infile(k(end)+1:end-12);
    end
    filename = strcat(infile,'sglcell_fft.txt');
    dlmwrite(filename, x, 'delimiter','\t');
    dlmwrite(filename, fft_toplot_array(:,1:(max/2)+1), '-append', 'delimiter','\t')

    
    
    % Calculate fft of inst vel x and y for whole files 
    figure(fft_fig);
    max = 0;
    for i = 1:length(all_data),
        if (max<size(all_data{i},1)),
            max = size(all_data{i},1);
        end
    end
    max = ceil(max/2)*2;
    count = 0;
    clear fft_array fft_toplot_array mean_fft_toplot_array;
    for i = 1:length(all_data),
        if (size(all_data{i},1) >= 10),
            count = count + 1;
            temp = all_data{i};
            vel_temp = sum(temp(2:end,11:12).^2,2)';
            fft_array(count, :) = zeros(1, max);
            fft_array(count, 1:length(vel_temp)) = detrend(vel_temp);
        end
    end
    for i = 1:size(fft_array,1),
        Y = fft(fft_array(i,:), max);
        Pyy = abs(Y)/max;
        fft_toplot_array(i,:) = Pyy;
    end
    mean_fft_toplot_array = mean(fft_toplot_array,1);
    
    
    % THIS ASSUMES YOUR CYCLE TIME IS 15s (so, a minute requires *4)  WRONG
    % ASSUMPTION NOW WRONG, TO BE FIXED IF THIS ANALYSIS BECOMES VALUABLE
    % (Francesco Marangoni 1/27/2011)
    x = (0:max/2)/max*4;
    plot(x,mean_fft_toplot_array(1:(max/2)+1), 'Color', colors{mod(cur_file-1, 7)+1}, 'LineWidth', 1.25);   
    xlabel('Frequency (cycles/min)');
    ylabel('Power');
    title('Average Power Spectrum of Velocity Change over Time');
    
    fclose(fid);
    
    [outfile, rem] = strtok(filelist{cur_file}, '.');
    outfile = strcat(outfile, '_fft.txt');
    fid = fopen(outfile, 'w');
    fprintf(fid, 'Freq (cyc/min)\tPower\n');
    for i = 1:length(x)
        fprintf(fid, '%5.2f\t%5.2f\n',x(i), mean_fft_toplot_array(i));
    end
    fclose(fid);
    
    waitbar(cur_file/num_files);
end
close(wh);

%distribution of 3DIV 
figure(distrib_fig);
x = 0:2:40;
set (gca, 'FontSize', 18);
bar(x, distrib_all', 1);
plot(x, distrib_all');
xlabel('Inst. 3D velocity (\mum/s)')
ylabel('Relative Frequency')
%title(sprintf('Distribution of 3D Instantaneous Velocities - %s', strrep(strtok(infile,'.'), '_', '\_')));
title('Distribution of 3D Instantaneous Velocities');
axis([-0.5 40.5 0 0.35])
%if num_files == 2,
%    p = ranksum(for_ttest{1}, for_ttest{2}, 0.05);
%    fprintf('Rank Sum Test for 3D Velocities -- p <= %d\n', p);
%end
legend(strrep(filelist,'_','\_'),'Location','SouthOutside');

figure(distrib_fig_two);
x = 0:2:40;
bar(x, distrib_2d', 1);
xlabel('Inst. 2D velocity (\mum/s)')
ylabel('Relative Frequency')
%title(sprintf('Distribution of 3D Instantaneous Velocities - %s', strrep(strtok(infile,'.'), '_', '\_')));
title('Distribution of 2D Instantaneous Velocities');
axis([-0.5 40.5 0 0.35])
%if num_files == 2,
%    p = ranksum(hist_data_2d{1}, hist_data_2d{2}, 0.05);
%    fprintf('Rank Sum Test for 2D Velocities -- p <= %d\n', p);
%end
legend(strrep(filelist,'_','\_'), 'Location','SouthOutside');

figure(scatter_fig);
legend(strrep(filelist,'_','\_'));

figure(fft_fig);
legend(strrep(filelist,'_','\_'));

figure(change_hist_fig);
x =-180:10:180;
bar(x, all_change', 1);
xlabel('Change in Cell Direction (\circ)');
ylabel('Relative Frequency');
%title(sprintf('Change in Cell Direction - %s', strrep(strtok(infile,'.'), '_', '\_')));
title('Change in Cell Direction');
axis([-180.5 180.5 0 0.1]);
legend(strrep(filelist,'_','\_'));

%figure(sample_fig);
%legend(filelist);


% %Migrating Cells

all_times = [];
for i = 1:length(all_data),
    temp = all_data{i};  % loads tracks one at a time
    all_times = [all_times; temp(2:end,22)]; % registers all timepoints modified Francesco 1/27/2011  to handle new column indexes
end
times = unique(all_times(find(~isnan(all_times))));  % gives back an array of timepoints without repetitions
pop_vel = cell(size(times));
motile = zeros(size(times));
nonmotile = zeros(size(times));
for curr_time = 1:length(times),  %one timepoint at a time
    for curr_cell = 1:length(all_data),  % one track at a time
        temp = all_data{curr_cell};
        temp_times = [NaN; temp(2:end,22)];   % modified Francesco 1/27/2011  to handle new column indexes - reads all timepoints for that track
        x = find(temp_times == times(curr_time));  % finds timepoints that are actually present in that track
        if ~isempty(x),
            if x == 2,
                cur_vel = temp(2, 19);  % modified Francesco 1/27/2011  to handle new column indexes
            else
                cur_vel = temp(x, 19);
            end
            if cur_vel > min_inst_vel,
                motile(curr_time) = motile(curr_time)+1;
            else
                nonmotile(curr_time) = nonmotile(curr_time)+1;
            end
            pop_vel{curr_time} = cat(1, pop_vel{curr_time}, cur_vel);
        end
        
    end
end
 
%change in cell motility over time
figure (cell_mot_chg_fig);
hold on;
plot(times, 100*(motile./(motile+nonmotile)),'b-')
title('Change in Population Cell Motility over Time');
xlabel('Time');
ylabel('Percent Motile Cells');
%create array to output
Motile_array = cat(2,times,(100*(motile./(motile+nonmotile))));
%output Motile array
%k = findstr(infile,'/');
%infile = infile(k(end)+1:end-12);
filename = strcat(infile,'motile_pop.txt');
dlmwrite(filename, Motile_array, '\t');
legend(strrep(filelist,'_','\_'));

% calculates and plots percentage of signaling cells over time
% FRANCESCO MARANGONI 1/27/2011
figure (sign_time_plot);
hold on;
sign_thresh = input('Insert the signaling index threshold above which the cell is signaling');
% pop_sign = cell(size(times));
signaling = zeros(size(times));
nonsignaling = zeros(size(times));
for curr_time = 1:length(times),  %one timepoint at a time
    for curr_cell = 1:length(all_data),  % one track at a time
        temp = all_data{curr_cell};
        temp_times = [NaN; temp(2:end,22)];   % modified Francesco 1/27/2011  to hanle new column indexes - reads all timepoints for that track
        x = find(temp_times == times(curr_time));  % finds timepoints that are actually present in that track
        if ~isempty(x),
            cur_sign = temp(x, 24);  
            macur_sign = movmean (cur_sign, (2*smooth_sig+1));
            if macur_sign > sign_thresh,
                signaling(curr_time) = signaling(curr_time)+1;
            else
                nonsignaling(curr_time) = nonsignaling(curr_time)+1;
            end
            % pop_sign{curr_time} = cat(1, pop_sign{curr_time}, cur_sign);
            % it would concatenate signaling values according to signaling
            % status
        end
        
    end
end
plot((times.*(cyc_time/60)), 100*(signaling./(signaling+nonsignaling)),'r-', 'Linewidth', 1.25)
title('Percentage of signaling cells over Time');
xlabel('Time');
ylabel('Percent of Signaling Cells');



%Construct and Plot Population inst. 2D velocity plots
pop_vel_mean = zeros(size(times));
pop_vel_error = zeros(size(times));
for curr_time = 1:length(times),
    pop_vel_mean(curr_time) = mean(pop_vel{curr_time});
    pop_vel_error(curr_time) = std(pop_vel{curr_time})/sqrt(length(pop_vel(curr_time)));
end
figure;
hold on;
plot(times, pop_vel_mean, 'b-');
plot(times, pop_vel_mean+pop_vel_error, 'r--');
plot(times, pop_vel_mean-pop_vel_error, 'r--');
title('Overall Population Mean +/- SEM Instantaneous 2D Velocity Over Time');
xlabel('Time');
ylabel('Instantaneous 2D Velocity');

%Plot the motile population's mean inst. 2D velocity
pop_vel_mean = zeros(size(times));
pop_vel_error = zeros(size(times));
for curr_time = 1:length(times),
    temp = pop_vel{curr_time};
    above_thresh_vel = temp(find(temp>=min_inst_vel));
    pop_vel_mean(curr_time) = mean(above_thresh_vel);
    pop_vel_error(curr_time) = std(above_thresh_vel)/sqrt(length(above_thresh_vel));
end
figure;
hold on;
plot(times, pop_vel_mean, 'b-');
plot(times, pop_vel_mean+pop_vel_error, 'r--');
plot(times, pop_vel_mean-pop_vel_error, 'r--');
title('Motile Population Mean +/- SEM Instantaneous 2D Velocity Over Time');
xlabel('Time');
ylabel('Instantaneous 2D Velocity');
%create array to output
TwoDIV_array = cat(2,times,pop_vel_mean,pop_vel_error);
%output 3DIV array
filename = strcat(infile,'2DIVmotpop.txt');
dlmwrite(filename, TwoDIV_array, '\t');


% Calculate and Plot 3D Inst Vel over Time

all_times = [];
for i = 1:length(all_data),
    temp = all_data{i};
    all_times = [all_times; temp(2:end,22)];    % modified Francesco 1/27/2011  to handle new column indexes
end
times = unique(all_times(find(~isnan(all_times))));
pop_vel = cell(size(times));
motile = zeros(size(times));
nonmotile = zeros(size(times));
for curr_time = 1:length(times),
    for curr_cell = 1:length(all_data),
        temp = all_data{curr_cell};
        temp_times = [NaN; temp(2:end,22)];  % modified Francesco 1/27/2011  to handle new column indexes
        x = find(temp_times == times(curr_time));
        if ~isempty(x),
            if x == 2,
                cur_vel = temp(2, 10);
            else
                cur_vel = temp(x, 10);
            end
            if cur_vel > min_inst_vel,
                motile(curr_time) = motile(curr_time)+1;
            else
                nonmotile(curr_time) = nonmotile(curr_time)+1;
            end
            pop_vel{curr_time} = cat(1, pop_vel{curr_time}, cur_vel);
        end
        
    end
end


%Construct and Plot Population inst. 3D velocity plots
pop_vel_mean = zeros(size(times));
pop_vel_error = zeros(size(times));
for curr_time = 1:length(times),
    pop_vel_mean(curr_time) = mean(pop_vel{curr_time});
    pop_vel_error(curr_time) = std(pop_vel{curr_time})/sqrt(length(pop_vel(curr_time)));
end
figure;
hold on;
plot(times, pop_vel_mean, 'b-');
plot(times, pop_vel_mean+pop_vel_error, 'r--');
plot(times, pop_vel_mean-pop_vel_error, 'r--');
title('Overall Population Mean +/- SEM Instantaneous 3D Velocity Over Time');
xlabel('Time');
ylabel('Instantaneous 3D Velocity');
%create array to output
ThreeDIV_array = cat(2,times,pop_vel_mean,pop_vel_error);
%output ThreeDIV array
filename = strcat(infile,'3DIVpop.txt');
dlmwrite(filename, ThreeDIV_array, '\t');


%Plot the motile population's mean inst. 3D velocity
pop_vel_mean = zeros(size(times));
pop_vel_error = zeros(size(times));
for curr_time = 1:length(times),
    temp = pop_vel{curr_time};
    above_thresh_vel = temp(find(temp>=min_inst_vel));
    pop_vel_mean(curr_time) = mean(above_thresh_vel);
    pop_vel_error(curr_time) = std(above_thresh_vel)/sqrt(length(above_thresh_vel));
end
figure;
hold on;
plot(times, pop_vel_mean, 'b-');
plot(times, pop_vel_mean+pop_vel_error, 'r--');
plot(times, pop_vel_mean-pop_vel_error, 'r--');
title('Motile Population Mean +/- SEM Instantaneous 3D Velocity Over Time');
xlabel('Time');
ylabel('Instantaneous 3D Velocity');
%create array to output
ThreeDIVmot_array = cat(2,times,pop_vel_mean,pop_vel_error);
%output Motile ThreeDIV array
filename = strcat(infile,'3DIVmotpop.txt');
dlmwrite(filename, ThreeDIVmot_array, '\t');

warning on all;