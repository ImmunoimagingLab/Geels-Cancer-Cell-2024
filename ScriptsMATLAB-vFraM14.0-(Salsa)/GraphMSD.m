% MODIFIED TO CONSIDER DISTANT JOINING FRANCESCO MARANGONI 2/5/2011
% would be better to pass cyc_time as variable Francesco Marangoni 2/5/2011
function [varargout] = GraphMSD(varargin)


cyc_time = 15;  % MODIFY HERE WHEN INTRODUCING CYC_TIME AS VARIABLE!!!!!
num_files = length(varargin);
filelist = varargin;

figh = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
wh = waitbar(0, 'Graphing file(s)...');
colors = {'k', 'r', 'g', 'b', 'c', 'm', 'y'};
%colors = {[0.0 0.502 0.251], 'r', 'b', 'c', 'm', 'k', 'y'};
figsingle = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');

for cur_file = 1:num_files,  % considers all open files one at a time
    clear infile;
    infile = filelist{cur_file};
    infile = strrep(infile, ' ', '');
    
    % initialize variables
    fid = fopen(infile);
    
    % init variables
    all_data = {};
    num_tracks = 1;
    numlines = 1;
    data = [];
    lengthmax = 0;
    
    % ignore first line
    tline = fgetl(fid);
    
    % read data out of file
    while ~feof(fid),
        % get the whole line
        tline = fgetl(fid);
        
        if (strcmp(upper(tline(1:3)), '---')),
            % so, we are at the end of a cell's track
            if length(data)>lengthmax,  %calculates max timepoint for tracks Francesco Marangoni 2/5/2011
                lengthmax = length(data);
            end
            all_data{num_tracks} = data;
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
                [token, tline] = strtok(tline, char(9));
                lineofdata{datacount} = token;
                datacount = datacount+1;
            end
            
            % transfer the chunks into numbers
            for i = 1:4,
                try
                    data(numlines, i) = str2num(lineofdata{i});
                catch
                    data(numlines, i) = nan;
                end
            end
            numlines = numlines + 1;
        end
    end
    
    %NEED TO GENERATE MEAN FOR MSD for all cells  MODIFIED FRANCESCO
    %MARANGONI 2/5/2011 to consider distant joining
    stop = 0;
    step = 1;
    avg_step_data = zeros(lengthmax, 2);
    step_data = zeros(lengthmax, length(all_data));  % added Francesco Marangoni 2/5/2011
    while ~stop,
        stop = 1;
        for i = 1:length(all_data),
            cur_data=all_data{i};  %  reads all tracks separately
            if size(cur_data,1) >= step,
                stop = 0;
                %step_data = [step_data cur_data(step, 3)];  % original
                %line
                act_step = cur_data(step,2)/cyc_time;  % Francesco Marangoni 2/5/2011
                step_data(act_step, i) = cur_data(step,3); % Francesco Marangoni 2/5/2011
            end
        end
        %cur_data
        %pause;
        warning off;
        %avg_step_data(step, 1) = mean(step_data);  % original line
        %avg_step_data(step, 2) = std(step_data)/sqrt(length(step_data));% original line
        warning on;
        step = step + 1;
    end
    max_step = step-1;
    for i = 1:lengthmax,  % computes mean and SEM of displacements at all timepoints  Francesco Marangoni 2/5/2011
        [r,c,v] = find (step_data(i,:));
        avg_step_data(i,1) = mean(v);
        avg_step_data(i,2) = std(v)/sqrt(length(v));
    end
    %avg_step_data(:,:)
    
    
    % Population MSD x/y plot with error bars of sem of mean displacement (D(t)) versus sqrt of time
    figure(figh);
    % x = sqrt([1:max_step]*cur_data(1,2)/60)'; % need the sqrt of time in minutes  %original line
    x = sqrt([1:lengthmax]*cyc_time/60)';  % proposed change Francesco Marangoni 2/5/2011  
    %errorbar(x, avg_step_data(:,1), avg_step_data(:,2));
    size(x);
    size(avg_step_data);
    errorbar([0; x], [0; avg_step_data(:,1)], [0; avg_step_data(:,2)], [0; avg_step_data(:,2)], colors{mod(cur_file-1, 7)+1}, 'LineWidth', 2);
    xlabel('Square root time (min^1^/^2)', 'FontSize', 18);
    ylabel('Mean displacement (\mum)', 'FontSize', 18);
    title('Mean Displacement versus Time', 'FontSize', 18);
    grid on;
    set (gca, 'FontSize', 18, 'LineWidth', 2);
    axis ([0 3 0 60]);
    axis square

    
    %MSD plot for each single cell
    figure(figsingle);
    xlabel('Square root time (min^1^/^2)', 'FontSize', 18);
    ylabel('Mean displacement (\mum)', 'FontSize', 18);
    title('Single Cell Mean Displacement versus Time', 'FontSize', 18);
    axis square
    grid off;
    set (gca, 'FontSize', 18, 'LineWidth', 2);
    axis ([0 3 0 50]);
    for i = 1:length(all_data),
        cur_data=all_data{i};
        plot([0; sqrt(cur_data(:,2)/60)], [0; cur_data(:,3)], colors{mod(cur_file-1, 7)+1}, 'LineWidth', 1);
    end

    
    
    fclose(fid);
    
    [outfile, rem] = strtok(filelist{cur_file}, '.');
    outfile = strcat(outfile, 'plot.txt');
    fid = fopen(outfile, 'w');
    fprintf(fid, 'Sqrt Time\tMean Disp.\tSEM\n');
    fprintf(fid, '%5.2f\t%5.2f\t%5.2f\n',0,0,0);
    for i = 1:length(x)
        fprintf(fid, '%5.2f\t%5.2f\t%5.2f\n',x(i),avg_step_data(i,1), avg_step_data(i,2));
    end
    fclose(fid);
    
    waitbar(cur_file/num_files, wh);

end

close(wh);
legend(strrep(filelist,'_','\_'));

figure(figh);
[legend_h,object_h,plot_h,text_strings] = legend(strrep(filelist,'_','\_'), 'Location', 'SouthOutside');
for i = 1:length(filelist),
    set(object_h(i), 'FontSize', 12);
end

figure(figsingle);
[legend_h,object_h,plot_h,text_strings] = legend(strrep(filelist,'_','\_'), 'Location', 'SouthOutside');
for i = 1:length(filelist),
    set(object_h(i), 'FontSize', 12);
end
for i = 1:length(plot_h),
    set(plot_h(i), 'Color', colors{mod(i-1, 7)+1});
end