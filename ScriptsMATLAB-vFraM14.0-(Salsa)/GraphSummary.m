function [varargout] = GraphSummary(varargin)

wh = waitbar(0, 'Graphing file(s)...');

num_files = length(varargin);
filelist = varargin;

%make plots
scatter_fig = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
contour_plot = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
surf_fig = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
distrib_fig = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
distrib_fig_two = figure; hold on;
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
set(gcf, 'InvertHardCopy', 'off');
colors = {'g', 'r', 'b', 'c', 'm', 'k', 'y'};
ccolors = [0 0 1; 1 1 0; 1 0 0; 0 1 1; 1 0 1; 0 0 0; 0 1 0];
colormap(ccolors);


all_data = [];

for cur_file = 1:num_files,
    clear infile;
    infile = filelist{cur_file};
    infile = strrep(infile, ' ', '');
    
    fid=fopen(infile);
    
    data = [];
    numlines = 1;
    
    %get first line, and ignore it
    tline = fgetl(fid);
    
    %get the actual data
    while ~feof(fid),
        tline = fgetl(fid);
        datacount = 1;
        lineofdata = {};
        
        while (length(tline) > 0),
            [token, tline] = strtok(tline, char(9));
            lineofdata{datacount} = token;
            datacount = datacount+1;
        end
        
        
        for i = 1:15,
            try
                data(numlines, i) = str2num(lineofdata{i});
            catch
                data(numlines, i) = nan;
            end
        end
        numlines = numlines + 1;
        
    end
    temp = data(find(~isnan(data(:,15))), 15);
    if (cur_file > 1),        
        if (length(temp) < size(all_data, 1)),
            temp(length(temp)+1:size(all_data, 1)) = nan;
        elseif (length(temp) > size(all_data, 1)),
            all_data(size(all_data, 1):length(temp), :) = nan;
        end
        all_data(:, cur_file) = temp;
    else
        all_data = temp;
    end
    
    vel_array{cur_file} = data(:,6:10);
    
    % 3D velocity histogram creation (not plotting)
    x = 0:2:40;
    n = hist(data(:,6),x);
    distrib_all(cur_file, :) = n/length(data(:,6));
    
    % 2D vel histogram creation (not plotting)
    hx = 0:2:40;
    hist_data_2d{cur_file} = sqrt(data(:,7).^2 + data(:,8).^2);
    n = hist(sqrt(data(:,7).^2 + (data(:,8)).^2), hx);
    distrib_2d(cur_file, :) = n/length(data(:,8));
    
    % Plot Scatterplot average velocity x versus average velocity y
    figure(scatter_fig);
    plot(data(:,7),data(:,8), '.', 'Color', colors{mod(cur_file-1,7)+1});
    xlabel('Average velocity in x (\mum/s)');
    ylabel('Average velocity in y (\mum/s)');
    %title(sprintf('Average velocity in x vs. Average velocity in y - %s', strrep(strtok(infile,'.'), '_', '\_')));
    title('Average velocity in x vs. Average velocity in y');
    axis([-40 40 -40 40]);
    grid on;
    
    % Average velocity contour plot
    xplot = -40:2:40;
    yplot = -40:2:40;
    for g = 1:length(xplot)-1,
        for h = 1:length(yplot)-1,
            tempx = logical(data(:,7) >= xplot(g)) & logical(data(:,7) <= xplot(g+1));
            tempy = logical(data(:,8) >= yplot(h)) & logical(data(:,8) <= yplot(h+1));
            zz(g,h) = sum(ismember(find(tempx), find(tempy)));
        end
    end
    figure(contour_plot);
    [tx, ty] = meshgrid(-39:2:39, -39:2:39);
    contour3(tx, ty, zz, 20, colors{mod(cur_file-1, 7)+1});
    xlabel('Average velocity in x (\mum/s)');
    ylabel('Average velocity in y (\mum/s)');
    %title(sprintf('Average velocity in x vs. Average velocity in y - %s', strrep(strtok(infile,'.'), '_', '\_')));
    title('Average velocity in x vs. Average velocity in y');
    axis([-40 40 -40 40]);
    grid on;
    
    figure(surf_fig);
    surfc(tx, ty, zz, 'FaceColor', colors{mod(cur_file-1,7)+1});
    xlabel('Average velocity in x (\mum/s)');
    ylabel('Average velocity in y (\mum/s)');
    %title(sprintf('Average velocity in x vs. Average velocity in y - %s', strrep(strtok(infile,'.'), '_', '\_')));
    title('Average velocity in x vs. Average velocity in y');
    axis([-40 40 -40 40]);
    
    waitbar(cur_file/num_files);
    fclose(fid);
end

% Chemotactic Index box plot
% figh = figure;
% figure(figh);
% boxplot(all_data, 0, '.', 1, 2);
% xlabel('Chemotactic Index');
% set(gca, 'XTickLabel', filelist);
% ylabel('');
% title('Chemotactic Index');

% Plot 3D histogram
figure(distrib_fig);
x = 0:2:40;
bar(x, distrib_all', 1);
colormap(ccolors);
xlabel('Average 3D velocity (\mum/s)')
ylabel('Relative Frequency')
%title(sprintf('Distribution of Average 3D Velocities - %s', strrep(strtok(infile,'.'), '_', '\_')));
title('Distribution of Average 3D Velocities');
axis([-0.5 40.5 0 0.35])
%if num_files == 2,
 %   p = ranksum(for_ttest{1}, for_ttest{2}, 0.05);
 %  fprintf('Rank Sum Test for 3D Velocities -- p <= %d\n', p);
 %end
legend(strrep(filelist,'_','\_'));

% Plot 2D histogram
figure(distrib_fig_two);
x = 0:2:40;
bar(x, distrib_2d', 1);
colormap(ccolors);
xlabel('Average 2D velocity (\mum/s)')
ylabel('Relative Frequency')
%title(sprintf('Distribution of 3D Average Velocities - %s', strrep(strtok(infile,'.'), '_', '\_')));
title('Distribution of Average 2D Velocities');
axis([-0.5 40.5 0 0.35])
% if num_files == 2,
%     p = ranksum(hist_data_2d{1}, hist_data_2d{2}, 0.05);
%     fprintf('Rank Sum Test for 2D Velocities -- p <= %d\n', p);
% end
legend(strrep(filelist,'_','\_'));

% close wait bar
close(wh);
