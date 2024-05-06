function output = my_nanmean(input)

output = mean(input(find(~isnan(input))));
