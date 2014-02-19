function count = hist3d(data, xrange, yrange, zrange)

for i=1:length(xrange)-1
    data((data(:,1)>xrange(i)) & (data(:,1)<=xrange(i+1)), 4) = i;
end

for i=1:length(yrange)-1
    data((data(:,2)>yrange(i)) & (data(:,2)<=yrange(i+1)), 5) = i;
end

for i=1:length(yrange)-1
    data((data(:,3)>zrange(i)) & (data(:,3)<=zrange(i+1)), 6) = i;
end

count = zeros(length(xrange)-1, length(yrange)-1, length(zrange)-1);

data = data(data(:,4)>0,:); % if a data point is out of the x range, throw it away
data = data(data(:,5)>0,:); % if a data point is out of the y range, throw it away
data = data(data(:,6)>0,:); % if a data point is out of the z range, throw it away

for i=1:size(data,1)
    count(data(i,4), data(i,5), data(i,6)) = count(data(i,4), data(i,5), data(i,6)) + 1;
end