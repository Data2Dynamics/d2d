%A Simple function that makes a 2D histgram developed by Sisi Ma (sisima[at]rci.rutgers.edu )

% Input:    data: two cols, x value; y value
%             xrange: range and bins for x value (edges)
%             yrange: range and bins for y value (edges)
%Output: Count of a specifice (x,y) bin combination; 
%       Suggested visualizing tool: I like to use imagesc; bar3 will work fine
%       too; have to change axis label though


function count=hist2d(data, xrange, yrange)



    for i=1:length(xrange)-1

        data((data(:,1)>xrange(i))&(data(:,1)<=xrange(i+1)),3)=i;
    end

    for i=1:length(yrange)-1

        data((data(:,2)>yrange(i))&(data(:,2)<=yrange(i+1)),4)=i;  

    end

    count=zeros(length(xrange)-1,length(yrange)-1);

    data=data(data(:,3)>0,:); % if a data point is out of the x range, throw it away
    data=data(data(:,4)>0,:);% if a data point is out of the y range, throw it away
    
    for i=1:size(data,1)
        count(data(i,3),data(i,4))=count(data(i,3),data(i,4))+1; 
    end



end

    
    