function [first,excluded,other,lat]=first_after(base_times,ref_times,timewindow)
%first = first base_times after ref_times (e.g. first mag after pellet (mag+))
%excluded = base_times within timewindow after each first (mags immediately after mag+, "exclude" because likely still a mag+)
%other = other base_times (mag-)
%lat = latency to first (for each ref_event)

excluded = [];
first=zeros(length(ref_times),1);
lat=zeros(length(ref_times),1);

for i = 1:length(ref_times)
    if ~isempty(min(base_times(ref_times(i) < base_times)))
        first(i) = min(base_times(ref_times(i) < base_times));
        lat(i) = first(i)-ref_times(i);
    else
        lat(i) = NaN;
    end
end
    
first = unique(first); %remove doubles
first = first(any(first,2),:);  %remove 0 rows

allothers = setdiff(base_times,first);

for f = 1:length(first)
    if any((allothers>first(f)) & (allothers-timewindow-first(f) < 0))
        exclude = allothers((allothers>first(f)) & (allothers-timewindow-first(f) < 0)).';
        excluded = [excluded; exclude'];
    end
end

excluded = unique(excluded);
other = setdiff(allothers, excluded);

end