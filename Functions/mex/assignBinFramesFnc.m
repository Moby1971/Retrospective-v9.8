function [card_assignments,resp_assignments] = assignBinFramesFnc(bin_times_card,bin_times_resp,resp_window,nr_klines,nr_card_frames,nr_resp_frames)

% Assignmnent of all k-line acquisitions to a specific bin
%
% INPUT : bin_times_card = array of all time-stamps (units of samples) of the whole acquisition
%
%                       >    <--- decreases ------ bin_times_card(j)
%                      ||            |
%                      ||  j-th bin  |
%                      ||            |
%                      ||            |
%       - increases -----> i-th time point / sample
%
% OUTPUT: For each time point / k-line sample assignment to cardiac bin (saw-tooth pattern)
%
%         /|      /|    /|        /|
%        / |     / |   / |      /  |
%       /  |   /   |  /  |    /    |  /     = PHASE BINNING, each heart-beat has equal number of bins
%      /   |  /    | /   |  /      | /                       despite differences in heart beat duration
%     /    |/      |/    |/        |/
%     12345 1 23 45 12345 1 2 3 4 5  <- bins
%


startpoint = round(bin_times_resp(2)+1);                        % first measurement to be considered
endpoint = round(bin_times_resp(length(bin_times_resp))-1);     % last measurement to be considered
loc_maxj_card = length(bin_times_card);                         % last bin border
loc_maxj_resp = length(bin_times_resp);

card_assignments = zeros(nr_klines,1);                          % zero = no assignment (breathing, begin/end of data)
parfor i=startpoint:endpoint                                    % start search to which heartbeat the measurement belongs

    j=loc_maxj_card;
    while j>1 && i<bin_times_card(j)  %#ok<*PFBNS>
        j=j-1;
    end
   
    card_assignments(i) = mod(j-1,nr_card_frames)+1;            % assign to bin frame number = j modulus nr_frames
    
    if nr_resp_frames==1 && resp_window(i)==1                   % if measurement is during respiration and only 1 resp state, put back to 0 to discard this k-line
        card_assignments(i) = 0;
    end

end

resp_assignments = zeros(nr_klines,1);                      % zero = no assignment (breathing, begin/end of data)
parfor i=startpoint:endpoint                                % start search to which respiration the measurement belongs
    j=loc_maxj_resp;
    while j>1 && i<bin_times_resp(j)
        j=j-1;
    end
    resp_assignments(i) = mod(j-1,nr_resp_frames)+1;        % assign to bin frame number = j modulus nr_frames
end