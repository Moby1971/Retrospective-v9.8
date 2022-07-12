function [sorted_kspace,sorted_averages] = sort2DdataFnc(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics,nr_reps,nrKsteps,traj,unsorted_kspace,card_bin_ass,includeWindow,resp_bin_ass,dyn_bin_ass)

sorted_kspace = complex(zeros(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics));
sorted_averages = zeros(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics);  

cnt = 0;

for slice=1:dimz                    % loop over slices
    for i=1:nr_reps                 % loop through all repetitions
        for j=1:nrKsteps            % loop through all the phase-encoding steps
            cnt = cnt + 1;
            if (card_bin_ass(cnt) > 0) && (includeWindow(cnt) == 1)         % if assigment = 0, this acquisition is discarded
                kline = traj(mod(cnt - 1,nrKsteps) + 1);             % the phase-encoding step using the trajectory info
                sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) = sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) + unsorted_kspace(1,i,slice,j,:,1);   % add the data to the correct k-position
                sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) = sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) + 1;        % increase the number of averages with 1
            end
        end
    end
end