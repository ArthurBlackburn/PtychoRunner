function [res_out, res] = particulate_FRC_res(r0, r1, params, nbins, res)
% Method used to determine the FRC resolution in manuscript linked to at 
% github.com/ArthurBlackburn/PtychoRunner.
%
% The method is based upon taking FRC measures from subregions and assessing 
% whether the correlating information in the patch is below a threshold, where 
% the threshold comes from Otsu's method appled to the mean FRC score.
% The distribution appears typically bi-modal: loosely some pathces are
% dominated by amorphous carbon, which is unlikely to form correlating
% features given the thickness of aC used in our work and the dynamic
% nature of aC at room temp, under high electron beam flux. Au particles in
% comparison are more robust against e-beam radiation and thin and so produce
% repeatable reconstructed information.
%
% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2023
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
% 
% *******************************************************************

% gets the FRC measures from all the subregions specified in params:
if nargin < 5
    res = CalcFRCMeas({r0, r1}, params{:}); 
end

% Get the mean FRC (Frc Mean = 'fm') for each sub-region   
res_shape = size(res);
fm = zeros(res_shape);

for ii = 1:res_shape(1)
    for jj = 1:res_shape(2)
        fm(ii, jj) = res(ii, jj).stat.fsc_mean;
    end 
end
% struct2cell then cell2mat can also do the above with some persuasion, 
% but these functions are more cumbersome  and clumsy than the above loop
% imo... e.g. would need to do something like the below, 
% fm = reshape(cell2mat({res.stat.fsc_mean}), size(res));
% but rework it as the above has error:
% "intermediate dot '.' indexing  produces a comma-separated list with n values,
% but it must produce a single value when followed by subsequent indexing operations".

% Analzyze the mean FSC/FRC scores for the subregions:
[counts, bin_edges] = histcounts(fm(:), 10); % 10 bins seems adequate.

% Threshold that minimizes the intra-class variance
thresh_norm = otsuthresh(counts);
% Convert relative position among bins to value of variable:
bin_centres = (bin_edges(1:(end-1)) + bin_edges(2:end))/2;
otsu_thresh = bin_centres(1) + thresh_norm*(bin_centres(end)-bin_centres(1));

% Exclude is 1 if we should exclude the data from consideration based on
% this threshold
exclude = fm < otsu_thresh;

% Now we need to do some kind of 'avaraging' of all the characteristics.
% Just summing the values in each frequency bin seems a fair way to go. We
% can also re-bin... 

% How many bins to use for the output summmed characteristic?
if nargin < 4 || isempty(nbins)
    num_freqs = zeros(1, numel(res));
    for jj = 1:numel(res)
        num_freqs(jj) = length(res(jj).freq);
    end
    nbins = min(num_freqs);
end

% Do the rebinning:
[new_freq_bins, rebinned_mean_FSC] = rebin_and_average(res, exclude, nbins, 'FSC');
[~, rebinned_mean_T] = rebin_and_average(res, exclude, nbins, 'T');

% Convert from normalized Nyquist frequency to spatial frequency.    
SF_per_A_from_Nfreq = @(freq) (freq./freq(end)) / (2*r0.recon_dat.recon_pixel_size*1e10);

res_out.spatial_freq = SF_per_A_from_Nfreq(new_freq_bins);
res_out.FRC = rebinned_mean_FSC;
res_out.thresh = rebinned_mean_T;

end

function [new_bin_centres, mean_pnts] = rebin_and_average(res, exclude, nbins, fname)
    
    % Rebin and average the values in each bin.
    % (Note: if all subregions are the same size, then all the below could be
    % made a whole lots simpler. The below covers the most general case which
    % we might investigate, where subregions are different sizes).

    % The new bin edges:
    new_freq_bin_edges = linspace(0,1,nbins+1);
    new_bin_centres = (new_freq_bin_edges(1:(end-1)) + new_freq_bin_edges(2:end))/2;
    
    % Setup comparison matrices:
    num_freqs = zeros(1, numel(res));
    for jj = 1:numel(res)
        num_freqs(jj) = length(res(jj).freq);
    end
    max_num_freqs = max(num_freqs);
    % Upper limits
    UL = new_freq_bin_edges(2:end);
    UL_mat = repmat(UL',[1 max_num_freqs]); 
        
    % Lower limits 
    LL = [new_freq_bin_edges(1:(end-1))]; 
    LL_mat = repmat(LL',[1 max_num_freqs]);
        
    FSC = zeros(numel(res), 0); % Second dimension will grow as appropriate in the loop
    bin_map = zeros(numel(res), 0);
    FRC_map = zeros(numel(res), 0);
        
    for ii = 1:numel(res)
	    % Each row contains the FSC for the result number.
        FSC(ii,1:(length(res(ii).FSC))) = res(ii).(fname);
        % The (normalized) frequencies for this result:
        normed_freq = (res(ii).freq)./max(res(ii).freq);
        % Sort and index into bins:
	    % put this into a matrix where the each row contains a repeat of the frequency points:
        compare_mat = repmat(normed_freq, [(nbins) 1]);
	    
	    % Which new frquency bin is the res.freq within?
	    % Frequencies above the lower edge of bin:
        above_lower = LL_mat(:, size(compare_mat,2)) <= compare_mat;
	    % Frequencies below the upper edge of bin?
        below_upper = compare_mat < UL_mat(:, size(compare_mat,2));
	    % Thus the freq is in this bin:
        in_bounds = above_lower & below_upper;
	    % convert to binary maps of indices:	
        inds = find(in_bounds);
        [a, b] = ind2sub(size(compare_mat),inds);
	    % a are the rows: row number gives the bin number;
	    % b are cols: col gives the index into the frequency output from the FSC program.
        bin_map(ii,1:length(a)) = a'; % bin numbers
        FRC_map(ii,1:length(b)) = b'; % related freqency point.
        % So finally we have a mapping of res.freq points (a) to new frequency bin numbers (b)
    end

    % Now let's take the mean number in each bin:
    mean_pnts = zeros(1,nbins);    
    
    for ki = 1:nbins
        scores_in_this_bin = [];
        for ii = 1:numel(res)
            if ~exclude(ii)
                scores_in_this_bin = [scores_in_this_bin, ...
                    FSC(ii, FRC_map(ii, (bin_map(ii,:) == ki) ) )]; %#ok<AGROW>                    
            end
        end
        mean_pnts(ki) = mean(scores_in_this_bin);
    end
        
end