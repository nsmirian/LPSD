function [image_roi,status, is_beam]=get_ROI(image, threshold_factor,bits, disp_choice)

if nargin < 4
    disp_choice = 0;
end

if nargin<3
    bits = 12;	
end

if nargin<2
    threshold_factor=2;
end

status = 'ok';	
beam_i_threshold = 0.122*2^(bits); % if bits=12,beam_i_threshold=500;

flag1 = 0;
flag2 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the mathematical definition of the matrix M:
% filter_size = 1.5;
% filter = fspecial('gaussian',round(filter_size*4+1),filter_size);

% M=zeros(7,7);
% for x=-3:3
%     for y=-3:3
%         M(x+4,y+4)=1/sqrt(2*pi*1.5^2)*exp(-(x^2+y^2)/2/1.5^2);
%     end
% end
% M=M/sum(sum(M));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 M = [0.0013    0.0041    0.0079    0.0099    0.0079    0.0041    0.0013
      0.0041    0.0124    0.0241    0.0301    0.0241    0.0124    0.0041
      0.0079    0.0241    0.0470    0.0587    0.0470    0.0241    0.0079
      0.0099    0.0301    0.0587    0.0733    0.0587    0.0301    0.0099
      0.0079    0.0241    0.0470    0.0587    0.0470    0.0241    0.0079
      0.0041    0.0124    0.0241    0.0301    0.0241    0.0124    0.0041
      0.0013    0.0041    0.0079    0.0099    0.0079    0.0041    0.0013];

image=double(image);

[mean_noise1,~,status_noise1] = estimate_noise(image);

if strcmp(status_noise1,'warning') 
    if disp_choice
    disp('Warning: Noise (1) estimation may failed!');
    end
    status = 'warning';
    flag1 = 1;
end

image = image - mean_noise1;

image_filt = conv2(image, M, 'same');

[mean_noise2,std_noise2,status_noise2] = estimate_noise(image_filt);

if strcmp(status_noise2,'warning')
    if disp_choice
        disp('Warning: Noise (2) estimation may failed!');
    end
    status = 'warning';
    flag2 = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the functionality of the following part of codes is still under test. It
% can be simply skipped.
if flag1*flag2
	length_image_filt = numel(image_filt);
	[y_data,x_data] = hist(image_filt(1:length_image_filt),length([min(image_filt(1:length_image_filt)):max(image_filt(1:length_image_filt))]));

	[par1, yFit1] = util_gaussFit(x_data,median_filter(y_data),0,1);
	
	k=0;
    fit_start = [];
    while isempty(fit_start);
        fit_start = find(y_data >= par1(1)*(1000-k)/1000, 1, 'last' ) + 2;
        k = k+1;
    end
	[par2, yFit2] = util_gaussFit(x_data(fit_start:end),median_filter(y_data(fit_start:end)),0,1);
	mean_noise2 = par2(2);
	std_noise2 = (par2(3)*par2(4)+par2(3));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

threshold = mean_noise2 + threshold_factor*std_noise2;

num_pix_x = size(image_filt,1); 
num_pix_y = size(image_filt,2); 

candidate_list = find(image_filt(1:(num_pix_x*num_pix_y))>=threshold);
num_candidates_start = length(candidate_list);

if num_candidates_start==0
	disp('Warning: ROI cannot be found!');
	status = 'warning';
end

roi = zeros(num_pix_x,num_pix_y);
roi(candidate_list) = 1;


image_filt_start=conv2(image_filt, ones(5,5)/5^2, 'same');
[~,start_roi_x] = max(max(image_filt_start));
[~,start_roi_y] = max(max(image_filt_start'));
roi=bwselect(logical(roi),start_roi_x,start_roi_y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% codes to be tested.

% conn_comp = bwconncomp(roi);
% conn_comp_pixel_index = conn_comp.PixelIdxList;
% for jj = 1:length(conn_comp_pixel_index)
%     area_size(jj) = length(conn_comp_pixel_index{1,jj});
% end
% [~,largest_area_idx] = max(area_size);
% largest_area = conn_comp_pixel_index{1, largest_area_idx};
% roi(~ismember(1:(num_pix_x*num_pix_y),largest_area))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_roi = roi.*image;

roi_thres_noise = image_roi<0;
image_roi(roi_thres_noise)=0;

% simple ideas for determing is_beam
temp = image_filt.*roi;
int_max = max(temp(temp > 0));
int_min = min(temp(temp > 0));
int_diff = int_max - int_min;
if int_diff > beam_i_threshold
    is_beam = 1;
else
    is_beam = 0;
end


return;









