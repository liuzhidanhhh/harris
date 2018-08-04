
function [locs,cornerness] = detHarrisCorners(img, corn_thresh, interactive)

%% parse the parameters
if(~exist('corn_thresh','var'))
    corn_thresh = 0.01;
end

if(~exist('interactive','var'))
    interactive = 0;
end

if(size(img,3)>1)
    img = rgb2gray(img);
end


%%
[nr,nc] = size(img);

%Filter for horizontal and vertical direction
dx = [1 0 -1];
dy = [1; 0; -1];


% Convolution of image with dx and dy
Ix = conv2(img, dx, 'same');
Iy = conv2(img, dy, 'same');

if(interactive)
    figure(1); 
    subplot(1,2,1); imagesc(Ix); axis equal image off; colormap('gray'); title('Ix');
    subplot(1,2,2); imagesc(Iy); axis equal image off; colormap('gray'); title('Iy');
    cdata = print('-RGBImage');
    imwrite(cdata, fullfile('corner', [name, '-grad.png']));
end

%% Raw Hessian Matrix
% Hessian Matrix Ixx, Iyy, Ixy
Ixx = Ix .* Ix;
Iyy = Iy .* Iy;
Ixy = Ix .* Iy;

if(interactive)
    figure(2); title('Before Smoothing');
    subplot(2,2,1); imagesc(Ixx); axis equal image off; colormap('gray'); title('Ixx');
    subplot(2,2,2); imagesc(Ixy); axis equal image off; colormap('gray'); title('Ixy');
    subplot(2,2,3); imagesc(Ixy); axis equal image off; colormap('gray'); title('Ixy');
    subplot(2,2,4); imagesc(Iyy); axis equal image off; colormap('gray'); title('Ixy');
    cdata = print('-RGBImage');
    imwrite(cdata, fullfile('corner', [name, '-rawHessian.png']));
end

%% Gaussian filter definition (https://en.wikipedia.org/wiki/Canny_edge_detector)
G = [2, 4, 5, 4, 2; 4, 9, 12, 9, 4;5, 12, 15, 12, 5;4, 9, 12, 9, 4;2, 4, 5, 4, 2];
G = 1/159.* G;

% Convolution with Gaussian filter
Ixx = conv2(Ixx, G, 'same');
Iyy = conv2(Iyy, G, 'same');
Ixy = conv2(Ixy, G, 'same');

if(interactive)
    figure(3); title('After Smoothing');
    subplot(2,2,1); imagesc(Ixx); axis equal image off; colormap('gray'); title('Ixx');
    subplot(2,2,2); imagesc(Ixy); axis equal image off; colormap('gray'); title('Ixy');
    subplot(2,2,3); imagesc(Ixy); axis equal image off; colormap('gray'); title('Ixy');
    subplot(2,2,4); imagesc(Iyy); axis equal image off; colormap('gray'); title('Ixy');
    cdata = print('-RGBImage');
    imwrite(cdata, fullfile('corner', [name, '-smoothHessian.png']));
end

%%
% calculate the lambda1 and lambda2
if(interactive)
    delta = (Ixx + Iyy).^2 - 4 * 1 * (Ixx.*Iyy - Ixy .* Ixy);
    lambda1 = 0.5 * ((Ixx + Iyy) + sqrt(delta));
    lambda2 = 0.5 * ((Ixx + Iyy) - sqrt(delta));
    figure(6); 
    subplot(1,2,1); imagesc(lambda1); axis equal image off; colormap('gray'); title('\lambda_1');
    subplot(1,2,2); imagesc(lambda2); axis equal image off; colormap('gray'); title('\lambda_2');
end

%%
% Calculate the corner responseling'luan
k = 0.04; % usually in the range[0.04 0.06]
R=zeros(nr,nc);
for h=1:nr  
    for w=1:nc  
        %计算M矩阵  
        M=[Ixx(h,w) Ixy(h,w);Ixy(h,w) Iyy(h,w)];  
        %计算R用于判断是否是边缘  
        R(h,w)=det(M) - k*(trace(M))^2;  
    end  
end  

rmax = max(R(:)); % get the maximum response of the cornerness
corner = zeros(nr,nc);
corner = (R>corn_thresh*rmax).*R;
%fun = @(x) max(x(:));
%Rlmax = nlfilter(R,[3 3], fun);
corner_peaks=imregionalmax(R);%区域最大值
if(interactive)
    figure(4); 
    subplot(1,2,1); imagesc(R); title('Before Non-Local Maximum Suppression');axis equal image off; colormap('gray');

    subplot(1,2,2); imagesc(R .* corner); title('After Non-Local Maximum Suppression');axis equal image off; colormap('gray');
    cdata = print('-RGBImage');
    imwrite(cdata, fullfile('corner', [name, '-cornerness.png']));
end

%[iLoc, jLoc] = find(R.*corner > corn_thresh * rmax);

%[iLoc,jLoc] = find(corner>0);
for h=1:nr  
    for w=1:nc  
        if(corner_peaks(h,w)==1)
            if(corner(h,w)<=0)
                corner_peaks(h,w)=0;
            end
        end
    end
end
[iLoc,jLoc] = find(corner_peaks==1);

if(interactive)
    figure(5); imagesc(img); hold on; axis equal image off; colormap('gray');
    for i = 1 : length(iLoc)
        plot(jLoc(i), iLoc(i), 'rx');
    end
    hold off;
    cdata = print('-RGBImage');size
    imwrite(cdata, fullfile('corner', [name, '-corners.png']));
end

% assign the output variables
locs = [jLoc iLoc];
cornerness = R(sub2ind(size(R),iLoc,jLoc));