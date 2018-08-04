function descps = extractNccFeature(img, Locs, halfsz)

% parse the input parameters
if(~exist('halfsz','var'))
   halfsz = [12,12];
else
    if(length(halfsz) <= 1)
        halfsz = [halfsz, halfsz];
    else
        halfsz = halfsz(1:2);
    end
end
halfsz = round(halfsz);
halfsz(halfsz<1) = 1;

%% ï¿½ï¿½È¡ï¿½ï¿½ï¿½ï¿½ï¿?

nc = size(img,3);
dim = prod(2*halfsz+1);
descps = zeros(size(Locs,1), nc * dim);  

img = double(img);
for i = 1 : size(Locs,1)
    x = Locs(i,1);
    y = Locs(i,2);
    
    xlo = max([1, x - halfsz(1)]);
    xhi = min([size(img,2), x + halfsz(1)]);
    ylo = max([1, y - halfsz(2)]);
    yhi = min([size(img,1), y + halfsz(2)]);
    j=1;
    for n = xlo : xhi
        for m = ylo : yhi
            if (j<nc * dim+1)
                descps(i,j) = descps(i,j)+ img(m,n);
                j=j+1;
            end
        end
    end

end

% do the normalization
descps = descps - repmat(mean(descps,2),[1 nc * dim]);
descps = descps ./ repmat(sqrt(sum(descps.^2,2)+1e-20),[1 nc*dim]);