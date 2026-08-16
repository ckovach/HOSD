function axes_handles = nimagesc(X,varargin)


title_font_size = 8;
stonset = [];
tickson = [];
clim = [min(min(min(X))) max(max(max(X)))];
fig = [];
varg = varargin;
i = 1;
flip = 0;
chlabels = [];
grid_size = [];
theta = 0;
title_height = 1;
noscale = 0;
while i <= length(varargin)
   switch lower(varg{i})
       case 'stonset'
           stonset = varargin{i+1};
           i = i+1;
       case 'tickson'
           tickson = varargin{i+1};
           i = i+1;
       case {'chlabels','labels'}
           chlabels = varargin{i+1};
           i = i+1;
       case 'caxis'
           clim = varargin{i+1};
           i = i+1;
       case 'figure'
           fig = varargin{i+1};
           i = i+1;
       case 'size'
           grid_size = varargin{i+1};
           i = i+1;
       case 'flip'
           flip = 1;
       case 'noscale'
           noscale = 1;
       case 'rotate'
           theta = varargin{i+1};
           i = i+1;           
       otherwise
           error([varargin{i},' is not a valid option.']);
   end         
   i = i+1;
end

if isempty(fig)
    fig = figure;
else
    figure(fig)
end

n = size(X,3);

if isempty(grid_size)
  grid_size = [ ceil(sqrt(n)), ceil(n/ceil(sqrt(n)))];
end



I = grid_size(1);
J = grid_size(2);

 index = reshape([1:prod(grid_size)]', grid_size(2), grid_size(1))';

 if flip == 1
    index = fliplr(index);
    index = reshape(index',prod(grid_size),1)';
    X = X(:,:,index);
    if ~isempty(chlabels)
        chlabels = chlabels(index);
    end        
end

if ~iscell(chlabels)
    chlabels = num2cell(chlabels);
end

plotsize = I*J;
axes_handles = [];
axes_pos = [];
for i = 1:n
    
    axes_handles(end+1) = subplot(I,J,i);
    axes_pos(end+1,:) = get(gca,'position');
    
    if noscale
         im = image(X(:,:,i));
     else
         im = imagesc(X(:,:,i));
     end
     
    if ~isempty(stonset)
      hold on
      for vl = 1:length(stonset)
        line([stonset(vl) stonset(vl)],[min(axis) max(axis)]);
      end    
    end
    
    if ~isempty(chlabels)
        
      t = title(num2str(chlabels{i}),'fontsize',title_font_size);
        
      set(t,'units','normalized')
      tpos = get(t,'position');
      set(t,'position',[tpos(1),title_height,tpos(3)])
    end
    
    if ~ismember(i, tickson)
        axis off
    end
            
    caxis(clim);
end
    

if theta ~= 0           %rotates subplots around the center
   thrad = pi*theta/180;
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

   ROT = [cos(thrad) sin(thrad); -sin(thrad) cos(thrad)];                               %Rotation matrix

   cent = mean(axes_pos(:,1:2));                                %Center of rotation is center of axes locations. Change to [.5 .5] to make it the center of the figure.
   
   tot_width = max(axes_pos(:,1) + axes_pos(:,3)) -min(axes_pos(:,1));
   tot_height = max(axes_pos(:,2) + axes_pos(:,4)) -min(axes_pos(:,2));
  
   
   rot_adj = [abs(cos(thrad)), abs(sin(thrad))];                                            % adjustment so that axes lie within the original box
   pos_adjustment = [(rot_adj*[tot_width, tot_height]')/tot_height,  (rot_adj*[tot_height, tot_width]')/tot_width];  %
   
   axes_pos_cent = axes_pos(:,1:2) - cent(ones(1,length(axes_pos)),:);
   axes_pos_cent_adj = axes_pos_cent./pos_adjustment(ones(length(axes_pos_cent),1),:);

   axes_pos_rotated = (ROT*axes_pos_cent_adj(:,1:2)')' + cent(ones(length(axes_pos),1),:);
   nusize = axes_pos(:,3:4)./pos_adjustment(ones(length(axes_pos),1),:);
   new_axes = [axes_pos_rotated, nusize];
   
   if mod(theta,90) == 0
       figure
       subplot(grid_size(2),grid_size(1),prod(grid_size))
       test_pos = get(gca,'position');
       close(gcf)
       new_axes(:,3:4) = ones(size(new_axes,1),1)*test_pos(3:4);
   end
   
   for i = 1:length(axes_handles)
      set(axes_handles(i),'position',new_axes(i,:))
      axes_pos(i,:) = get(axes_handles(i),'position');    
  end
end
