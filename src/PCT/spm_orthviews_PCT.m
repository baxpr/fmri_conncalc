function varargout = spm_orthviews(action,varargin)
% Display Orthogonal Views of a Normalized Image
% FORMAT H = spm_orthviews('Image',filename[,position])
% filename - name of image to display
% area     - position of image
%            -  area(1) - position x
%            -  area(2) - position y
%            -  area(3) - size x
%            -  area(4) - size y
% H        - handle for ortho sections
% FORMAT spm_orthviews('BB',bb)
% bb       - bounding box
%            [loX loY loZ
%             hiX hiY hiZ]
%
% FORMAT spm_orthviews('Redraw')
% Redraws the images
%
% FORMAT spm_orthviews('Reposition',centre)
% centre   - X, Y & Z coordinates of centre voxel
%
% FORMAT spm_orthviews('Space'[,handle])
% handle   - the view to define the space by
% with no arguments - puts things into mm space
%
% FORMAT spm_orthviews('MaxBB')
% sets the bounding box big enough display the whole of all images
%
% FORMAT spm_orthviews('Resolution',res)
% res      - resolution (mm)
%
% FORMAT spm_orthviews('Delete', handle)
% handle   - image number to delete
%
% FORMAT spm_orthviews('Reset')
% clears the orthogonal views
%
% FORMAT spm_orthviews('Pos')
% returns the co-ordinate of the crosshairs in millimetres in the
% standard space.
%
% FORMAT spm_orthviews('Pos', i)
% returns the voxel co-ordinate of the crosshairs in the image in the
% ith orthogonal section.
%
% FORMAT spm_orthviews('Xhairs','off') OR spm_orthviews('Xhairs')
% disables the cross-hairs on the display.
%
% FORMAT spm_orthviews('Xhairs','on')
% enables the cross-hairs.
%
% FORMAT spm_orthviews('Interp',hld)
% sets the hold value to hld (see spm_slice_vol).
%
% FORMAT spm_orthviews('AddBlobs',handle,XYZ,Z,mat)
% Adds blobs from a pointlist to the image specified by the handle(s).
% handle   - image number to add blobs to
% XYZ      - blob voxel locations (currently in millimeters)
% Z        - blob voxel intensities
% mat      - matrix from millimeters to voxels of blob.
% This method only adds one set of blobs, and displays them using a
% split colour table.
%
% spm_orthviews('AddColouredBlobs',handle,XYZ,Z,mat,colour)
% Adds blobs from a pointlist to the image specified by the handle(s).
% handle   - image number to add blobs to
% XYZ      - blob voxel locations (currently in millimeters)
% Z        - blob voxel intensities
% mat      - matrix from millimeters to voxels of blob.
% colour   - the 3 vector containing the colour that the blobs should be
% Several sets of blobs can be added in this way, and it uses full colour.
% Although it may not be particularly attractive on the screen, the colour
% blobs print well.
%
% FORMAT spm_orthviews('window',hld,[w1 w2])
% Set the intensity window for specified handles
%
% FORMAT spm_orthviews('autowindow',hld)
% Automagically set the window
%
% FORMAT spm_orthviews('addcolorbar',hld)
% Add color bar to specified handles
%
% FORMAT spm_orthviews('colorbar2',hld,x)
% Supply a scale factor x to define right hand axis of colorbar.
%
% FORMAT spm_orthviews('Register',hReg)
% See spm_XYZreg for more information.
%
% FORMAT spm_orthviews('RemoveBlobs',handle)
% Removes all blobs from the image specified by the handle(s).
%
% FORMAT spm_orthviews('Register',hReg)
% hReg      - Handle of HandleGraphics object to build registry in.
% See spm_XYZreg for more information.
%
%_______________________________________________________________________
% @(#)spm_orthviews.m	2.25 John Ashburner & Matthew Brett 01/03/19
% UM Bios mods
%   Changed default interpolation to nearest neighbor (0)
%   Added gap to center of cross hairs
%   Added support for MIPs
%   Based on diag_orthviews v1.19, 01/12/31
% @(#)spm_orthviews.m	1.4 TE Nichols 02/12/17


% The basic fields of st are:
%         n        - the number of images currently being displayed
%         vols     - a cell array containing the data on each of the
%                    displayed images.
%         Space    - a mapping between the displayed images and the
%                    mm space of each image.
%         bb       - the bounding box of the displayed images.
%         centre   - the current centre of the orthogonal views
%         callback - a callback to be evaluated on a button-click.
%         xhairs   - crosshairs off/on
%         hld      - the interpolation method
%         fig      - the figure that everything is displayed in
%         mode     - the position/orientation of the sagittal view.
%                    - currently always 1
% 
%         st.registry.hReg \_ See spm_XYZreg for documentation
%         st.registry.hMe  /
% 
% For each of the displayed images, there is a non-empty entry in the
% vols cell array.  Handles returned by "spm_orthviews('Image',.....)"
% indicate the position in the cell array of the newly created ortho-view.
% Operations on each ortho-view require the handle to be passed.
% 
% When a new image is displayed, the cell entry contains the information
% returned by spm_vol (type help spm_vol for more info).  In addition,
% there are a few other fields, some of which I will document here:
% 
%         premul - a matrix to premultiply the .mat field by.  Useful
%                  for re-orienting images.
%         window - either 'auto' or an intensity range to display the
%                  image with.
% 
%         area   - Position of region containing three Ax,Co,Sag axes
%
%         ax     - a cell array containing an element for the three
%                  views.  The fields of each element are handles for
%                  the axis, image and crosshairs:
%                               ax  - Axis handle for one image  
%			         d  - image handle               
%			        lx  - horizontal cross-hair line 
%			        ly  - vertical cross-hair line   
% 
%         blobs  - optional.  Is there for using to superimpose blobs.
%                  vol     - 3D array of image data
%                  mat     - a mapping from vox-to-mm (see spm_vol, or
%                            help on image formats).
%                  max     - maximum intensity for scaling to.  If it
%                            does not exist, then images are auto-scaled.
% 		 
%                  There are two colouring modes: full colour, and split
%                  colour.  When using full colour, there should be a
%                  'colour' field for each cell element.  When using
%                  split colourscale, there is a handle for the colorbar
%                  axis.
% 		 
%                  colour  - if it exists it contains the
%                            red,green,blue that the blobs should be
%                            displayed in.
%                  cbar    - handle for colorbar (for split colourscale).
%
%         in     - Intensity
%
%

global st;

if isempty(st), reset_st; end;

spm('Pointer','watch');

if nargin == 0, action = ''; end;
action = lower(action);

switch lower(action),
case {'image','mip'},
	H = specify_image(varargin{1});
	if ~isempty(H)
		if length(varargin)>=2, st.vols{H}.area = varargin{2}; end;
		if isempty(st.bb), st.bb = maxbb; end;
		bbox;

		if strcmp(action,'mip')
		  st.vols{H}.mip = struct('bbOld',[],'SpaceOld',[],'hldOld',[],'premulOld',[],...
					  'imgt',[],'imgc',[],'imgs',[]);
		end

		redraw(H);
	end;
	varargout{1} = H;

case 'mipify',
	for i=valid_handles(varargin{1})
	  if isfield(st.vols{i},'mip')
	    st.vols{i} = rmfield(st.vols{i},'mip');
	  else
	    st.vols{i}.mip = struct('bbOld',[],'SpaceOld',[],'hldOld',[],'premulOld',[],...
				    'imgt',[],'imgc',[],'imgs',[]);
	  end
	  redraw(i);
	end

case 'submip',
        handle = varargin{1};

	for i=valid_handles(handle),
	  if st.mode == 1,
	    axpos = get(st.vols{i}.ax{1}.ax,'Position');
	    tmp = get(st.vols{i}.ax{3}.ax,'Position');
	    offx = tmp(1) - sum(axpos([1 3]));
	    H = spm_orthviews('mip',st.vols{i}.fname,[(axpos(1)+axpos(3)+offx) (axpos(2)+0.005) axpos(3) axpos(4)]);
	  end
	end
	varargout{1} = H;

case 'addcolorbar',
	if nargin<2 
	  addcolorbar(valid_handles(1:24));
	elseif nargin<3
	  addcolorbar(valid_handles(varargin{1}));
	elseif nargin<3
	  addcolorbar(valid_handles(varargin{1}),varargin{2});
	end

case 'colorbartitl',
	if nargin<2 
	  setcolorbartitl(valid_handles(1:24),varargin{1:end});
	  updatecolorbar(valid_handles(varargin{1:24}));
	else
	  setcolorbartitl(valid_handles(varargin{1}),varargin{2:end});
	  updatecolorbar(valid_handles(varargin{1}));
	end

case 'addtitle',
	if nargin<2 
	  addtitle(valid_handles(1:24));
	else
	  addtitle(valid_handles(varargin{1}));
	end

case 'rmtitle',
	if nargin<2 
	  rmtitle(valid_handles(1:24));
	else
	  rmtitle(valid_handles(varargin{1}));
	end

case 'clickcolorbar',
	hs = varargin{1};
	M  = varargin{2};
	[c,h] = findcbarclick;
	for i = valid_handles(hs)
	  MnMx = [NaN NaN];
	  if strcmp(M,'max')
	    MnMx(2) = c;
	  else  % min
	    MnMx(1) = c;
	  end
	  updatecolorbar(i,MnMx)
	  redraw(i)
	end
	
case 'colorbar2',
	setscale2(valid_handles(varargin{1}),varargin{2});
	updatecolorbar(valid_handles(varargin{1}));
	show_inten(valid_handles(varargin{1}));

case 'bb',
	if length(varargin)> 0 & all(size(varargin{1})==[2 3]), st.bb = varargin{1}; end;
	bbox;
	redraw_all;

case 'redraw',
	redraw_all;
	eval(st.callback);
	if isfield(st,'registry'),
		spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
	end;

case 'reposition',
	if length(varargin)<1, tmp = findcent;
	else, tmp = varargin{1}; end;
	if length(tmp)==3, st.centre = tmp(:); end;
	redraw_all;
	eval(st.callback);
	if isfield(st,'registry'),
		spm_XYZreg('SetCoords',st.centre,st.registry.hReg,st.registry.hMe);
	end;

case 'setcoords',
	st.centre = varargin{1};
	st.centre = st.centre(:);
	redraw_all;
	eval(st.callback);

case 'space',
	if length(varargin)<1,
		st.Space = eye(4);
		st.bb = maxbb;
		redraw_all;
	else,
		space(varargin{1});
		redraw_all;
	end;
	bbox;

case 'maxbb',
	st.bb = maxbb;
	bbox;
	redraw_all;

case 'resolution',
	resolution(varargin{1});
	bbox;
	redraw_all;

case 'window',
	if length(varargin)<2,
		win = 'auto';
	elseif strcmp(varargin{2},'ask'), 
	  win = spm_input('Range','+1','e','',2);
	elseif length(varargin{2})==2
	  win = varargin{2};
	end
	for i=valid_handles(varargin{1}),
		st.vols{i}.window = win;
	end;
	redraw(varargin{1});

case 'autowindow',
	if length(varargin)<1, hdls = 1:24; else, hdls = varargin{1}; end
	for i=valid_handles(hdls),
		if length(varargin)>1
		  st.vols{i}.window = autowindow(i,varargin{2});
		else
		  st.vols{i}.window = autowindow(i);
		end
	end;
	redraw(hdls);

case 'delete',
	my_delete(varargin{1});

case 'move',
	move(varargin{1},varargin{2});
	% redraw_all;

case 'reset',
	my_reset;

case 'pos',
	if isempty(varargin),
		H = st.centre(:);
	else,
		H = pos(varargin{1});
	end;
	varargout{1} = H;

case 'interp',
	st.hld = varargin{1};
	redraw_all;

case 'xhairs',
	xhairs(varargin{1});

case 'xhairsgap',
	xhairsgap(varargin{1});

case 'register',
	register(varargin{1});

case 'addblobs',
	addblobs(varargin{1}, varargin{2},varargin{3},varargin{4});
	redraw(varargin{1});

case 'addcolouredblobs',
	addcolouredblobs(varargin{1}, varargin{2},varargin{3},varargin{4},varargin{5});

case 'addimage',
	addimage(varargin{1}, varargin{2});
	redraw(varargin{1});

case 'addcolouredimage',
	addcolouredimage(varargin{1}, varargin{2},varargin{3});
 
case 'addtruecolourimage',
 % spm_orthviews('Addtruecolourimage',handle,filename,colourmap,prop,mx,mn)
 % Adds blobs from an image in true colour
 % handle   - image number to add blobs to [default 1]
 % filename of image containing blob data [default - request via GUI]
 % colourmap - colormap to display blobs in [GUI input]
 % prop - intensity proportion of activation cf grayscale [0.4]
 % mx   - maximum intensity to scale to [maximum value in activation image]
 % mn   - minimum intensity to scale to [minimum value in activation image]
 %
 if nargin < 2
   varargin(1) = {1};
 end
 if nargin < 3
   varargin(2) = {spm_get(1, 'img', 'Image with activation signal')};
 end
 if nargin < 4
   actc = [];
   while isempty(actc)
     actc = getcmap(spm_input('Colourmap for activation image', '+1','s'));
   end
   varargin(3) = {actc};
 end
 if nargin < 5
   varargin(4) = {0.4};
 end
 if nargin < 6
   actv = spm_vol(varargin{2});
   varargin(5) = {max([eps maxval(actv)])};
 end
 if nargin < 7
   varargin(6) = {min([0 minval(actv)])};
 end

 addtruecolourimage(varargin{1}, varargin{2},varargin{3}, varargin{4}, ...
		    varargin{5}, varargin{6});
 redraw(varargin{1});

case 'rmblobs',
	rmblobs(varargin{1});
	redraw(varargin{1});

otherwise,
	warning('Unknown action string')
end;

spm('Pointer');
return;


%_______________________________________________________________________
%_______________________________________________________________________
function addtitle(handle)
global st

for i=valid_handles(handle)
	str = spm_str_manip(st.vols{i}.fname,'k25');
	
	if isfield(st.vols{i},'titl_ax')
	  set(st.vols{i}.titl_h,'String',str)
	else
	  area = st.vols{i}.area;
	  ys =  get(st.vols{i}.ax{2}.ax,'Position');
	  ys = ys(2)+ys(4);
	  ax = axes('Parent',st.fig,...
		    'Position',[area(1)+area(3)*0  ys+area(4)*0.015 ...
		                        area(3)*1     area(4)*0.04],...
		    'visible','off');
	  st.vols{i}.titl_ax = ax;
	  st.vols{i}.titl_h  = text(0.01,0,str,'HorizontalAlignment','left','Parent',ax);
	end
end

%_______________________________________________________________________
%_______________________________________________________________________
function rmtitle(handle)
global st
for i = valid_handles(handle)
	if isfield(st.vols{i},'titl_ax')
	  delete(st.vols{i}.titl_ax);
	  st.vols{i} = rmfield(st.vols{i},{'titl_ax','titl_h'});
	end
end


%_______________________________________________________________________
%_______________________________________________________________________
function addcolorbar(handle)
global st
for i=valid_handles(handle)
	if st.mode == 0,
		axpos = get(st.vols{i}.ax{2}.ax,'Position');
	else,
		axpos = get(st.vols{i}.ax{1}.ax,'Position');
	end;
	ax = axes('Parent',st.fig,...
		  'Position',[(axpos(1)+axpos(3)+0.05) (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
		  'Box','on');
	mx = maxval(st.vols{i});
	mn = minval(st.vols{i});
	st.vols{i}.window = [mn mx];
	h_im = image([0 1],[mn mx],[1:64]','Parent',ax);
	set(ax,'YDir','normal','XTickLabel',[],'XTick',[]);
	if isfield(st.vols{i},'cb1titl'),
	  set(get(ax,'ylabel'),'String',st.vols{i}.cb1titl); 
	end
	sc2 = getscale2(i);
	if ~isnan(sc2)
	  ax2 = axes('Parent',st.fig,...
		     'Position',get(ax,'Position'),...
		     'Xlim',[0 1],'Ylim',[mn mx]*sc2,...
		     'YAxisLocation','right',...
		     'XAxisLocation','top',...
		     'color','none','Xtick',[],...
		     'HitTest','off',...
		     'UIContextMenu',get(UIC,'UIContextMenu'));
	  set(ax,'box','off','UserData',ax2);
	  if isfield(st.vols{i},'cb2titl'),
	    set(get(ax2,'ylabel'),'String',st.vols{i}.cb2titl); 
	  end
	end
	st.vols{i}.cbar = ax;

	tmp = sprintf('Cbar-%d-%d',st.fig,i);
	h = uicontextmenu('Tag',tmp,'UserData',ax,'Parent',st.fig);
	tmp = num2str(i);
	uimenu(h,'Label','Adjust intensity window...')
	uimenu(h,'Separator','on','Label','Set max',...
	       'CallBack',['spm_orthviews(''clickcolorbar'',' tmp ',''max'')'],...
	       'Interruptible','off','BusyAction','Cancel');
	uimenu(h,'Separator','off','Label','Set min',...
	       'CallBack',['spm_orthviews(''clickcolorbar'',' tmp ',''min'')'],...
	       'Interruptible','off','BusyAction','Cancel');
	uimenu(h,'Separator','off','Label','Manual',...
	       'CallBack',['spm_orthviews(''window'',' tmp ',''ask'')'],...
	       'Interruptible','off','BusyAction','Cancel');
	uimenu(h,'Separator','off','Label','Reset',...
	       'CallBack',['spm_orthviews(''window'',' tmp ')'],...
	       'Interruptible','off','BusyAction','Cancel');
	uimenu(h,'Separator','on','Label','Set max - all images',...
	       'CallBack',['spm_orthviews(''clickcolorbar'',1:24,''max'')'],...
	       'Interruptible','off','BusyAction','Cancel');
	uimenu(h,'Separator','off','Label','Set min - all images',...
	       'CallBack',['spm_orthviews(''clickcolorbar'',1:24,''min'')'],...
	       'Interruptible','off','BusyAction','Cancel');
	uimenu(h,'Separator','off','Label','Manual - all images',...
	       'CallBack',['spm_orthviews(''window'',1:24,''ask'')'],...
	       'Interruptible','off','BusyAction','Cancel');
	uimenu(h,'Separator','off','Label','Reset - all images',...
	       'CallBack',['spm_orthviews(''window'',1:24)'],...
	       'Interruptible','off','BusyAction','Cancel');
	set(h_im,'UIContextMenu',h)


end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function updatecolorbar(handle,MnMx)
% Optional 2-vector MnMx sets min and max; either may be NaN, indicating
% it should be obtained from the image data.
global st
if nargin<2, MnMx = []; end
for i=valid_handles(handle)
	if isfield(st.vols{i},'cbar')
	  if isempty(MnMx) | any(isnan(MnMx))
	    if strcmp(st.vols{i}.window,'auto'),
	      mx = maxval(st.vols{i});   % This should never happen
	      mn = minval(st.vols{i});
	    else
	      mx = st.vols{i}.window(2);
	      mn = st.vols{i}.window(1);
	    end
	    if ~isempty(MnMx)
	      if ~isnan(MnMx(1)), mn = MnMx(1); end
	      if ~isnan(MnMx(2)), mx = MnMx(2); end
	    end
	  else
	    mn = MnMx(1); 
	    mx = MnMx(2);
	  end
	  st.vols{i}.window = [mn mx];
	  ax = st.vols{i}.cbar;
	  set(ax,'Ylim',[mn mx]);
	  if isfield(st.vols{i},'cb1titl'),
	    set(get(ax,'ylabel'),'String',st.vols{i}.cb1titl); 
	  end
	  set(get(ax,'child'),'Ydata',[mn mx]);
	  sc2 = getscale2(i);
	  ax2 = NaN;
	  if ~isempty(get(ax,'UserData'))
	    % 2nd colorbar axes exists
	    ax2 = get(ax,'UserData');
	    set(ax2,'Ylim',[mn mx]*sc2);
	  elseif ~isnan(sc2)
	    % Create 2nd colorbar axes
	    UIC = get(ax,'Children');
	    ax2 = axes('Parent',st.fig,...
		       'Position',get(ax,'Position'),...
		       'Xlim',[0 1],'Ylim',[mn mx]*sc2,...
		       'YAxisLocation','right',...
		       'XAxisLocation','top',...
		       'color','none','Xtick',[],...
		       'HitTest','off',...
		       'UIContextMenu',get(UIC,'UIContextMenu'));
	    set(ax,'box','off','UserData',ax2);
	  end
	  if ~isnan(ax2) & isfield(st.vols{i},'cb2titl'),
	    set(get(ax2,'ylabel'),'String',st.vols{i}.cb2titl); 
	  end
	end
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addblobs(handle, xyz, t, mat)
global st
for i=valid_handles(handle),
	if ~isempty(xyz),
		rcp      = round(xyz);
		dim      = max(rcp,[],2)';
		off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
		vol      = zeros(dim)+NaN;
		vol(off) = t;
		vol      = reshape(vol,dim);
		st.vols{i}.blobs=cell(1,1);
		if st.mode == 0,
			axpos = get(st.vols{i}.ax{2}.ax,'Position');
		else,
			axpos = get(st.vols{i}.ax{1}.ax,'Position');
		end;
		ax = axes('Parent',st.fig,'Position',[(axpos(1)+axpos(3)+0.25) (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
			'Box','on');
		mx = max([eps max(t)]);
		mn = min([0 min(t)]);
		image([0 1],[mn mx],[1:64]' + 64,'Parent',ax);
		set(ax,'YDir','normal','XTickLabel',[]);
		st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'cbar',ax,'max',mx, 'min',mn);
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addimage(handle, fname)
global st
for i=valid_handles(handle),
	vol = spm_vol(fname);
	mat = vol.mat;
	st.vols{i}.blobs=cell(1,1);
	if st.mode == 0,
		axpos = get(st.vols{i}.ax{2}.ax,'Position');
	else,
		axpos = get(st.vols{i}.ax{1}.ax,'Position');
	end;
	ax = axes('Parent',st.fig,'Position',[(axpos(1)+axpos(3)+0.15) (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
		'Box','on');
	mx = max([eps maxval(vol)]);
	mn = min([0 minval(vol)]);
	image([0 1],[mn mx],[1:64]' + 64,'Parent',ax);
	set(ax,'YDir','normal','XTickLabel',[]);
	st.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'cbar',ax,'max',mx);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addcolouredblobs(handle, xyz, t, mat,colour)
global st
for i=valid_handles(handle),
	if ~isempty(xyz),
		rcp      = round(xyz);
		dim      = max(rcp,[],2)';
		off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
		vol      = zeros(dim)+NaN;
		vol(off) = t;
		vol      = reshape(vol,dim);
		if ~isfield(st.vols{i},'blobs'),
			st.vols{i}.blobs=cell(1,1);
			bset = 1;
		else,
			bset = length(st.vols{i}.blobs)+1;
		end;
		axpos = get(st.vols{i}.ax{2}.ax,'Position');
		mx = max([eps maxval(vol)]);
		mn = min([0 minval(vol)]);
		st.vols{i}.blobs{bset} = struct('vol',vol,'mat',mat,'max',mx,'min',mn,'colour',colour);
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addcolouredimage(handle, fname,colour)
global st
for i=valid_handles(handle),

	vol = spm_vol(fname);
	mat = vol.mat;
	if ~isfield(st.vols{i},'blobs'),
		st.vols{i}.blobs=cell(1,1);
		bset = 1;
	else,
		bset = length(st.vols{i}.blobs)+1;
	end;
	axpos = get(st.vols{i}.ax{2}.ax,'Position');
	mx = max([eps maxval(vol)]);
	mn = min([0 minval(vol)]);
	st.vols{i}.blobs{bset} = struct('vol',vol,'mat',mat,'max',mx,'min',mn,'colour',colour);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function addtruecolourimage(handle,fname,colourmap,prop,mx,mn)
% adds true colour image to current displayed image  
global st
for i=valid_handles(handle),
  vol = spm_vol(fname);
  mat = vol.mat;
  if ~isfield(st.vols{i},'blobs'),
    st.vols{i}.blobs=cell(1,1);
    bset = 1;
  else,
    bset = length(st.vols{i}.blobs)+1;
  end;
% axpos = get(st.vols{i}.ax{2}.ax,'Position');
  c = struct('cmap', colourmap,'prop',prop);
  st.vols{i}.blobs{bset} = struct('vol',vol,'mat',mat,'max',mx,'min',mn,'colour',c);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function rmblobs(handle)
global st
for i=valid_handles(handle),
	if isfield(st.vols{i},'blobs'),
		for j=1:length(st.vols{i}.blobs),
			if isfield(st.vols{i}.blobs{j},'cbar') & ishandle(st.vols{i}.blobs{j}.cbar),
				delete(st.vols{i}.blobs{j}.cbar);
			end;
		end;
		st.vols{i} = rmfield(st.vols{i},'blobs');
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function register(hreg)
global st
tmp = uicontrol('Position',[0 0 1 1],'Visible','off','Parent',st.fig);
h   = valid_handles(1:24);
if ~isempty(h),
	tmp = st.vols{h(1)}.ax{1}.ax;
	st.registry = struct('hReg',hreg,'hMe', tmp);
	spm_XYZreg('Add2Reg',st.registry.hReg,st.registry.hMe, 'spm_orthviews');
else,
	warning('Nothing to register with');
end;
st.centre = spm_XYZreg('GetCoords',st.registry.hReg);
st.centre = st.centre(:);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function xhairs(arg1),
global st
st.xhairs = 0;
opt = 'on';
if ~strcmp(arg1,'on'),
	opt = 'off';
else,
	st.xhairs = 1;
end;
for i=valid_handles(1:24),
	for j=1:3,
		set(st.vols{i}.ax{j}.lx,'Visible',opt); 
		set(st.vols{i}.ax{j}.ly,'Visible',opt);  
	end; 
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function xhairsgap(arg1),
global st
st.xhairsgap = arg1;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function H = pos(arg1)
global st
H = [];
for arg1=valid_handles(arg1),
	is = inv(st.vols{arg1}.premul*st.vols{arg1}.mat);
	H = is(1:3,1:3)*st.centre(:) + is(1:3,4);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function my_reset
global st
if ~isempty(st) & isfield(st,'registry') & ishandle(st.registry.hMe),
	delete(st.registry.hMe); st = rmfield(st,'registry');
end;
my_delete(1:24);
reset_st;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function my_delete(arg1)
global st
for i=valid_handles(arg1),
	kids = get(st.fig,'Children');
	for j=1:3,
		if any(kids == st.vols{i}.ax{j}.ax),
			set(get(st.vols{i}.ax{j}.ax,'Children'),'DeleteFcn','');
			delete(st.vols{i}.ax{j}.ax);
		end;
	end;
	st.vols{i} = [];
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function resolution(arg1)
global st
res      = arg1/mean(svd(st.Space(1:3,1:3)));
Mat      = diag([res res res 1]);
st.Space = st.Space*Mat;
st.bb    = st.bb/res;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function move(handle,pos)
global st
for handle = valid_handles(handle),
	st.vols{handle}.area = pos;
end;
bbox;
% redraw(valid_handles(handle));
return;
%_______________________________________________________________________
%_______________________________________________________________________
function bb = maxbb
global st
mn = [Inf Inf Inf];
mx = -mn;
for i=valid_handles(1:24),
	bb = [[1 1 1];st.vols{i}.dim(1:3)];
	c = [	bb(1,1) bb(1,2) bb(1,3) 1
		bb(1,1) bb(1,2) bb(2,3) 1
		bb(1,1) bb(2,2) bb(1,3) 1
		bb(1,1) bb(2,2) bb(2,3) 1
		bb(2,1) bb(1,2) bb(1,3) 1
		bb(2,1) bb(1,2) bb(2,3) 1
		bb(2,1) bb(2,2) bb(1,3) 1
		bb(2,1) bb(2,2) bb(2,3) 1]';
	tc = st.Space\(st.vols{i}.premul*st.vols{i}.mat)*c;
	tc = tc(1:3,:)';
	mx = max([tc ; mx]);
	mn = min([tc ; mn]);
end;
bb = [mn ; mx];
return;
%_______________________________________________________________________
%_______________________________________________________________________
function space(arg1)
global st
if ~isempty(st.vols{arg1})
	num = arg1;
	Mat = st.vols{num}.premul(1:3,1:3)*st.vols{num}.mat(1:3,1:3);
	Mat = diag([sqrt(sum(Mat.^2)) 1]);
	Space = (st.vols{num}.mat)/Mat;
	bb = [1 1 1;st.vols{num}.dim(1:3)];
	bb = [bb [1;1]];
	bb=bb*Mat';
	bb=bb(:,1:3);
	st.Space  = Space;
	st.bb = bb;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function H = specify_image(arg1, arg2)
global st
H=[];
ok = 1;
eval('V = spm_vol(arg1);','ok=0;');
if ok == 0,
	fprintf('Can not use image "%s"\n', arg1);
	return;
end;

ii = 1;
while ~isempty(st.vols{ii}), ii = ii + 1; end;

DeleteFcn = ['spm_orthviews(''Delete'',' num2str(ii) ');'];
V.ax = cell(3,1);
for i=1:3,
	ax = axes('Visible','off','DrawMode','fast','Parent',st.fig,'DeleteFcn',DeleteFcn,...
		'YDir','normal');
	d  = image(0,'Tag','Transverse','Parent',ax,...
		'DeleteFcn',DeleteFcn);
	set(ax,'Ydir','normal');
	lx = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
	ly = line(0,0,'Parent',ax,'DeleteFcn',DeleteFcn);
	if ~st.xhairs,
		set(lx,'Visible','off');
		set(ly,'Visible','off');
	end;
	V.ax{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
end;
V.premul    = eye(4);
V.window    = 'auto';
st.vols{ii} = V;

H = ii;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function bbox
global st
Dims = diff(st.bb)'+1;

TD = Dims([1 2])';
CD = Dims([1 3])';
if st.mode == 0, SD = Dims([3 2])'; else, SD = Dims([2 3])'; end;

un    = get(st.fig,'Units');set(st.fig,'Units','Pixels');sz=get(st.fig,'Position');set(st.fig,'Units',un);
sz    = sz(3:4);
sz(2) = sz(2)-40;

for i=valid_handles(1:24),
	area = st.vols{i}.area(:);
	area = [area(1)*sz(1) area(2)*sz(2) area(3)*sz(1) area(4)*sz(2)];
	if st.mode == 0,
		sx   = area(3)/(Dims(1)+Dims(3))/1.02;
	else,
		sx   = area(3)/(Dims(1)+Dims(2))/1.02;
	end;
	sy   = area(4)/(Dims(2)+Dims(3))/1.02;
	s    = min([sx sy]);

	offy = (area(4)-(Dims(2)+Dims(3))*1.02*s)/2 + area(2);
	sky = s*(Dims(2)+Dims(3))*0.02;
	if st.mode == 0,
		offx = (area(3)-(Dims(1)+Dims(3))*1.02*s)/2 + area(1);
		skx = s*(Dims(1)+Dims(3))*0.02;
	else,
		offx = (area(3)-(Dims(1)+Dims(2))*1.02*s)/2 + area(1);
		skx = s*(Dims(1)+Dims(2))*0.02;
	end;

	DeleteFcn = ['spm_orthviews(''Delete'',' num2str(i) ');'];

	% Transverse
	set(st.vols{i}.ax{1}.ax,'Units','pixels', ...
		'Position',[offx offy s*Dims(1) s*Dims(2)],...
		'Units','normalized','Xlim',[0 TD(1)]+0.5,'Ylim',[0 TD(2)]+0.5,...
		'Visible','on','XTick',[],'YTick',[]);

	% Coronal
	set(st.vols{i}.ax{2}.ax,'Units','Pixels',...
		'Position',[offx offy+s*Dims(2)+sky s*Dims(1) s*Dims(3)],...
		'Units','normalized','Xlim',[0 CD(1)]+0.5,'Ylim',[0 CD(2)]+0.5,...
		'Visible','on','XTick',[],'YTick',[]);

	% Sagittal
	if st.mode == 0,
		set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
			'Position',[offx+s*Dims(1)+skx offy s*Dims(3) s*Dims(2)],...
			'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
			'Visible','on','XTick',[],'YTick',[]);
	else,
		set(st.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
			'Position',[offx+s*Dims(1)+skx offy+s*Dims(2)+sky s*Dims(2) s*Dims(3)],...
			'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
			'Visible','on','XTick',[],'YTick',[]);
	end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function redraw_all
global st
redraw(1:24);
return;
%_______________________________________________________________________
function mx = maxval(vol)
if isstruct(vol),
	mx = -Inf;
	for i=1:vol.dim(3),
		tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
		imx = max(tmp(find(finite(tmp))));
		if ~isempty(imx),mx = max(mx,imx);end
	end;
else,
	mx = max(vol(find(finite(vol))));
end;
%_______________________________________________________________________
function mn = minval(vol)
if isstruct(vol),
        mn = Inf;
        for i=1:vol.dim(3),
                tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
		imn = min(tmp(find(finite(tmp))));
		if ~isempty(imn),mn = min(mn,imn);end
        end;
else,
        mn = min(vol(find(finite(vol))));
end;

%_______________________________________________________________________
function mo = modeval(vol,T)

mo = spm_mode(vol);


%_______________________________________________________________________
function qs = quantvals(vol,q1,q2)

if isstruct(vol),
	[IQR Qs n] = spm_iqr(vol);
	binw       = 2*IQR/n^(1/3);
	[n x]      = spm_histvol(vol,-binw);
else,
	IQR   = iqr(vol(finite(vol)));
	n     = sum(finite(vol));
	binw  = 2*IQR/n^(1/3);
        [n x] = hist(vol(:),range(vol)/binw);
end;
nbin  = length(n);
cn    = cumsum(n)/sum(n);
qs    = [max([max([find(cn <= q1)]) 1])  ...
	 min([max([find(cn <= q2)]) nbin])];
if diff(qs)~=0
  qs   = x(qs);
else
  qs   = x(qs(1))+(x(2)-x(1))/2*[-1 1];
end

%_______________________________________________________________________
%_______________________________________________________________________
function redraw(arg1)
global st
bb   = st.bb;
Dims = diff(bb)'+1;
is   = inv(st.Space);
cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);

for i = valid_handles(arg1),

	if isfield(st.vols{i},'mip')

	  TD = Dims([1 2]);
	  CD = Dims([1 3]);
	  if st.mode==0
	    SD = Dims([3 2]);
	  else
	    SD = Dims([2 3]);
	  end
	  
	  mip = st.vols{i}.mip;
	  
	  % Do we need to recreate mip?
	  if (isempty(mip.bbOld) | any(any(mip.bbOld~=bb)) | ...
	      isempty(mip.SpaceOld) | any(any(mip.SpaceOld~=st.Space)) | ...
	      isempty(mip.hldOld) | mip.hldOld~=st.hld | ...
	      isempty(mip.premulOld) | any(any(mip.premulOld~=st.vols{i}.premul)))
	      
	    M = st.vols{i}.premul*st.vols{i}.mat;
	    TM0 = [	1 0 0 -bb(1,1)+1
			0 1 0 -bb(1,2)+1
			0 0 1 -bb(1,3)+1
			0 0 0 1];
	    TM = inv(TM0*(st.Space\M));

	    ok=1;
	    eval('[imgt,imgc,imgs] = spm_project_vol(st.vols{i},TM,Dims,st.hld);','ok=0');
	    if (ok)
	      imgt = imgt';
	      imgc = imgc';
	      if st.mode == 1
		imgs = fliplr(imgs');
	      end

	      st.vols{i}.mip.imgt = imgt;
	      st.vols{i}.mip.imgc = imgc;
	      st.vols{i}.mip.imgs = imgs;

	      st.vols{i}.mip.bbOld = bb;
	      st.vols{i}.mip.SpaceOld = st.Space;
	      st.vols{i}.mip.hldOld = st.hld; 
	      st.vols{i}.mip.premulOld = st.vols{i}.premul;
	      
	    end
	  else % used saved mip

	    imgt = mip.imgt;
	    imgc = mip.imgc;
	    imgs = mip.imgs;
	    
	    ok=1;
	    
	  end % need to recreate mip?

	  hide_inten(i);

	else % Not mip:
	  
	  M = st.vols{i}.premul*st.vols{i}.mat;
	  TM0 = [	1 0 0 -bb(1,1)+1
			0 1 0 -bb(1,2)+1
			0 0 1 -cent(3)
			0 0 0 1];
	  TM = inv(TM0*(st.Space\M));
	  TD = Dims([1 2]);
	  
	  CM0 = [	1 0 0 -bb(1,1)+1
			0 0 1 -bb(1,3)+1
			0 1 0 -cent(2)
			0 0 0 1];
	  CM = inv(CM0*(st.Space\M));
	  CD = Dims([1 3]);
	  
	  if st.mode ==0,
	    SM0 = [	0 0 1 -bb(1,3)+1
			0 1 0 -bb(1,2)+1
			1 0 0 -cent(1)
			0 0 0 1];
	    SM = inv(SM0*(st.Space\M)); SD = Dims([3 2]);
	  else,
	    SM0 = [	0  1 0 -bb(1,2)+1
			0  0 1 -bb(1,3)+1
			1  0 0 -cent(1)
			0  0 0 1];
	    SM0 = [	0 -1 0 +bb(2,2)+1
			0  0 1 -bb(1,3)+1
			1  0 0 -cent(1)
			0  0 0 1];
	    SM = inv(SM0*(st.Space\M));
	    SD = Dims([2 3]);
	  end;

	  ok=1;
	  eval('imgt  = (spm_slice_vol(st.vols{i},TM,TD,st.hld))'';','ok=0;');
	  eval('imgc  = (spm_slice_vol(st.vols{i},CM,CD,st.hld))'';','ok=0;');
	  eval('imgs  = (spm_slice_vol(st.vols{i},SM,SD,st.hld))'';','ok=0;');
	  
	  TMO  = [	1 0 0 -cent(1)
			0 1 0 -cent(2)
			0 0 1 -cent(3)
			0 0 0 1];
	  TM   = inv(TMO*(st.Space\M)); TM = TM(1:3,4);
	  st.vols{i}.in = spm_sample_vol(st.vols{i},TM(1),TM(2),TM(3),st.hld);
	  show_inten(i);
	  addtitle(i);
	  
	end % mip yes/no

	if (ok==0), fprintf('Image "%s" can not be resampled\n', st.vols{i}.fname);
	else,
		if strcmp(st.vols{i}.window,'auto'),
			mx = -Inf; mn = Inf;
			if ~isempty(imgt),
				mx = max([mx max(max(imgt))]);
				mn = min([mn min(min(imgt))]);
			end;
			if ~isempty(imgc),
				mx = max([mx max(max(imgc))]);
				mn = min([mn min(min(imgc))]);
			end;
			if ~isempty(imgs),
				mx = max([mx max(max(imgs))]);
				mn = min([mn min(min(imgs))]);
			end;
			if mx==mn, mx=mn+eps; end;
		else,
			mx = st.vols{i}.window(2);
			mn = st.vols{i}.window(1);
			r=min([mn mx]);imgt = max(imgt,r); r=max([mn mx]);imgt = min(imgt,r);
			r=min([mn mx]);imgc = max(imgc,r); r=max([mn mx]);imgc = min(imgc,r);
			r=min([mn mx]);imgs = max(imgs,r); r=max([mn mx]);imgs = min(imgs,r);
		end;

		updatecolorbar(i);
		
	
		if isfield(st.vols{i},'blobs'),
			if ~isfield(st.vols{i}.blobs{1},'colour'),
				% Add blobs for display using the split colourmap
				scal = 64/(mx-mn);
				dcoff = -mn*scal;
				imgt = imgt*scal+dcoff;
				imgc = imgc*scal+dcoff;
				imgs = imgs*scal+dcoff;

				vol = st.vols{i}.blobs{1}.vol;
				mat = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				if isfield(st.vols{i}.blobs{1},'max'),
					mx = st.vols{i}.blobs{1}.max;
				else,
					mx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
					st.vols{i}.blobs{1}.max = mx;
				end;
				if isfield(st.vols{i}.blobs{1},'min'),
					mn = st.vols{i}.blobs{1}.min;
				else,
					mn = min([0 minval(st.vols{i}.blobs{1}.vol)]);
					st.vols{i}.blobs{1}.min = mn;
				end;

				vol  = st.vols{i}.blobs{1}.vol;
				M    = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				tmpt = spm_slice_vol(vol,inv(TM0*(st.Space\M)),TD,[0 NaN])';
				tmpc = spm_slice_vol(vol,inv(CM0*(st.Space\M)),CD,[0 NaN])';
				tmps = spm_slice_vol(vol,inv(SM0*(st.Space\M)),SD,[0 NaN])';

				sc   = 64/(mx-mn);
				off  = 65.51-mn*sc;
				msk  = find(finite(tmpt)); imgt(msk) = off+tmpt(msk)*sc;
				msk  = find(finite(tmpc)); imgc(msk) = off+tmpc(msk)*sc;
				msk  = find(finite(tmps)); imgs(msk) = off+tmps(msk)*sc;

				cmap = get(st.fig,'Colormap');
				if size(cmap,1)~=128
					figure(st.fig)
					spm_figure('Colormap','gray-hot')
				end;
			elseif isstruct(st.vols{i}.blobs{1}.colour),
				% Add blobs for display using a defined
                                % colourmap

				% colourmaps
				gryc = [0:63]'*ones(1,3)/63;
				actc = ...
				    st.vols{1}.blobs{1}.colour.cmap;
				actp = ...
				    st.vols{1}.blobs{1}.colour.prop;
				
				% scale grayscale image, not finite -> black
				imgt = scaletocmap(imgt,mn,mx,gryc,65);
				imgc = scaletocmap(imgc,mn,mx,gryc,65);
				imgs = scaletocmap(imgs,mn,mx,gryc,65);
				gryc = [gryc; 0 0 0];
				
				% get max for blob image
				vol = st.vols{i}.blobs{1}.vol;
				mat = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				if isfield(st.vols{i}.blobs{1},'max'),
					cmx = st.vols{i}.blobs{1}.max;
				else,
					cmx = max([eps maxval(st.vols{i}.blobs{1}.vol)]);
				end;

				% get blob data
				vol  = st.vols{i}.blobs{1}.vol;
				M    = st.vols{i}.premul*st.vols{i}.blobs{1}.mat;
				tmpt = spm_slice_vol(vol,inv(TM0*(st.Space\M)),TD,[0 NaN])';
				tmpc = spm_slice_vol(vol,inv(CM0*(st.Space\M)),CD,[0 NaN])';
				tmps = spm_slice_vol(vol,inv(SM0*(st.Space\M)),SD,[0 NaN])';
				
				% actimg scaled round 0, black NaNs
                                topc = size(actc,1)+1;
				tmpt = scaletocmap(tmpt,-cmx,cmx,actc,topc);
				tmpc = scaletocmap(tmpc,-cmx,cmx,actc,topc);
				tmps = scaletocmap(tmps,-cmx,cmx,actc,topc);
				actc = [actc; 0 0 0];
				
				% combine gray and blob data to
                                % truecolour
				imgt = reshape(actc(tmpt(:),:)*actp+ ...
					       gryc(imgt(:),:)*(1-actp), ...
					       [size(imgt) 3]);
				imgc = reshape(actc(tmpc(:),:)*actp+ ...
					       gryc(imgc(:),:)*(1-actp), ...
					       [size(imgc) 3]);
				imgs = reshape(actc(tmps(:),:)*actp+ ...
					       gryc(imgs(:),:)*(1-actp), ...
					       [size(imgs) 3]);
				
				
			else,
				% Add full colour blobs - several sets at once
				scal = 1/(mx-mn);
				dcoff = -mn*scal;
				imgt = repmat(imgt*scal+dcoff,[1,1,3]);
				imgc = repmat(imgc*scal+dcoff,[1,1,3]);
				imgs = repmat(imgs*scal+dcoff,[1,1,3]);

				for j=1:length(st.vols{i}.blobs),
					vol = st.vols{i}.blobs{j}.vol;
					mat = st.vols{i}.premul*st.vols{i}.blobs{j}.mat;
					if isfield(st.vols{i}.blobs{j},'colour'),
						colour = st.vols{i}.blobs{j}.colour;
					else,
						colour = [1 0 0];
					end;
					if isfield(st.vols{i}.blobs{j},'max'),
						mx = st.vols{i}.blobs{j}.max;
					else,
						mx = max([eps max(st.vols{i}.blobs{j}.vol(:))]);
						st.vols{i}.blobs{j}.max = mx;
					end;
					if isfield(st.vols{i}.blobs{j},'min'),
						mn = st.vols{i}.blobs{j}.min;
					else,
						mn = min([0 min(st.vols{i}.blobs{j}.vol(:))]);
						st.vols{i}.blobs{j}.min = mn;
					end;

					vol  = st.vols{i}.blobs{j}.vol;
					M    = st.vols{i}.premul*st.vols{i}.blobs{j}.mat;
					tmpt = spm_slice_vol(vol,inv(TM0*(st.Space\M)),TD,[0 NaN])'/(mx-mn)+mn;
					tmpc = spm_slice_vol(vol,inv(CM0*(st.Space\M)),CD,[0 NaN])'/(mx-mn)+mn;
					tmps = spm_slice_vol(vol,inv(SM0*(st.Space\M)),SD,[0 NaN])'/(mx-mn)+mn;
					tmpt(find(~finite(tmpt))) = 0;
					tmpc(find(~finite(tmpc))) = 0;
					tmps(find(~finite(tmps))) = 0;

					tmp  = cat(3,tmpt*colour(1),tmpt*colour(2),tmpt*colour(3));
					imgt = (repmat(1-tmpt,[1 1 3]).*imgt+tmp);
					tmp = find(imgt<0); imgt(tmp)=0; tmp = find(imgt>1); imgt(tmp)=1;

					tmp  = cat(3,tmpc*colour(1),tmpc*colour(2),tmpc*colour(3));
					imgc = (repmat(1-tmpc,[1 1 3]).*imgc+tmp);
					tmp = find(imgc<0); imgc(tmp)=0; tmp = find(imgc>1); imgc(tmp)=1;

					tmp  = cat(3,tmps*colour(1),tmps*colour(2),tmps*colour(3));
					imgs = (repmat(1-tmps,[1 1 3]).*imgs+tmp);
					tmp = find(imgs<0); imgs(tmp)=0; tmp = find(imgs>1); imgs(tmp)=1;
				end;
			end;
		else,
			scal = 64/(mx-mn);
			dcoff = -mn*scal;
			imgt = imgt*scal+dcoff;
			imgc = imgc*scal+dcoff;
			imgs = imgs*scal+dcoff;
		end;

		callback = 'spm_orthviews(''Reposition'');';

		% GapFrac controls size of gap, specifying the fraction of the median dimension
		GapFrac = st.xhairsgap;
		GapSp   = (diff(bb)'+1)*[-1 NaN 1]*GapFrac;
		Gap     = (cent-bb(1,:)'+1)*[1 NaN 1] + [1 1 1]'*median(GapSp) - 0.5;
		GapPt   = [1 NaN 1]; % opposing points
		
		set(st.vols{i}.ax{1}.d,'ButtonDownFcn',callback, 'Cdata',imgt);
		set(st.vols{i}.ax{1}.lx,'ButtonDownFcn',callback,...
				  'Xdata',[0 Gap(1,:) TD(1)]+0.5,...
				  'Ydata',[1 GapPt 1]*(cent(2)-bb(1,2)+1));
		set(st.vols{i}.ax{1}.ly,'ButtonDownFcn',callback,...
				  'Ydata',[0 Gap(2,:) TD(2)]+0.5,...
				  'Xdata',[1 GapPt 1]*(cent(1)-bb(1,1)+1));

		set(st.vols{i}.ax{2}.d,'ButtonDownFcn',callback, 'Cdata',imgc);
		set(st.vols{i}.ax{2}.lx,'ButtonDownFcn',callback,...
				  'Xdata',[0 Gap(1,:) CD(1)]+0.5,...
				  'Ydata',[1 GapPt 1]*(cent(3)-bb(1,3)+1));
		set(st.vols{i}.ax{2}.ly,'ButtonDownFcn',callback,...
				  'Ydata',[0 Gap(3,:) CD(2)]+0.5,...
				  'Xdata',[1 GapPt 1]*(cent(1)-bb(1,1)+1));

		set(st.vols{i}.ax{3}.d,'ButtonDownFcn',callback,'Cdata',imgs);
		if st.mode ==0,
			set(st.vols{i}.ax{3}.lx,'ButtonDownFcn',callback,...
					  'Xdata',[0 Gap(3,:) SD(1)]+0.5,...
					  'Ydata',[1 GapPt 1]*(cent(2)-bb(1,2)+1));
			set(st.vols{i}.ax{3}.ly,'ButtonDownFcn',callback,...
					  'Ydata',[0 Gap(2,:) SD(2)]+0.5,...
					  'Xdata',[1 GapPt 1]*(cent(3)-bb(1,3)+1));
		else,
			Gap(2,:) = (bb(2,2)+1-cent(2))*[1 NaN 1] + median(GapSp);

			set(st.vols{i}.ax{3}.lx,'ButtonDownFcn',callback,...
					  'Xdata',[0 Gap(2,:) SD(1)]+0.5,...
					  'Ydata',[1 GapPt 1]*(cent(3)-bb(1,3)+1));
			set(st.vols{i}.ax{3}.ly,'ButtonDownFcn',callback,...
					  'Ydata',[0 Gap(3,:) SD(2)]+0.5,...
					  'Xdata',[1 GapPt 1]*(bb(2,2)+1-cent(2)));
		end;
	end;
end;
drawnow;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function [c,h] = findcbarclick
global st
obj    = get(st.fig,'CurrentObject');
c      = [];
h      = [];
if isempty(obj), return, end
for i=valid_handles(1:24)
  if isfield(st.vols{i},'cbar')
    Cb = st.vols{i}.cbar;
    if (Cb == obj) | (Cb == get(obj,'Parent'))
      c = get(Cb,'CurrentPoint');
      h  = i;
      break
    end
  end
end
if ~isempty(c)
  c = c(1,2);
end
return;
%_______________________________________________________________________
%_______________________________________________________________________
function centre = findcent
global st
obj    = get(st.fig,'CurrentObject');
centre = [];
cent   = [];
cp     = [];
for i=valid_handles(1:24),
	for j=1:3,
		if ~isempty(obj),
			if any([st.vols{i}.ax{j}.d  ...
				st.vols{i}.ax{j}.lx ...
				st.vols{i}.ax{j}.ly]== obj)
				cp = get(get(obj,'Parent'),'CurrentPoint');
			elseif (st.vols{i}.ax{j}.ax == obj),
				cp = get(obj,'CurrentPoint');
			end;
		end;
		if ~isempty(cp),
			cp   = cp(1,1:2);
			is   = inv(st.Space);
			cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);
			switch j,
				case 1,
				cent([1 2])=[cp(1)+st.bb(1,1)-1 cp(2)+st.bb(1,2)-1];
				case 2,
				cent([1 3])=[cp(1)+st.bb(1,1)-1 cp(2)+st.bb(1,3)-1];
				case 3,
				if st.mode ==0,
					cent([3 2])=[cp(1)+st.bb(1,3)-1 cp(2)+st.bb(1,2)-1];
				else,
					cent([2 3])=[st.bb(2,2)+1-cp(1) cp(2)+st.bb(1,3)-1];
				end;
			end;
			break;
		end;
	end;
	if ~isempty(cent), break; end;
end;
if ~isempty(cent), centre = st.Space(1:3,1:3)*cent(:) + st.Space(1:3,4); end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function handles = valid_handles(handles)
global st;
handles = handles(:)';
handles = handles(find(handles<=24 & handles>=1 & ~rem(handles,1)));
for h=handles,
	if isempty(st.vols{h}), handles(find(handles==h))=[]; end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function reset_st
global st
fig     = spm_figure('FindWin','Graphics');
bb      = [ [-78 78]' [-112 76]' [-50 85]' ];
st      = struct('n', 0, 'vols',[], 'bb',bb,'Space',eye(4),'centre',[0 0 0],'callback',';','xhairs',1,'xhairsgap',0.05,'hld',0,'fig',fig,'mode',1);
st.vols = cell(24,1);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function img = scaletocmap(inpimg,mn,mx,cmap,miscol)
if nargin < 5, miscol=1;end
cml = size(cmap,1);
scf = (cml-1)/(mx-mn);
img = round((inpimg-mn)*scf)+1;
img(find(img<1))=1; 
img(find(img>cml))=cml;
img(~finite(img)) = miscol;
%_______________________________________________________________________
%_______________________________________________________________________
function cmap = getcmap(acmapname)
% get colormap of name acmapname
if ~isempty(acmapname)
  cmap = evalin('base',acmapname,'[]');
  if isempty(cmap) % not a matrix, is .mat file?
    [p f e] = fileparts(acmapname);
    acmat = fullfile(p, [f '.mat']);
    if exist(acmat, 'file')
      s = struct2cell(load(acmat));
      cmap = s{1};
    end
  end
end
if size(cmap, 2)~=3
  warning('Colormap was not an N by 3 matrix')
  cmap = [];
end
return

% Set coloar bar label(s)
%_______________________________________________________________________
function setcolorbartitl(hndl,titl1,titl2)
global st
if nargin<2, titl1='';end
if nargin<3, titl2='';end
for handle = valid_handles(hndl)
	st.vols{handle}.cb1titl = titl1;
	st.vols{handle}.cb2titl = titl2;
end

% Set an alternate scalefactor for image; for two-sided color bar
%_______________________________________________________________________
function setscale2(hndl,sc)
global st
for handle = valid_handles(hndl)
	st.vols{handle}.scale2 = sc;
end

%-Get alternate scalefactor for image; for two-sided color bar
%_______________________________________________________________________
function sc = getscale2(hndl)
global st
for handle = valid_handles(hndl)
	if isfield(st.vols{handle},'scale2')
	  sc = st.vols{handle}.scale2;
	else
	  sc = NaN;
	end
end

%-Automagically window
% Centers window about the mode
% Optionally takes quantile cut offs
%_______________________________________________________________________
function win = autowindow(hndl,qt)
global st
if nargin<2, qt = [0.1 0.9]; end
for handle = valid_handles(hndl)
	mn = minval(st.vols{handle});
	mx = maxval(st.vols{handle});
	qs = quantvals(st.vols{handle},qt(1),qt(2));
	mo = modeval(st.vols{handle},qs);
	if (mn > -(mx-mn)/64)
	  % If min is virtualy positive, use mode trick
	  win = [0 2*mo];
	else
	  d   = max(abs(mo-qs));
	  win = [mo-d mo+d];
	end
end


%-Display intensity under cursor
%_______________________________________________________________________
function show_inten(hndl)
global st
for i = valid_handles(hndl)
	str = sprintf('%g',st.vols{i}.in);
	sc  = getscale2(hndl);
	if isnan(sc)
	  str2 = ''; 
	else
	  str2 = sprintf('[%g,  %g]',st.vols{i}.in,st.vols{i}.in*sc);
	end

	Chng = 0;
	if isfield(st.vols{i},'in_ax')
	  % Have we added a 2nd colorscale since last time?
	  Chng = xor(length(st.vols{i}.in_h)==1,isnan(sc)); 
	end
	if isfield(st.vols{i},'in_ax') & ~Chng
	  set(st.vols{i}.in_h(1),'String',str)
	  if ~isempty(str2)
	    set(st.vols{i}.in_h(2),'String',str2)
	  end
	else
	  if isfield(st.vols{i},'in_ax')
	    delete(st.vols{i}.in_ax)
	  end
	  area = st.vols{i}.area;
	  if st.mode == 0,
	    axpos = get(st.vols{i}.ax{2}.ax,'Position');
	  else,
	    axpos = get(st.vols{i}.ax{1}.ax,'Position');
	  end;

%	  axpos1 = get(st.vols{i}.ax{1}.ax,'Position');
%	  axpos3 = get(st.vols{i}.ax{3}.ax,'Position');
%	  ax = axes('Parent',st.fig,...
%		    'Position',[axpos3(1)     axpos1(2)...
%		                axpos3(3)     area(4)*0.08],...
%		    'visible','off');

	  ys =  get(st.vols{i}.ax{2}.ax,'Position');
	  ys = ys(2)+ys(4);
	  ax = axes('Parent',st.fig,...
		    'Position',[area(1)+area(3)*0  ys+area(4)*0.005 ...
		                        area(3)*1     area(4)*0.04],...
		    'visible','off');
	  st.vols{i}.in_ax = ax;
	  if isnan(sc)
	    h = text(1,0.25,str,'HorizontalAlignment', 'right','Parent',ax);
	  else
	    h = [text(1,0.25,str2,'HorizontalAlignment', 'right','Parent',ax)];
	  end
	  st.vols{i}.in_h = h;
	end
end


%-Hide intensity
%_______________________________________________________________________
function hide_inten(hndl)
global st
for i = valid_handles(hndl)
	if isfield(st.vols{i},'in_ax')
	  delete(st.vols{i}.in_ax);
	  st.vols{i} = rmfield(st.vols{i},{'in_ax','in_h'});
	end
end

