function [ mcmap, cdata_mapped ] = multicmap( h, cmaps, clims )
% MULTICMAP Apply multiple colormaps to image objects
%
%   [MCMAP,CDATA]=MULTICMAP(H,CMAPS,CLIMS) given handles to 
%   image objects plotted in the current axes, along with 
%   specification of individual colormaps for the image 
%   objects, this function combines the individual colormaps 
%   into a single multi-colormap matrix, and computes re-mapped 
%   CData matrices for individual image objects, so that they
%   point to their respective colormaps within the multi-colormap.
%
%   This function assumes that CData of each image object 
%   is a 2D matrix of indexes to the figures colormap 
%   (which is the default behavior in MATLAB), as opposed
%   to CData being a 3D matrix of true color RGB image values.
%
%   Not that this routine is a generalization of examples 
%   presented in [1].
%
%   Known limitations: currently multiple subplots in a single
%                      figure are not supported.
%
%   Inputs
%           H numeric array (or structure) containing handles
%             to image objects plotted in the current axes
%
%           CMAPS column cell array (or structure) of colormaps 
%                 for the image objects specified in H. Each 
%                 colormap should conform to MATLAB's colormap
%                 specifications, i.e., it should be an m-by-3 
%                 matrix of real numbers between 0.0 and 1.0. 
%                 Each row of a given colormap matrix is an RGB 
%                 vector that defines one color. See "doc colormap" 
%                 for further information regarding MATLAB's colormaps.
%
%           CLIMS m-by-2 matrix, or structure of two element 
%                 row vectors, each containing lower and upper 
%                 limits for each image object's CData indexes.
%                 See "doc imagesc" for detailed diagram and 
%                 explanation of CLIMS.
%
%   Outputs 
%           MCMAP m-by-3 colormap matrix composed of individual 
%                 colormaps supplied in CMAPS.
%
%           CDATA re-mapped CData matrices as members of a cell array
%                 for individual image objects specified in H.
%
%   References
%           [1] MathWorks: Product Support, 
%               "1215 - Using Multiple Colormaps in a Single Figure"
%               url: http://www.mathworks.com/support/tech-notes/1200/1215.html
%
%   Example
%            clear all; close all; clc; 
%        
%            % Read the sample sample images 
%            data.earth = load( 'earth' );
%            data.penny = load( 'penny' );
%        
%            % Get indexed images into variables
%            photo.earth = data.earth.X;
%            photo.penny = data.penny.P;
%        
%            % Get dimensions of each image
%            dims.earth = size( photo.earth );
%            dims.penny = size( photo.penny );
%        
%            % Define x- and y-points for each image
%            x.earth = [ 1:dims.earth(1) ];
%            x.penny = [ 1:dims.penny(1) ] + round((dims.earth(1)-dims.penny(1))/2);
%            y.earth = 1:dims.earth(2);
%            y.penny = [ 1:dims.penny(2) ] + round((dims.earth(2)-dims.penny(2))/2);
%        
%            % Define colormaps for each image like so:
%            maps.earth = data.earth.map;
%            maps.penny = hot(64);
%        
%            % Define transparency for visible parts of each image
%            alpha.earth = 1;        % fully opaque
%            alpha.penny = 0.7;      % 30% transparency
%        
%            % Define which parts of the second image are to be fully transparent (invisible)
%            AlphaData.penny = ones( size( photo.penny ) ) * alpha.penny;
%            AlphaData.penny( photo.penny<5*min(photo.penny(:)) ) = 0;
%        
%            % Plot images with their respective colormaps
%            hfig = figure( 'Position', [ 400 10 600 600 ], 'PaperPositionMode', 'auto', 'color', 'w' );
%        
%            % Make sure both images get retained in the current axes
%            hold on;
%        
%            % Plot images and retain handles to the image objects
%            h.earth = image( x.earth, y.earth, photo.earth );
%            h.penny = image( x.penny, y.penny, photo.penny );
%        
%            % Apply transparency settings
%            set( h.earth, 'AlphaData', alpha.earth );
%            set( h.penny, 'AlphaData', AlphaData.penny );
%        
%            % Apply axes limits 
%            xlim( [ min(struct2array(x)) max(struct2array(x)) ] );
%            ylim( [ min(struct2array(y)) max(struct2array(y)) ] );
%        
%            % Make sure our images are not up-side-down
%            axis ij square off
%        
%            % Apply colormaps to their respective image objects
%            multicmap( h, maps );
%
%   See also EXAMPLE_TWO_IMAGES, EXAMPLE_THREE_IMAGES, EXAMPLE_SPEECH.

%   Author: Kamil Wojcicki, UTD, February 2012.


    % check for correct number of input arguments
    if nargin<2 || nargin>3
        error( sprintf('Incorrect number of input arguments.\nType "help %s" for usage help.\n', mfilename) ); 
    end

    % if image object handles were passed-in as values for members 
    % of a structure, then convert the structure to a numeric array
    if isstruct( h ) 
        h = struct2array( h );
    end 

    % determine number of image objects
    N = length( h );

    % if individual colormaps were passed-in as matrices for members
    % of a structure, then convert the structure to a cell array
    if isstruct( cmaps )
        cmaps = struct2cell( cmaps );
    end 

    % combine colormaps for individual image objects into
    % a single multi-colormap (as a m-by-3 matrix)
    mcmap = cell2mat( cmaps );


    % determine if data limits were specified
    if nargin==3

        % if the data limits were specified as 1-by-2 arrays
        % for member of a structure, then convert the structure
        % to a cell array of data limits
        if isstruct( clims )
            clims = struct2cell( clims );
        end

        % if the data limits are as 1-by-2 arrays, members of 
        % a cell array, then convert the cell array into m-by-2 matrix
        if iscell( clims )
            clims = cell2mat( clims );
        end

        % determine the dimensionality of the data limits matrix
        [ rows, columns ] = size( clims );

        % for each image object there should be lower
        % and upper limit for the data (i.e., exactly 2)
        assert( columns==2 );

        % if only one array of data limit values was supplied
        % then use it for all image objects by replicating 
        % it for each
        if rows==1
            clims = repmat( clims, N, 1 );

        % otherwise, make sure that the number of supplied
        % data limits matched the number of image objects
        else
            assert( rows==N );
        end

    end


    % for each image object compute CData re-mapped 
    % to index to the new multi-colormap.
    for n = 1:N

        % get the image data (assumed to be indices
        % into the figure's colormap)
        cdata = get( h(n), 'CData' );

        % if climits were supplied, apply them
        if nargin==3
            cdata( cdata<clims(n,1) ) = clims(n,1);
            cdata( cdata>clims(n,2) ) = clims(n,2);
        end

        % figure-out the maximum, minimum and range
        % of the current image object's CData
        cdata_max = max( cdata(:) );
        cdata_min = min( cdata(:) );
        cdata_range = cdata_max - cdata_min;

        % get the size of the colormap for the 
        % current image object
        cmap_size = size( cmaps{n}, 1 );

        % re-map the CData for the current image object to its colormap
        cdata_mapped{n} = min( cmap_size, round((cmap_size-1)*(cdata-cdata_min)/cdata_range)+1 );

        % re-map the CData for the current image object to its 
        % corresponding colormap within the multi-colormap
        if n>1, cdata_mapped{n} = cdata_mapped{n} + max( cdata_mapped{n-1}(:) ); end

    end


    % if no outputs have been requested, then update data for
    % for the individual image objects to the re-mapped data,
    % also set data limits (for plotting) and apply the new colormap
    if nargout==0

        % for each image object update CData to re-mapped data
        for n = 1:N
            set( h(n), 'CData', cdata_mapped{n} );
        end

        % set CData axis limits (for plotting)
        caxis( [ min(cdata_mapped{1}(:)) max(cdata_mapped{end}(:)) ] );

        % apply the new colormap
        colormap( mcmap );

    end


% EOF