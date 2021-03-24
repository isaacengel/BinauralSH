function violins = violinplot(data, cats, varargin)
% EDITED from Bastian Bechtold's violinplot function to support 
% subcategories and a few other little things
% Last edit: 8th May 2020, Isaac Engel - isaac.engel(at)imperial.ac.uk
%
%Violinplots plots violin plots of some data and categories
%   VIOLINPLOT(DATA) plots a violin of a double vector DATA
%
%   VIOLINPLOT(DATAMATRIX) plots violins for each column in
%   DATAMATRIX.
%
%   VIOLINPLOT(TABLE), VIOLINPLOT(STRUCT), VIOLINPLOT(DATASET)
%   plots violins for each column in TABLE, each field in STRUCT, and
%   each variable in DATASET. The violins are labeled according to
%   the table/dataset variable name or the struct field name.
%
%   VIOLINPLOT(TABLE/STRUCT/DATASET, CATEGORIES)
%   same as the above, but plot one violin per subcategory, as indicated 
%   in the vector CATEGORIES. The violins are labeled according to
%   the table/dataset variable name or the struct field name, and
%   color-coded per subcategory.
%
%   VIOLINPLOT(DATAMATRIX, CATEGORYNAMES) plots violins for each
%   column in DATAMATRIX and labels them according to the names in the 
%   cell-of-strings CATEGORYNAMES. If CATEGORYNAMES is a 2-element
%   cell array, the second element indicates the subcategories, which are 
%   plotted with separate color-coded violins.
%
%   VIOLINPLOT(DATA, CATEGORIES) where double vector DATA and vector
%   CATEGORIES are of equal length; plots violins for each category in
%   DATA. If CATEGORIES is a 2-cell array, plot separate color-coded
%   violins for each subcategory, as indicated in the second cell of
%   CATEGORIES.
%
%
%   violins = VIOLINPLOT(...) returns an object array of
%   <a href="matlab:help('Violin')">Violin</a> objects.
%
%   VIOLINPLOT(..., 'PARAM1', val1, 'PARAM2', val2, ...)
%   specifies optional name/value pairs for all violins:
%     'Width'        Width of the violin in axis space.
%                    Defaults to 0.3
%     'Bandwidth'    Bandwidth of the kernel density estimate.
%                    Should be between 10% and 40% of the data range.
%     'ViolinColor'  Fill color of the violin area and data points.
%                    Defaults to the next default color cycle.
%     'ViolinAlpha'  Transparency of the violin area and data points.
%                    Defaults to 0.3.
%     'EdgeColor'    Color of the violin area outline.
%                    Defaults to [0.5 0.5 0.5]
%     'BoxColor'     Color of the box, whiskers, and the outlines of
%                    the median point and the notch indicators.
%                    Defaults to [0.5 0.5 0.5]
%     'MedianColor'  Fill color of the median and notch indicators.
%                    Defaults to [1 1 1]
%     'ShowData'     Whether to show data points.
%                    Defaults to true
%     'ShowNotches'  Whether to show notch indicators.
%                    Defaults to false
%     'ShowMean'     Whether to show mean indicator
%                    Defaults to false
%     'GroupOrder'   Cell of category names in order to be plotted.
%                    Defaults to alphabetical ordering
%     'GroupOrder2'  Cell of subcategory names in order to be plotted.
%                    Defaults to alphabetical ordering
%     'Colors'       N x 3 matrix to set individual colors for each violin.
%                    If there are subcategories, this defines their color
%                    code. If a single value is used, it works like
%                    'ViolinColor'.
%                    With subcategories, it defaults to Parula color map.
%                    Without subcategories, it defaults like ViolinColor.
%     'DataSize'     Size of the data points
%                    Defaults to 36
%     'MedianSize'   Size of the median indicator
%                    Defaults to 36
%     'LineWidth'    Width of whiskers and violin edges
%                    Defaults to 0.5 points
%     'MeanWidth'    Width of the mean line
%                    Defaults to 1 point

% Copyright (c) 2016, Bastian Bechtold
% This code is released under the terms of the BSD 3-clause license

    hascategories = exist('cats','var') && not(isempty(cats));
    hassubcategories = hascategories && iscell(cats) && ~iscellstr(cats);
    if hascategories
        if ~isvector(cats) && ~iscell(cats)
            error('''cats'' must be a vector or a cell array of vectors.')
        end
        if hassubcategories
            if numel(cats) > 2
                warning('Only the first 2 cells of ''cats'' are used. The rest are ignored.')
            end
            cats2 = cats{2}; % 2nd element is the subcategories
            cats = cats{1}; % 1st element is the main categories
        end
        if isa(data, 'dataset') || isstruct(data) || istable(data)
            if hassubcategories
                warning('Only the first cell of ''cats'' is used. The rest are ignored.')
            end
            if isvector(cats) && numel(cats) ~= size(data,1)
                error('The vector of categories must have as many elements as rows in the dataset')
            end
        end
        % If the cat vectors are columns from tables, convert to vector
        if istable(cats)
            cats = cats{:,:};
        end
        if hassubcategories && istable(cats2)
            cats2 = cats2{:,:};
        end
    end
    
%     % If the data is a single column of a table, convert to vector
%     if istabledata)
    
    %parse the optional grouporder argument 
    %if it exists parse the categories order 
    % but also delete it from the arguments passed to Violin
    grouporder = {};
    idx=find(strcmpi(varargin, 'GroupOrder'));
    if ~isempty(idx) && numel(varargin)>idx
        if iscell(varargin{idx+1})
            grouporder = varargin{idx+1};
            varargin(idx:idx+1)=[];
        else
            error('Second argument of ''GroupOrder'' optional arg must be a cell of category names')
        end
    end
    
    %parse the optional grouporder2 argument 
    %if it exists parse the subcategories order 
    % but also delete it from the arguments passed to Violin
    grouporder2 = {};
    idx=find(strcmpi(varargin, 'GroupOrder2'));
    if ~isempty(idx) && numel(varargin)>idx
        if iscell(varargin{idx+1})
            grouporder2 = varargin{idx+1};
            varargin(idx:idx+1)=[];
        else
            error('Second argument of ''GroupOrder2'' optional arg must be a cell of category names')
        end
    end
    
    %parse the optional colors argument 
    %if it exists parse the color values
    % but also delete it from the arguments passed to Violin
    colors = [];
    idx=find(strcmpi(varargin, 'Colors'));
    if ~isempty(idx) && numel(varargin)>idx
        if ismatrix(varargin{idx+1}) && size(varargin{idx+1},2) == 3
            colors = varargin{idx+1};
            varargin(idx:idx+1)=[];
        else
            error('Second argument of ''Colors'' optional arg must be a Nx3 matrix')
        end
    end
    
    %% tabular data, with or without subcategories
    if isa(data, 'dataset') || isstruct(data) || istable(data)
        if isa(data, 'dataset')
            colnames = data.Properties.VarNames;
        elseif istable(data)
            colnames = data.Properties.VariableNames;
        elseif isstruct(data)
            colnames = fieldnames(data);
        end
        if hascategories
            if isempty(grouporder2)
                seccats = categorical(cats);
            else
                seccats = categorical(cats, grouporder2);
            end
            seccatnames = categories(seccats); 
            nseccats = length(seccatnames);
        end
        catnames = {};
        for n=1:length(colnames)
            if isnumeric(data.(colnames{n}))
                catnames = [catnames colnames{n}];
            end
        end
        
        if ~istable(data) || ~hascategories
            for n=1:length(catnames)
                thisData = data.(catnames{n});
                if ~isempty(colors)
                    if size(colors,1) == 1 % if only one color
                        varargin_ = ['ViolinColor',colors(1,:),varargin];
                    else
                        varargin_ = ['ViolinColor',colors(n,:),varargin];
                    end
                else
                    varargin_ = varargin;
                end
                violins(n) = Violin(thisData, n, varargin_{:});
            end
            xlim([0.5 (length(catnames))+0.5]);
            set(gca, 'xtick', 1:length(catnames), 'xticklabels', catnames);  
        else
            % Set color palette
            if ~isempty(colors)
                c = colors;
            else
                c = parula(nseccats+1);
            end
            % Prepare legend
            dummy_lines = [];
            for m=1:nseccats
            	dummy_lines(end+1) = plot([NaN,NaN],'color',[c(m,:),0.5],'LineWidth',8); hold on
            end
            % Draw violin plots
            for n=1:length(catnames)
                thisData = data.(catnames{n});
                mainpos(n) = 1 + (nseccats+0.5) .* (n-1) + (nseccats-1)/2;
                for m=1:nseccats
                    ind = (n-1)*nseccats + m;
                    rows = seccats == seccatnames(m);
                    thisData_ = thisData(rows,:);
                    secpos = 1+(n-1)*(nseccats+0.5) + (m-1);
                    varargin_ = ['ViolinColor',c(m,:),varargin];
                    violins(ind) = Violin(thisData_, secpos, varargin_{:});
                end
            end
            xlim([0.25 secpos+0.75]);
            set(gca, 'xtick', mainpos, 'xticklabels', catnames);
            % Show legend
            [hh,icons,plots,txt]=legend(dummy_lines,seccatnames,'Location','best');
            for m=1:nseccats
                pos = icons(m).Position;
                icons(m).Position = [0.23 pos(2) 0];
                icons(2*m+nseccats).XData = [0.05 0.2];
                icons(2*m+nseccats-1).XData = [0.05 0.2];
            end
            set(hh,'box','off','Color',[1 1 1 1])
        end

    % 1D data, one category and subcategory for each data point
    elseif hassubcategories && numel(data) == size(cats,1)
        if isempty(grouporder)
            cats = categorical(cats);
        else
            cats = categorical(cats, grouporder);
        end
        if isempty(grouporder2)
            seccats = categorical(cats2);
        else
            seccats = categorical(cats2, grouporder2);
        end
        catnames = categories(cats); % categories() may not work here
        seccatnames = categories(seccats); 
        nseccats = length(seccatnames);
        
        % Set color palette
        if ~isempty(colors)
            c = colors;
        else
            c = parula(nseccats+1);
        end
        
        % Prepare legend
        dummy_lines = [];
        for m=1:nseccats
            dummy_lines(end+1) = plot([NaN,NaN],'color',[c(m,:),0.5],'LineWidth',8); hold on
        end
            
        % Draw violin plots
        for n=1:length(catnames)
            thisCat = catnames(n);
            thisData = data(cats == thisCat);
            thisSeccats = seccats(cats == thisCat);
            mainpos(n) = 1 + (nseccats+1) .* (n-1) + (nseccats-1)/2;
            for m=1:nseccats
                ind = (n-1)*nseccats + m;
                rows = thisSeccats == seccatnames(m);
                thisData_ = thisData(rows,:);
                secpos = (n-1)*(nseccats+1) + m;
                varargin_ = ['ViolinColor',c(m,:),varargin];
                violins(ind) = Violin(thisData_, secpos, varargin_{:});
            end
        end
        xlim([0 (length(catnames))*(nseccats+1)]);
        set(gca, 'xtick', mainpos, 'xticklabels', catnames);
        % Show legend
        [hh,icons,plots,txt]=legend(dummy_lines,seccatnames,'Location','best');
        for m=1:nseccats
            pos = icons(m).Position;
            icons(m).Position = [0.23 pos(2) 0];
            icons(2*m+nseccats).XData = [0.05 0.2];
            icons(2*m+nseccats-1).XData = [0.05 0.2];
        end
        set(hh,'box','off','Color',[1 1 1 1])  
        
    % 1D data, one category for each data point
    elseif hascategories && numel(data) == numel(cats)
        if isempty(grouporder)
            cats = categorical(cats);
        else
            cats = categorical(cats, grouporder);
        end

        catnames = categories(cats);
        for n=1:length(catnames)
            thisCat = catnames{n};
            thisData = data(cats == thisCat);
            violins(n) = Violin(thisData, n, varargin{:});
        end
        xlim([0.5 (length(catnames))+0.5]);
        set(gca, 'xtick', 1:length(catnames), 'xticklabels', catnames);

    % 1D data, no categories
    elseif not(hascategories) && isvector(data)
        violins = Violin(data, 1, varargin{:});
        set(gca, 'xtick', 1);

    % 2D data with categories and subcategories
    elseif ismatrix(data) && hassubcategories
        if isempty(grouporder)
            cats = categorical(cats);
        else
            cats = categorical(cats, grouporder);
        end
        if isempty(grouporder2)
            seccats = categorical(cats2);
        else
            seccats = categorical(cats2, grouporder2);
        end
        catnames = categories(cats); % categories() may not work here
        seccatnames = categories(seccats); 
        nseccats = length(seccatnames);
        
        % Set color palette
        if ~isempty(colors)
            c = colors;
        else
            c = parula(nseccats+1);
        end
        
        % Prepare legend
        dummy_lines = [];
        for m=1:nseccats
            dummy_lines(end+1) = plot([NaN,NaN],'color',[c(m,:),0.5],'LineWidth',8); hold on
        end
        
        % Draw violin plots
        for n=1:length(catnames)
            thisCat = catnames(n);
            thisData = data(:, n);
            mainpos(n) = 1 + (nseccats+1) .* (n-1) + (nseccats-1)/2;
            for m=1:nseccats
                ind = (n-1)*nseccats + m;
                rows = seccats == seccatnames(m);
                thisData_ = thisData(rows,:);
                secpos = (n-1)*(nseccats+1) + m;
                varargin_ = ['ViolinColor',c(m,:),varargin];
                violins(ind) = Violin(thisData_, secpos, varargin_{:});
            end
        end
        xlim([0 (length(catnames))*(nseccats+1)]);
        set(gca, 'xtick', mainpos, 'xticklabels', catnames);
        % Show legend
        [hh,icons,plots,txt]=legend(dummy_lines,seccatnames,'Location','best');
        for m=1:nseccats
            pos = icons(m).Position;
            icons(m).Position = [0.23 pos(2) 0];
            icons(2*m+nseccats).XData = [0.05 0.2];
            icons(2*m+nseccats-1).XData = [0.05 0.2];
        end
        set(hh,'box','off','Color',[1 1 1 1])        
        
    % 2D data with or without categories, without subcategories
    elseif ismatrix(data)
        for n=1:size(data, 2)
            thisData = data(:, n);
            violins(n) = Violin(thisData, n, varargin{:});
        end
        xlim([0.5 size(data, 2)+0.5]);
        set(gca, 'xtick', 1:size(data, 2));
        if hascategories && length(cats) == size(data, 2)
            set(gca, 'xticklabels', cats);
        end

    end

end
