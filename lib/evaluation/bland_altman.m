%########################################################################
%
%	- PPGI Toolbox - 
%   A MATLAB toolbox for Photoplethysmography Imaging (PPGI)
%
% Author   : Christian S. Pilz
% Company  : The Nature of Space of Time
% Date     : 07.05.2019
%
% Contact  : cpi@partofthestars.com
% Web Page : www.partofthestars.com
%
% Version  : Alpha RA 1.0
%
%########################################################################
%
%	bland_altman.m:
%
% Description:
%
%   A static MATLAB class to draw a Blant-Altman and correlation graph 
%   for two datasets.
%
% Usage:
%
% bland_altman(data1, data2)
% bland_altman(data1, data2,label) - Names of data sets. Formats can be
%   - {'Name1'}
%   - {'Name1, 'Name2'}
%   - {'Name1, 'Name2', 'Units'}
% bland_altman(data1, data2,label,tit,gnames)
% bland_altman(data1, data2,label,tit,gnames,corrinfo) - specifies what
% information to display on the correlation chart as a cell of string in 
% order of top to bottom. The following codes are available:
%  - 'eq'     slope and intercept equation
%  - 'int%'   intercept as % of mean values
%  - 'r'      pearson r-value
%  - 'r2'     pearson r-value squared
%  - 'rho'    Spearman rho value
%  - 'SSE'    sum of squared error
%  - 'n'      number of data points used
%  * if not specified or empty, default is: {'eq';'r2';'SSE';'n'}
% bland_altman(data1, data2,label,tit,gnames,corrinfo,BAinfo) - specifies what
% information to display on the Bland-Altman chart as for corrinfo, but
% with the following codes:
%  - 'RPC'    reproducibility coefficient (1.96*SD)
%  - 'RPC(%)' reproducibility coefficient and % of mean values
%  - 'CV'     coefficient of variation (SD of mean values in %)
%  * if not specified or empty, default is: {'RPC(%)';'CV'}
% bland_altman(data1, data2,label,tit,gnames,corrinfo,BAinfo,limits) -
% specifies the axes limits:
%  - scalar - lower limit (eg. 0)
%  - [min max] - specifies minimum and maximum
%  - 'tight' - minimum and maximum of data.
%  - 'auto' - plot default. {default option}
% bland_altman(data1, data2,label,tit,gnames,corrinfo,BAinfo,limits,colors) -
% specify the order of group colors. (eg. 'brg' for blue, red, green) or
% RGB columns.
% bland_altman(data1, data2,label,tit,gnames,corrinfo,BAinfo,limits,colors,symbols)
% - specify the order of symbols. (eg. 'sod' for squares, circles, dots).
% Alternatively can be set to 'Num' to display the subject number.
% bland_altman(fig, ...) - specify a figure handle in which to
% display
% figure in which the Bland-Altman and correlation will be displyed 
% bland_altman(ah, ...) - specify an axes which will be replaced by the 
% Bland-Altman and correlation exes.
% cr = BlandAltman(...) - return the coefficient of reproducibility
% (1.96*sd)
% [cr fig] = bland_altman(...) - also return the figure handles
% [cr fig sstruct] = bland_altman(...) - also return the structure of
% statistics for the analysis
%

classdef bland_altman
    
   methods(Static)
       function [cr, fig, sstruct] = draw(varargin)
            % Added feature to control y-axis scaling on BA figure. Follow code below
            % for details. May become a parameter with demand.
            fixylim = true;

            if isscalar(varargin{1}) && isequal(size(varargin{1}),[1 1]) && ishandle(varargin{1})
                shift = 1;
                fig = varargin{1};
            else
                shift = 0;
                fig = [];
            end
            data1 = varargin{shift+1};
            data2 = varargin{shift+2};
            if nargin>=shift+3
                label = varargin{shift+3};
            else
                label = '';
            end
            if nargin>=shift+4
                tit = varargin{shift+4};
            else
                tit = '';
            end
            if nargin>=shift+5
                gnames = varargin{shift+5};
            else
                gnames = '';
            end
            if nargin>=shift+6 && ~isempty(varargin{shift+6})
                if ischar(varargin{shift+6})
                    corrinfo = varargin(shift+6);
                else
                    corrinfo = varargin{shift+6};
                end
            else
                corrinfo = {'eq';'r2';'SSE';'n'};
            end
            if nargin>=shift+7 && ~isempty(varargin{shift+7})
                if ischar(varargin{shift+7})
                    BAinfo = varargin(shift+7);
                else
                    BAinfo = varargin{shift+7};
                end
            else
                BAinfo = {'RPC(%)';'CV'};
            end

            if nargin>=shift+8
                axesLimits = varargin{shift+8};
            else
                axesLimits = 'auto';
            end

            if nargin>=shift+9 && ~isempty(varargin{shift+9})
                colors = varargin{shift+9};
            else
                colors = 'rbgmcky';
            end
            if ~ischar(colors)
                if size(colors,2)~=3
                    if size(colors,1)==3
                        colors = colors';
                    else
                        error('Colors must be specified in either character codes or RGB');
                    end
                end
            elseif size(colors,1)==1
                colors = colors';
            end

            if nargin>=shift+10 && ~isempty(varargin{shift+10})
                symb = varargin{shift+10};
            else
                symb = 'sodp^v';
            end

            markersize = 4;
            xdatamode = 1; % use the (1) mean of data1 and data2 or (2) data1 for x position on bland altman

            units = '';
            if iscell(label)
                if length(label)==1
                    labelx = label{1};
                    labely = label{1};
                    labelm = label{1};
                    labeld = ['\Delta ' label{1}];
                elseif length(label)==2
                    labelx = label{1};
                    labely = label{2};
                    labelm = ['Mean ' label{1} ' & ' label{2}];
                    labeld = [label{2} ' - ' label{1}];
                else % units also provided
                    units = label{3};
                    labelx = [label{1} ' (' units ')'];
                    labely = [label{2} ' (' units ')'];
                    labelm = ['Mean ' label{1} ' & ' label{2} ' (' units ')'];
                    labeld = [label{2} ' - ' label{1} ' (' units ')'];
                end	
            else
                labelx = label;
                labely = label;
                labelm = label;
                labeld = ['\Delta ' label];
            end
            if isempty(units)
                unitsstr = '';
            else
                unitsstr = [' ' units];
            end

            s = size(data1);
            if ~isequal(s,size(data2));
                error('data1 and data2 must have the same size');
            end

            switch length(s)
                case 1
                    s = [s 1 1];
                case 2
                    s = [s 1];
                case 3
                otherwise
                    error('Data have too many dimension');
            end
            n = s(1); % number of elements in each group
            groups = numel(data1)/n;
            if size(colors,1)<s(3)
                error('More groups than colors specified. Use the colors input variable to specify colors for each group.');
            end
            if ~strcmpi(symb,'Num') && length(symb)<s(2)
                error('More subgroups than symbolss specified. Use the symbols input variable to specify symbols for each subgroup, or use the ''Num'' option.');
            end

            data1 = reshape(data1, [numel(data1),1]);
            data2 = reshape(data2, [numel(data2),1]);
            mask = isfinite(data1) & isnumeric(data1) & isfinite(data2) & isnumeric(data2);

            if isempty(fig)
                fig = figure;
                set(fig,'units','centimeters','position',[3 3 20 10],'color','w');
                cah = subplot(121); set(gca,'FontSize',8)
                dah = subplot(122); set(gca,'FontSize',8)
            elseif strcmpi(get(fig,'type'),'figure')
                cah = subplot(121);
                dah = subplot(122);
            elseif strcmpi(get(fig,'type'),'axes')
                ah = fig;
                pos = get(ah,'position');
                fig = get(ah,'parent');
                delete(ah);
                cah = axes('parent',fig,'position',[pos(1) pos(2) pos(3)/2 pos(4)]);
                dah = axes('parent',fig,'position',[pos(1)+pos(3)/2 pos(2) pos(3)/2 pos(4)]);
            else
                error('What in tarnations is the handle that was passed to Bland-Altman????')
            end
            set(cah,'tag','Correlation Plot');
            set(dah,'tag','Bland Altman Plot');

            %% Correlation
            hold(cah,'on');
            for groupi=1:groups
                if strcmpi(symb,'Num')
                    for i=1:n
                        text(data1((groupi-1)*n+i),data2((groupi-1)*n+i),num2str(i),'parent',cah,'fontsize',markersize,'color',colors(floor((groupi-1)/s(2))+1,:), 'HorizontalAlignment','Center', 'VerticalAlignment','Middle');
                    end
                else
                    if s(3)==1
                        marker = symb(1);
                        color = colors(groupi,:);
                    else
                        marker = symb(rem(groupi-1,s(2))+1);
                        color = colors(floor((groupi-1)/s(2))+1,:);
                    end
                    ph=plot(cah,data1((groupi-1)*n+(1:n)),data2((groupi-1)*n+(1:n)),marker,'color',color);
                    set(ph,'markersize',markersize);
                end
            end
            % Linear regression
            [polyCoefs, S] = polyfit(data1(mask),data2(mask),1);
            r = corrcoef(data1(mask),data2(mask)); r=r(1,2);
            rho = corr(data1(mask),data2(mask),'type','Spearman');
            N = sum(mask);
            SSE = sqrt(sum((polyval(polyCoefs,data1(mask))-data2(mask)).^2)/(N-2));

            if ischar(axesLimits)
                if strcmpi(axesLimits,'Auto')
                    % Workaround - Add invisible minimum and maximum point to fix Auto axes limits (text
                    % does not count for axis('auto')
                    if strcmpi(symb,'Num')
                        mindata = min( min(data1(mask)), min(data2(mask)) );
                        maxdata = max( max(data1(mask)), max(data2(mask)) );
                        ph = plot(cah, [mindata maxdata], [mindata maxdata], '.', 'Visible','on');
                    end
                    axesLimits = axis(cah); 
                    axesLimits(1) = min(axesLimits(1),axesLimits(3));
                    axesLimits(2) = max(axesLimits(2),axesLimits(4));
                    if strcmpi(symb,'Num')
                        delete(ph);
                    end
                elseif strcmpi(axesLimits,'Tight')
                    axesLimits(1) = min( min(data1(mask)), min(data2(mask)) );
                    axesLimits(2) = max( max(data1(mask)), max(data2(mask)) );
                else
                    error(['Unknown axes limit option (' axesLimits ') detected.']);
                end
            else
                if length(axesLimits)==1
                    a = axis(cah);
                    axesLimits(2) = max(a(2),a(4));
                else
                    % Do nothing
                end
            end
            axesLimits(3) = axesLimits(1);
            axesLimits(4) = axesLimits(2);

            axis(cah,axesLimits); axis(cah,'square');
            plot(cah,axesLimits(1:2), polyval(polyCoefs,axesLimits(1:2)),'-k');
            h = plot(cah,axesLimits(1:2),axesLimits(1:2),':'); set(h,'color',[0.6 0.6 0.6]);
            if 0 % Add 95% CI lines
                xfit = axesLimits(1):(axesLimits(2)-axesLimits(1))/100:axesLimits(2);
                [yfit, delta] = polyconf(polyCoefs,xfit,S);
                h = [plot(cah,xfit,yfit+delta);...
                    plot(cah,xfit,yfit-delta)];
                set(h,'color',[0.6 0.6 0.6],'linestyle','-');
            end
            corrtext = {};
            for i=1:length(corrinfo)
                switch lower(corrinfo{i})
                    case 'eq'
                        if polyCoefs(2)>0
                            corrtext = [corrtext; ['y=' num2str(polyCoefs(1),3) 'x+' num2str(polyCoefs(2),3)]];
                        else
                            corrtext = [corrtext; ['y=' num2str(polyCoefs(1),3) 'x' num2str(polyCoefs(2),3)]];
                        end
                    case 'int%', corrtext = [corrtext; ['intercept=' num2str(polyCoefs(2)/mean(data1+data2)*2*100,3) '%']];
                    case 'r2', corrtext = [corrtext; ['r^2=' num2str(r^2,4)]];
                    case 'r', corrtext = [corrtext; ['r=' num2str(r,4)]];
                    case 'rho', corrtext = [corrtext; ['rho=' num2str(rho,4)]];
                    case 'sse', corrtext = [corrtext; ['RMSE=' num2str(SSE,2) unitsstr]];
                    case 'n', corrtext = [corrtext; ['n=' num2str(N)]];
                end
            end
            text(axesLimits(1)+0.01*(axesLimits(2)-axesLimits(1)),axesLimits(1)+0.9*(axesLimits(2)-axesLimits(1)),corrtext,'parent',cah);
            xlabel(cah,labelx); ylabel(cah,labely);

            %% Differences
            set(dah,'units','normalized');
            hold(dah,'on');
            for groupi=1:groups
                d1 = data1((groupi-1)*n+(1:n));
                d2 = data2((groupi-1)*n+(1:n));
                dif = d2-d1;
                if strcmpi(symb,'Num')
                    for i=1:n
                        if xdatamode==1
                            text(mean([d1(i),d2(i)]), dif(i), num2str(i), 'parent',dah,'fontsize',markersize,'color',colors(floor((groupi-1)/s(2))+1,:));
                        else
                            text(d1(i), dif(i), num2str(i), 'parent',dah,'fontsize',markersize,'color',colors(floor((groupi-1)/s(2))+1,:));
                        end
                    end
                else
                    if s(3)==1
                        marker = symb(1);
                        color = colors(groupi,:);
                    else
                        marker = symb(rem(groupi-1,s(2))+1);
                        color = colors(floor((groupi-1)/s(2))+1,:);
                    end

                    if xdatamode==1
                        ph = plot(dah,mean([d1,d2],2),dif,marker,'color',color);
                    else
                        ph = plot(dah,d1,dif,marker,'color',color);
                    end
                    set(ph,'markersize',markersize);
                end
            end
            axis(dah,'square')
            xlabel(dah,labelm); ylabel(dah,labeld);

            % add std-dev lines
            st = std(data2(mask)-data1(mask));
            mn = mean(data2(mask)-data1(mask));
            [h, p] = ttest(data2(mask)-data1(mask),0);
            cr = 1.96*st;

            % fix limits to +/- data limit
            if fixylim 
                a = [axesLimits(1:2) [-1 1]*abs(axesLimits(2)-axesLimits(1))/2];
                axis(dah, a);
            end

            plot(a(1:2),mn+[0 0],'k')
            % plot(a(1:2),mn+st*[1 1],'k')
            % plot(a(1:2),mn-st*[1 1],'k')
            plot(a(1:2),mn+cr*[1 1],':k')
            plot(a(1:2),mn-cr*[1 1],':k')
            a = axis(dah);
            if fixylim
                fontsize = 6;
                text(a(2),mn+cr, [num2str(mn+cr,2) ' (+1.96SD)'],'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize);
                text(a(2),mn,[num2str(mn,2) ' [p=' num2str(p,2) ']'],'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize);
                text(a(2),mn-cr, [num2str(mn-cr,2) ' (-1.96SD)'],'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize);
            else
                fontsize = 8;
                text(a(2),mn+cr,{'+1.96SD',num2str(mn+cr,2)},'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize);
                text(a(2),mn,{num2str(mn,2),['p=' num2str(p,2)]},'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize);
                text(a(2),mn-cr,{num2str(mn-cr,2),'-1.96SD'},'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize);
            end
            BAtext = {};

            for i=1:length(BAinfo)
                switch lower(BAinfo{i})
                    case 'rpc', BAtext = [BAtext; ['{\bfRPC: ' num2str(cr,2) unitsstr '}']];
                    case 'rpc(%)', BAtext = [BAtext; ['{\bfRPC: ' num2str(cr,2) unitsstr '} (' num2str(100*cr/mean((data2(mask)+data1(mask))/2),2) '%)']];
                    case 'cv', BAtext = [BAtext; ['CV: ' num2str(100*st/mean((data1(mask)+data2(mask))/2),2) '%']];
                    case 'p', BAtext = [BAtext; ['p-value: ' num2str(p,4) ]]; 
                    case 'ks' % Kolmogorov-Smirnov test that difference-data is Gaussian
                        ddata = data1(:)-data2(:);
                        [h, p] = kstest((ddata-mean(ddata))/std(ddata));
                        BAtext = [BAtext; ['KS p-value: ' num2str(p)]];
                    case 'kurtosis' % Kolmogorov-Smirnov test that difference-data is Gaussian
                        ddata = data1(:)-data2(:);
                        BAtext = [BAtext; ['kurtosis: ' num2str(kurtosis(ddata))]];
                end
            end
            text(a(2),a(4),BAtext,'interpreter','tex','HorizontalAlignment','right','VerticalAlignment','top');

            if ~isempty(tit)
                h = suptitle(tit);
                set(h,'interpreter','tex');
            end

            % Add legend
            if ~strcmpi(symb,'Num') && ~isempty(gnames)
                lh = legend('show');
                if iscell(gnames)
                    if length(gnames)==2 
                        if iscell(gnames{1}) 
                            temp = cell(1,groups);
                            for groupi=1:length(gnames{1})
                                for j=1:length(gnames{2})
                                    temp{groupi+(j-1)*length(gnames{1})} = [gnames{1}{groupi} '-' gnames{2}{j}];
                                end
                            end	
                            gnames = temp;
                        elseif iscell(gnames{2})
                            gnames = strcat(gnames{1}, '-', gnames{2});
                        end
                    end
                end
                cpos = get(cah,'Position');
                dpos = get(dah,'Position');
                set(cah,'Position',cpos+[0 0.07 0 0]);
                set(dah,'Position',dpos+[0 0.07 0 0]);
                set(lh,'string',gnames,'orientation','horizontal');
                drawnow;
                set(lh,'units','normalized');
                pos = get(lh,'position'); pos = min(pos(3),0.9);
                set(lh,'position',[(1-pos)/2 0.02 pos 0.05]);
            else
            % 	set(lh,'visible','off');
            end

            if nargout>2
                sstruct = struct('N',N,...
                    'CR', cr,...
                    'r',r,...
                    'r2',r^2,...
                    'SSE',SSE,...
                    'rho',rho,...
                    'Slope',polyCoefs(1),...
                    'Intercept',polyCoefs(2));
            end
       end
   end
end

