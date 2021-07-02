function progressbar(varargin)
    %   set interval between one upgrade and the other in seconds
    drawInterval = 0.5;
    %   horizontal padding - horizontal distance between the bar and other
    %   bar/borders in pixels.
    hpad = 10;
    %   vertical padding - horizontal distance between the bar and other
    %   bar/borders in pixels.
    vpad = 10;
    %   bar width and height
    barWidth = 400;
    barHeight = 30;
    %   progressbar(Name, Tasks) //
    %       - If there is a Name progress bar, add Tasks Task Completed to Name ProgressBar
    %       - Else creates a 'Name' progress bar of Tasks TotalTasks
    %   progressbar(Name, "Delete") //
    %       - If there is a Name progress bar Delete it
    %   progressbar(Name) //
    %       - If there is a Name progress bar, add 1 Task Completed to Name ProgressBar
    %   progressbar(Tasks) //
    %       - Add Tasks to last created Progressbar Task Completed
    %   progressbar() //
    %       - Add Tasks to last created Progressbar Task Completed
    % Upgrade data Struct
    persistent progressbarData lastBar lastUpgrade progFigure upgradeNow;
    
    % If progressbarData/lastCPUtime doesn't exist, create it
    if (isempty(progressbarData))
        progressbarData = containers.Map;
    end
    if (isempty(progFigure))
        progFigure = " ";
    end
    if (isempty(lastUpgrade))
        lastUpgrade = 0;
    end
    if (isempty(upgradeNow))
        upgradeNow=true;
    end
    barName = '';
    tasksDone = 1;
    removeBar = false;
    reset=false;
    if nargin == 2
        % First must be chararray/string, second num or 'remove'
        if(ischar(varargin{1}) || isstring(varargin{1}))
            barName = varargin{1};
        else
            return
        end
        if(isnumeric(varargin{2}))
            tasksDone = varargin{2};
        elseif(ischar(varargin{2}) || isstring(varargin{2}))
            if(varargin{2}=="delete")
                removeBar=true;
            elseif(varargin{2}=="reset")
                reset=true;
            else
                return
            end
        else
            return
        end
    elseif nargin == 1
        % First must be chararray/string or num
        if(isnumeric(varargin{1}))
            tasksDone = varargin{1};
        elseif(ischar(varargin{1}) || isstring(varargin{1}))
            barName = varargin{1};
        else
            return
        end
    elseif nargin == 0
        % already set
    else
        return
    end

    %If bar not assigned, use last Created bar
    if(barName=="")
        barName=lastBar;
    end

    %If no bar, do nothing

    if(isempty(barName))
        return
    end

    %Check if barName exist and create a new bar with taskDone if it
    %doesn't.
    %Delete bar id removeBar is true
    if isKey(progressbarData, barName)
        if (removeBar)
            
            remove(progressbarData, barName);
            upgradeNow=true;
        else
            barStruct=progressbarData(barName);
            barStruct.tasksDone = barStruct.tasksDone + tasksDone;
            if(barStruct.tasksDone<barStruct.totalTasks)
                progressbarData(barName)=barStruct;
            else
                validDelete(barStruct, 'label')
                validDelete(barStruct, 'patch')
                validDelete(barStruct, 'axis')
                remove(progressbarData, barName);
                upgradeNow=true;
            end
        end
    else
        lastBar = barName;
        barStruct = struct;
        barStruct.totalTasks = tasksDone;
        barStruct.tasksDone = 0;
        progressbarData(barName)=barStruct;
        upgradeNow=true;
    end
    nBar = length(progressbarData);
    if(nBar<1 || reset )
        if ishandle(progFigure)
            delete(progFigure) % Close progress bar
        end
        clear progressbarData lastBar lastUpgrade progFigure upgradeNow
        return
    end
    %Avoid upgrading too often
    CPUTime = cputime;
    if(CPUTime > lastUpgrade+drawInterval)
        % Upgrade graphic representation if more than drawintervalPassed
        figureWidth = barWidth+2*hpad;
        figureHeight = barHeight*nBar+vpad*(2*nBar+1);
        if(~isgraphics(progFigure))
            upgradeNow=true;
            progFigure = figure('Resize', 'off', 'MenuBar', 'none', ...
                'Position', [figureWidth figureHeight figureWidth figureHeight]);
        end
        progFigure.Position(3) = figureWidth;
        progFigure.Position(4) = figureHeight;
        keys = progressbarData.keys;
        if upgradeNow
            figure(progFigure)
            %redraw axis if one bar was removed or if first call
            for i=1:nBar
                barStruct = progressbarData(keys{i});
                currentBarHeight = figureHeight-(i)*(2*vpad+barHeight);
                if(isfield(barStruct,'axis'))
                    if(isgraphics(barStruct.axis))
                        validDelete (barStruct, 'axis');
                        validDelete (barStruct, 'patch');
                    end
                end
                barStruct.axis = axes( 'Units', 'pixels', 'XLim', [0 1], 'YLim', [0 1], 'Box', 'on', 'ytick', [], 'xtick', [],...
                    'Position', [ hpad currentBarHeight barWidth barHeight]);
                barStruct.patch = patch( 'XData', [0 0 0 0], 'YData', [0 0 1 1] );
                barStruct.label = text(0, 1.25,keys{i});
                progressbarData(keys{i}) = barStruct;
            end
            upgradeNow=false;
        end
        for i=1:nBar
            barStruct = progressbarData(keys{i});
            barStruct.patch.XData(2:3) = barStruct.tasksDone / barStruct.totalTasks;
            barStruct.label.String = sprintf('%s - %0.0f / %0.0f',...
                keys{i},barStruct.tasksDone,barStruct.totalTasks);
            progressbarData(keys{i}) = barStruct;
        end
        lastUpgrade = CPUTime;
        drawnow
    end
end

function validDelete(struct, name)
    if(isfield(struct, name))
        delete(struct.(name))
    end
end