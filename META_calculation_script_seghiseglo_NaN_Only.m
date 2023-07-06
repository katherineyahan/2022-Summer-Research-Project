
%This macro predicts the timing of META, allows the user to confirm/adjust
% this prediction of NaN META according to the graph or projection, then predicts the
% timing of DECOND

clearvars -except Celloutput Emboutput 

% User is prompted to select the source folder from ImageJ macros then sets
% this as the parent folder and uses the
% "Adj_ROIs_and_Projections" folder within it.
folder = uigetdir;
parent=folder;
if ismac
    slash = '/';
elseif ispc
    slash = '\';
end
folder = horzcat(folder, slash, 'Adj_ROIs_and_Projections');
fileList = getAllFiles(folder);

% Find all the '.tif' and all of the 'combproj' files in 'fileList'
boo = strfind(fileList,'tif');
foo = strfind(fileList, 'combproj');
% cellfun Apply function to each cell in cell array
% Finds indices of fileList correpsonding to both a '.tif' and 'combproj' file
moo = find(~cellfun('isempty', boo)&(~cellfun('isempty', foo)));
% Creates a new fileList consisting of only the files at the entries found
% above
if ~isempty(moo)
   combproj_fileList = fileList(moo, 1);
end

%Gets the index of the final slash in the combined projection file names and "_" (just before
%H2Bcombproj_adj.tif) to extract the ID of each projection (embryo name &
%cell ID)
idxStart = strfind(combproj_fileList,slash);
[n,m] = cellfun(@size,idxStart);
first = cellfun(@(x) x(end),idxStart,'UniformOutput',false);
first=cell2mat(first);
first=first+1;
idxLast = strfind(combproj_fileList,'H2Bcombproj_adj.tif');
last = cell2mat(idxLast);
last=last-1;

combproj_fileList = string(combproj_fileList);
combproj_IDs = extractBetween(combproj_fileList, first, last);

combproj_fileList = cellstr(combproj_fileList);
combproj_IDs = cellstr(combproj_IDs);

% Get the name of the genotype corresponding to these files. 
% Find index of last and second last slash character in folder path. 
% The strain name will be between these two indices 
% (ex. "Y:\Leah\OD421\Adj_ROIs_and_Projections")

ind = find(folder==slash);
l = length(ind);
start = (ind(length(ind)-1))+1;
finish = (ind(length(ind)))-1;
strain = folder(start:finish);

finish = (ind(l-1));

% Generates string corresponding to destination of csv output file
% Output file will be saved in selected directory (now parent) as
% "strain_Celloutput.mat" and "strain_Emboutput.mat"

CelloutputPath = strcat(convertCharsToStrings(parent), convertCharsToStrings(slash), convertCharsToStrings(strain), '_Celloutput.mat');
EmboutputPath = strcat(convertCharsToStrings(parent), convertCharsToStrings(slash), convertCharsToStrings(strain), '_Emboutput.mat');

% clearvars -except Celloutput Emboutput CelloutputPath EmboutputPath combproj_fileList combproj_IDs

% Checks if the Celloutput struct exists. If it does not, it checks whether the output files above exist.
% If these files exist, we load them into the workspace.
A = exist('Celloutput');
if A~=1
    if isfile(CelloutputPath)
        load(CelloutputPath, 'Celloutput');
    end
    if isfile(EmboutputPath)
        load(EmboutputPath, 'Emboutput');
    end
    if ~isfile(CelloutputPath) | ~isfile(EmboutputPath)
        error('No Celloutput variable in the work space and file could not be loaded.');
    end
end

[~,kk] = size(Celloutput);
[~, ccc] = size(Emboutput);

% Checks if the 'scored' field exists in Celloutput. If it does not, then
% this field is created and filled with an empty array at every row. If it
% does exist, it checks whether there is something in this field and
% inititalizes the for loop at the first empty row

if ~isfield(Celloutput, 'METAscored')
    for j = 1:1:kk
        Celloutput(j).METAscored = [];
    end
end

% Iterate through all of the entries (all of the cells) in the Cellouput
% struct and predict META for all P cells

for j=1:1:kk
    
    pCells={'EMS','P1','P2','E','MS','P3','C','P4','D'};
    cell=Celloutput(j).cell;
    if any(strcmp(pCells,cell)) && isempty(Celloutput(j).METAscored)
         
         % Import H2B channel measurement values from the Celloutput struct
         H2B_table = Celloutput(j).H2BTable;
         frames = H2B_table.Frame;
        % Skip cells where 3 frames of fewer measured
        if numel(frames) > 3
             fr_Rt = H2B_table.FrameRate_sec_;
             fr_Rt = fr_Rt(1);
             X1 = frames;
        
             %Convert frames to time in seconds, with frame 1 = time 0
             X2 = (X1).*fr_Rt;%%
              
             %Calculate the seghi as a percent of seglo area
             Y8 = H2B_table.seghiArea./H2B_table.segloArea;
             
             %Normalize the Y-axis data
             N8 = Y8./nanmax(Y8);
               
             % To calculate first derivative: calculate the change in N8 over
             % time in frames. To keep the arrays the same size, we fill the
             % first row with an NaN and each subsequent row will be the first
             % derivative (i.e. slope) for the values in N8 at that row and one row earlier.
              
             X2b = NaN;
             X2b(2:max(size(X2)),1) = X2(1:max(size(X2))-1,1);
             N8b = NaN;
             N8b(2:max(size(N8)),1) = N8(1:max(size(N8))-1,1);
             d8 = (N8-N8b)./(X2-X2b);
            
             % Normalize first derivatives by dividing by maximum absolute value:
             max_d8 = max(abs(d8));
             D8 = d8./max_d8;
           
             %Find the frame of the highest local max; it usually occurs
             %before META
             max_D8=max(D8);
             
             %Find the corresponding frame and index
             idx=find(D8==max_D8);
             max_fr=X1(idx);
             [s,~]=size(D8);
             
             if ~isnan(Celloutput(j).NEBDscrd)
                 NEBDfr=Celloutput(j).NEBDscrd/fr_Rt;
                 Celloutput(j).NEBDfr=int8(NEBDfr);
             end
                 
             %If the NEBD is not NaN
             if ~isnan(Celloutput(j).NEBDfr)
               idx_NEBD=find(X1==Celloutput(j).NEBDfr);
           
              %find the local min after the peak;eliminate those who has less
              %than 10 frames and the frames after NEBD are less than 7
                if s>10 && s-idx_NEBD>=7
                   local_mins=islocalmin(D8);
                   upper_bound=find(X1==Celloutput(j).NEBDfr+4);
                   local_maxs=islocalmax(D8); 
                   for o=find(X1==Celloutput(j).NEBDfr+1):1:upper_bound
                     if local_mins(o)==1
                        NEBDscr=D8(o);
                        break;
                     end
                   end
             
                 %if the value is less than the drop after NEBD, then take the frame before it as the metaphase  
                 %The peak has to be after NEBD in order to apply this
                 %criteria
                   if max_fr>Celloutput(j).NEBDfr
                     for k=idx:1:s
                        if local_mins(k)==1
                          if D8(k)<NEBDscr
                           pred=k-1;
                           METApred=k-1;
                           break;
                          end
                        end
                     end
                   else
                        %Try use segloCirc
                        Y1=H2B_table.segloCirc;
                        N1 = Y1./nanmax(Y1);
                        %Find the lowest local min after the NEBD local min
                        min_segloCirc=100;
                        [s,~]=size(Y1);
                        local_mins_segloCirc=islocalmin(Y1);
                        for k=o+1:1:s
                          if local_mins_segloCirc(k)==1
                            if Y1(k)<min_segloCirc
                            min_segloCirc=Y1(k);
                            end
                          end
                        end
                      i=find(Y1==min_segloCirc);
                      METApred=i-1;
                 
                   end
                   %META usually falls between three and seven frames after
                   %NEBD
                   if METApred > idx_NEBD+2 && METApred < idx_NEBD+8
                       x2=X2(METApred);%%
                       Celloutput(j).METApred=x2;
                   else
                       Celloutput(j).METApred=NaN;
                   end
               else
                  Celloutput(j).METApred=NaN;
                end
                 
                  
            else
               Celloutput(j).METApred=NaN;
             end
        end
    end
end

%Score the META that does not have a predicted value
for j=1:1:kk

    pCells={'EMS','P1','P2','E','MS','P3','C','P4','D'};
    cell=Celloutput(j).cell;
    if any(strcmp(pCells,cell)) && isempty(Celloutput(j).METAscored) && isnan(Celloutput(j).METApred)
        H2B_table = Celloutput(j).H2BTable;
        frames = H2B_table.Frame;
        
         fr_Rt = H2B_table.FrameRate_sec_;
         fr_Rt = fr_Rt(1);
         X1 = frames;
        
         %Convert frames to time in seconds, with frame 1 = time 0
          X2 = (X1).*fr_Rt;%%
          
          %SegloCirc
           Y1=H2B_table.segloCirc;
           N1= Y1./nanmax(Y1);
              
          %Calculate the seghi as a percent of seglo area
          Y8 = H2B_table.seghiArea./H2B_table.segloArea;
             
          %Normalize the Y-axis data
           N8 = Y8./nanmax(Y8);
               
           % To calculate first derivative: calculate the change in N8 over
           % time in frames. To keep the arrays the same size, we fill the
           % first row with an NaN and each subsequent row will be the first
           % derivative (i.e. slope) for the values in N8 at that row and one row earlier.
              
           X2b = NaN;
           X2b(2:max(size(X2)),1) = X2(1:max(size(X2))-1,1);
           N8b = NaN;
           N8b(2:max(size(N8)),1) = N8(1:max(size(N8))-1,1);
           d8 = (N8-N8b)./(X2-X2b);
            
           % Normalize first derivatives by dividing by maximum absolute value:
           max_d8 = max(abs(d8));
           D8 = d8./max_d8;
           
           %Find the frame of the highest local max; it usually occurs
           %before META
           max_D8=max(D8);
             
           %Find the corresponding frame and index
           idx=find(D8==max_D8);
        
           x1=Celloutput(j).NEBDscrd;
        
           emb = convertCharsToStrings(Celloutput(j).embName);
           Celloutput(j).embName = emb;
           cell = convertCharsToStrings(Celloutput(j).cell);
           Celloutput(j).cell = cell;
        
           T = struct2table(Emboutput);
        
          % Finds the number of embryos in our dataset and the index of
          % the current embryo. These values will be displayed on the
          % graph so that the user can tell how far along they are.
          allEmbNames = T{:,1};
          a = strfind(allEmbNames, emb);
          b = find(~cellfun('isempty', a));
          [nembs, ~] = size(allEmbNames);
        
          % Finds the number of cells for the current embryo and the
          % index of the current cell. These values will be displayed on
          % the graph so that the user can tell how far along they are.
          allCells = T.cells(b);
          allCells = allCells{1, 1};
          c = strcmp(allCells, cell);
          d = find(c==1);
          ncells = T.numbCells(b);
        
         %Center the graphs on the scored frame for NEBD, but if this
         %calculation has failed then use the frame corresponding to the max
         %value for seghi per seglo
         if ~isnan(x1)
             midFrame = round(x1);
         else
             % To centre the graphs on the peak RID sum/max proj
             midFrame = round(X2(idx));
         end
        
         % Ideally the x-axis would be roughly the same time scale for all graphs
         % (i.e. the same number of time points, since spreading different
         % numbers of time points over the same length x-axis can distort the curves
         % and bias/complicate the scoring). Shows 30 minute range.
         xstart_adj = midFrame - 900;
         xstop_adj = midFrame + 900;
        
         %Plot normalized experimental data and derivative
         %Label graph with predicted NEBD
         figure1 = figure('units','normalized','outerposition',[0 0 1 1]);
         % Create axes
         % sets the x-axis tick values
         axes1 = axes('Parent',figure1,'XTick',[xstart_adj:60:xstop_adj],'XGrid','on');
         box(axes1,'on');
         hold(axes1,'all');
         
         curr1 = plot(X2,N1,'DisplayName', ' segloCirc', 'Marker','.','Color','#52127E'); %purple
         curr2 = plot(X2,D8,'DisplayName','dy/dx seghiperseglo','Marker','.','Color','b');%blue
        
         title({Celloutput(j).embName;Celloutput(j).cell},'Interpreter','none');
         legend('Location', 'southwest');
         xlabel('Time (sec)');
         axis([xstart_adj xstop_adj -1 1])
         
         % Annotate graph with embryo and cell indices
        annotation(figure1,'textbox',...
            [0.15 0.8 0.1 0.04],...
            'String',['Celloutput index ' num2str(d) '/' num2str(ncells)],...
            'FontSize',12,...
            'FitBoxToText','on');
        
        annotation(figure1,'textbox',...
            [0.15 0.85 0.1 0.04],...
            'String',['Emboutput index ' num2str(b) '/' num2str(nembs)],...
            'FontSize',12,...
            'FitBoxToText','on');
        
         xL = get(gca, 'XLim');
         yL = get(gca, 'YLim');
         % Plot a dashed line to mark y=0
         plot(xL, [0 0], 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', '--','DisplayName', 'y=0');
         
          % Annotate graph with scored NEBD if it is not NaN
          if ~isnan(x1)
            t1 = text(x1,0.8,'NEBD');
            NEBD_line = plot([x1 x1],yL, 'Color', 'c', 'LineWidth',1, 'DisplayName', 'NEBD');
          end
          %Annotate the peak of seghi per seglo
            x2=X2(idx);
            t2=text(x2,0.8,'META');
            META_line=plot([x2 x2],yL,'Color', 'r', 'Linewidth', 1, 'DisplayName', 'META');
           
         % Show most recent graph window
         shg;
         
         answer = questdlg({'Press enter to confirm and move to next cell. ';
            'To correct the scoring, press "s". ';
            'To save progress and quit, press "q". ';
          },'Correct NEBD?','OK','OK');
      
          w = NaN;
          w = waitforbuttonpress;
          while w==0
              w = waitforbuttonpress;
          end
          % w=1 if enter, s or q pressed
          value = NaN;
          value = double(get(gcf, 'CurrentCharacter'));
          % value=13 if enter, value=115 if s, value=113 if q
        
          % Add clause to deal with scenario where user pushes button other
          % than enter, s or q
        
         while value~=13 && value~=115 && value~=113
             w = NaN;
             w = waitforbuttonpress;
           while w==0
               w = waitforbuttonpress;
           end
           % w=1 if enter, s or q pressed
           value = NaN;
           value = double(get(gcf, 'CurrentCharacter'));
         end
        
         while value~=13
           if value==113
                % User clicked "q", save our structs and kill script. 
                close(figure1);
                save(CelloutputPath, 'Celloutput');
                save(EmboutputPath, 'Emboutput');
                return
           elseif value==115
                % User clicked "s", rescore plot
                answer = questdlg({'To adjust score using graph, click on plot area. ';
                    'To adjust score using movie of nucleus, press space bar. ';
                    },'Correct NEBD?','OK','OK');
                w = NaN;
                w = waitforbuttonpress;
                % If user clicks on plot, w=0
                while w==0
                    % replace predicted line with user timepoint and let
                    % user confirm scoring
         
                     [x2,y2] = ginput(1);
                    %text(x1,y1,'NEBD');
                    delete(META_line);
                    delete(t2);
                    if (x2<=xstop_adj)&&(x2>=xstart_adj)
                        META_line = plot([x2 x2],yL, 'Color', 'r', 'LineWidth',1, 'DisplayName', 'META');
                        t2 = text(x2,0.8,'META');
                    end
                    x2 = round(x2);
                
                    if (x2 == 5000) || (x2>xstop_adj) || (x2<xstart_adj)
                        x2 = NaN;
                    end
                    w = NaN;
                    w = waitforbuttonpress;
                    % If user presses return, w=1 and while loop should be
                    % exited; if user clicks again on plot, steps should be
                    % repeated allowing for rescore. If user hits space
                    % bar, w=1 and while loop also exited.
                end
                % If user pressed return, value=13 and while loop should be
                % exited; if user pressed space bar, value=32
                value = NaN;
                value = double(get(gcf, 'CurrentCharacter'));
                while value==32
                    % User pressed space bar, opens projection to score
                    % NEBD visually from movie. If NEBD cannot be scored
                    % visually, leave input as 'NaN' and press 'Ok'.
                    id = strcat(emb,'_',cell,"_");
                    boo = strcmp(combproj_IDs,id);
                    file = combproj_fileList{boo, 1};
                    [stack, img_read] = tiffread2(file); % stack is a structure with the px values for each tiff in stack.data. img_read = number of images in stack
                    % Convert stack.date into multidimensional array
                    [~, n] = size(stack);
                    % aa and bb are image dimensions (xy)
                    [aa, bb] = size(stack(1).data);
                    stack_dbl = NaN(aa,bb,n);
                    for m = 1:1:n
                        stack_dbl(:,:,m) = double(stack(m).data);
                    end
                    outstack = uint8(stack_dbl);
                    Hh = implay(outstack,4);
                    
                     % User inputs frame value corresponding to congression start
                    % and stop, as assessed from image file. If not possible to
                    % score --> enter NaN;
                    answer = questdlg({'Manually enter frame for each event. ';
                    'Enter NaN if event is not scorable. ';
                    'Leave blank to keep score on plot. ';
                    },'Correct?','OK','OK');
                    prompt = {'Enter frame for NEBD:'};
                    dlgtitle = id;
                    dims = [1 35];
                    definput = {''};
                    opts.WindowStyle = 'normal';
                    answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
                    figure(Hh.Parent) % Brings implay window to front
                    META = str2num(cell2mat(answer(1,1)));
                    
                    close(Hh);
                    figure(figure1)

                    delete(META_line);
                    delete(t2);
                    % Easier to enter values in frames rather than seconds,
                    % so we need to convert this value back into seconds
                    META = (META)*fr_Rt;%%
                    x2 = META;
                    if ~isnan(META)
                        META_line = plot([META META],yL, 'Color', 'r', 'LineWidth',1, 'DisplayName', 'META');
                        t2 = text(META,0.8,'META');
                    end
                    w = NaN;
                    w = waitforbuttonpress;
                    if w==1
                        value = NaN;
                        value = double(get(gcf, 'CurrentCharacter'));
                        % If user is satisfied and presses enter, w=1 and value=13. 
                        % If user is not satisfied and presses space bar again, while loop should repeat.
                    elseif w==0
                        % If w=0, user clicked on plot
                        while w==0
                            % replace predicted line with user timepoint and let
                            % user confirm scoring
                            [x2,y2] = ginput(1);
                            %text(x1,y1,'NEBD');
                            delete(META_line);
                            delete(t2);
                            if (x2<=xstop_adj)&&(x2>=xstart_adj)
                                META_line = plot([x2 x2],yL, 'Color', 'r', 'LineWidth',1, 'DisplayName', 'META');
                                t2 = text(x2,0.8,'META');
                            end
                            x2 = round(x2);

                            if (x2 == 5000) || (x2>xstop_adj) || (x2<xstart_adj)
                                x2 = NaN;
                            end
                            w = NaN;
                            w = waitforbuttonpress;
                            % If user presses return, w=1 and while loop should be
                            % exited; if user clicks again on plot, steps should be
                            % repeated allowing for rescore. If user hits space
                            % bar, w=1 and while loop also exited.
                        end
                        % If user pressed return, value=13 and while loop should be
                        % exited; if user pressed space bar, value=32
                        value = NaN;
                        value = double(get(gcf, 'CurrentCharacter'));
                    end
                end
           end
         end
         % User pressed enter key, so they are satisfied with NEBD timing.
        % Close the graph.
        close(figure1);
        
        Celloutput(j).METAscrd = x2;
        Celloutput(j).METAscored = datetime;
    end
end

%In the scoring step we are plotting time in seconds on the x-axis and when
%a user corrects NEBD by clicking on the plot, we round the value of x.
%This works when we are using frames as it will round NEBD to the nearest
%integar which should be the frame closest to where we lciked. However, for
%time, this rounding leads to times of NEBD that are not represented in the
%time sampling of the movie (e.g. 908 seconds when the time step was 30
%seconds). Therefore we need to correct the scored NEBD values so that they
%are equal to the nearest sampled timepoint.

for j = 1:1:kk
    METAscrd = Celloutput(j).METAscrd;
    if ~isnan(METAscrd)
        FrRt = Celloutput(j).H2BTable.FrameRate_sec_(1);
        Time_frames = Celloutput(j).H2BTable.Frame(:);
        Time_sec = Time_frames.*FrRt;
        METAscrd = repmat(METAscrd,length(Time_sec),1);
        boo = abs(Time_sec-METAscrd);
        METAcorr = Time_sec(boo==min(boo));
        Celloutput(j).METAscrd = METAcorr;
    end
end

% Save the structs
save(CelloutputPath, 'Celloutput');
save(EmboutputPath, 'Emboutput');

    
    
            
            