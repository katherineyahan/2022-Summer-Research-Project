% This macro predicts the timing of DECOND, allows the user to confirm/adjust
% this prediction according to the graph or projection, then predicts the
% timing of Prometaphase, Metaphase, Anaphase, and Exit
%This script needs to be run after S4
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
CompiledPath = strcat(convertCharsToStrings(parent), convertCharsToStrings(slash), convertCharsToStrings(strain), '_Compiled.mat');

% Checks if the Celloutput struct exists. If it does not, it checks whether the output files above exist.
% If these files exist, we load them into the workspace.
A = exist('Celloutput');
if A~=1
    if isfile(CompiledPath)
        load(CompiledPath);
    end
    if isfile(EmboutputPath)
        load(EmboutputPath, 'Emboutput');
    end
    
    if ~isfile(CompiledPath) | ~isfile(EmboutputPath)
        error('No Celloutput variable in the work space and file could not be loaded.');
    end
end

[~,kk] = size(Celloutput);
[~, ccc] = size(Emboutput);


% Checks if the 'scored' field exists in Celloutput. If it does not, then
% this field is created and filled with an empty array at every row. If it
% does exist, it checks whether there is something in this field and
% inititalizes the for loop at the first empty row

if ~isfield(Celloutput, 'DECONDscored')
    for j = 1:1:kk
        Celloutput(j).DECONDscored = [];
    end
end

% Iterate through all of the entries (all of the cells) in the Cellouput
% struct and predict DECOND for all cells based on a combination of RID, seghiperseglo, segloArea

for j=1:1:kk
        H2B_table = Celloutput(j).H2BTable;
        frames = H2B_table.Frame;
        % Skip cells where 3 frames of fewer measured
        if numel(frames) > 3 && ~isnan(Celloutput(j).NEBDscrd)
           
           
             
            fr_Rt = H2B_table.FrameRate_sec_;
            fr_Rt = fr_Rt(1);
             X1 = frames;
             
            %Calculate the frame that NEBD occurs
            NEBDfr=Celloutput(j).NEBDscrd/fr_Rt;
            Celloutput(j).NEBDfr=int8(NEBDfr);
        
             %Convert frames to time in seconds, with frame 1 = time 0
              X2 = (X1).*fr_Rt;
              X2b = NaN;
              X2b(2:max(size(X2)),1) = X2(1:max(size(X2))-1,1);
              
              
              %Find the range that DECOND could possibly occur
              Lower_bound=Celloutput(j).NEBDfr+5;
              Higher_bound=Celloutput(j).NEBDfr+13;
              
             
          
              %First, find the lowest point in H2B RID after NEBD
              Y3=Celloutput(j).H2Bmeas_bc.nucRID;
              N3=Y3./nanmax(Y3);
              [a,~]=size(N3);
              NEBD=find(X1==Celloutput(j).NEBDfr)+5;
              %Check if there are enough frames after NEBD
             if X1(a)>= Higher_bound-5                 
                 
                 %Calculate the seglo area
                  Y9 = H2B_table.segloArea;
                  N9 = Y9./max(Y9);
                  
                  %Calculate the derivative of seglo area
                  N9b = NaN;
                  N9b(2:max(size(N9)),1) = N9(1:max(size(N9))-1,1);             
                  d9 = (N9-N9b)./(X2-X2b);
                  D9 = d9./max(abs(d9));
                  
                  %Smooth the curve
                  sD9 = movmean(D9,3);
                  
                  [ss,~]=size(sD9);
                  hasDECOND=false;
                  pos_sD9=0; %the point that sD9 starts to increase above 0
                  
                 %Check if the seglo area starts to increase for 3 consecutive frames, if not, that
                 %means the imaging did not capture DECOND
                  for k=NEBD+3:1:ss
                      if sD9(k)>0 
                          if k+1<ss && k+2<=ss
                              if sD9(k+1)>0 && sD9(k+2)>0 
                                  pos_sD9=k;
                                  hasDECOND=true;
                                break
                              end
                          end
                      end
                  end
                  
                  if hasDECOND==true  
                      
                     %Find the lowest local min in H2B RID after NEBD
                     b=N3(NEBD:a,:);
                     min_N3=min(b(:));
                     %Find its corresponding frame and index
                     min_idx=find(N3==min_N3);
                     min_fr=X1(min_idx);
            
             
                    %If the frame falls between the range, and the seglo
                    %area starts to increase
                    if min_fr >= Lower_bound && min_fr <= Higher_bound 
                        if sD9(min_idx)>0
                            %Check if it's too late; has two consistent
                            %points which the sD9 is above 0 before
                            if sD9(min_idx-1)>0 && sD9(min_idx-2)>0
                                %If it has the sD9 of two frames before it
                                %are above 0, take the point that sD9
                                %starts to increase above 0 as the DECOND
                                Celloutput(j).DECONDpred=X2(pos_sD9);
                            else
                               Celloutput(j).DECONDpred=X2(min_idx);
                            end
                        end
                    %If not, take the frame before the seghiperseglo drop below
                    %0.1
                    else
                      %Calculate the seghi as a percent of seglo area
                      Y8 = H2B_table.seghiArea./H2B_table.segloArea;
                      N8 = Y8./nanmax(Y8);
                      [s,~]=size(N8);
                      NEBD_idx=find(X1==Celloutput(j).NEBDfr)+6;
                      %Find the peak of seghiperseglo after NEBD; it
                      %is usually near META
                      META_idx=find(X1==Celloutput(j).NEBDfr)+3;
                      c=N8(META_idx:META_idx+7,:);
                      max_N8=max(c(:));
                      peak=find(N8==max_N8);
                      before_idx=0;
                      %Find the point after the peak
                      for i=peak+1:1:s
                         if N8(i)<0.1
                             before_idx=i-1;
                             if sD9(before_idx)>0 %the seglo area starts to increase
                                 if sD9(before_idx-1)>0 && sD9(before_idx-2)>0
                                     Celloutput(j).DECONDpred=X2(pos_sD9);
                                 else
                                     Celloutput(j).DECONDpred=X2(before_idx);
                                 end
                             end
                            break
                         end
                      end
                   
                    end
                
                  else
                      Celloutput(j).DECONDpred=NaN;
                  end
             else
              Celloutput(j).DECONDpred=NaN;
             end
                 
        else
            Celloutput(j).DECONDpred=NaN;
              
        end
end


% Iterate through all of the entries (all of the cells) in the Cellouput
% struct and check/confirm the predicted NEBD
for j = 1:1:kk
    if isempty(Celloutput(j).DECONDscored)
        % Import H2B channel measurement values from the Celloutput struct
        H2B_table = Celloutput(j).H2BTable;
        frames = H2B_table.Frame;
        fr_Rt = H2B_table.FrameRate_sec_;
        fr_Rt = fr_Rt(1);
        
        X1 = frames;
        
        %Convert frames to time in seconds, with frame 1 = time 0
         X2 = (X1).*fr_Rt;
         X2b = NaN;
         X2b(2:max(size(X2)),1) = X2(1:max(size(X2))-1,1);
        
        Y3=Celloutput(j).H2Bmeas_bc.nucRID;
        N3=Y3./nanmax(Y3);
        
        %Calculate the seghi as a percent of seglo area
         Y8 = H2B_table.seghiArea./H2B_table.segloArea;
         N8 = Y8./nanmax(Y8);
        
        %Calculate the seglo area
         Y9 = H2B_table.segloArea;
         N9 = Y9./max(Y9);
                  
        %Calculate the derivative of seglo area
         N9b = NaN;
         N9b(2:max(size(N9)),1) = N9(1:max(size(N9))-1,1);             
         d9 = (N9-N9b)./(X2-X2b);
         D9 = d9./max(abs(d9));
                  
         %Smooth the curve
         sD9 = movmean(D9,3);         
         [ss,~]=size(sD9);
         hasDECOND=false;
         pos_sD9=0; %the point that sD9 starts to increase above 0
                  
         %Check if the seglo area starts to increase for 3 consecutive frames, if not, that
         %means the imaging did not capture DECOND
          for k=NEBD+3:1:ss
               if sD9(k)>0 
                    if k+1<ss && k+2<=ss
                          if sD9(k+1)>0 && sD9(k+2)>0 
                               pos_sD9=k;
                               hasDECOND=true;
                               break
                           end
                     end
                end
          end
           
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
             midFrame = round(X2(pos_sD9));
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
         
         curr8 = plot(X2,sD9,'DisplayName','segloArea','Marker','.','Color','b');
         curr7 = plot(X2,N8,'DisplayName', 'seghiperseglo', 'Marker','.','Color','#52127E'); %purple
         curr10=plot(X2,N3,'DisplayName','nucRID','Marker','.','Color','m');
         
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
          
          %x2=Celloutput(j).METAscrd;
          %if ~isnan(x2)
           % t2=text(x2,0.8,'META');
           % META_line=plot([x2 x2],yL,'Color', 'r', 'Linewidth', 1, 'DisplayName', 'META');
         % end
          
          x3=Celloutput(j).DECONDpred;
          if ~isnan(x3)
              t3=text(x3,0.8,'DECOND');
              DECOND_line=plot([x3 x3],yL,'Color', 'b', 'Linewidth', 1, 'DisplayName', 'DECOND');
          else 
              x3=X2(pos_sD9);
              t3=text(x3,0.8,'DECOND');
              DECOND_line=plot([x3 x3],yL,'Color', 'b', 'Linewidth', 1, 'DisplayName', 'DECOND');
          end
          
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
                save(CompiledPath, 'Celloutput');
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
         
                     [x3,y3] = ginput(1);
                    %text(x1,y1,'NEBD');
                    delete(DECOND_line);
                    delete(t3);
                    if (x3<=xstop_adj)&&(x3>=xstart_adj)
                        DECOND_line = plot([x3 x3],yL, 'Color', 'b', 'LineWidth',1, 'DisplayName', 'DECOND');
                        t3 = text(x3,0.8,'DECOND');
                    end
                    x3 = round(x3);
                
                    if (x3 == 5000) || (x3>xstop_adj) || (x3<xstart_adj)
                        x3 = NaN;
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
                    DECOND = str2num(cell2mat(answer(1,1)));
                    
                    close(Hh);
                    figure(figure1)

                    delete(DECOND_line);
                    delete(t3);
                    % Easier to enter values in frames rather than seconds,
                    % so we need to convert this value back into seconds
                    DECOND = (DECOND)*fr_Rt;%%
                    x3 = DECOND;
                    if ~isnan(DECOND)
                        DECOND_line = plot([DECOND DECOND],yL, 'Color', 'b', 'LineWidth',1, 'DisplayName', 'DECOND');
                        t3 = text(DECOND,0.8,'DECOND');
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
                            [x3,y3] = ginput(1);
                            %text(x1,y1,'NEBD');
                            delete(DECOND_line);
                            delete(t3);
                            if (x3<=xstop_adj)&&(x3>=xstart_adj)
                                DECOND_line = plot([x3 x3],yL, 'Color', 'b', 'LineWidth',1, 'DisplayName', 'DECOND');
                                t3 = text(x3,0.8,'DECOND');
                            end
                            x3 = round(x3);

                            if (x3 == 5000) || (x3>xstop_adj) || (x3<xstart_adj)
                                x3 = NaN;
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
        
        Celloutput(j).DECONDscrd = x3;
        Celloutput(j).DECONDscored = datetime;
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
    DECONDscrd = Celloutput(j).DECONDscrd;
    if ~isnan(DECONDscrd)
        FrRt = Celloutput(j).H2BTable.FrameRate_sec_(1);
        Time_frames = Celloutput(j).H2BTable.Frame(:);
        Time_sec = Time_frames.*FrRt;
        DECONDscrd = repmat(DECONDscrd,length(Time_sec),1);
        boo = abs(Time_sec-DECONDscrd);
        DECONDcorr = Time_sec(boo==min(boo));
        Celloutput(j).DECONDscrd = DECONDcorr;
    end
end

% Save the structs
save(CompiledPath, 'Celloutput');
save(EmboutputPath, 'Emboutput');

        
