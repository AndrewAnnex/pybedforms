% NEWDUNES - A script for simulating bedforms and crossbedding.  This
% script outputs image files for use as individual frames of Quicktime
% movies or printing.
% File written by Dave Rubin and modified by Carissa Carter
% Last modified in 4/2005

clear all
clf


%%******** NOVICE: INPUT FILE NAMES PROMPT ******************************%%
HowManyFiles=input('How many files do you want to run?  :');
for i=1:HowManyFiles
    FileList(i,:)=cellstr(input('Name of input parameter file? (Example: for fig16.m enter fig16)  :','s'))  ;
end
FilesToRun=char(FileList);
%%***********************************************************************%%


%%******** ADVANCED: MANUALLY SET INPUT FILES, SEPARATED BY COMMAS ******%%
% FilesToRun=char('fig5','fig55','fig46n')
%%***********************************************************************%%

colormode=input('Enter "color" for movies/tiffs or "bw" for single image postscripts :','s');

moviemaker;
FilesToRun
NumberOfFiles=length(FilesToRun(:,1));

for i=1:NumberOfFiles;
    CurrentFig=FilesToRun(i,:);
for FrameNumber = 1:NumberOfFrames   % Finish loop once for each movie frame.


%%******** ADVANCED: USE TO START A MOVIE IN THE MIDDLE  *****************%    
% if FrameNumber == 34  %Set starting framenumber and uncomment this line*
%                        %You must also uncomment line # 124              *
%%*************************************************************************


DuneInit
  ZHORIZ = ZHORIZ + dZHO(FrameNumber);
  TRENDF = TRENDF + dTrend(FrameNumber);
  TRENDS = TRENDS + dTrend(FrameNumber);
  TRENDT = TRENDT + dTrend(FrameNumber);
  
  for N = 1 : NBEDSH
    TIME = 1 + dT(FrameNumber) - N;
    if MORVRT == false; break; end
    if FirstRun
      DuneTopo

    sandsurface;   

      if max(max(z)) > zref
         MORHRZ = true;
      end
      
      zcont = z;
      FirstRun = false;

    end
    if MORHRZ
      x = XSurf;
      y = YSurf;
      DuneTopo;
      zcont = min(z,zcont);
      if mod(TIME,INTXBD) == 0;
        [c,h]=contour3(x,y,zcont,[zref zref]);                  %surface contours
        set(h,'linewidth', 1.25,'edgecolor',[0.81 0.63 0.42])    %surface contour colors  - change linewidth here!! 1.25 is thick line
      end
      pause (0.1)
      if max(max(zcont)) <= zref; MORHRZ = false; end
    end
    if MORVRT
      x = XEdge;
      y = YEdge;
      DuneTopo;

      if mod(TIME,INTXBD) == 0;
       
%%******** FOR MOVIE FRAME / COLOR OUTPUT  ******************************%%
switch colormode
    case {'color','c','Color','COLOR'}
         plot3(XEdge(300:399),YEdge(300:399),ZBED(300:399),'Color', [0.72 .63 .42],'Linewidth',1.25) %adjust colors and linewidth
         plot3(XEdge(200:299),YEdge(200:299),ZBED(200:299),'Color', [0.72 .63 .42],'Linewidth',1.25) %adjust colors and linewidth
         plot3(XEdge(100:199),YEdge(100:199),ZBED(100:199),'Color', [0.72 .63 .42],'Linewidth',1.25) %darker face contours and linewidth
         plot3(XEdge(1:100),YEdge(1:100),ZBED(1:100),'Color', [0.870588 0.721569 0.529412],'Linewidth',1.25) %lighter face contours and linewidth
%%***********************************************************************%%

%%******** FOR POSTSCRIPT B/W OUTPUT  ***********************************%%
    case {'bw','BW','Bw','postscript','black and white'}
        plot3(XEdge(100:199),YEdge(100:199),ZBED(100:199),'Color', [0.1 .1 .1])
        plot3(XEdge(1:100),YEdge(1:100),ZBED(1:100),'Color', [0.1 .1 .1])       
%%***********************************************************************%%
end
      end
      
      if max(ZBED) <= -30; MORVRT = false; end
   end

  end  % End of loop to finish each image. 
  hold off  
  
  ImageName = sprintf('%s%05.0f', FILENM, FrameNumber)
  
%%******** FOR MOVIE FRAME / COLOR OUTPUT *******************************%%
switch colormode
    case {'color','c','Color','COLOR'}
        print('-dtiff','-r250',ImageName)  
        clf                                
%%***********************************************************************%%

%%******** FOR POSTSCRIPT B/W OUTPUT  ***********************************%%
    case {'bw','BW','Bw','postscript','black and white'}
        set(gca,'drawmode','fast')  
        saveas(gcf,ImageName,'ai')  
        if NumberOfFiles > 1
            clf
        end
%%***********************************************************************%%
end

%  end   %% IF STARTING A MOVIE IN THE MIDDLE UNCOMMENT THIS LINE 

end  % End of loop to make movie sequence.
end  % End of loop to run multiple files.

