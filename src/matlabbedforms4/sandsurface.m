%--------------SAND SURFACE----  
% SANDSURFACE - A subroutine of NEWDUNES that adds a sandy texture to the dune models, or a shaded BW surface if selected by the user.
% File written by Carissa Carter
% Last modified in 3/2005

switch colormode
    case {'color','c','Color','COLOR'}
  set(gcf, 'Visible','on');    %turns on/off the on-screen print of each figure
  set(axes,'Position',[.05 .02685 .9 .9463]);  %minimize white space around figure
  axis vis3d
  set(gcf, 'Renderer', 'zbuffer')

%%plot surface with underlying colormap
  Bedformplot = surf(x,y,min(z,zref));
  axis([0-CenterShift 100-CenterShift 0-CenterShift 100-CenterShift -30 60]);
  colormap(copper);
  brighten(.5);
  material dull;
  shading interp
  axis off
  hold on
  
%%adds a sand grain texture to above surface   
  load sandpic5  
  sand5=double(sandpic5)/255;
  Bedformplot = surf(x,y,min(z,zref),sand5,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
  
%%adds a bottom to the figure  
  xbottom=[1-CenterShift 100-CenterShift 100-CenterShift 1-CenterShift];
  ybottom=[1-CenterShift 1-CenterShift 100-CenterShift 100-CenterShift];
  zbottom=[-30 -30 -30 -30];
  fill3(xbottom,ybottom,zbottom,[.6 .3 0])
  
%%brighten the area
  shiny=camlight(270,30,'infinite');
  set(shiny,'Color',[0.933333 0.866667 0.509804])
   
%%individual boxes for each face, brown color
  xbox1=[1-CenterShift XEdge(1:100) 1-CenterShift];
  ybox1 = [100-CenterShift YEdge(1:100) 1-CenterShift];
  zbox1=[-30 ZBED(1:100) -30];
  fill3(xbox1,ybox1,zbox1,[.6 .3 .1],'edgecolor',[0.870588 0.721569 0.529412],'linewidth',.85)
  
  xbox2=[1-CenterShift XEdge(100:199) 100-CenterShift];
  ybox2=[1-CenterShift YEdge(100:199) 1-CenterShift];
  zbox2=[-30 ZBED(100:199) -30];
  fill3(xbox2,ybox2,zbox2,[.6 .3 0],'edgecolor',[0.72 .63 .42],'linewidth',.85)
  
  xbox3=[100-CenterShift XEdge(200:299) 100-CenterShift];
  ybox3=[1-CenterShift YEdge(200:299)  100-CenterShift];
  zbox3=[-30 ZBED(200:299) -30];
  fill3(xbox3,ybox3,zbox3,[.6 .3 0],'edgecolor',[0.870588 0.721569 0.529412],'linewidth',.85) 
  
  xbox4=[100-CenterShift XEdge(300:399) 1-CenterShift];
  ybox4=[100-CenterShift YEdge(300:399) 100-CenterShift];
  zbox4=[-30 ZBED(300:399) -30];
  fill3(xbox4,ybox4,zbox4,[.6 .3 0],'edgecolor',[0.72 .63 .42],'linewidth',.85)  
  
%%add yellowish lighting 
  shine=light('Position',[0-CenterShift,50,0],'Color',[.621569 0.42549 0.0331373],'Style','local');
  
  hold on

  
  
    case {'bw','BW','Bw','postscript','black and white'}
%SANDSURFACE_BW - A subroutine of NEWDUNES that allows for postscript printing of the
%dune models in black and white
%Last modified in 4/2005

  set(gcf, 'Visible','on');    %turns on/off the on-screen print of each figure
  set(axes,'Position',[.05 .02685 .9 .9463]);  %minimize white space around figure

  axis vis3d
  set(gcf, 'Renderer', 'painters')
 
%%plot surface with pseudo lighting
   bplotX=x;
   bplotY=y;
   bplotZ=min(z,zref);
   Bedformplot = surfl(x,y,min(z,zref),[60,100],[.1,.4,.3,100]);
   axis([0-CenterShift 100-CenterShift 0-CenterShift 100-CenterShift -30 60]);
   colormap(gray);
   shading interp
   material dull;
   axis off
   hold on
  
 %individual boxes for each face, gray color
  xbox1=[1-CenterShift XEdge(1:100) 1-CenterShift];
  ybox1 = [100-CenterShift YEdge(1:100) 1-CenterShift];
  zbox1=[-30 ZBED(1:100) -30];
  fill3(xbox1,ybox1,zbox1,[.9 .9 .9])
  
  xbox2=[1-CenterShift XEdge(100:199) 100-CenterShift];
  ybox2=[1-CenterShift YEdge(100:199) 1-CenterShift];
  zbox2=[-30 ZBED(100:199) -30];
  fill3(xbox2,ybox2,zbox2,[.7 .7 .7])
  
  hold on
  
end

