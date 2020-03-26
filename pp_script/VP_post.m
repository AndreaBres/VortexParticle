%% Creazione Video Vortex Particles
close all; clearvars
fileID = fopen('vortex_particle.in');

% upload di dati che ci interessano
n_blades = cell2mat(textscan(fileID, '%*s %*s %d',1, 'HeaderLines', 2));
n_elements = cell2mat(textscan(fileID, '%*s %*s %d',1, 'HeaderLines', 1));
n_panels = cell2mat(textscan(fileID, '%*s %*s %d',1, 'HeaderLines', 1));

%numero di nuove particelle ad ogni frame
new_p = n_blades*(n_elements+1);
% numero di frame
nf = 100;

% DA LEGGERE SOLAMENTE SE SI HA VOGLIA DI DIVENTARE SCEMI ALTRIMENTI VAI A
% LINEA 18
% create a default color map ranging from blue to light white
blue = [0, 0, 1];
white = [255, 255, 255]/255;
colors_p = [linspace(blue(1),white(1),nf)', linspace(blue(2),white(2),nf)', linspace(blue(3),white(3),nf)'];
 
%distanza in x della particella più lontana dal rotore DA MODIFICARE (O
%MEGLIO DA VERIFICARE)
last_particle = 2.14


fclose(fileID);
%%

% Creo "videoVP":
videoVP = VideoWriter('00_VP.avi');

% Apro "videoVP":
open(videoVP);

% se proprio vuoi saperlo questo serve per rendere il grafico più simile al
% vero, altrimenti sembra che la pala si deforma
figure('Position',[0 0 1080 2/(last_particle+1)*1080]);

% Ciclo sui vari Frame
for i = 1:nf
 % devo fare degli if perchè all'inizio ho solo la pala, poi i pannelli e poi le particelle   
    if i == 1
       eval(['BC = importdata(''./solution/blade_coord/coord', num2str(i), '.dat'');'])

       frame = plot3(BC(:,1), BC(:,2), BC(:,3),'k');
       axis([-1 last_particle -1 1 -1 1]);
       hold on
%        quiver(0, 0, 0, 1, 0, 0, 'r')
%        quiver(0, 0, 0, 0, 1, 0, 'b')
%        quiver(0, 0, 0, 0, 0, 1, 'g')

    elseif i == 2 && nf == 3
       eval(['BC = importdata(''./solution/blade_coord/coord', num2str(i), '.dat'');'])
       eval(['WP = importdata(''./solution/wake_panel_coord/coord', num2str(i), '.dat'');'])

       plot3(BC(:,1), BC(:,2), BC(:,3),'k');
       axis([-1 last_particle -1 1 -1 1]);
       hold on
       frame = plot3(WP(:,1), WP(:,2), WP(:,3),'b');
%        quiver(0, 0, 0, 1, 0, 0, 'r')
%        quiver(0, 0, 0, 0, 1, 0, 'b')
%        quiver(0, 0, 0, 0, 0, 1, 'g')
       
       
    elseif(i == n_panels + 2)
       % importo i dati del frame attuale
       eval(['BC = importdata(''./solution/blade_coord/coord', num2str(i), '.dat'');'])
       eval(['WP = importdata(''./solution/wake_panel_coord/coord', num2str(i), '.dat'');']) 
       eval(['P = importdata(''./solution/particle_coord/coord', num2str(i), '.dat'');'])
       
       %pale vengono plottate nere
       plot3(BC(:,1), BC(:,2), BC(:,3),'k');
       % l'asse x è settato scegliendo qual è la distanza della particella
       % più lontana dal rotore
       axis([-1 last_particle -1 1 -1 1]);
       hold on
       % pannelli vengono plottati blu
       plot3(WP(:,1), WP(:,2), WP(:,3),'b');
       
       % prurito mentale di me e lello, commentare per il quieto vivere
       % il prurito mentale ha avuto il sopravvento: le particelle vengono
       % plottate con un colore diverso in base alla loro eta'. Si va da
       % blu a bianco
       for j = 1:i - (n_panels+2)
          scatter3(P(1+(j-1)*new_p:j*new_p,1),P(1+(j-1)*new_p:j*new_p,2),P(1+(j-1)*new_p:j*new_p,3))
       end  
       
       % ultimo scatter3 in cui imposto il frame tanto per andare sul
       % sicuro
       frame = scatter3(P(end,1),P(end,2),P(end,3),'w');
          
       % disegnamo il sistema di riferimento ground   
%        quiver(0, 0, 0, 1, 0, 0, 'r')
%        quiver(0, 0, 0, 0, 1, 0, 'b')
%        quiver(0, 0, 0, 0, 0, 1, 'g')
    end
       
    currFrame = getframe(gcf);
    
    writeVideo(videoVP, currFrame);
    
    delete(frame)
    hold off
end

close(videoVP);


