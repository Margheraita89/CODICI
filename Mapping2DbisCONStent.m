clear all
close all
clc

% Lumen
path = '/Users/macbook/Documents/CODICI/SectionsTotali/';
sections_number = 49;
num_pti_sez = 102;

sections = LoadSections(path,sections_number,num_pti_sez);

%Outer
path = '/Users/macbook/Documents/CODICI/SectionsTotaliOuter/';
sections_number = 49;
num_pti_sez = 102;

sections_out = LoadSections(path,sections_number,num_pti_sez);

% Stent sottoforma di centerline
filename = 'FFCent2.dat';
delimiterIn = ' ';
headerlinesIn = 1;
cd('/Users/macbook/Documents/CODICI/')
A = importdata(filename,delimiterIn,headerlinesIn);
stent = [A.data(:,1),A.data(:,2),A.data(:,3)];

%sections_out = LoadSections(path,sections_number,num_pti_sez);
%counter = 1;

cd(path);
str = '-*';
str_stent = 'o';
counter = 1

for i = 1:size(sections,3)
    plot_3D(sections(:,:,i), str, counter, 'black')
    hold on;
    data = sections(:,:,i);
    text(data(1,1), data(1,2), data(1,3),num2str(i),'FontSize',14)
    plot_3D(stent, str_stent, counter, 'red')
end
counter = counter + 1;

% Da modificare all'occorrenza dato lo stent
% idx_stent = [10,38];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

knots = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14];
coefs =  [ -120.7  -50.0 0 1.0;
       -50.0 -120.7 0 1.0;
       50.0 -120.7 0 1.0;
       120.7  -50.0 0 1.0;
       120.7   50.0 0 1.0;
       50.0  120.7 0 1.0;
       -50.0  120.7 0 1.0;
       -120.7   50.0 0 1.0;
       -120.7  -50.0 0 1.0;
       -50.0 -120.7 0 1.0;
       50.0 -120.7 0 1.0];

curve  = nrbmak(coefs',knots);

figure(counter)
nrbplot(curve,100);
hold on;
nrbctrlplot(curve);
hold on;
for i = 1:curve.number(1)
   s = num2str(i);
   text(curve.coefs(1,i),curve.coefs(2,i),curve.coefs(3,i),s,'FontSize',14);
end
counter = counter + 1;

curve_cl = FromUntoCl(curve);
figure(counter)
nrbplot(curve_cl,100);
hold on;
nrbctrlplot(curve_cl);
hold on;
for i = 1:curve_cl.number(1)
   s = num2str(i);
   text(curve_cl.coefs(1,i),curve_cl.coefs(2,i),curve_cl.coefs(3,i),s,'FontSize',14);
end
counter = counter + 1;

height = 210;
% Estrudo la curva in direzione longitudinale per il lumen
srf = nrbextrude(curve,[0.0 0.0 height]);
% Raffino il template in 3 diverse 

nx = 6;
ny_1 = 8;
ny_2 = 27;
ny_3 = 6;
ref = 'k';
e = [0,1];
srf_1 = nurb_refinement_surf(srf,e,[nx,ny_1],ref,1);
srf_2 = nurb_refinement_surf(srf,e,[nx,ny_2],ref,1);
srf_3 = nurb_refinement_surf(srf,e,[nx,ny_3],ref,1);

% CreateFEAPInputFile(nurbs)

[srf_1_in,counter] = DoMappingUnclamped(srf_1,sections(:,:,1:10),srf_1,counter);
[srf_2_in,counter] = DoMappingUnclamped(srf_2,sections(:,:,10:38),srf_2,counter);
[srf_3_in,counter] = DoMappingUnclamped(srf_3,sections(:,:,38:49),srf_3,counter);

[CPs1,CPs_w1,Ws1]=GetCPs(srf_1_in);
[CPs2,CPs_w2,Ws2]=GetCPs(srf_2_in);
[CPs3,CPs_w3,Ws3]=GetCPs(srf_3_in);

counter = counter+1;

figure(counter)
str = '-*';
for i = 1:size(CPs1,3)
    plot_3D(CPs1(:,:,i), str, counter, 'red');
    text(CPs1(1,1,i),CPs1(1,2,i),CPs1(1,3,i),num2str(i),'FontSize',14);
    hold on;
end
for i = 1:size(CPs2,3)
    plot_3D(CPs2(:,:,i), str, counter, 'blue');
    text(CPs2(1,1,i),CPs2(1,2,i),CPs2(1,3,i),num2str(i),'FontSize',14);
    hold on;
end
for i = 1:size(CPs3,3)
    plot_3D(CPs3(:,:,i), str, counter, 'black');
    text(CPs3(1,1,i),CPs3(1,2,i),CPs3(1,3,i),num2str(i),'FontSize',14);
end
counter = counter + 1;

srf_in{1} = srf_1_in;
srf_in{2} = srf_2_in;
srf_in{3} = srf_3_in;

figure(counter)
nrbplot(srf_1_in,[100,100])
hold on;
nrbplot(srf_2_in,[100,100])
hold on;
nrbplot(srf_3_in,[100,100])
hold on;
counter = counter + 1;

% Esporto in paraview
path_paraview = '/Users/macbook/Documents/CODICI/';
ExportParaview(srf_1_in, 'Srf_1', counter, path, path_paraview)
ExportParaview(srf_2_in, 'Srf_2', counter, path, path_paraview)
ExportParaview(srf_3_in, 'Srf_3', counter, path, path_paraview)

[srf_1_out,counter] = DoMappingUnclamped(srf_1_in,sections_out(:,:,1:10),srf_1,counter);
[srf_2_out,counter] = DoMappingUnclamped(srf_2_in,sections_out(:,:,10:38),srf_2,counter);
[srf_3_out,counter] = DoMappingUnclamped(srf_3_in,sections_out(:,:,38:49),srf_3,counter);

[CPs1o,CPs_w1,Ws1]=GetCPs(srf_1_out);
[CPs2o,CPs_w2,Ws2]=GetCPs(srf_2_out);
[CPs3o,CPs_w3,Ws3]=GetCPs(srf_3_out);

counter = counter+1;

figure(counter)
str = '-*';
for i = 1:size(CPs1o,3)
    plot_3D(CPs1o(:,:,i), str, counter, 'red');
    text(CPs1o(1,1,i),CPs1o(1,2,i),CPs1o(1,3,i),num2str(i),'FontSize',14);
    hold on;
end
for i = 1:size(CPs2o,3)
    plot_3D(CPs2o(:,:,i), str, counter, 'blue');
    text(CPs2o(1,1,i),CPs2o(1,2,i),CPs2o(1,3,i),num2str(i),'FontSize',14);
    hold on;
end
for i = 1:size(CPs3o,3)
    plot_3D(CPs3o(:,:,i), str, counter, 'black');
    text(CPs3o(1,1,i),CPs3o(1,2,i),CPs3o(1,3,i),num2str(i),'FontSize',14);
end
counter = counter + 1;

figure(counter)
nrbplot(srf_1_out,[100,100])
hold on;
nrbplot(srf_2_out,[100,100])
hold on;
nrbplot(srf_3_out,[100,100])
hold on;
counter = counter + 1;

%Esporto in paraview
path_paraview = '/Users/macbook/Documents/CODICI/';
ExportParaview(srf_1_out, 'Srf_1_out', counter, path, path_paraview)
ExportParaview(srf_2_out, 'Srf_2_out', counter, path, path_paraview)
ExportParaview(srf_3_out, 'Srf_3_out', counter, path, path_paraview)

srf_out{1} = srf_1_out;
srf_out{2} = srf_2_out;
srf_out{3} = srf_3_out;

% Creo il volume
VolOut{1} = CreoVolume(srf_1_in, srf_1_out);
VolOut{2} = CreoVolume(srf_2_in, srf_2_out);
VolOut{3} = CreoVolume(srf_3_in, srf_3_out);

%% Export the volume in Paraview
% r = nrbeval(VolOut, {linspace(VolOut.knots{1}(VolOut.order(1)),VolOut.knots{1}(end-VolOut.order(1)+1),100) ...
%                      linspace(VolOut.knots{2}(VolOut.order(2)),VolOut.knots{2}(end-VolOut.order(2)+1),100)...
%                      linspace(VolOut.knots{3}(VolOut.order(3)),VolOut.knots{3}(end-VolOut.order(3)+1),20)});
% u=zeros(100000000,1);
% u=zeros(100000000,1);
% mkdir(path_paraview,'RISULTATI');
% path1 = strcat(path_paraview,'RISULTATI');
% cd(path1);
% msh_to_vtk_mod(r,u,'Vol_UNCLAMPED_ROSSI','u');
% cd(strcat(path,'..'))
% counter = counter + 1;
% cd(strcat(path,'..'));
                 
% MainPrestressUnclamped(VolOut, srf_in);

% [C,Cxi,Cet,Cze] = bezierExtraction3D(VolOut.knots{1},VolOut.knots{2},VolOut.knots{3},VolOut.order(1)-1,VolOut.order(2)-1,VolOut.order(3)-1,1,0,0);

cd('/Users/macbook/Documents/CODICI/');
feap_name = CreateFEAPInputFile(srf_in,VolOut);

% RIcavo la thickness del vaso CPs per CPs
% fileID = fopen('thickness.txt','w');

% cd('/home/margherita/Documents/UNCLAMPED/VORP/GRAZIOSI/');
% 
% for i = 1:size(srf_in.coefs,3)
%     sin = srf_in.coefs(:,:,i);
%     sout = srf_out.coefs(:,:,i);
%     sin = sin';
%     sout = sout';
%     for j = 1:size(sin,1)
%         dist(i,j) = sqrt((sout(j,1) - sin(j,1))^2 + (sout(j,2) - sin(j,2))^2 + (sout(j,3) - sin(j,3))^2);
%         fprintf(fileID,'%6.2f %6.2f %6.2f %6.2f\n',sin(i,1),sin(i,2),sin(i,3),dist(i,j));
%     end
%     fileID = RicavoThickness(fileID,sin,sout);
% end

function [approx,counter] = DoMappingUnclamped(nurbs,sections,srf,counter)
    % Mapping 2D
    ncp_x = size(nurbs.coefs,2);
    ncp_y = size(nurbs.coefs,3);
    ncp = ncp_x*ncp_y;
    npt_x = size(sections,1);
    npt_y = size(sections,3);

    X_1 = [];
    Y_1 = [];
    Z_1 = [];
    idx = 0;

    for j = 1:npt_x
        for i = 1:npt_y
            X_1 = [X_1; sections(j,1,i)];
            Y_1 = [Y_1; sections(j,2,i)];
            Z_1 = [Z_1; sections(j,3,i)];
            idx = idx + 1;
        end
    end

    X = reshape(X_1,npt_y,npt_x);
    Y = reshape(Y_1,npt_y,npt_x);
    Z = reshape(Z_1,npt_y,npt_x);

    N = npt_x*npt_y;

    p = nurbs.order(1);
    q = nurbs.order(2);
    uknots = nurbs.knots{1};
    vknots = nurbs.knots{2};
    unique_v = unique(vknots);

    % unique_u = unique(uknots);
    % Nu = spcol(uknots, p, linspace(unique_u(1),unique_u(end),npt_x));
    Nu = spcol(uknots, p, linspace(nurbs.knots{1}(nurbs.order(1)),nurbs.knots{1}(length(nurbs.knots{1})-nurbs.order(1)+1),npt_x));
    Nv = spcol(vknots, q, linspace(unique_v(1),unique_v(end),npt_y))';
    % Compute basis functions using another function to check them
    % u = linspace (unique_u(1),unique_u(end),npt_x);  
    % u = linspace(nurbs.knots(nurbs.order),nurbs.knots(length(nurbs.knots)-nurbs.order+1),npt_x)
    % s = findspan (ncp-1, p-1, u, uknots);  
    % B = basisfun (s, u, p-1, uknots);

    v = zeros(1,ncp);
    C = zeros(N,ncp);
    A = zeros(ncp_x,ncp_y);

    k = 0;
    k_stent_iniziale = 0;
    k_stent_finale = 0;
    C1 = [];
    for j = 1:npt_y
        for i = 1:npt_x
            k = k + 1;
            A = Nu(i,:)'*Nv(:,j)';
            v = reshape(A,1,ncp);
            C(k,:) = v;
        end
    end

    x = reshape(X',N,1);
    y = reshape(Y',N,1);
    z = reshape(Z',N,1);

    id = 1:ncp_x:size(C,2);
    id_f = id(2:end)-1;
    id_f = [id_f,size(C,2)];

    % Prima sezione
    % C1 = C(:,1:ncp_x);
    % % Seconda sezione
    % C2 = C(:,ncp_x+1:end);
    % % Sommo le colonne i cui CPs coincidono
    % C1(:,1:3) = C1(:,1:3) + C1(:,9:11);
    % C2(:,1:3) = C2(:,1:3) + C2(:,9:11);
    % % Shrink delle matrici
    % C1 = C1(:,1:ncp_x-nurbs.order(1)+1);
    % C2 = C2(:,1:ncp_x-nurbs.order(1)+1);
    % C3 = [C1,C2];

    for i = 1:size(id,2)
        tmp = C(:,id(i):id_f(i));
        tmp(:,1:3) = tmp(:,1:3) + tmp(:,end-nurbs.order(1)+2:end);
        tmp = tmp(:,1:ncp_x-nurbs.order(1)+1);
        C1 = [C1, tmp];
    end

    B = C1\[x,y,z]; 

    %idx_coeff_stent = find(not((C1(k_stent_iniziale:k_stent_finale,:)\[x(k_stent_iniziale:k_stent_finale,:),y(k_stent_iniziale:k_stent_finale,:),z(k_stent_iniziale:k_stent_finale,:)])==0));
    
    Bx = reshape(B(:,1),[],srf.number(2));
    Bx = [Bx; Bx(1:nurbs.order(1)-1,:)]';

    By = reshape(B(:,2),[],srf.number(2));
    By = [By; By(1:nurbs.order(1)-1,:)]';

    Bz = reshape(B(:,3),[],srf.number(2));
    Bz = [Bz; Bz(1:nurbs.order(1)-1,:)]';

    approx.form = 'B-NURBS';
    approx.number = [srf.number(1) srf.number(2)];
    for i = 1:srf.number(2)
        approx.coefs(1,:,i) = Bx(i,:);
        approx.coefs(2,:,i) = By(i,:);
        approx.coefs(3,:,i) = Bz(i,:);
        approx.coefs(4,:,i) = ones(1,size(Bx,2));
    end
    approx.knots = {srf.knots{1},srf.knots{2}};
    approx.order = [srf.order(1) srf.order(2)];

    figure(counter)
    nrbplot(approx,[100,100]);
    counter = counter + 1;

    figure(counter)
    nrbctrlplot(approx)
    counter = counter +1;
    
    % Compute error: the difference between the observed response value yi and the fitted response value ŷi, 
    % and is identified as the error associated with the data.
    Rx = x - C1*B(:,1);
    [valx,idx] = max(Rx);
    Ry = y - C1*B(:,2);
    [valy,idy] = max(Ry);
    Rz = z - C1*B(:,3);
    [valz,idz] = max(Rz);
    
    figure(counter)
    quiver3(x,y,z,Rx,Ry,Rz,'LineWidth',2,'Color',[0 0 0]);
    counter = counter + 1;
end

function feap_name = CreateFEAPInputFile(NurbsIn, VolOut)
    % Elimina file vecchia ottimizzazione
    s = system('rm *.vts');
    % s = system('rm Imesh*');
    s = system('rm Omesh*');
    s = system('rm Lmesh*');
    s = system('rm Mmesh*');
    s = system('rm IXmesh*');
    s = system('rm feapname');
    %s = system('rm TIED.txt');

    %% Creo file FEAP
    mat_name='Isotropic';
    mat_prop=[3 0.33];
    mat_num = [1,2,1];

    mat_name_stent='Neohooke';
    mat_prop = [3 0.3]
    
    quadrature = [2,2,2];
    % quadrature = [4,3,3];
    feap_name=strcat('IPSD');

    fid = fopen(feap_name,'wt');
    feap_header(fid,mat_name,mat_prop,quadrature);
    
    fid1 = fopen('Coordinates.txt','wt');

    counter_knold=0;
    counter_cpold=1;
    counter_sideold=1;
    cp_end=0;
    counter_bc=1;
    idx_stent = [10,38];

    for i=1:size(VolOut,2)
        
        % KNOT VECTORS
        counter_kn=feap_knotsun(fid,VolOut{i},counter_knold);
        counter_knold=counter_kn;

        % CONTROL POINTS
        counter_cp=feap_cps(fid,VolOut{i},counter_cpold);
        counter_cpold=counter_cp;

        % SIDES
        %[csd,counter_blocks]=feap_sides_stent(fid,VolOut(i),counter_kn,counter_sideold,cp_end,idx_stent);
        csd=feap_sides(fid,VolOut{i},counter_kn,counter_sideold,cp_end);
        counter_sideold=csd(end)+1;
        n_thick=size(VolOut{i}.coefs,4);
        cp_end=(counter_sideold-1)*n_thick;
        
        % BLOCK
        feap_block(fid,counter_kn,csd,mat_num(i));  
        counter_bc=counter_cp+1;  
    end

    for i = 1:size(NurbsIn,2)
        % MESH di SUPERFICIE INTERNA
        % KNOT VECTORS
        counter_kn=feap_knotssun(fid,NurbsIn{i},counter_knold);
        counter_knold=counter_kn;

        %CONTROL POINTS
        counter_cp=feap_cpss(fid,NurbsIn{i},counter_cpold);
        counter_cpold = counter_cp;

        %SIDES
        csd=feap_sidess(fid,NurbsIn{i},counter_knold,counter_sideold,cp_end);
        counter_sideold=csd(end)+1;
        cp_end=counter_cpold-1;

        %BLOCK
        feap_blocks(fid,csd);
        counter_sideold=counter_sideold+1;
        counter_bc=counter_cp+1; 
    end

    feap_closure(fid)

end

function [fileID] = RicavoThickness(fileID, sin, sout)

%     h = figure;
%     w = figure;
%     y = figure;
%     z = figure;
    
    sin = sin';
    sout = sout';
    
    sin = sin(:,1:3);
    sout = sout(:,1:3);
    
    % Cerco il baricentro delle due sezioni
    ptin = zeros(size(sin,1),3);
    ptout = zeros(size(sout,1),3);

    bin = trova_bar(sin);
    bout = trova_bar(sout);
    
    % Traslo i punti nel baricentro
    ptin(:,1) = sin(:,1)-bin(1);
    ptin(:,2) = sin(:,2)-bin(2);
    ptin(:,3) = sin(:,3)-bin(3);
    
    ptout(:,1) = sout(:,1)-bout(1);
    ptout(:,2) = sout(:,2)-bout(2);
    ptout(:,3) = sout(:,3)-bout(3);
         
%     figure(h)
%     plot3(squeeze(sin(:,1)),squeeze(sin(:,2)),squeeze(sin(:,3)),'*');
%     hold on;
%     plot3(ptin(:,1),ptin(:,2),ptin(:,3),'o');
%     hold on;
%     title('PUNTI TRASLATI NEL BARICENTRO (S\_IN)');
%     grid on;
%         
%     figure(w)
%     plot3(squeeze(sout(:,1)),squeeze(sout(:,2)),squeeze(sout(:,3)),'*');
%     hold on;
%     plot3(ptout(:,1),ptout(:,2),ptout(:,3),'o');
%     hold on;
%     title('PUNTI TRASLATI NEL BARICENTRO (S\_OUT)');
%     grid on;
    
    % Mi metto in un sistema di riferimento locale
    ptin_1 = [ptin(1,1), ptin(1,2), ptin(1,3)];
    ptin_2 = [ptin(5,1), ptin(5,2), ptin(5,3)];
    ptin_1=ptin_1./norm(ptin_1);
    ptin_2=ptin_2./norm(ptin_2);
    
    ptout_1 = [ptout(1,1), ptout(1,2), ptout(1,3)];
    ptout_2 = [ptout(5,1), ptout(5,2), ptout(5,3)];
    ptout_1=ptout_1./norm(ptout_1);
    ptout_2=ptout_2./norm(ptout_2);

    % Plotto i punti che prendo
%     figure(y)
%     plot3(ptin_1(:,1),ptin_1(:,2),ptin_1(:,3),'ro');
%     plot3(ptin_2(:,1),ptin_2(:,2),ptin_2(:,3),'bo');
%     hold on;
%     plot3(ptout_1(:,1),ptout_1(:,2),ptout_1(:,3),'ro');
%     plot3(ptout_2(:,1),ptout_2(:,2),ptout_2(:,3),'bo');

    % Calcolo la normale dei punti sul piano (srf_in)
    n = cross(ptin_1, ptin_2);
    % Calcolo i versori nella base originale ed in quella traslata
    E1 = [1,0,0];
    E2 = [0,1,0];
    E3 = [0,0,1];
    e1 = ptin_1;
    e3 = n./norm(n);
    e2 = cross(e1,e3);
    % Riempio la matrice di rotazione
    Rot1 = [dot(E1,e1), dot(E2,e1), dot(E3,e1);
         dot(E1,e2), dot(E2,e2), dot(E3,e2);
         dot(E1,e3), dot(E2,e3), dot(E3,e3)];
    
    % Ruoto tutti i punti
    for j = 1:size(ptin,1)
       ptinr(j,:) = Rot1*ptin(j,:)';
    end
    
    % Calcolo la normale dei punti sul piano (srf_out)
    n = cross(ptout_1, ptout_2);
    % Calcolo i versori nella base originale ed in quella traslata
    E1 = [1,0,0];
    E2 = [0,1,0];
    E3 = [0,0,1];
    e1 = ptout_1;
    e3 = n./norm(n);
    e2 = cross(e1,e3);
    % Riempio la matrice di rotazione
    Rot2 = [dot(E1,e1), dot(E2,e1), dot(E3,e1);
         dot(E1,e2), dot(E2,e2), dot(E3,e2);
         dot(E1,e3), dot(E2,e3), dot(E3,e3)];
    
    % Ruoto tutti i punti
    for j = 1:size(ptout,1)
       ptoutr(j,:) = Rot2*ptout(j,:)';
    end
    
    for i = 1:size(ptinr,1)
        dist(i) = sqrt((ptoutr(i,1) - ptinr(i,1))^2 + (ptoutr(i,2) - ptinr(i,2))^2 + (ptoutr(i,3) - ptinr(i,3))^2);
        fprintf(fileID,'%6.2f %6.2f %6.2f %6.2f\n',sin(i,1),sin(i,2),sin(i,3),dist(i));
    end
    
end
    
function [] = ExportParaview(srf, name, counter, path, path_paraview)
    figure(counter)
    q=nrbplot_mod_unclamped(srf,[500,500]);
    u=zeros(1000000,1);
    u=zeros(1000000,1);
    mkdir(path_paraview,'RISULTATI');
    path1 = strcat(path_paraview,'RISULTATI');
    cd(path1);
    msh_to_vtk_mod(q,u,name,'u');
    cd(strcat(path,'..'))
    counter = counter + 1;
    mkdir(path_paraview,'RISULTATI_RHINO');
    %path_rhino = strcat(path_paraview,'RISULTATI_RHINO');
    %cd(path_rhino)
    %igesout(srf_in,'PsdSrf');
    cd(strcat(path,'..'));
end 

function [VolOut] = CreoVolume(srf_in,srf_out)
    VolOut.form='B-NURBS';
    VolOut.dim=4;
    VolOut.number=[srf_in.number(1) srf_in.number(2) 2];
    VolOut.coefs(:,:,:,1)=srf_in.coefs;
    VolOut.coefs(:,:,:,2)=srf_out.coefs;
    % Creo lo spessore: metto una lineare tra le due superfici (nel caso non usi vettori unclamped)
    VolOut.knots={srf_in.knots{1} srf_in.knots{2} [0 0 1 1]};
    % Order nel caso in cui metta una lineare e poi faccio la degree elevation
    VolOut.order=[srf_in.order(1) srf_in.order(2) 2];
    ref = 'k';
    e = [0 0 1];
    sub = [0 0 1];
    VolOut = nurb_refinement_vol(VolOut,e,sub,ref);
end
