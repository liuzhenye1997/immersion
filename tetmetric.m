%初始化网格
function [points,faces]=tetmetric(dataname)
    addpath('./tetemtric_code')
    
    %% input: specify which model and which time to compare
    %dataname = list(2*obj-1).name(1:end-12); % folder that contains the input meshes and results
    interpmethod ='SQP';  %["ARAP" "FFMP" "GE" "ABF" "SQP"];

    basedir = '.';
    %% load input
    basedir = string(sprintf('%s/%s/', basedir,"tetinterp_meshes"));
    source_name=basedir +  dataname+"_source.mesh";
    [x1, tets] = readMESH( source_name );
    [x2,~]=readMESH( basedir + dataname+"_target.mesh" );
    % x2 = load( basedir + "target.txt" );

    %% global alignment based on anchors
    anchorId = tets(1,1:3)';
    fAlign = @(z) alignToAnchors(z, [anchorId x2(anchorId,:)]);

    %% show iterpolation results
    switch interpmethod
        case 'SQP'
            I = SQP_interp(x1, x2, tets, 1);
        case 'ARAP'
            I = arap_interp(x1, x2, tets);
        case 'GE'
            I = ShapeSpace_interp(x1, x2, tets);
        case 'FFMP'
            I = FFMP_interp(x1, x2, tets);
        otherwise
            clear I;
    end

    %% show 2 interpolation results
    n=100;
    for i=0:n
        t=i/n;
        if strcmp(interpmethod, 'ABF')
            xInterp = ABF_interp(x1, x2, tets, t);
        else
            xInterp = I.interp(t);
        end

    %     subplot(2,2,i*2); 
        %h = fDrawMesh(x1, tets); axis off; axis equal; set(h, 'cdata', genCData(x1), 'facecolor', 'flat'); colormap jet; caxis([0 3]);
        %set(h, 'vertices', fAlign(xInterp), 'cdata', genCData(xInterp), 'facecolor', 'flat'); title(sprintf('t=%g', t));

        name="./SQP/"+dataname+sprintf("_%f.mesh",t);
        save_tet(xInterp,tets,name);
        [points,faces]=Extract_surface(xInterp,tets);
        objname=sprintf("./SQP/obj/%s_%f.obj",dataname,t);
        save_obj(points,faces,objname);
    end
  
end