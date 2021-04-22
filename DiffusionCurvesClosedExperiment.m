domains = {'flower', 'airfoil'};
modes = {'biharmonic', 'odeco', 'odeco-rotpi'};

for domaini = 1:numel(domains)
    domain = domains{domaini};
    
    switch domain
        case 'flower'
            [V,F,b,bc] = flower_closed();
        case 'airfoil'
            [V,F,b,bc] = airfoil();
    end
    for modei = 1:numel(modes)
        mode = modes{modei};
        
        rotfield = false;
        if strcmp(mode,'odeco-rotpi')
            rotfield = true;
        end
        
        switch mode
            case 'biharmonic'
                dcm = 'biharmonic';
            case {'odeco', 'odeco-rotpi'}
                dcm = 'odeco';
        end
        [us,zs] = DiffusionCurvesClosed(V, F, b, bc, dcm, rotfield);
        
        figure;
        hold on;
        for i = 1:numel(V)
            tsurf(F{i}, V{i}, 'FaceVertexCData',us{i});
            outl = outline(F{i});
            lX = [V{i}(outl(:,1),1)'; V{i}(outl(:,2),1)'];
            lY = [V{i}(outl(:,1),2)'; V{i}(outl(:,2),2)'];
            plot(lX, lY, '-k');
        end
        hold off;
        axis equal;
        shading interp;
        saveas(gcf, ['diffusioncurves-' domain '-' mode '.png']);
        
        if numel(zs)>0
            figure;
            hold on;
            for i = 1:numel(V)
                VisualizeFrameField2D(ProcessMesh2D([V{i} zeros(size(V{i},1),1)], F{i}), zs{i}, ceil(sqrt(sum(doublearea(V{i},F{i})))*50));
            end
            hold off;
            saveas(gcf, ['diffusioncurves-field-' domain '-' mode '.png']);
        end
    end
end