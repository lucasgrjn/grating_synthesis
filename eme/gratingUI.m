function gratingUI(obj)

    
    dx = obj.domain.discretization(1);
    dz = obj.domain.discretization(2);
    zf = obj.domain.domain(2);
    xf = obj.domain.domain(1);
    
    lambda1 = obj.domain.wavelengthSpectrum(1);
    lambda2 = obj.domain.wavelengthSpectrum(2);
    dlambda = obj.domain.wavelengthSpectrum(3);
    
    centerLambdaIndex = round(length(obj.domain.k0)/2);
    
    % Create a figure and axes
    h = figure('Position', [300 -250 850 650]);
    
    set(gcf, 'Color', 'w');
    
    
    imagesc(dz:dz:zf,dx:dx:xf,abs(obj.diel)); set(gca,'Ydir', 'normal'); colorbar;
    xlabel('z (\mum)','FontSize',14); ylabel('x (\mum)','FontSize',14); set(gca,'FontSize',12);
    title('Index Distribution', 'FontSize', 14);
    
    
    
   % Create push button
    popup = uicontrol('Style', 'popup',...
           'String', {'Index Distribution','Coupling vs Z-offset and Angle', 'Field (real) at Center Wavelength',...
           'Coupling Vs Wavelength', 'Reflection Vs Wavelength', 'Field (amplitude) at Center Wavelength'},...
           'FontSize', 12, ...
           'Position', [20 600 250 50],...
           'Callback', @what2plot);       

    
    function what2plot(source,~)
        if source.Value == 1
           imagesc(dz:dz:zf,dx:dx:xf,abs(obj.diel)); set(gca,'Ydir', 'normal'); colorbar;
           xlabel('z (\mum)','FontSize',14); ylabel('x (\mum)','FontSize',14); set(gca,'FontSize',12);
           title('Index Distribution', 'FontSize', 14);
           colormap('default');
        end
        
        if source.Value == 2
            imagesc(obj.fiberCoup.angleVec,obj.fiberCoup.zOffset,obj.fiberCoup.coup(:,:,centerLambdaIndex)); colorbar;
            xlabel('Angle (degrees)','FontSize',14); ylabel('z offset','FontSize',14); set(gca,'FontSize',12);
            title('Coupling Efficiency', 'FontSize', 14);
            colormap('default');
        end
        
        if source.Value == 3
            if obj.domain.polarization == 0 
                imagesc(dz:dz:zf, dx:dx:xf, real(obj.fullFields.Ey(:,:,centerLambdaIndex))); 
                set(gca, 'Ydir', 'normal'); colormap('redbluedark'); 
                xlabel('z (\mum)','FontSize',14); ylabel('x (\mum)','FontSize',14); set(gca,'FontSize',12);
                title('Field (real) at Center Wavelength', 'FontSize', 14);
                maxx = max(max(abs(real(obj.fullFields.Ey(:,:,1)))));
                caxis([-maxx maxx]); 
            else
                imagesc(dz:dz:zf, dx:dx:xf, real(obj.fullFields.Hy(:,:,centerLambdaIndex))); 
                set(gca, 'Ydir', 'normal'); colormap('redbluedark'); 
                maxx = max(max(abs(real(obj.fullFields.Hy(:,:,1)))));
                caxis([-maxx maxx]); 
            end
        end
        
        if source.Value == 4
            
            if lambda1 ~= lambda2
            [tempmax, I1] = max(obj.fiberCoup.coup(:,:,centerLambdaIndex)); [maxCoup, I2] = max(tempmax);
            optZOffset = obj.fiberCoup.zOffset(I1(I2));
            optAngle = obj.fiberCoup.angleVec(I2);

            %Get the efficiency at that z-offset and angle for each wavelength
            for jj = 1:length(lambda1:dlambda:lambda2)
                coupVLambda(jj)= obj.fiberCoup.coup(I1(I2),I2,jj);
            end

            plot(lambda1:dlambda:lambda2,coupVLambda, 'LineWidth',1.5); xlabel('Wavelength (\mum)','FontSize', 14);
            ylabel('Coupling Efficiency','FontSize', 14);  axis([lambda1 lambda2 0 1]);
            set(gcf,'Color', 'w'); set(gca,'FontSize',12); title('Grating Response','FontSize',14); grid on;
            else
                
            fprintf('\n You only simulated one wavelength \n');
            end
        end
        
         if source.Value == 5
            
            if lambda1 ~= lambda2
            [tempmax, I1] = max(obj.fiberCoup.coup(:,:,centerLambdaIndex)); [maxCoup, I2] = max(tempmax);
            optZOffset = obj.fiberCoup.zOffset(I1(I2));
            optAngle = obj.fiberCoup.angleVec(I2);

            %Get the efficiency at that z-offset and angle for each wavelength
            for jj = 1:length(lambda1:dlambda:lambda2)
                reflVLambda(jj)= obj.scatterProperties.PowerRefl(jj,1);
            end

            plot(lambda1:dlambda:lambda2,reflVLambda*10^2, 'LineWidth',1.5); xlabel('Wavelength (\mum)','FontSize', 14);
            ylabel('Reflection % into Waveguide','FontSize', 14);  axis([lambda1 lambda2 0 ceil(max(reflVLambda*10^2))]);
            set(gcf,'Color', 'w'); set(gca,'FontSize',12); title('Reflection','FontSize',14); grid on;
            else
                
            fprintf('\n You only simulated one wavelength \n');
            end
         end
        
         if source.Value == 6
            if obj.domain.polarization == 0 
                imagesc(dz:dz:zf, dx:dx:xf, abs(obj.fullFields.Ey(:,:,centerLambdaIndex))); 
                set(gca, 'Ydir', 'normal'); colormap('redbluedark'); 
                xlabel('z (\mum)','FontSize',14); ylabel('x (\mum)','FontSize',14); set(gca,'FontSize',12);
                title('Field (amplitude) at Center Wavelength', 'FontSize', 14);
                maxx = max(max(abs(real(obj.fullFields.Ey(:,:,1)))));
                caxis([-maxx maxx]); 
            else
                imagesc(dz:dz:zf, dx:dx:xf, abs(obj.fullFields.Hy(:,:,centerLambdaIndex))); 
                set(gca, 'Ydir', 'normal'); colormap('redbluedark'); 
                maxx = max(max(abs(real(obj.fullFields.Hy(:,:,1)))));
                caxis([-maxx maxx]); 
            end
        end
        
    end
    


    if verLessThan('matlab', '8.4')
        
        fprintf('\n GUI does not work with versions of matlab before 2014b \n');
        figure;
        imagesc(dz:dz:zf,dx:dx:xf,obj.diel); set(gca,'Ydir', 'normal'); colorbar;
           xlabel('z (\mum)','FontSize',14); ylabel('x (\mum)','FontSize',14); set(gca,'FontSize',12);
           title('Index Distribution', 'FontSize', 14);
           colormap('default');
        figure;
        imagesc(obj.fiberCoup.angleVec,obj.fiberCoup.zOffset,obj.fiberCoup.coup(:,:,centerLambdaIndex)); colorbar;
            xlabel('Angle (degrees)','FontSize',14); ylabel('z offset','FontSize',14); set(gca,'FontSize',12);
            title('Coupling Efficiency', 'FontSize', 14);
            colormap('default');
        figure;
        if obj.domain.polarization == 0 
                imagesc(dz:dz:zf, dx:dx:xf, real(obj.fullFields.Ey(:,:,centerLambdaIndex))); 
                set(gca, 'Ydir', 'normal'); colormap('redbluedark'); 
                xlabel('z (\mum)','FontSize',14); ylabel('x (\mum)','FontSize',14); set(gca,'FontSize',12);
                title('Field at Center Wavelength', 'FontSize', 14);
                maxx = max(max(abs(real(obj.fullFields.Ey(:,:,1)))));
                caxis([-maxx maxx]); 
            else
                imagesc(dz:dz:zf, dx:dx:xf, real(obj.fullFields.Hy(:,:,centerLambdaIndex))); 
                set(gca, 'Ydir', 'normal'); colormap('redbluedark'); 
                maxx = max(max(abs(real(obj.fullFields.Hy(:,:,1)))));
                caxis([-maxx maxx]); 
         end
        figure;
         if lambda1 ~= lambda2
            [tempmax, I1] = max(obj.fiberCoup.coup(:,:,centerLambdaIndex)); [maxCoup, I2] = max(tempmax);
            optZOffset = obj.fiberCoup.zOffset(I1(I2));
            optAngle = obj.fiberCoup.angleVec(I2);

            %Get the efficiency at that z-offset and angle for each wavelength
            for jj = 1:length(lambda1:dlambda:lambda2)
                coupVLambda(jj)= obj.fiberCoup.coup(I1(I2),I2,jj);
            end

            plot(lambda1:dlambda:lambda2,coupVLambda, 'LineWidth',1.5); xlabel('Wavelength (\mum)','FontSize', 14);
            ylabel('Coupling Efficiency','FontSize', 14);  axis([lambda1 lambda2 0 1]);
            set(gcf,'Color', 'w'); set(gca,'FontSize',12); title('Grating Response','FontSize',14); grid on;
            else
                
            fprintf('\n You only simulated one wavelength \n');
          end
        figure;
         if lambda1 ~= lambda2
            [tempmax, I1] = max(obj.fiberCoup.coup(:,:,centerLambdaIndex)); [maxCoup, I2] = max(tempmax);
            optZOffset = obj.fiberCoup.zOffset(I1(I2));
            optAngle = obj.fiberCoup.angleVec(I2);

            %Get the efficiency at that z-offset and angle for each wavelength
            for jj = 1:length(lambda1:dlambda:lambda2)
                reflVLambda(jj)= obj.scatterProperties.PowerRefl(jj,1);
            end

            plot(lambda1:dlambda:lambda2,reflVLambda*10^2, 'LineWidth',1.5); xlabel('Wavelength (\mum)','FontSize', 14);
            ylabel('Reflection % into Waveguide','FontSize', 14);  axis([lambda1 lambda2 0 ceil(max(reflVLambda*10^2))]);
            set(gcf,'Color', 'w'); set(gca,'FontSize',12); title('Reflection','FontSize',14); grid on;
            else
                
            fprintf('\n You only simulated one wavelength \n');
         end
    end
        
end