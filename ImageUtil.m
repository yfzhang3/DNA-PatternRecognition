classdef ImageUtil
    %ImageUtil image processing functions that can be used and tested
    %outside of the UI app.
    %  Constant defines parameters used for cell recognization
    
    properties (Constant)
        SIZE_AOPEN =100;
        SIZE_AOPEN_LARGE=3000;
        RADIUS =55
        SIZE_OPEN =3
        NUM_AC_ITER=6
        PROPERTY_CENTROID='Centroid'
        LINE_WIDTH=1
        HUGE_BW_LEVEL=0.99
        
        % parameters to control how intensity is caculcuated
        USE_EQ_INTENSITY=false
        REMOVE_LINE4BGM=true
        USE_LINE_REGION=false
        USE_DEFAULT_BGM=true
        LINE_REGION_RADIUS=2
    end
    
    
    methods (Static)
        %line detection for individual cell
        function [x1, y1, x2, y2] = lineDetection(image, debug)
            I = image;
            I = im2gray(I(:,:,1));
            
            level = ImageUtil.HUGE_BW_LEVEL;
            BW = imbinarize(I,level);
            [H,theta,rho] = hough(BW);
            P = houghpeaks(H,1,'threshold',ceil(0.3*max(H(:))));
            x = theta(P(:,2));
            y = rho(P(:,1));
            lines = houghlines(BW,theta,rho,P,'FillGap',90,'MinLength',2.5);
            if (~isempty(lines))
                xy = [lines(1).point1; lines(1).point2];
                
                x1 = xy(1,1);
                y1 = xy(1,2);
                x2 = xy(2,1);
                y2 = xy(2,2);
            else
                x1 = 0;
                y1 = 0;
                x2 = 0;
                y2 = 0;
            end
            
            if(debug)
                label ='grayscale image and BW';
                figure('Name',label','NumberTitle','off');
                imshowpair(I, BW, 'montage');
                
                figure('Name','Line detection','NumberTitle','off');
                imshow(I), 
                hold on
                plot(x,y,'s','color','black');
                if ~isempty(lines)
                plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
                plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
                plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
                end
                hold off
            end
            
        end
        
        
        %generate a
        function [bw4, thresh] = getCleanedBW( image, debug, ftitle, thresh)
            I_eq = adapthisteq(image);
            if(isnan(thresh))
                thresh= graythresh(I_eq);
            end
            bw = imbinarize(I_eq, thresh);
            %bw = imbinarize(I_eq, 'adaptive');  //not working well for
            %cropped cell images
            
            %             I=imadjust(image);
            %             I_eq = adapthisteq(I);
            %             bw = imbinarize(I_eq, graythresh(I_eq));
            
            %             bw2 = imfill(bw,'holes');
            %             %          bw3 = imopen(bw2, strel('disk',ImageUtil.SIZE_OPEN));
            %             bw3 = imopen(bw2,  ones(ImageUtil.SIZE_OPEN,ImageUtil.SIZE_OPEN));
            
            
            %          bw3 = imopen(bw2, strel('disk',ImageUtil.SIZE_OPEN))
            bw2 = imfill(bw,'holes');
            bw3 = imopen(bw2,  ones(ImageUtil.SIZE_OPEN,ImageUtil.SIZE_OPEN));
            
            
            bw4 = bwareaopen(bw3, ImageUtil.SIZE_AOPEN);
            bw41=bw4;
            mask = bwconvhull(bw4,'objects');
            
            mask = bwareaopen(mask, ImageUtil.SIZE_AOPEN_LARGE);
            
            convex = activecontour(I_eq, mask, ImageUtil.NUM_AC_ITER, 'edge');
            bw4=convex;
            
            
            if(debug)
                bw4_p = bwperim(bw4,8);
                overlay1 = imoverlay(image, bw4_p, [.3 1 .3]);
                images={image, I_eq, bw, bw2, bw3, bw41, convex, bw4_p, overlay1};
                labels= {'Original', 'contrast enhanced', 'black&white', 'filled hole', 'opened', ...
                    'area opened',  'active contoured', 'perimeter', 'overlay' };
                
                
                ImageUtil.debugImages(ftitle, 2, 5,images, labels);
            end
        end
        
        
        
        %finds all cell centroids in image
        function cellCentroid = locateCentroids(image, debug, label)
            [bwmask, ~]= ImageUtil.getCleanedBW(image, false, label, NaN);
            
            s = regionprops(bwmask,image,{'Centroid','WeightedCentroid'});
            
            numObj= numel(s);
            weightedCellCentroid=zeros(numObj, 2);
            cellCentroid=zeros(numObj, 2);
            for k = 1 : numObj
                weightedCellCentroid(k,:) = round(s(k).WeightedCentroid);
            %    weightedCellCentroid(k,:) = [round(s(k).WeightedCentroid(1)), round(s(k).WeightedCentroid(2))];
                cellCentroid(k,:) = [round(s(k).Centroid(1)), round(s(k).Centroid(2))];
            end
            
            if(debug)
                figure('Name', ['Locate Centroids: ' label], 'NumberTitle','off');
                
                %show image with background masked
                %maskedI=image;
                %maskedI(bwmask == 0) = 0;
                %imshow(maskedI)
                
                %show image with identified region perimeter
                bw_p = bwperim(bwmask,8);
                overlay = imoverlay(image, bw_p, [.3 1 .3]);
                imshow(overlay);
                
                title('Weighted (red) and Unweighted (blue) Centroids');
                hold on
                for k = 1 : numObj
                    plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*')
                    plot(cellCentroid(k,1), cellCentroid(k,2), 'bo')
                    text(cellCentroid(k,1)-5, cellCentroid(k,2)-5, num2str(k));
                end
                hold off
            end
        end
        
        
        % find lines in image, prevCentroids: the original set of the centroids
        % fileName are only for showing detail info
        %return set of centroids of cells with lines detected in it.
        % the value of returned line position is relative to the centroids
        function [centroids, X1, X2, Y1, Y2]= findLines( image, prevCentroids, fileName)
            [numCentroid, ~] = size(prevCentroids);
            index=1;
            for ci = 1:numCentroid
                currC = prevCentroids(ci,:);
                cx=currC(1);
                cy=currC(2);
                
                currCell = ImageUtil.cropImage(image, cx, cy, ImageUtil.RADIUS);
                [x1, y1, x2, y2] = ImageUtil.lineDetection(currCell, false);
                if(x1==x2 && y1==y2)
                    fprintf(2, 'No line found in %i th cell of %s\n', ci, fileName );
                    continue;
                end
                
                centroids(index,:) = currC;
                [rCX, rCY]=ImageUtil.absoluteToRelative(cx, cy, cx, cy, ImageUtil.RADIUS);
                
                X1(index) = x1-rCX;
                Y1(index) = y1-rCY;
                X2(index) = x2-rCX;
                Y2(index) = y2-rCY;
                index = index+1;
            end
        end
        
        function [bwmask, x,y] = adjustACentroid( image, cx, cy,ci, debug)
            range=25;
            radius=ImageUtil.RADIUS+range;
            
            cImage=ImageUtil.cropImage(image, cx, cy, radius);
            [bwmask, ~]= ImageUtil.getCleanedBW(cImage, debug, sprintf('adjust centroid for %i', ci), NaN);
            
            s = regionprops(bwmask,cImage,{'Centroid', 'PixelIdxList'});
            numObj = numel(s);
            
            [rcx, rcy] = ImageUtil.absoluteToRelative(cx, cy,cx,cy, radius);
            found = false;
            for k = 1 : numObj
                newCx= round(s(k).Centroid(1));
                newCy= round(s(k).Centroid(2));
                if((numObj ==1) ||((abs(rcx-newCx) < range) && (abs(rcy-newCy)<range)))
                    [x, y] =ImageUtil.relativeToAbsolute(newCx,newCy, cx,cy, radius);
                    found=true;
                else
                    bwmask(s(k).PixelIdxList(:))=0;
                end
            end
            if(~found)
                fprintf(2, "Did not find centroid for %ith cell, please make sure the image is adjusted\n", ci);
                [x, y]= deal(cx, cy);
                debug=true;
            end
            
            
            %              bwmask=ImageUtil.cropLogicImage(bw4, rcx, rcx, ImageUtil.RADIUS);
            %              mask = bwareaopen(bwmask, ImageUtil.SIZE_AOPEN_LARGE);
            if(debug)
                %                 show cell image and old & new centeroid
                %                ImageUtil.debugImages('adjust one Centroid', 1, 2,{bw4, cImage}, {'mask','Original (red) & adjusted (blue) Centroids'});
                ImageUtil.debugImages(sprintf('adjust %i cell Centroid',ci), 1, 2,{bwmask, cImage}, {'mask','Original (red) & adjusted (blue) Centroids'});
                hold on
                plot(rcx, rcy, 'r*')
                if(found)
                    plot(newCx,newCy, 'bo')
                end
                hold off
            end
            
        end
        %baseCs is array of all centroids
        function [bwmask, adjustedCs] = adjustCentroids(image, baseCs)
            adjustedCs = zeros(size(baseCs));
            range=25;
            radius=ImageUtil.RADIUS+range;
            for ci=1: length(baseCs)
                cx= baseCs(ci, 1);
                cy= baseCs(ci, 2);
                
                [masks{ci}, adjustedCs(ci,1), adjustedCs(ci,2)] = ImageUtil.adjustACentroid(image, cx, cy, ci, false);
                
                
            end
            bwmask=ImageUtil.overlayMasks(size(image), baseCs, masks, radius);
            
        end
        
        
        function [subI, xMin, yMin] = cropImage(image, CX, CY, radius)
            
            [sizeImage1,sizeImage2] = size(image);
            xMin = max(1, CX - radius);
            xMax = min(sizeImage1, CX + radius);
            yMin = max(1, CY - radius);
            yMax = min(sizeImage2, CY + radius);
            
            subI = image(yMin:yMax, xMin:xMax);
        end
        
        
        function [subI, xMin, yMin] = cropLogicImage(bwmask, CX, CY, radius)
            
            [sizeImage1,sizeImage2] = size(image);
            xMin = max(1, CX - radius);
            xMax = min(sizeImage1, CX + radius);
            yMin = max(1, CY - radius);
            yMax = min(sizeImage2, CY + radius);
            
            subI = bwmask( xMin:xMax, yMin:yMax);
        end
        
        
        %translate relative position in smaller images to absolute position in the larger image,
        %the smaller image is cropped centered at (CX, CY) with given radius
        function [AX, AY] =relativeToAbsolute(X, Y, CX, CY, radius)
            xMin = max(1, CX - radius);
            yMin = max(1, CY - radius);
            AX= xMin+ X-1;
            AY= yMin+ Y-1;
        end
        
        %translate absolute position in larger images to relative position in the smaller image,
        %the smaller image is cropped centered at (CX, CY) with given radius
        function [RX, RY] =absoluteToRelative(X, Y, CX, CY, radius)
            xMin = max(1, CX - radius);
            yMin = max(1, CY - radius);
            RX= X-xMin+1;
            RY= Y-yMin+1;
        end
        
        %display images in m row and n columns.
        function debugImages(ftitle, m, n, images, labels)
            %displays ftitle
            figure('Name',ftitle','NumberTitle','off');
            %matlab property with displaying figures
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
            
            % Get rid of tool bar and pulldown menus that are along top of figure.
            set(gcf, 'Toolbar', 'none', 'Menu', 'none');
            
            for ii = 1:length(images)
                %organizes imgs into grid, iith img
                subplot(m, n, ii);
                imshow(images{ii});
                text(.5,0, labels{ii}, 'horiz','center',...
                    'vert','top', 'FontSize',10, 'units','normalized');
            end
        end
        
        
        
        
        %generate image with original image and identified centroids and
        %Lines
        function syncImage(image, fname, destFile, centroids, LX1, LX2,LY1, LY2, mask,closeFig)
            currFig = figure('Name', sprintf('Images with detected lines: %s', fname), 'NumberTitle','off');
            if(isnan(mask))
                %original and changed image are same 
                image_m = image;
            else
                %only cell perimeter is white, all else black, 8 is for
                %pixel connectivty
                bw4_p = bwperim(mask,8);
                %overlays green perim onto cell 
                image_m = imoverlay(image, bw4_p, [.2 0.8 .2]);
            end
            
            imshowpair(image_m, image,'montage');
            %           imshow(image);
            %ensures new plots do not delete existing plots
            hold on
            for ci = 1:length(LX1)
                currC = round(centroids(ci,:));
                %fprintf(1,'plot %ith centroid, location: %i:%i\n', ci, currC(1), currC(2));
                %gets abs line endpt coords
                ax1 = LX1(ci)+currC(1);
                ax2 = LX2(ci)+currC(1);
                ay1 = LY1(ci)+currC(2);
                ay2 = LY2(ci)+currC(2);
                
                plot(currC(1), currC(2), 'o', 'Color','red');
                text(currC(1)-5, currC(2)-5, num2str(ci));
                plot([ax1, ax2], [ay1, ay2],'LineWidth',2,'Color','green');
            end
            hold off
            if(~isnan(destFile))
                saveas(gcf, destFile);
            end
            if(closeFig)
                close(currFig);
            end
        end
        
        %graphs intensity of lines over course of images
        function [] = graphPlot(intensity, vname, includes)
            
            
            figure('Name', vname, 'NumberTitle','off');
            
            ax = axes;
            %         ax.ColorOrder = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
            ax.LineStyleOrder = {'-','--',  ':', '-.'};
            hold on
            title(sprintf('%s vs. Time',vname));
            xlabel('Image Index')
            ylabel(vname);
            
            
            maxX=size(intensity,1);
            
            num=size(intensity,2);
            index=0;
            for k=1:num
                if(isempty(includes) || any(includes(:) == k))
                    
                    index=index+1;
                    %style = styles(mod(index, 4)+1);
                    p=plot(intensity(:,k), 'DisplayName',sprintf('cell%d',k));
                    p.LineWidth = 1.5;
                    
                    
                    y=intensity(maxX,k);
                    text(maxX,y, sprintf('cell%d',k)); %draw at line end
                    %                 y=intensity(1,k);
                    %                 text(0.8,y, sprintf('c%d',k));
                end
            end
            
            str=ImageUtil.printSettings();
            dim = [.2 .5 .3 .3];
            %annotation('textbox',dim,'String',str,'FitBoxToText','on');
            
            hold off
            legend show
        end
        
        function value=printSettings()
            str1='';
            if ImageUtil.USE_EQ_INTENSITY
                str1='EQed,';
            end
            str2='';
            if ImageUtil.REMOVE_LINE4BGM
                str2 =sprintf('BGMIgnoreLineSurround(%i),',ImageUtil.LINE_REGION_RADIUS);
            end
            str3='';
            if  ImageUtil.USE_LINE_REGION
                str3= sprintf('UseLineSurroundArea(%i),',ImageUtil.LINE_REGION_RADIUS);
            end
            str4='';
            if  ImageUtil.USE_DEFAULT_BGM
                str4= sprintf('DefaultBGM');
            end
            value =strcat(str1, str2, str3,str4);
        end
        
        function [] = write(intensity, fname)
            for i = 1:size(intensity, 1)
                time(i, 1) = i;
            end
            for i = 1:size(intensity, 2)
                cellIndex(1, i) = i;
            end
            writecell({'time'}, fname, 'Range', 'A2');
            writecell({'cellIndex'}, fname, 'Range', 'B1');
            writematrix(time, fname, 'Range', 'A3');
            writematrix(cellIndex, fname, 'Range', 'C1');
            writematrix(intensity, fname, 'Range', 'C3');
            %      winopen('output.xlsx');
        end
        
        function [intensities,bgs]=linesIntensities(fulimage, centroids, X1, X2,Y1,Y2, masks)
            radius=ImageUtil.RADIUS;
            intensities=zeros(length(X1),1);
            numel = length(X1);
            bgs=zeros(numel,1);
            for k = 1:numel
                currC = centroids(k,:);
                cx = currC(1); cy=currC(2);
                currCell= ImageUtil.cropImage(fulimage, cx, cy, radius);
                [rCX, rCY]=ImageUtil.absoluteToRelative(cx, cy, cx, cy, radius);
                
                rx1 = X1(k) +rCX;
                ry1 = Y1(k) +rCY;
                rx2 = X2(k) +rCX;
                ry2 = Y2(k) +rCY;
                
                if(isnan(masks))
                    mask= NaN;
                else
                    mask = masks{k};
                end
                [intensities(k), bgs(k)]= ImageUtil.lineIntensity(currCell,rx1,rx2,ry1,ry2,k,mask, false);
            end
            
        end
        %finds intensity of the line in a cell
        %I=cropped cell img
        function [intensity,meanBackgroundIntensity] = lineIntensity(I, x1,x2,y1,y2,ci, mask, debug)
            
            if ImageUtil.USE_DEFAULT_BGM
                %only uses mid 50% background values
                meanBackgroundIntensity= ImageUtil.calcBgMeanIntensity(I, x1,x2,y1,y2,ci, mask, debug);
            else
                %uses all background values
                meanBackgroundIntensity= ImageUtil.calcBgMeanIntensity2(I, x1,x2,y1,y2,ci, mask, debug);
            end
            %increases contrast in images using CLAHE
            if ImageUtil.USE_EQ_INTENSITY
                I= adapthisteq(I);
            end
            
            if ImageUtil.USE_LINE_REGION
                %calcs intensity in rect region surrounding line
                lintensity = ImageUtil.calcIntensity2(I, x1, x2, y1, y2) ;
            else
                lintensity = mean(improfile(I, [x1 x2], [y1 y2], 200));
            end
            %take diff of two intensities because only taking intensity of
            %line can mess w/ data if line has bright background
            %taking diff allows to compare more effectively with other
            %lines
            intensity = lintensity- meanBackgroundIntensity ;
            
        end
        
        function value=calcBgMeanIntensity2(I, x1,x2,y1,y2,ci, mask, debug)
            if(isnan(mask))
                I_adj = imadjust(I);
                [mask, ~] = ImageUtil.getCleanedBW(I_adj, debug, sprintf( 'get BW mask for %ith cell', ci), NaN);
            end
            %makes line area in mask all black so bgrnd intensity doesnt
            %include line intensity
            if(ImageUtil.REMOVE_LINE4BGM)
                [sizex, sizey]=size(mask);
                bw=ImageUtil.mask4Line(sizex, sizey, x1, x2, y1, y2);
                mask(bw)=0;
            end
            
            if(ImageUtil.USE_EQ_INTENSITY)
                I= adapthisteq(I);
            end
            s=regionprops(mask,I,'MeanIntensity');
            value =s.MeanIntensity;
        end
        
        function value=calcBgMeanIntensity(I, x1,x2,y1,y2,ci, mask, debug)
            if(isnan(mask))
                I_adj = imadjust(I);
                [mask1, ~] = ImageUtil.getCleanedBW(I_adj, debug, sprintf( 'get BW mask for %ith cell', ci), NaN);
            end
            
            if(ImageUtil.REMOVE_LINE4BGM)
                [sizex, sizey]=size(mask1);
                bw=ImageUtil.mask4Line(sizex, sizey, x1, x2, y1, y2);
                mask=mask1;
                mask(bw==1)= 0;
                if debug
                    imshowpair(mask, I, 'Montage');
                end
            else
                mask=mask1;
            end
            
            %use original image, not the adjusted for intensity calculation
            if(ImageUtil.USE_EQ_INTENSITY)
                I_eq= adapthisteq(I);
            else
                I_eq=I;
            end
            
            numPixels = 0;
            background = [];
            for i = 1 : size(I_eq, 1)
                for j = 1 : size(I_eq, 2)
                    if mask(i, j) == 1
                        numPixels = numPixels + 1;
                        background(numPixels) = int32(I_eq(i, j));
                    end
                end
            end
            
            
            meanBackgroundIntensity =0;
            if(numPixels >0)
                backgroundIntensity = 0;
                background = sort(background);
                IQR = background(round(numPixels*3/4))-background(round(numPixels/4));
                count = 0;
                for i = 1 : numPixels%removes outlier values in cell, counts remaining and averages them
                    if ((background(round(numPixels/2))-(IQR*3/2) <= background(i)) && (background(round(numPixels/2))+(IQR*3/2) >= background(i)))
                        backgroundIntensity = backgroundIntensity + background(i);
                        count = count + 1;
                    end
                end
                meanBackgroundIntensity=backgroundIntensity/count;
            end
            value=meanBackgroundIntensity;
        end
        
        function value=calcIntensity2(I, x1, x2, y1, y2)
            [sizex, sizey]=size(I);
            bw=ImageUtil.mask4Line(sizex, sizey, x1, x2, y1, y2);
            s=regionprops(bw,I,'MeanIntensity');
            value =s.MeanIntensity;
        end
        
        %makes bw mask of area surrounding line
        function bw=mask4Line(sizex, sizey, x1, x2, y1, y2)
            margin=ImageUtil.LINE_REGION_RADIUS/2;
            
            len = sqrt((x2-x1)^2 +(y2-y1)^2);
            a=margin*(x2-x1)/len;
            b=margin*(y2-y1)/len;
            
            
            X=  round([x1-a-b x1-a+b  x2+a+b x2+a-b ]);
            Y = round([y1+a-b y1-a-b  y2-a+b y2+a+b]);
            
            bw = poly2mask(X,Y,sizex,sizey);
        end
        
        function recreateFolder(folder)
            try
                rmdir(folder, 's');
            catch
            end
            mkdir(folder);
        end
        
        function [bwmask, X1, X2, Y1, Y2] = adjustLines (I, centroids,X1, X2, Y1, Y2, fname, debug)
            [numCentroid, ~] = size(centroids);
            masks=cell(numCentroid);
            for ci = 1:numCentroid
                currC = centroids(ci,:);
                cx=currC(1);
                cy=currC(2);
                
                [rcx, rcy]=ImageUtil.absoluteToRelative(cx, cy, cx, cy, ImageUtil.RADIUS);
                [ox1, ox2, oy1, oy2] = deal(X1(ci)+rcx, X2(ci)+rcx, Y1(ci)+rcy, Y2(ci)+rcy);
                
                currCell = ImageUtil.cropImage(I, cx, cy, ImageUtil.RADIUS);
                [masks{ci}, ~] = ImageUtil.getCleanedBW(currCell, false, sprintf( 'get BW mask for %ith cell', ci), NaN);

                [x1, y1, x2, y2] = ImageUtil.adjustALine(currCell, ox1, ox2, oy1, oy2, sprintf('%s-C%i',fname, ci), masks{ci}, debug);
                if(x1==x2 && y1==y2)
                    fprintf(2, 'No adjusted line found in %i th cell of %s\n', ci, fname );
                    continue;
                end
                
                X1(ci) = x1-rcx;
                Y1(ci) = y1-rcy;
                X2(ci) = x2-rcx;
                Y2(ci) = y2-rcy;
                
                %  I = imadjust(currCell);
            end
            
            bwmask=ImageUtil.overlayMasks(size(I), centroids, masks, ImageUtil.RADIUS);
            
        end
        
        function [x1, y1, x2, y2] = adjustALine(I, ox1, ox2, oy1, oy2, cname, mask, debug)
            if ~isnan(mask)
               I(~mask)=0;
            end
            
            [nx1, ny1, nx2, ny2] = ImageUtil.lineDetection(I, false );
            
            if(ox1 >ox2)
                [ox1, ox2, oy1, oy2]= deal(ox2, ox1, oy2, oy1);
            end
            if(nx1 >nx2)
                [nx1, nx2, ny1, ny2]= deal(nx2, nx1, ny2, ny1);
            end
            [dx1, dy1] = deal(ox2-ox1, oy2-oy1);
            len1= norm([dx1, dy1],2);
            
            if(nx1==nx2 && ny1==ny2) %no line detected, do nothing
                [x1, x2, y1, y2]= deal(ox1, ox2, oy1, oy2);
                if debug
                    fprintf( 1, 'AdjustLine, no line detected for %s: old (%d,%d, %d, %d)L=%d \n', ...
                        cname, ox1, ox2, oy1, oy2, len1);
                end
                
                return;
            end
            
            dif=[ox1, ox2, oy1, oy2] -[nx1, nx2, ny1, ny2];
            
            [dx2, dy2] = deal( nx2-nx1, ny2-ny1);
            
            len2=norm([dx2, dy2],2);
            
            difnorm = norm(dif, 2);
            if ( difnorm <3)  %small movement,just use the new line
                [x1,  x2,y1, y2]= deal(nx1, nx2, ny1, ny2);
                if debug
                    fprintf( 1, 'AdjustLine, very small movement, use new line for %s: old (%d,%d, %d, %d)L=%d, New (%d,%d, %d, %d)L=%d\n', ...
                        cname, ox1, ox2, oy1, oy2, len1, nx1, nx2, ny1, ny2, len2);
                end
                return;
            end
            
            
            %find line cross.
            %      [xcs,ycs] = polyxpoly([ox1,ox2], [oy1, oy2], [nx1,nx2],[ny1, ny2]);
            
            warning off MATLAB:singularMatrix
            
            A=[-dy1, dx1; -dy2, dx2];
            B=[dx1*oy1-dy1*ox1; dx2*ny1-dy2*nx1];
            X=linsolve(A, B);
            xc = X(1);
            yc = X(2);
            
            m1= dy1/dx1;
            m2= dy2/dx2;
            
            c1= oy1-m1*ox1;
            c2= ny1-m2*nx1;
            
            if (isnan(xc) || isnan(yc) || abs(xc)==Inf  || abs(yc)==Inf) %parallel
                %                 [x1, y1, x2, y2]= deal(nx1, ny1, nx2, ny2);
                %
                %
                %
                %
                %
                %             if debug
                %                 fprintf( 2, 'Equation has no solution for %s, use new line: old (%d,%d, %d, %d)L=%d, Detected (%d,%d, %d, %d)L=%d\n', ...
                %                     cname, ox1, ox2, oy1, oy2, len1, nx1, nx2, ny1, ny2, len2);
                %             end
                %                 %y = mx+c, distance = (c2-c1)/sqrt(1+m^2), m2 and m1 should
                

                
                ldist= abs(c2-c1)/sqrt(1+m2^2);
                
                
                if(abs(ldist) >5)  % probably a different line, just ignore and use old line.
                    [x1,  x2,y1, y2]= deal(ox1,   ox2, oy1, oy2);
                    if debug
                        fprintf( 2, 'AdjustLine, too big parallel difference , use old line for %s: old (%d,%d, %d, %d)L=%d, New (%d,%d, %d, %d)L=%d\n', ...
                            cname, ox1, ox2, oy1, oy2, len1, nx1, nx2, ny1, ny2, len2);
                    end
                else
                    lcx= ldist*dy2/len2;
                    lcy = ldist*dx2/len2;
                    x1= ox1+lcx;
                    x2= ox2+lcx;
                    y1= oy1+lcy;
                    y2= oy2+lcy;
                    if debug
                        fprintf( 2, 'AdjustLine, parallel movement for %s: old (%d,%d, %d, %d)L=%d, New (%d,%d, %d, %d)L=%d\n', ...
                            cname, ox1, ox2, oy1, oy2, len1, nx1, nx2, ny1, ny2, len2);
                    end
                    
                    
                end
                return;
                
            end
            
            if difnorm >20 && abs(m1-m2) >0.2 || len2< 0.3*len1
                [x1,  x2,y1, y2]= deal(ox1, ox2, oy1, oy2);
                if debug
                    fprintf( 1, 'Too big difference %s, use old line: old (%d,%d, %d, %d)L=%d, Detected (%d,%d, %d, %d)L=%d\n', ...
                        cname, ox1, ox2, oy1, oy2, len1, nx1, nx2, ny1, ny2, len2);
                end
                return;
            end
            
            A2= [dx1/len1,dy1/len1]'*[dx2/len2, dy2/len2];
            
            p1=round([ox1-xc, oy1-yc]*A2 +[xc, yc]);
            p2=round([ox2-xc, oy2-yc] *A2+[xc, yc]);
            [x1, x2, y1, y2]= deal(p1(1), p2(1), p1(2), p2(2) );
            
            %
            %             p1=[ox1-xc, oy1-yc]*A2 +[xc, yc];
            %             p2=[ox2-xc, oy2-yc] *A2+[xc, yc];
            
            %             [x1, x2]=deal(floor(p1(1)), ceil(p2(1)));
            %
            %             if ny2 >ny1
            %               [y1, y2]= deal(floor(p1(2)), ceil(p2(2)));
            %             else
            %               [y1, y2]= deal(ceil(p1(2)), floor(p2(2)));
            %             end
            %
            
            
            
            if(debug)
                figure('Name', sprintf('Adjust Line %s New (Red) vs ol line (Green) vs Detected(yellow)', cname), 'NumberTitle','off');
                imshow(I);
                hold on
                plot(xc,yc,'s','color','black');
                
                plot([ox1, ox2],[oy1, oy2], 'LineWidth',2,'Color','green');
                plot([nx1, nx2],[ny1, ny2], 'LineWidth',2,'Color','yellow');
                plot([x1, x2],[y1, y2], 'LineWidth',2,'Color','red');
 
                hold off
                newNorm = sqrt((x2-x1)^2 +(y2-y1)^2);
                
                if(newNorm < len1-2)
                    fprintf(2, 'Line is too short %s: old (%d,%d, %d, %d)L=%d, Detected (%d,%d, %d, %d)L=%d, New (%d,%d, %d, %d)L=%d\n', ...
                        cname, ox1, ox2, oy1, oy2, len1, nx1, nx2, ny1, ny2, len2, x1, x2, y1, y2, newNorm);
                else
                    if(newNorm > len1+2)
                        fprintf(2, 'Line is too long %s: old (%d,%d, %d, %d)L=%d, Detected (%d,%d, %d, %d)L=%d, New (%d,%d, %d, %d)L=%d\n', ...
                            cname, ox1, ox2, oy1, oy2, len1, nx1, nx2, ny1, ny2, len2, x1, x2, y1, y2, newNorm);
                        
                    else
                        fprintf( 1, 'Line changed %s: old (%d,%d, %d, %d)L=%d, Detected (%d,%d, %d, %d)L=%d, New (%d,%d, %d, %d)L=%d\n', ...
                            cname, ox1, ox2, oy1, oy2, len1, nx1, nx2, ny1, ny2, len2, x1, x2, y1, y2, newNorm);
                    end
                    
                end
            end
        end
        
        function bigMask = overlayMasks(sizeI, centroids, bwmasks, radius)
            bigMask = zeros(sizeI);
            numel = length(bwmasks);
            %     figure;
            for k=1:numel
                
                xMin = max(1, centroids(k,1) - radius);
                yMin = max(1, centroids(k,2) - radius);
                [sizex, sizey]= size(bwmasks{k});
                for i= 1:sizex
                    for j=1:sizey
                        if bwmasks{k}(i,j)
                            [xb, yb]= deal(xMin+j-1,yMin+i-1);
                            bigMask(yb,xb ) = 1;
                        end
                    end
                end
                
                %           imshowpair(bwmasks{k}, bigMask, 'Montage');
            end
            
            %          imshow(bigMask);
            
        end
        
        
        
    end
end

