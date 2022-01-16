classdef Test
    %IMAGEUTILTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        %imageFilePath
        debug=true;
    end
    
    methods (Static)
        
        function  testLineDetection()
            cells=[1, 4, 5, 9];
            %cells= 1:12;
            for k=cells
                filePath = sprintf('%s/testdata/cell-u-15-%i.tif', pwd,k);
                I = imread(filePath);
                [x1, y1, x2, y2] = ImageUtil.lineDetection(I,Test.debug);
                fprintf(1,'For %ithcell, length of detected line is: %i, (%d,%d,%d,%d)\n', k, round(sqrt((x2-x1)^2 +(y2-y1)^2)), x1, y1, x2, y2);
                
                [x1, y1, x2, y2] = ImageUtil.adjustALine(I, x1-5, x2-5, y1-5, y2-5, sprintf('cell%i', k), NaN, true);
                
                [x1, y1, x2, y2] = ImageUtil.adjustALine(I, x1-5, x2+5, y1+5, y2-5,sprintf('cell%i', k), NaN, true);
                fprintf(1,'length of adjusted line is: %i,  (%d,%d,%d,%d)\n', round(sqrt((x2-x1)^2 +(y2-y1)^2)),  x1, y1, x2, y2);
                
            end
            
        end
        function  testLineIntensity()
            
            cells=[4, 5, 6, 7, 9, 12];
            cells=[4];
            % cells = [1:12];
            for k=cells
                filePath = sprintf('%s/testdata/cell-%i.tif', pwd,k);
                I = imread(filePath);
                [intensity, bg] = ImageUtil.lineIntensity(I, 40, 80, 40, 80, k,  NaN, Test.debug);
                fprintf(1,'Intensity of line is: %d, mean background intensity: %d\n', intensity, bg);
            end
            
        end
        
        
        function  testGenerateBWMask()
            fname='updated-15.tif';
            %           fname='updated-16.tif';  %this image has cells together
            %            fname='cell-9-80.tif';
            filePath = sprintf('%s/testdata/%s', pwd, fname);
            I = imread(filePath);
            
            %             gmag=imgradient(I);
            %             BW1 = edge(I,'Canny');
            %
            %             BW2 = edge(I,'Prewitt');
            %
            %             BW3 = edge(I,'Sobel');
            %             BW4 = edge(I,'Roberts');
            %             BW5 = edge(I,'log');
            %             BW6 = edge(I,'zerocross');
            %
            %
            %             ImageUtil.debugImages(filePath,2, 5, {I, gmag, BW1, BW2, BW3,BW4,BW5,BW6}, ...
            %                 {'updated image', 'gradient magnitude', 'Çanny Edge', 'Prewitt Edge','Sobel','Roberts','log','zerocross'});
            
            
            [mask, thresh] = ImageUtil.getCleanedBW(I, Test.debug, fname, NaN);
            
        end
        
        function  testLocateCenters()
            fname='updated-15.tif';
            filePath = sprintf('%s/testdata/%s', pwd, fname);
            I = imread(filePath);
            centroids= ImageUtil.locateCentroids(I, Test.debug, fname);
        end
        
        function testFindLines()
            fname='updated-15.tif';
            filePath = sprintf('%s/testdata/%s', pwd, fname);
            I = imread(filePath);
            centroids= ImageUtil.locateCentroids(I, false, fname);
            [centroids, X1, X2, Y1, Y2]= ImageUtil.findLines(I, centroids, fname);
            %          ImageUtil.syncImage(I, NaN, centroids, NaN, NaN,NaN, NaN, false);
            
            [mask, centroids1] = ImageUtil.adjustCentroids(I, centroids);
            [centroids12, X1, X2, Y1, Y2] =ImageUtil.findLines(I, centroids1,fname);
            ImageUtil.syncImage(I, fname,NaN, centroids12, X1, X2,Y1, Y2, mask, false);
            
            
            fname='updated-1.tif';
            filePath = sprintf('%s/testdata/%s', pwd, fname);
            I = imread(filePath);
            [mask, centroids2]= ImageUtil.adjustCentroids(I, centroids12);
            ImageUtil.syncImage(I, fname, NaN,centroids2, X1, X2,Y1, Y2, mask, false);
            
        end
        
        function testadjustACentroidOriginal()
        end
        
        function testadjustACentroid()
            fname='updated-17.tif';
            filePath = sprintf('%s/testdata/%s', pwd, fname);
            I = imread(filePath);
            centroids= ImageUtil.locateCentroids(I, false, fname);
            cells=[8,6,5,7,9,12,13];
            %cells=[5];
   %         cells=[13];
            I_eq = adapthisteq(I);
            thresh= graythresh(I_eq);
            
            for k = cells
                cx=centroids(k,1);
                cy=centroids(k,2);
                
                [bwmask, acx, acy]= ImageUtil.adjustACentroid(I, cx, cy, k ,Test.debug);
                fprintf(1, '%ith cell is moved from (%i,%i) to (%i,%i\n)', k, cx,cy,acx, acy);
            end
            
        end
        
        function testLinesIntensity()
            fname='original-17.tif';
            I = imread(sprintf('%s/testdata/%s', pwd, fname));
            
            I_ad=imadjust(I);
            centroids= ImageUtil.locateCentroids(I_ad, true, fname);
            
            [mask, centroids1] = ImageUtil.adjustCentroids(I_ad, centroids);
            [centroids12, X1, X2, Y1, Y2] =ImageUtil.findLines(I_ad, centroids1,fname);
            ImageUtil.syncImage(I_ad, fname,NaN, centroids12, X1, X2,Y1, Y2, mask, false);
            
            [intensities, bg]=ImageUtil.linesIntensities(I, centroids12, X1, X2,Y1,Y2, NaN);
            value=zeros(1, length(intensities));
            value(1, :)=intensities(:);
            
        end
        
        
        
        function  cropCellsToFile()
            folder='testdata/';
            try
                rmdir(folder, 's');
            catch
            end
            mkdir(folder);
            mkdir('testdata/original');
            
            rootFolder = '/Users/sijiepeng/Downloads/rm/hrm';
            imageType = '.tif';
            rootFolderList = dir([rootFolder '/*' imageType]);
            rootFolderList = rootFolderList(1:end);
            
            for fi = 1:length(rootFolderList)
                fileName = rootFolderList(fi).name;
                I = imread([rootFolder '/' fileName] );
                
                name=sprintf('original-%02i.tif',fi);
                imwrite(I, [folder 'original/' name]);
                
                I_adj=imadjust(I);
                name=sprintf('updated-%02i.tif',fi);
                imwrite(I_adj, [folder name]);
            end
            
            
            fname='updated-17.tif';
            filePath = sprintf('%s/%s/%s', pwd, folder, fname);
            I = imread(filePath);
            centroids= ImageUtil.locateCentroids(I, false, fname);
            [~, adjustedCs] = ImageUtil.adjustCentroids(I, centroids);
            
            
            for fi = 1:length(rootFolderList)
                fileName = rootFolderList(fi).name;
                
                I = imread([rootFolder '/' fileName] );
                I_adj=imadjust(I);
                
                numel = size(centroids, 1);
                % [adjustedCs, X1, X2, Y1, Y2] =ImageUtil.findLines(I, centroids,fname);
                for ci = 1:numel
                    CX=centroids(ci,1);
                    CY=centroids(ci,2);
                    [subI,~,~] = ImageUtil.cropImage(I, CX, CY, 80);
                    name=sprintf('cell-o-%i-%i-80.tif',fi,ci);
                    imwrite(subI, [folder name]);
                    
                    [subI,~,~] = ImageUtil.cropImage(I_adj, CX, CY, 80);
                    name=sprintf('cell-u-%i-%i-80.tif',fi,ci);
                    imwrite(subI, [folder name]);
                    %               imfinfo([folder name]);
                end
                 numel = size(adjustedCs, 1);
                for ci = 1:numel
                    CX=adjustedCs(ci,1);
                    CY=adjustedCs(ci,2);
                    [subI,~,~]= ImageUtil.cropImage(I, CX, CY, ImageUtil.RADIUS);
                    name=sprintf('cell-o-%i-%i.tif',fi,ci);
                    imwrite(subI, [folder name]);
                    
                    [subI,~,~]= ImageUtil.cropImage(I_adj, CX, CY, ImageUtil.RADIUS);
                    name=sprintf('cell-u-%i-%i.tif',fi,ci);
                    imwrite(subI, [folder name]);
                end
            end
        end
        
        function calcIntensity()
            %preparation
            updatedfolder='testdata/updated/';
            ImageUtil.recreateFolder(updatedfolder);
            syncFolder='Pattern Recognition/tempFold/';
            ImageUtil.recreateFolder(syncFolder);
            
            fname='updated-15.tif';
            filePath = sprintf('%s/testdata/%s', pwd, fname);
            I = imread(filePath);
            
            centroids= ImageUtil.locateCentroids(I, false, fname);
            numel = length(centroids);
            if(numel>30)
                errorMessage = sprintf('Too many cell regions (%i) detected in %s', numel, fullFileName);
                uiwait(warndlg(errorMessage));
            end
            [baseCentroids, X1, X2, Y1, Y2] =ImageUtil.findLines(I, centroids,fname);
            [~, baseCentroids] = ImageUtil.adjustCentroids(I, baseCentroids);
            [baseCentroids, X1, X2, Y1, Y2] =ImageUtil.findLines(I, baseCentroids,fname);
            ImageUtil.syncImage(I, fname,NaN, baseCentroids, X1, X2,Y1, Y2, NaN, false);
            
            
            
           % rootFolder = '/Users/sijiepeng/Downloads/rm/hrm';
           rootFolder = 'testdata/original/';
            imageType = '.tif';
            rootFolderList = dir([rootFolder '/*' imageType]);
            rootFolderList = rootFolderList(1:end);
            
            intensity=zeros(length(rootFolderList), length(baseCentroids));
            bg=zeros(length(rootFolderList), length(baseCentroids));
            for fi = 1:length(rootFolderList)
                fileName = rootFolderList(fi).name;
                currImage = imread([rootFolder '/' fileName] );
                I_adj=imadjust(currImage);
                
                %                imwrite(I, [updatedfolder sprintf('updated-%i.tif',fi);]);
                
                destPath=[syncFolder fileName];
                
                [bwmask1, centroids]  = ImageUtil.adjustCentroids( I_adj, baseCentroids);
                [bwmask, aX1, aX2, aY1, aY2] = ImageUtil.adjustLines (I_adj, centroids,X1, X2, Y1, Y2, sprintf('f%i',fi), true);
                
                if Test.debug
                    ImageUtil.syncImage(I_adj, ['adjusted--' fileName], destPath, centroids,aX1, aX2, aY1, aY2, bwmask, false);
                end
                
                [intensity(fi, :),bg(fi,:) ]=ImageUtil.linesIntensities(currImage, centroids, X1, X2, Y1, Y2, NaN);
            end
            
            cell=[1, 2, 3];
            cell =[];
            ImageUtil.graphPlot(intensity,'intensity', cell);
            ImageUtil.graphPlot(bg,'mean bg intensity',cell);
            measured = intensity+bg;
            ImageUtil.graphPlot(measured,'Measured intensity',cell);
            comp = measured./bg;
            ImageUtil.graphPlot(comp,'measured line intensity/background',cell);
            %             saveas(gcf, 'Pattern Recognition/tempFold/lineIntensityGraph');
            ImageUtil.write(intensity, 'intensity.xlsx');
            ImageUtil.write(measured, 'measured-intensity.xlsx');
            ImageUtil.write(bg, 'bgintensity.xlsx');
            ImageUtil.write(comp, 'intensity-ratio.xlsx');
        end
        
        function testRegionBasedIntensity
            x1= 30, x2=50,y1=30, y2=80;
            hw=2;
            k=(y2-y1)/(x2-x1);
            
            
            x = [x1-hw x1+hw  x2+hw x2-hw ];
            y = [y1+hw y1-hw  y2-hw y2+hw];
            
            bw = poly2mask(x,y,111,111);
            figure
            imshow(bw);
            hold on
            plot([x1,x2], [y1,y2], 'b','LineWidth',2)
            hold off
            
            
            
            L = bwlabel(bw);
            
            
            %            imshow(label2rgb(L, @jet, [.7 .7 .7]))
            %
            %            s = regionprops(L, 'PixelIdxList', 'PixelList');
            
            filePath = sprintf('%s/testdata/cell-4.tif', pwd);
            I = imread(filePath);
            
            s = regionprops(L,I,'MeanIntensity')
            
            s2 = regionprops(bw,I,'MeanIntensity')
        end
        function  testMaskLine()
            fname='cell-9-80.tif';
            filePath = sprintf('%s/testdata/%s', pwd, fname);
            I = imread(filePath);
            [x1, x2, y1, y2]=deal(10, 140, 60, 20);
            bw=ImageUtil.mask4Line(100, 200, x1, x2, y1, y2);
            figure;
            imshow(bw);
            hold on
            plot([x1,x2], [y1,y2], 'y','LineWidth',2)
            hold off
        end
        
        
        function testAdjustLine()
            fname='cell-u-9-80.tif';
            
   %         [fname, ox1, ox2, oy1, oy2]= deal('cell-u-1-3.tif', 60,81, 76, 56);
            
            
            
            %Line is too short f17-C9: old (42,79, 40, 73)L=4.957822e+01, Detected (53,75, 52, 73)L=3.041381e+01, New (52,75, 51, 73)L=3.182766e+01
            [fname, ox1, ox2, oy1, oy2]= deal('cell-u-17-9.tif', 42,79, 40, 73);
            
            
      %      [fname, ox1, ox2, oy1, oy2]= deal('cell-u-6-1.tif', 51,84, 23, 62);
            
          %  [fname, ox1, ox2, oy1, oy2]= deal('cell-u-9-9.tif', 42,79, 41, 74);
          
           
       %   [fname, ox1, ox2, oy1, oy2]= deal('cell-u-16-6.tif', 38,83, 40, 29);
          
       %    [fname, ox1, ox2, oy1, oy2]= deal('cell-u-11-3.tif', 60,81, 76, 56);
           
         %  [fname, ox1, ox2, oy1, oy2]= deal('cell-u-7-6.tif', 38,83, 40, 29);

            filePath = sprintf('%s/testdata/%s', pwd, fname);
            I = imread(filePath);
            [x1, y1, x2, y2] = ImageUtil.adjustALine(I, ox1, ox2, oy1, oy2, fname, NaN, true);
            
        end
        
        
        
    end
    
    
    
end

