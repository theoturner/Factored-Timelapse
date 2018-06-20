function factoredTimeLapseVideo
    
    % Due to the way this algorithm samples data (as in Lawrence et al.,
    % 2006), it can on occasion construct a slightly discoloured image. If
    % the result is unsatisfactory, please run the program again.
    
    
    close all;
    
    
    % The percentage of the time that a given pixel is in the shade. This
    % varies depending on the timelapse used, adjust it accordingly. 
    percentageInShadow = 50;
    
    
    
    % KEY STEP 1: READ TIMELAPSE AND DISTINGUISH GROUD PIXELS FROM SKY
    % PIXELS.
    
    % Mask the sky, leaving only ground pixels. Variable name pun not
    % intended.
    groundTruthImage = imread('timelapse/groundTruth.png');
    [height, width] = size(groundTruthImage);
    frameDimension = height * width;

    groundTruth = find(groundTruthImage)';
    numberIlluminated = size(groundTruth, 2);

    
    % Read in timelapse frames.
    timelapseVideo = VideoReader('timelapse/video.mp4');
    length = round(timelapseVideo.Duration * timelapseVideo.FrameRate);
    
    timelapseFrames = zeros(length, frameDimension, 3);
    
    % ReadFrame doesn't allow specific frame access without manual frame
    % index computation from the elapsed time and framerate, hence we have
    % to use an equivalent while loop.
    count = 1;
    while(hasFrame(timelapseVideo))
        
        thisFrame = readFrame(timelapseVideo);
        timelapseFrames(count, :, :) = reshape(thisFrame, [frameDimension, 3]);
        
        count = count + 1;
        
    end
    
    
    % The spatio-temporal volume F(t) = Isky(t) + Ssun(t) * Isun(t). This
    % is the volume we can factor into its original components.
    spatioTemporalVolume = timelapseFrames(:, groundTruth, :);

    clear timelapseFrames;
    
    
    
    % KEY STEPS 2-4: FACTORISE LIGHT - SEE factorComponents() FOR EACH
    % STEP. THESE STEPS ARE DONE INDEPENDENTY ON EACH COLOUR CHANNEL, I.E.
    % KEY STEPS 2-4 ARE REPEATED FOR EACH OF THE RED, GREEN AND BLUE
    % CHANNELS.
    
    % All of the per-channel shift maps, skylight and sunlight basis curves
    % are unnecessary for our visualisation. If you'd like to observe these
    % factors of the spatio-temporal volume, just replace the '~'
    % characters with your chosen variable names. Shift maps, skylight
    % basis curves and sunlight basis curves are the last three items in
    % the array of outputs (in order).
    [redShadowVolume, redSkyVolume, redSunVolume, redSkylight, redSunlight, ... 
        redShift, ~, ~] = factorComponents(spatioTemporalVolume(:, :, 1), percentageInShadow, length, numberIlluminated);
    
    [greenShadowVolume, greenSkyVolume, greenSunVolume, greenSkylight, greenSunlight, ...
        ~, ~, ~] = factorComponents(spatioTemporalVolume(:, :, 2), percentageInShadow, length, numberIlluminated);
    
    [blueShadowVolume, blueSkyVolume, blueSunVolume, blueSkylight, blueSunlight, ...
        ~, ~, ~] = factorComponents(spatioTemporalVolume(:, :, 3), percentageInShadow, length, numberIlluminated);
    
    
    % KEY STEP 5: VISUALISE RESULTS.
    
    % Simply plot the above decomposed factors over the course of the
    % frames if you want to observe their behaviour over time:
    % plot(1 : length, [factor]);
    
    
    % Plot shift map. Note the values are normalised to better distinguish
    % shifts by eye.
    redShift = imcomplement(redShift);
    [smallestShift, largestShift] = deal(min(redShift(:)), max(redShift(:)));
    baselineShift = redShift - smallestShift;
    redShift = baselineShift ./ (largestShift - smallestShift);
    shiftMap = zeros(height, width);
    shiftMap(groundTruth') = redShift(:);
      
    figure;
    imshow(shiftMap);
    title('Shift map (normalised)');
    
    clear shift shiftMap;
    
    
    % Take shadow map, apply a bilateral filter and visualise.
    inverseShadowFunction = zeros(height, width, length);
    inverseShadowFrame = groundTruthImage;
    
    for count = 1 : length
        
        % Only need to use a single shadow volume for visualisation.
        inverseShadowFrame(groundTruth') = redShadowVolume(count, :);
        inverseShadowFunction(:, :, count) = bilateralFilter(inverseShadowFrame, height, width, 1);
        
    end
    
    implay(inverseShadowFunction);
    
    clear inverseShadowFunction inverseShadowFrame;
    
    
    skylightImage = reconstructImage(groundTruth', redSkylight, greenSkylight, blueSkylight, height, width);   
    figure;
    imshow(skylightImage);
    title('Skylight image');
    
    sunlightImage = reconstructImage(groundTruth', redSunlight, greenSunlight, blueSunlight, height, width);
    figure;
    imshow(sunlightImage);
    title('Sunlight image');
    
    
    % Reconstructed volumes F(r) = Isky(r) + Ssun(r) * Isun(r).
    redVolume = redSkyVolume + redShadowVolume .* redSunVolume;
    greenVolume = greenSkyVolume + greenShadowVolume .* greenSunVolume;
    blueVolume = blueSkyVolume + blueShadowVolume .* blueSunVolume;
    
    fullVolume = zeros(length, numberIlluminated, 3);
    [fullVolume(:, :, 1), fullVolume(:, :, 2), fullVolume(:, :, 3)] = deal(redVolume, greenVolume, blueVolume);
    
    clear redVolume greenVolume blueVolume;
    
    reconstructedSpatioTemporalVolume = reconstructVolume(fullVolume, groundTruthImage, groundTruth, height, width, length, frameDimension);
    
    clear fullVolume;
    
    
    % Compute and visualise image reconstruction error over the course of
    % the timelapse
    fullSpatioTemporalVolume = zeros(length, frameDimension, 3);
    fullSpatioTemporalVolume(:, (frameDimension - numberIlluminated + 1) : frameDimension, :) = spatioTemporalVolume;
    
    frameRMSE = zeros(length, 1);
    for count = 1 : length
        
        pixelDifferences = abs(reconstructedSpatioTemporalVolume(count, :, 1) - fullSpatioTemporalVolume(count, :, 1)).^2;
        frameRMSE(count,1) = sqrt(sum(pixelDifferences(:)) / frameDimension);
        
    end
    
    figure;
    plot(1 : length, frameRMSE, 'r');
    xlabel('time');
    ylabel('RMS error');
    title('Image reconstruction error over time');
    
    clear reconstructedSpatioTemporalVolume spatioTemporalVolume;

end



function [shadowVolume, skyVolume, sunVolume, skylight, sunlight, ...
    shift, skylightBasis, sunlightBasis] = factorComponents(spatioTemporalVolume, percentageInShadow, length, numberIlluminated)

    % KEY STEP 2: DETECT SHADOWS.
    shadowVolume = zeros(length, numberIlluminated);

    % Stop the shadow map jumping around
    stableVolume = zeros(length, numberIlluminated);   
    for count = 1 : numberIlluminated
        
        thisVolume = spatioTemporalVolume(:, count);
        stableVolume(:, count) = smooth(thisVolume, 20);
        
    end
    
    % The shadow map is 0 at each intensity less than the minimum deemed
    % not 'in the shade' and 1 otherwise (minimum as per section 4 of the
    % Sunkavalli et al. 2007 publication)
    fullShadowBinary = sort(spatioTemporalVolume, 1);
    shadowMap = fullShadowBinary(1 : round((percentageInShadow / 100) * length), :);
    perPixelMaximumIntensity = median(shadowMap, 1) .* (3 / 2); 
    for count = 1 : length
        
        for illuminatedCount = 1 : numberIlluminated
            
            if stableVolume(count, illuminatedCount) > perPixelMaximumIntensity(illuminatedCount)
            
                shadowVolume(count, illuminatedCount) = 1;
                
            end
            
        end
        
    end
    
    
    % KEY STEP 3: EXTRACT SKYLIGHT AND SKYLIGHT BASIS CURVES USING
    % ALTERNATING CONSTRAINED LEAST SQUARES ON GROUPS IN THE FULL
    % SPATIO-TEMPORAL LIGHT VOLUME.
    [lightSlice, basisSlice] = generate1DSlices(spatioTemporalVolume, length, numberIlluminated);
    
    [skylight, skylightBasis, ~, ...
        ~] = leastSquaresSky(spatioTemporalVolume, 1 - shadowVolume, lightSlice, basisSlice, length, numberIlluminated);
    
    skyVolume = skylight * skylightBasis;
    skyVolume = skyVolume';


    % KEY STEP 4: EXTRACT SUNLIGHT AND SUNLIGHT BASIS CURVES USING
    % ALTERNATING CONSTRAINED LEAST SQUARES ON GROUPS IN THE SUNLIGHT
    % COMPONENT OF THE SPATIO-TEMPORAL LIGHT VOLUME AND APPLYING PER-PIXEL
    % TIME SHIFTS.
    sunVolume = spatioTemporalVolume - skyVolume;
    for count = 1 : numel(sunVolume)
        
        if sunVolume(count) < 0
            
            sunVolume(count) = 0;
            
        end
        
    end
    
    [sunLightSlice, sunBasisSlice] = generate1DSlices(sunVolume, length, numberIlluminated);
    
    [sunlight, sunlightBasis, shift, ...
        shiftedBasis] = leastSquaresSun(sunVolume, shadowVolume, sunLightSlice, sunBasisSlice, length, numberIlluminated);
    
    % Recover the sun volume from the sunlight factor and shifted basis
    % curve.
    for count = 1 : numberIlluminated
    
        sunVolume(:, count) = sunlight(count, :) * shiftedBasis(:, count);
    
    end

end


% Compute 1D slices from the given light volume by grouping the data and
% sampling the required number of group means (numberIlluminated for light
% and length for basis) to form the vector (Lawrence et al, 2006). Note the
% differing orientations of the light and basis slices.
function [lightSlice, basisSlice] = generate1DSlices(lightVolume, length, numberIlluminated)

    lightVolume = lightVolume';

    [~, groupMeans] = kmeans(lightVolume, 25);
    groupLength = numel(groupMeans);
    groupMeans = reshape(groupMeans, [1, groupLength]);

    lightSlice = randsample(groupMeans, numberIlluminated, 1);
    lightSlice = reshape(lightSlice, [numberIlluminated, 1]);
    
    basisSlice = randsample(groupMeans, length, 1);
    basisSlice = reshape(basisSlice, [1, length]);
    
end


% Compute least squares for the skylight component.
function [light, basis, shift, shiftedBasis] = leastSquaresSky(lightVolume, confidence, light, basis, length, numberIlluminated)
    
    lightVolume = lightVolume';
    confidence = confidence';
    
    shift = zeros(numberIlluminated, 1);
    shiftedBasis = basis';
    
    light = leastSquaresLight(lightVolume, light, shiftedBasis, confidence, numberIlluminated);

    % Alternating constrained least squares for the basis component
    % (Lawrence et al., 2006). See leastSquaresLight() for explanation.
    for count = 1 : length
        
        thisConfidence = confidence(:, count);
        thisLightSlice = lightVolume(:, count);
        predictedBasis = double(light .* thisConfidence);
        predictedLight = double(thisLightSlice .* thisConfidence);
        basis(:, count) = lsqlin(predictedBasis, predictedLight);
        
    end
    
end


% Compute least squares for the sunlight component.
function [light, basis, shift, shiftedBasis] = leastSquaresSun(lightVolume, confidence, light, basis, length, numberIlluminated)

    lightVolume = lightVolume';
    confidence = confidence';
    
    % Compute the shifts: for each row of the sun volume (single pixel in
    % the groundTruth pixels over all frames), find the greatest difference
    % in cross-correlation between it and the basis curve. This minimises
    % the per-pixel error.
    shift = zeros(numberIlluminated, 1);
    for count = 1 : numberIlluminated
        
        [crossCorrelation, differenceMarkers] = xcorr(lightVolume(count, :), basis);
        
        largestDifference = abs(crossCorrelation(1));
        largestDifferenceLocation = 1;
        
        crossCorrelationSize = numel(crossCorrelation);
        
        if crossCorrelationSize > 1
            
            for correlationCount = 2 : crossCorrelationSize
                
                thisDifference = abs(crossCorrelation(correlationCount));
                
                if thisDifference > largestDifference
                    
                    largestDifference = thisDifference;
                    largestDifferenceLocation = correlationCount;
                    
                end
                
            end
            
        end
        
        shift(count) = differenceMarkers(largestDifferenceLocation);
        
    end
    
    
    % Limit the shifts and get locations in the light volume
    shiftedLocations = length * (1 : numberIlluminated) + shift';
    for illuminatedCount = 1 : numberIlluminated
        
        thisLocation = shiftedLocations(illuminatedCount);
        largestShiftLocation = max(thisLocation, (length * (illuminatedCount - 1) + 1));
        shiftedLocations(illuminatedCount) = min(largestShiftLocation, length * illuminatedCount);
        
    end
    
    thisShift = shift(count, :);
    shiftedBasis = zeros(1, length);
    
    if thisShift <= 0
        
        shiftLocations = 1 : thisShift + length;
        
    else
        
        shiftLocations = thisShift + 1 : length;
        
    end
    
    thisLength = numel(shiftLocations);
    shiftedBasis(shiftLocations) = basis(:, 1 : thisLength);
    shiftedBasis = shiftedBasis';
    
    
    light = leastSquaresLight(lightVolume, light, shiftedBasis, confidence, numberIlluminated);
    
    % Alternating constrained least squares for the basis component
    % (Lawrence et al., 2006). See leastSquaresLight() for explanation.
    thisLightSlice = lightVolume(shiftedLocations)';
    for count = 1 : length
        
        thisConfidence = confidence(:, count);
        predictedBasis = double(light .* thisConfidence);
        predictedLight = double(thisLightSlice .* thisConfidence);
        basis(:, count) = lsqlin(predictedBasis, predictedLight);
        
    end
    
    
    % Apply shifts to basis curve
    length = numel(basis);
    numberIlluminated = numel(shift);
    shiftedBasis = zeros(length, numberIlluminated);
    
    for count = 1 : numberIlluminated
        
        thisShift = shift(count, :);
        thisShiftedBasis = zeros(1, length);
        
        if thisShift >= 0

            shiftLocations = thisShift + 1 : length;

        else

            shiftLocations = 1 : thisShift + length;

        end
        
        thisLength = numel(shiftLocations);
        thisShiftedBasis(shiftLocations) = basis(:, 1 : thisLength);
        shiftedBasis(:, count) = thisShiftedBasis;
        
    end
    
end


% Alternating constrained least squares for the light component (Lawrence
% et al., 2006). We decompose with the second dimension of the weight
% matrix being 1, so this is reduced to a linear constrained least squares
% problem for each row of the matrix.
function light = leastSquaresLight(lightVolume, light, shiftedBasis, confidence, numberIlluminated)
    
    for count = 1 : numberIlluminated
        
        % Apply confidencence matrix: this is the shadow volume for the
        % sunlight component and the complement of it for the skylight
        % component. To apply, simply pre-multiply each component by the
        % confidence.
        thisConfidence = confidence(count, :)';
        predictedBasis = double(shiftedBasis .* thisConfidence);
        predictedLight = double(lightVolume(count, :)' .* thisConfidence);
        
        % Compute least squares. Here, we find a vector v such that the
        % norm of (1/2)*(predictedBasis * v - predictedLight) is a minimum.
        light(count, :) = lsqlin(predictedBasis, predictedLight)';
        
    end

end


% Reconstruct a given volume after factorisation by combining channels for
% pixels covered by the ground truth.
function reconstructedVolume = reconstructVolume(fullVolume, groundTruthImage, groundTruth, height, width, length, frameDimension)

    reconstructedVolume = zeros(height, width, length, 3);

    for channelCount = 1 : 3
        
        singleChannelVolume = zeros(height, width, length);
        
        for frameCount = 1 : length
            
            groundTruthImage(groundTruth') = fullVolume(frameCount, :, channelCount);
            singleChannelVolume(:, :, frameCount) = groundTruthImage;

        end
        
        reconstructedVolume(:, :, :, channelCount) = singleChannelVolume;
        
    end
    
    reconstructedVolume = reshape(reconstructedVolume, [frameDimension, length, 3]);
    reconstructedVolume = permute(reconstructedVolume, [2, 1, 3]);
    
end


% Reconstruct images in the same way as volumes.
function reconstructedImage = reconstructImage(groundTruth, redVolumeComponent, greenVolumeComponent, blueVolumeComponent, height, width)
    
    reconstructedImage = zeros(height, width, 3);
    [redChannel, greenChannel, blueChannel] = deal(zeros(height, width));
    
    redChannel(groundTruth) = redVolumeComponent(:);
    greenChannel(groundTruth) = greenVolumeComponent(:);
    blueChannel(groundTruth) = blueVolumeComponent(:);
    
    [reconstructedImage(:, :, 1), reconstructedImage(:, :, 2), reconstructedImage(:, :, 3)] = deal(redChannel, greenChannel, blueChannel);
    
end


% A simple bilateral filter for preserving the edges of the shadow map.
% This is adapted from Tomasi and Manduchi, 1998 and Banterle et al., 2012.
function smoothedImage = bilateralFilter(originalImage, height, width, diameter)

    smoothedImage = zeros(size(originalImage));
    windowPixelIntensities = zeros(4, 1);
    
    for heightCount = diameter + 1 : (height - diameter)
        
        for widthCount = diameter + 1: (width - diameter)
            
            originalPixel = originalImage(heightCount, widthCount);
            filteredIntensity = 0;
            
            for filterWindowCount = 1 : diameter
                
                % Get intensities around the window.
                windowPixelIntensities(1) = originalImage(heightCount - filterWindowCount, widthCount);                
                windowPixelIntensities(2) = originalImage(heightCount, widthCount - filterWindowCount);
                windowPixelIntensities(3) = originalImage(heightCount, widthCount + filterWindowCount);               
                windowPixelIntensities(4) = originalImage(heightCount + filterWindowCount, widthCount);
                             
                % Find intensity difference, square, add to the original
                % intensity and scale.
                squareDifferences = sum((windowPixelIntensities - double(originalPixel)).^2);
                filteredIntensity = filteredIntensity + 4 * (squareDifferences + originalImage(heightCount, widthCount)) / diameter;
                
            end
            
            smoothedImage(heightCount, widthCount) = filteredIntensity;
            
        end
        
    end
    
end