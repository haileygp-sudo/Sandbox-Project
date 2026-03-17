[sampleStruct, probStruct, Comments] = scfread('sample.scf');
%figure
%hold on
%plot(sampleStruct.A);
%plot(sampleStruct.C);
%plot(sampleStruct.G);
%plot(sampleStruct.T);
%legend('A','C','G','T');

newSampleStructA = removeNoise(sampleStruct.A);
newSampleStructC = removeNoise(sampleStruct.C);
newSampleStructG = removeNoise(sampleStruct.G);
newSampleStructT = removeNoise(sampleStruct.T);
figure;
hold on;
plot(newSampleStructA);
plot(newSampleStructC);
plot(newSampleStructG);
plot(newSampleStructT);
legend('A (Filtered)', 'C (Filtered)', 'G (Filtered)', 'T (Filtered)');

[prob,chain] = calculateProb(newSampleStructA, newSampleStructC, newSampleStructG, newSampleStructT);
disp(prob)
chainStr = "";
for (i = 1:length(chain))
    chainStr = strcat(chainStr, chain(i));
end
disp(chainStr);
disp(strlength(chainStr));
nCounter = 0;
for (i = 1:742)
    if (strcmp(probStruct.base(i), 'N') == 1)
        nCounter = nCounter + 1;
    end
end
disp(nCounter);
    


function newSampleStruct = removeNoise(sampleData)
    startPoint = 1; %creates the start point to cut the data from
    trimThreshold = mean(sampleData) * 0.1; %creates a threshold for "valuable" data

    for i = 1:length(sampleData) 
        if sampleData(i) > trimThreshold %if data is above the threshold to be considered "usefUl"
            startPoint = i; %we then update the starting point of the data
            break; %stop iterating because we have the starting point we want
        end
    end
    
    sampleData = sampleData(startPoint:end); %cut off unnecessary data at the beginning

    x = 1:length(sampleData); %creates our x values
    x = x'; %turns it into a row array

    thresholdValue = mean(sampleData) * 0.7; %creates a threshold for what values are considered "noise"
    noiseIndices = find(sampleData < thresholdValue); %finds the indices on each sampleData to find when the value is just noise
    noiseValues = sampleData(noiseIndices); %finds the actual value of the signal strength for that base pair

    p = polyfit(noiseIndices, noiseValues, 5); %creates parameters for the noise values of the sample (excludes all useful values)
    noiseThreshold = polyval(p, x); %creates a rough curve of what the noise is using parameters

    newSampleStruct = sampleData - noiseThreshold; %subtracts the "noise" from the actual data

    newSampleStruct(newSampleStruct < 0) = 0; %any negative values turn positive
end

function [baseProb, baseChain] = calculateProb(structA, structC, structG, structT)
    % Assume that the chromatogram data (A, C, G, T) and probabilities are given

% Downsample the chromatogram to the base sequence length
lengthA = length(structA);
lengthC = length(structC);
lengthG = length(structG);
lengthT = length(structT);
minLength = lengthA;
if (lengthC < minLength)
    minLength = lengthC;
end
if (lengthG < minLength)
    minLength = lengthG;
end
if (lengthT < minLength)
    minLength = lengthT;
end
newStructA = structA(1:minLength);
newStructC = structC(1:minLength);
newStructG = structG(1:minLength);
newStructT = structT(1:minLength);


baseSequenceLength = 677; % Example base sequence length
downsamplingFactor = minLength / baseSequenceLength;

% Initialize base sequence and probabilities
baseChain = cell(1, baseSequenceLength);
baseProb = zeros(baseSequenceLength, 4); % Columns for A, C, G, T probabilities

for i = 1:baseSequenceLength
    % Get the corresponding chromatogram indices for this base
    startIdx = round((i - 1) * downsamplingFactor) + 1;
    endIdx = round(i * downsamplingFactor);
    
    % Extract the intensities for each base in this window
    intensities = [sum(newStructA(startIdx:endIdx)), ...
                   sum(newStructC(startIdx:endIdx)), ...
                   sum(newStructG(startIdx:endIdx)), ...
                   sum(newStructT(startIdx:endIdx))];
    
    % Normalize the intensities to get the probabilities
    totalIntensity = sum(intensities);
    probabilities = intensities / totalIntensity;
    
    % Store the probabilities
    baseProb(i, :) = probabilities;
    
    % Determine the base with the highest probability
    [~, baseIdx] = max(probabilities);
    
    % Assign the base to the sequence
    if baseIdx == 1
        baseChain{i} = 'A';
    elseif baseIdx == 2
        baseChain{i} = 'C';
    elseif baseIdx == 3
        baseChain{i} = 'G';
    else
        baseChain{i} = 'T';
    end
end
end

% Display the result
